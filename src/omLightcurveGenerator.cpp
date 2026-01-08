#include "lightcurveGenerator.h"
#include "backupGenerator.h"
#include "astroFns.h"
#include "VBMicrolensingLibrary.h"
#include "singleLens.h"
#include "ephem.h"
#include "coords.h"
#include "argsort.h"
#include<time.h>
#include<vector>

#include<fstream>
#include<sstream>
#include<iomanip>
#include<numeric>
#include<sys/stat.h>
#define DEBUGVAR 1
// This assumes Event->vbm has already been initialized and configured


//General lightcurve generator that incorporates heirarchical orbital motion for multiple lenses and sources

void lightcurveGenerator(struct filekeywords* Paramfile, struct event *Event, struct obsfilekeywords World[], struct slcat *Sources, struct slcat *Lenses, ofstream& logfile_ptr)
{

  if(Paramfile->verbosity>=3)
    cout << "At lightcurveGenerator start, Tol=" 
	 << Event->vbm->Tol 
	 << ", RelTol=" 
	 << Event->vbm->RelTol 
	 << std::endl;


  //if the event is saturated in each band, no need to calculate the lightcurve
  if(Event->nepochs==0 || Event->allsat)
    {
      return;
    }
  vector<int> idxshift;
  int shiftedidx;
  const double timeout = Paramfile->lc_timeout;
  const bool enforce_timeout = (timeout > 0.0);
  time_t starttime = time(NULL);
  bool timed_out = false;
  int obsidx;
  coords c;
  double rho, u;
  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    idxshift.push_back(Event->nepochsvec[obsidx]);

  Event->vbm->SetMethod(VBMicrolensing::Method::Multipoly);

  Event->xsrc.clear();
  Event->ysrc.clear();
  Event->mu_src.clear();
  
  Event->xsrc.resize(Event->nsrc);
  Event->ysrc.resize(Event->nsrc);
  Event->mu_src.resize(Event->nsrc);
  for(int i=0; i<Event->nsrc; i++)
    {
      Event->xsrc[i].resize(Event->nepochs);
      Event->ysrc[i].resize(Event->nepochs);
      Event->mu_src[i].resize(Event->nepochs);
    }
  
  Event->xlens.resize(Event->nlens);
  Event->ylens.resize(Event->nlens);
  for(int i=0; i<Event->nlens; i++)
    {
      Event->xlens[i].resize(Event->nepochs);
      Event->ylens[i].resize(Event->nepochs);
    }
  
  //Conventions:
  //Track the apparent motion of the centers of mass of the source and lens, then compute offsets from them
  //for each component

  double antipode_ra = c.fold(Event->ra + PI,0,twoPi); //Used for computing orbits
  double antipode_dec = c.fold(-Event->dec,-PI,PI);
      
  vector<double> s_delta(3,0.0); //Offset of the chosen lens from its center of mass at tref
  vector<double> l_delta(3,0.0); //Offset of the chosen source from its center of mass at tref

  int sn = Event->source;
  int ln = Event->lens;

  double rEsrc = Event->rE * Sources->data[sn][Sources->DIST]/Lenses->data[ln][Lenses->DIST];

  double cosa = cos(Event->alpha);
  double sina = sin(Event->alpha);
  
  //Setup orbits for the source(s)

  vector<vector<orbitalElements> > s_elements; //Orbital elements classes for the sources

  int nsrc = 1;
  if(Paramfile->multiple_sources && Event->scompanions.size()>0)
    {
      nsrc = 1 + Event->scompanions.size();

      //Ultimately, Sort the source companions by orbit size, compute orbits heirarchically

      //For now, we are only going to deal with binary sources
      nsrc=2;

      if(Paramfile->verbosity>=1) cout << "Multiple source event, nsrc=" << nsrc << endl;

      s_elements.resize(nsrc);
      for(int i=0;i<nsrc;i++)
	{
	  s_elements[i].resize(nsrc-1);
	}

      int sc = Event->scompanions[0];
	  
      //Binary case, orbit about a barycenter
      //semimajor axis (AU), eccentricity,
      //double a, a0, da; //semimajor axis           (AU)
      //double e, e0, de; //eccentricity             (rad? not unitless?)
      //double I, I0, dI; //inclination              (deg)
      //double L, L0, dL; //mean longitude           (deg)
      //double w, w0, dw; //longitude of perihelion  (deg) \varomega
      //double O, O0, dO; //longitude of asc node    (deg) \Omega

      //Work in ecliptic coordinates so we can use all the tools of the class
      //double a1,e,I1,L1,w1,O1,dL1;
      //double a2,e,I2,L2,w2,O2,dL2;
      double acomb = (1.0+Event->scomp_q[0])*Event->scomp_a[0];
      double a1 = acomb-Event->scomp_a[0];

      //orbitalElements(double a, double e, double I, double L, double w, double O, double dL_, double epoch_=J2000)
      s_elements[1][0] = orbitalElements(Event->scomp_a[0], Event->scomp_e[0], Event->scomp_I[0], Event->scomp_L0[0], Event->scomp_w[0], Event->scomp_O[0], Event->scomp_dL[0]);
      s_elements[0][0] = orbitalElements(-a1, Event->scomp_e[0], Event->scomp_I[0], Event->scomp_L0[0], Event->scomp_w[0], Event->scomp_O[0], Event->scomp_dL[0]);

      

      //Compute the origin shift relative to the center of mass of the lens - normalization by rEsrc will be done later
      vector<double> xp;      
      for(int j=0;j<int(s_elements[0].size());j++)
	{
	  s_elements[0][j].viewfrom(Event->tref,antipode_ra,antipode_dec,&xp);
	  s_delta[0] += xp[0]; s_delta[1] += xp[1]; s_delta[2] += xp[2];
	}
      
    }
  //Hold the magnifications
  vector<double> mu(nsrc,0.0);

  //Setup orbits for the lens(es)

  vector<vector<orbitalElements> > l_elements; //Orbital elements classes for the lenses
  int nlens = 1 + Event->lcompanions.size() + Event->p_a.size();
  l_elements.clear();
  l_elements.resize(nlens);

  if(nlens==2)
    {
      if(Paramfile->multiple_lenses && Event->lcompanions.size()>0)
	{

	  if(Paramfile->verbosity>=1) cout << "Lens is a binary star, no planets" << endl;
	  
	  //the binary lens is a binary star
	  l_elements[0].resize(1);
	  l_elements[1].resize(1);
	  double acomb = (1.0+Event->lcomp_q[0])*Event->lcomp_a[0];
	  double a1 = acomb-Event->lcomp_a[0];

	  //orbitalElements(double a, double e, double I, double L, double w, double O, double dL_, double epoch_=J2000)
	  l_elements[1][0] = orbitalElements(Event->lcomp_a[0], Event->lcomp_e[0],
					     Event->lcomp_I[0], Event->lcomp_L0[0],
					     Event->lcomp_w[0], Event->lcomp_O[0],
					     Event->lcomp_dL[0]);
	  l_elements[0][0] = orbitalElements(-a1, Event->lcomp_e[0],
					     Event->lcomp_I[0], Event->lcomp_L0[0],
					     Event->lcomp_w[0], Event->lcomp_O[0],
					     Event->lcomp_dL[0]);

	  //Compute the origin shift relative to the center of mass of the lens - I don't think this is needed
	  //vector<double> xp;      
	  //for(int j=0;j<int(s_elements[0].size());j++)
	  //  {
	  //    s_elements[0][j].viewfrom(Event->tref,antipode_ra,antipode_dec,&xp);
	  //    s_delta[0] += xp[0]; s_delta[1] += xp[1]; s_delta[2] += xp[2];
	  //  }
	}
      else 
	{
	  if(Paramfile->verbosity>=1) cout << "Lens is a single star and one planet" << endl;
	  //we have a planet and a star
	  //the binary lens is a single star and planet
	  l_elements[0].resize(1); // the star
	  l_elements[1].resize(1); //the planet
	  double acomb = (1.0+Event->p_q[0])*Event->p_a[0];
	  double a1 = acomb-Event->p_a[0];

	  //orbitalElements(double a, double e, double I, double L, double w, double O, double dL_, double epoch_=J2000)
	  l_elements[1][0] = orbitalElements(Event->p_a[0], Event->p_e[0], Event->p_I[0], Event->p_L0[0], Event->p_w[0], Event->p_O[0], Event->p_dL[0]);
	  l_elements[0][0] = orbitalElements(-a1, Event->p_e[0], Event->p_I[0], Event->p_L0[0], Event->p_w[0], Event->p_O[0], Event->p_dL[0]);

	  //Compute the origin shift relative to the center of mass of the lens
	  //vector<double> xp;      
	  //for(int j=0;j<int(s_elements[0].size());j++)
	  //  {
	  //    s_elements[0][j].viewfrom(Event->tref,antipode_ra,antipode_dec,&xp);
	  //    s_delta[0] += xp[0]; s_delta[1] += xp[1]; s_delta[2] += xp[2];
	  //  }
	}
    }
  else
    {
      //We have to sort through multiple scenarios
      if(Paramfile->multiple_lenses && Event->lcompanions.size()>0)
	{
	  //We have planets/moons and a binary star
	  //Add the binary star to the planet orbit array, then we can sort it
	  for(int i=0;i<Event->lcompanions.size();i++)
	    {
	      Event->p_a.push_back(Event->lcomp_a[i]);
	      Event->p_e.push_back(Event->lcomp_e[i]);
	      Event->p_I.push_back(Event->lcomp_I[i]);
	      Event->p_L0.push_back(Event->lcomp_L0[i]);
	      Event->p_w.push_back(Event->lcomp_w[i]);
	      Event->p_O.push_back(Event->lcomp_O[i]);
	      Event->p_dL.push_back(Event->lcomp_dL[i]);
	      Event->p_q.push_back(Event->lcomp_q[i]);
	      Event->p_orbtype.push_back(-1); //to represent a binary star
	    }
	}
      else
	{
	  //We have a single star and planets/moons
	  //Nothing needed here
	  
	}

      //Determine the heirarchy of orbits
      auto orbsize_order = argsort(Event->p_a);
      if(Paramfile->verbosity>1)
	{
	  cout << "Orbit order:" << endl;
	  for(auto os : orbsize_order) cout << os << " " << Event->p_a[os] << endl;
	  
	}
      int moons = 0;
      int barycenters = 0;
      int circumbinary=0;
      int distantbinary=0;
      int mixedbinary=0;
      for(auto oc : Event->p_orbtype)
	{
	  if(oc==3) moons++;
	}

      //The possible systems
      //Binary star + planet(s)
      //   //P or S type
      //Binary star + planet(s) + moon(s)
      //   //P or S type
      //Single star + planet(s)
      //Single star + planet(s) + moon(s)

      if(Paramfile->multiple_lenses && Event->lcompanions.size()>0)
	{
	  //Figure out the type of binary star system we have
	  //The binary star companion has an orbtype of -1, if the last in the orbsize_order list is it, its a distant companion
	  if(Event->p_orbtype[orbsize_order.back()]==-1)
	    {
	      distantbinary=1;
	      if(Paramfile->verbosity>=1) cout << "Lens involves a distant binary star and at least one planet, nlens=" << Event->nlens << " nplanets=" << Event->nplanets << endl;
	    }
	  else
	    {
	      //The closest non-moon orbit is the stellar companion --> circumbinary
	      for(auto oso : orbsize_order)
		{
		  if(Event->p_orbtype[oso]==3) continue;
		  if(Event->p_orbtype[oso]==-1) circumbinary=1;
		  break;
		}
		  
	      if(circumbinary==1)
		{
		  if(Paramfile->verbosity>=1) cout << "Lens involves a close binary star and at least one circumbinary planet, nlens " << Event->nlens << " nplanets=" << Event->nplanets << endl;		  
		}
	      else
		{
		  mixedbinary=1; //the binary star is between two planets, this is quite possibly extremely unphysical
		  if(Paramfile->verbosity>=1) cout << "Lens is classified as a mixed binary that might be unphysical/unstable, nlens " << Event->nlens << " nplanets=" << Event->nplanets << endl;
		}
	    }	 	  
	}

      //Now construct orbits
      
      //For the moon-planet system we need a barycenter that will orbit the star,
      //then the moon and planet will orbit the barycenter

      double qsys;

      if(moons>0)
	{
	  double mbary = Event->p_mass[0]; //ratio of the planet+moons system relative to the planet (for now)
	  double qbary = 1.0;
	  double q;
	  
	  if(Paramfile->verbosity>=1) cout << "Lens involves a planet with moons, nlens=" << nlens << " nmoons=" << moons << endl;
	  for(auto idx : orbsize_order)
	    {
	      //Only one planet allowed with moons, the planet will be the first in the planet list [0], the second in the orbit list [1]
	      //Work through moon orbits in order of orbit size so the more distant ones deal with a central barycenter and total inner mass
	      
	      if(Event->p_orbtype[idx]!=3) continue; //skip non-moons
	      
	      q = Event->p_mass[idx]/mbary; //ratio of moon mass to all internal mass
	      double acomb = (qbary+q)*Event->p_a[idx];
	      double a1 = acomb-Event->p_a[idx];
	      double pfix = sqrt(qbary+q);
	      Event->p_period[idx] /= pfix;
	      l_elements[idx+1].push_back(orbitalElements(Event->p_a[idx], Event->p_e[idx], Event->p_I[idx], Event->p_L0[idx], Event->p_w[idx], Event->p_O[idx], Event->p_dL[idx]*pfix));
	      //reflex orbit of the planet due to the moons
	      l_elements[1].push_back(orbitalElements(-a1, Event->p_e[idx], Event->p_I[idx], Event->p_L0[idx], Event->p_w[idx], Event->p_O[idx], Event->p_dL[idx]*pfix));
	      for(auto jdx : orbsize_order)
		{
		  //reflex orbit of any inner moons due to the current one
		  if(jdx==idx) break;
		  if(Event->p_orbtype[jdx]==3)
		    {
		      l_elements[jdx+1].push_back(orbitalElements(-a1, Event->p_e[idx], Event->p_I[idx], Event->p_L0[idx], Event->p_w[idx], Event->p_O[idx], Event->p_dL[idx]*pfix));
		    }
		}
	      mbary += Event->p_mass[idx];
	      qbary += q;
	      
	    } //end loop over moons
	  qsys = qbary*Event->p_q[0]; //mass ratio of the moon system to the first star
	} //end if moons
      else qsys=1.0;

      //We've added the moons, now work through the other bodies, treating the planet as the combined mass of it with its moons

      //Start with the first star
      double qbary=1;
      double q;
      int planet_yet=0;
      
      for(auto idx : orbsize_order)
	{
	  if(Event->p_orbtype[idx]==3) continue; //ignore the moons, they are orbiting the planet
	  
	  if(moons>0 && Event->p_orbtype[idx]!=-1)
	    {
	      //if there are moons, there is only one planet
	      q = qsys; //use the mass of the planet moon system
	      planet_yet=1; //modify the moon system with the reflex orbit
	    }
	  else q = Event->p_q[idx];
		  
	  double acomb = (qbary+q)*Event->p_a[idx];
	  double a1 = acomb-Event->p_a[idx];
	  double pfix = sqrt(qbary+q);
	  Event->p_period[idx] /= pfix;
	  l_elements[idx+1].push_back(orbitalElements(Event->p_a[idx], Event->p_e[idx], Event->p_I[idx], Event->p_L0[idx], Event->p_w[idx], Event->p_O[idx], Event->p_dL[idx]*pfix));
	  //reflex motion of the main star
	  l_elements[0].push_back(orbitalElements(-a1, Event->p_e[idx], Event->p_I[idx], Event->p_L0[idx], Event->p_w[idx], Event->p_O[idx], Event->p_dL[idx]*pfix));
	  int skip_unless_moon=0;
	  //add the reflex motion to any other bodies inside this one's orbit
	  for(auto jdx : orbsize_order)
	    {
	      if(jdx==idx) //we've reached this object, the rest can be skipped unless they are a moon
		{
		  skip_unless_moon=1;
		  continue;
		}
	      //skip if it is a moon, unless the planet is inside this idx object
	      if(Event->p_orbtype[jdx]!=3 || planet_yet==1)
		{
		  //this works because the above continue will skip it for the first planet_yet==1 which is the planet itself
		  l_elements[jdx+1].push_back(orbitalElements(-a1, Event->p_e[idx], Event->p_I[idx], Event->p_L0[idx], Event->p_w[idx], Event->p_O[idx], Event->p_dL[idx]*pfix));
		}
	    }
	  qbary += q;
		      
		
	} //end loop over other bodies

    } //end if >=2 lenses


  double* lens_parameters = new double[nlens*3];
  lens_parameters[2] = 1.0;
  for(int i=1;i<nlens;i++)
    {
      //cout << nlens << " " << i << " " << 3*i+2 << " " << Event->p_q.size() << endl;
      lens_parameters[3*i+2] = Event->p_q[i-1]; //Event->lcomp_q[i];
    }

  //Check the orbital elements arrays are right
  if(Paramfile->verbosity>=2)
    {
      cout << "l_elements.size()=" << l_elements.size();
      for(auto lel : l_elements)
	{
	  cout << " " << lel.size();
	}
      cout << endl;
      for(auto lorb : l_elements)
	{
	  for(auto lel : lorb) lel.print_elements();
	}
      cout << "s_elements.size()=" << s_elements.size();
      for(auto sel : s_elements)
	{
	  cout << " " << sel.size();
	}
      cout << endl;
      for(auto sorb : s_elements)
	{
	  for(auto sel : sorb) sel.print_elements();
	}
    }

  //Calculate the lightcurve
  for(int idx=0; idx<Event->nepochs; idx++)
    {
      if(enforce_timeout)
	{
	  time_t now = time(NULL);
	  if (difftime(now, starttime) > timeout)
	    {
	      timed_out = true;
	      break;
	    }
	}
      obsidx = Event->obsidx[idx];
      shiftedidx=idx-idxshift[obsidx];
      int filt = World[obsidx].filter;
      double amp = 0.0;

      //Compute parallax shift of the source barycenter
      double tt = (Event->epoch[idx] - Event->t0) / Event->tE_r;
      double uu = Event->u0;
      
      if (Paramfile->pllxMultiplyer)
	{
	  tt += Event->pllx[obsidx].tshift[shiftedidx];
	  uu += Event->pllx[obsidx].ushift[shiftedidx];
	}
      Event->umin=min(Event->umin,qAdd(tt,uu));

      //Compute the location of the reference point of the source relative to the reference point of the lens
      double xs0 = tt * cosa - uu * sina;
      double ys0 = tt * sina + uu * cosa;

      vector<double> xs(nsrc,0.0); //source position in the plane of the sky, ecliptic sky coordinates in AU
      vector<double> ys(nsrc,0.0);
      vector<double> ds(nsrc,0.0); //source distance

      if(nsrc>1)
	{
	  for(int i=0;i<nsrc;i++)
	    {
	      vector<double> xp;
	      for(int j=0;j<int(s_elements[i].size());j++)
		{
		  s_elements[i][j].viewfrom(Event->jdtimes[obsidx][idx],antipode_ra,antipode_dec,&xp);
		  xs[i] += xp[0]; ys[i] += xp[1]; ds[i] += xp[2];
		}
	      xs[i] -= s_delta[0]; ys[i] -= s_delta[1]; ds[i] -= s_delta[2];
	      xs[i] /= rEsrc; ys[i] /= rEsrc; ds[i] /= rEsrc;
	      //Rotation needed here?
	      
	      xs[i] += xs0; ys[i] += ys0;
	    
	      Event->xsrc[i][idx] = xs[i];
	      Event->ysrc[i][idx] = ys[i];
	    }
	}
      else
	{
	  xs[0] = xs0; ys[0] = ys0;
	  Event->xsrc[0][idx] = xs0; Event->ysrc[0][idx] = ys0; 
	}


      vector<double> xl(nlens,0.0); //lens position in the plane of the sky, ecliptic sky coordinates in AU
      vector<double> yl(nlens,0.0);
      vector<double> dl(nlens,0.0); //lens distance 

      if(nlens>1)
	{
	  for(int i=0;i<nlens;i++)
	    {
	      vector<double> xp;      
	      for(int j=0;j<int(l_elements[i].size());j++)
		{
		  if(Paramfile->verbosity>=3) l_elements[i][j].print_elements();
		  l_elements[i][j].viewfrom(Event->jdtimes[obsidx][idx],antipode_ra,antipode_dec,&xp);
		  if(Paramfile->verbosity>=3) cout << setprecision(16) << i << " " << j << " " << Event->jdtimes[obsidx][idx] << " " << xp[0] << " " << xp[1] << " " << xp[2] << endl;
		  xl[i] += xp[0]; yl[i] += xp[1]; dl[i] += xp[2];
		}
	      xl[i] -= l_delta[0]; yl[i] -= l_delta[1]; dl[i] -= l_delta[2];


	      xl[i] /= Event->rE; yl[i] /= Event->rE; dl[i] /= Event->rE;

	      //rotations needed here?
	      
	      lens_parameters[3*i+0] = xl[i];
	      lens_parameters[3*i+1] = yl[i];
	      Event->xlens[i][idx] = xl[i];
	      Event->ylens[i][idx] = yl[i];
	    }

	  Event->vbm->SetLensGeometry(nlens,lens_parameters);
	}
      else
	{
	  xl[0] = 0.0; yl[0] = 0.0;
	  Event->xlens[0][idx] = 0.0; Event->ylens[0][idx] = 0.0; 
	}

      

      //Finally ready to compute magnifications
      
      if(nlens==1)
	{
	  for(int is=0;is<nsrc;is++)
	    {
      	      if(is==0) rho = Event->rs;	
	      else rho = Event->scomp_rs[is-1];
	      u = qAdd(xs[is],ys[is]);
	      mu[is] = Event->vbm->ESPLMag2(u, rho);
	      //handle astrometry
	    }
	}
      else if(nlens==2)
	{
	  double s = qAdd(xl[1]-xl[0],yl[1]-yl[0]);
	  double q = Event->p_q[0];
	  double rot = atan2(yl[1],xl[1]);
	  double cr = cos(-rot); double sr = sin(-rot);

	  for(int is=0;is<nsrc;is++)
	    {
	      if(is==0) rho = Event->rs;	
	      else rho = Event->scomp_rs[is-1];
	      //rotate coordintates to binary axis
	      //Binary mag works from the center of mass, so translate source to CoM, then rotate
	      double xs_com_ecl = xs[is] + l_delta[0];
	      double ys_com_ecl = ys[is] + l_delta[1];
	      double xsi = cr*xs_com_ecl - sr*ys_com_ecl;
	      double ysi = sr*xs_com_ecl + cr*ys_com_ecl;
	      
	      mu[is] = Event->vbm->BinaryMag2(s,q,xsi, ysi, rho);
	      //handle astrometry
	      //rotate astrometry back
	    }
	}
      else
	{	 
	  for(int is=0;is<nsrc;is++)
	    {
	      if(is==0) rho = Event->rs;
	      else rho = Event->scomp_rs[is-1];
	      mu[is] = Event->vbm->MultiMag2(xs[is], ys[is], rho);
	      //handle astrometry
	    }
	}

      //Here Atrue is magnification, but later it gets converted into fractional flux
      Event->Atrue[idx] = mu[0];
      for(int is=1;is<nsrc;is++)
	{
	  Event->Atrue[idx] += Event->scomp_fsofs1[is-1][filt] * (mu[is]-1.0);
	}


      // Keep track of highest magnification
      if (Event->Atrue[idx] > Event->Amax)
	{
	  Event->Amax = amp;
	  Event->peakpoint = idx;
	}
      
    } //End for epoch

  delete [] lens_parameters;
  
  if(timed_out)
    {
      Event->lcerror = LCGEN_TIMEOUT_ERR;
      Event->detected = 0;
      Event->deterror = 0;
      if(Paramfile->verbosity >= 1)
	{
	  cout << "Lightcurve generation timed out after " << timeout << " seconds" << endl;
	}
      if(logfile_ptr.good())
	{
	  logfile_ptr << "Lightcurve generation timed out after " << timeout << " seconds" << std::endl;
	}
      return;
    }


      
}
