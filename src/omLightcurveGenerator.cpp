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
#include<cmath>
//#include<thread>
//#include<chrono>


#define DEBUGVAR 1
// This assumes Event->vbm has already been initialized and configured

//General lightcurve generator that incorporates heirarchical orbital motion for multiple lenses and sources

//void sleep_thread(int n)
//{
//  this_thread::sleep_for(chrono::milliseconds(n));
//}

void lightcurveGenerator(struct filekeywords* Paramfile, struct event *Event, struct obsfilekeywords World[], struct slcat *Sources, struct slcat *Lenses, ofstream& logfile_ptr)
{

  if(Paramfile->verbosity>=1) cout << "lcgen Event->nlens: " << Event->nlens << endl;
  cout << "Skip magnification = " << Paramfile->skip_magnification << endl;
  
  if(Paramfile->verbosity>=3)
    cout << "At lightcurveGenerator start, Tol=" 
	 << Event->vbm->Tol 
	 << ", RelTol=" 
	 << Event->vbm->RelTol 
	 << endl;

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

  //if(Event->nlens<4) Event->vbm->SetMethod(VBMicrolensing::Method::Multipoly);
  //else Event->vbm->SetMethod(VBMicrolensing::Method::Nopoly);
  Event->vbm->SetMethod(VBMicrolensing::Method::Nopoly);
  Event->vbm->a1 = Event->gamma;

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
  int sc = -1;
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

      sc = Event->scompanions[0];
	  
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
      //double acomb = (1.0+Event->scomp_q[0])*Event->scomp_a[0];
      //double a1 = acomb-Event->scomp_a[0];
      double acomb = Event->scomp_a[0];
      double a1 = 1.0/(1.0+Event->scomp_q[0])*Event->scomp_a[0];
      double a2 = Event->scomp_q[0]/(1.0+Event->scomp_q[0])*Event->scomp_a[0];

      //orbitalElements(double a, double e, double I, double L, double w, double O, double dL_, double epoch_=J2000)
      s_elements[1][0] = orbitalElements(a1, Event->scomp_e[0], Event->scomp_I[0], Event->scomp_L0[0], Event->scomp_w[0], Event->scomp_O[0], Event->scomp_dL[0]);
      s_elements[0][0] = orbitalElements(-a2, Event->scomp_e[0], Event->scomp_I[0], Event->scomp_L0[0], Event->scomp_w[0], Event->scomp_O[0], Event->scomp_dL[0]);

      

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
  int nlens = Event->nlens; //1 + Event->lcompanions.size() + Event->p_a.size();
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
	  //double acomb = (1.0+Event->lcomp_q[0])*Event->lcomp_a[0];
	  //double a1 = acomb-Event->lcomp_a[0];
	  double acomb = Event->lcomp_a[0];
	  double a1 = 1.0/(1.0+Event->lcomp_q[0])*Event->lcomp_a[0];
	  double a2 = Event->lcomp_q[0]/(1.0+Event->lcomp_q[0])*Event->lcomp_a[0];


	  //orbitalElements(double a, double e, double I, double L, double w, double O, double dL_, double epoch_=J2000)
	  l_elements[1][0] = orbitalElements(a1, Event->lcomp_e[0],
					     Event->lcomp_I[0], Event->lcomp_L0[0],
					     Event->lcomp_w[0], Event->lcomp_O[0],
					     Event->lcomp_dL[0]);
	  l_elements[0][0] = orbitalElements(-a2, Event->lcomp_e[0],
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
	  //double acomb = (1.0+Event->p_q[0])*Event->p_a[0];
	  //double a1 = acomb-Event->p_a[0];
	  double mbary=Lenses->data[ln][Lenses->MASS];
	  double m = Event->p_mass[0];
	  double q = m/mbary;
	  double acomb = Event->p_a[0];
	  double a1 = 1.0/(1.0+q)*Event->p_a[0];
	  double a2 = q/(1.0+q)*Event->p_a[0];
	  mbary+=m;
	  Event->p_period[0] = sqrt(cube(acomb)/mbary);
	  Event->p_dL[0] = 360.0/Event->p_period[0];


	  //orbitalElements(double a, double e, double I, double L, double w, double O, double dL_, double epoch_=J2000)
	  l_elements[1][0] = orbitalElements(a1, Event->p_e[0], Event->p_I[0], Event->p_L0[0], Event->p_w[0], Event->p_O[0], Event->p_dL[0]);
	  l_elements[0][0] = orbitalElements(-a2, Event->p_e[0], Event->p_I[0], Event->p_L0[0], Event->p_w[0], Event->p_O[0], Event->p_dL[0]);

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
	  if(Paramfile->verbosity>=2) cout << "lcompanions.size()=" << Event->lcompanions.size() << endl;
	  for(int i=0;i<Event->lcompanions.size();i++)
	    {
	      Event->p_mass.push_back(Event->lcomp_mass[i]);
	      Event->p_period.push_back(Event->lcomp_period[i]);
	      Event->p_a.push_back(Event->lcomp_a[i]);
	      Event->p_e.push_back(Event->lcomp_e[i]);
	      Event->p_I.push_back(Event->lcomp_I[i]);
	      Event->p_L0.push_back(Event->lcomp_L0[i]);
	      Event->p_w.push_back(Event->lcomp_w[i]);
	      Event->p_O.push_back(Event->lcomp_O[i]);
	      Event->p_dL.push_back(Event->lcomp_dL[i]);
	      Event->p_q.push_back(Event->lcomp_q[i]);
	      Event->p_s0.push_back(0.0);
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
      Event->moons = 0;
      int barycenters = 0;
      Event->circumbinary=0;
      Event->distantbinary=0;
      Event->mixedbinary=0;
      for(auto oc : Event->p_orbtype)
	{
	  if(oc==3) Event->moons++;
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
	      Event->distantbinary=1;
	      if(Paramfile->verbosity>=1) cout << "Lens involves a distant binary star and at least one planet, nlens=" << Event->nlens << " nplanets=" << Event->nplanets << endl;
	    }
	  else
	    {
	      //The closest non-moon orbit is the stellar companion --> circumbinary
	      for(auto oso : orbsize_order)
		{
		  if(Event->p_orbtype[oso]==3) continue;
		  if(Event->p_orbtype[oso]==-1) Event->circumbinary=1;
		  break;
		}
		  
	      if(Event->circumbinary==1)
		{
		  if(Paramfile->verbosity>=1) cout << "Lens involves a close binary star and at least one circumbinary planet, nlens " << Event->nlens << " nplanets=" << Event->nplanets << endl;		  
		}
	      else
		{
		  Event->mixedbinary=1; //the binary star is between two planets, this is quite possibly extremely unphysical
		  if(Paramfile->verbosity>=1) cout << "Lens is classified as a mixed binary that might be unphysical/unstable, nlens " << Event->nlens << " nplanets=" << Event->nplanets << endl;
		}
	    }	 	  
	}

      //Now construct orbits
      
      //For the moon-planet system we need a barycenter that will orbit the star,
      //then the moon and planet will orbit the barycenter

      double msysmoons;

      if(Event->moons>0)
	{
	  double mbary = Event->p_mass[0]; //ratio of the planet+moons system relative to the planet (for now)
	  double q;
	  
	  if(Paramfile->verbosity>=1) cout << "Lens involves a planet with moons, nlens=" << nlens << " nmoons=" << Event->moons << endl;
	  for(auto idx : orbsize_order)
	    {
	      //Only one planet allowed with moons, the planet will be the first in the planet list [0], the second in the orbit list [1]
	      //Work through moon orbits in order of orbit size so the more distant ones deal with a central barycenter and total inner mass
	      
	      if(Event->p_orbtype[idx]!=3) continue; //skip non-moons

	      if(Paramfile->verbosity>=1) cout << "Adding moon idx=" << idx << " to planet 1 " << endl;
	      
	      q = Event->p_mass[idx]/mbary; //ratio of moon mass to all internal mass
	      double acomb = Event->p_a[idx];
	      double a1 = 1.0/(1.0+q)*acomb;
	      double a2 = q/(1.0+q)*acomb;

	      mbary += Event->p_mass[idx];
	      Event->p_period[idx] = sqrt(cube(acomb)/mbary);
	      Event->p_dL[idx] = 360.0/Event->p_period[idx];
	      if(Paramfile->verbosity>=2) cout << "Orbit moon idx=" << idx << " a1=" << a1 << " acomb=" << acomb << endl;
	      l_elements[idx+1].push_back(orbitalElements(a1, Event->p_e[idx], Event->p_I[idx], Event->p_L0[idx], Event->p_w[idx], Event->p_O[idx], Event->p_dL[idx]));
	      //reflex orbit of the planet due to the moons
	      if(Paramfile->verbosity>=2) cout << "Reflex planet jdx=" << 1 << " object idx=" << idx << " a2=" << a2 << endl;
	      l_elements[1].push_back(orbitalElements(-a2, Event->p_e[idx], Event->p_I[idx], Event->p_L0[idx], Event->p_w[idx], Event->p_O[idx], Event->p_dL[idx]));
	      for(auto jdx : orbsize_order)
		{
		  //reflex orbit of any inner moons due to the current one
		  if(jdx==idx) break;
		  if(Event->p_orbtype[jdx]==3)
		    {
		      if(Paramfile->verbosity>=2) cout << "Reflex moon jdx=" << jdx << " object idx=" << idx << " a2=" << a2 << endl;
		      l_elements[jdx+1].push_back(orbitalElements(-a2, Event->p_e[idx], Event->p_I[idx], Event->p_L0[idx], Event->p_w[idx], Event->p_O[idx], Event->p_dL[idx]));
		    }
		}
	      
	    } //end loop over moons
	  msysmoons = mbary;
	} //end if moons


      //We've added the moons, now work through the other bodies, treating the planet as the combined mass of it with its moons

      //Start with the first star
      double mbary=Lenses->data[ln][Lenses->MASS];
      double q,m;
      int planet_moon=0;
      
      for(auto idx : orbsize_order)
	{
	  if(Event->p_orbtype[idx]==3) continue; //ignore the moons, they are orbiting the planet and motion is accounted for in planet barycenter, we'll shift them when we shift the planet
	  
	  if(Event->moons>0 && Event->p_orbtype[idx]!=-1)
	    {
	      //if there are moons, there is only one planet
	      m = msysmoons; //use the mass of the planet moon system
	      planet_moon=1; //modify the moon system with the reflex orbit
	    }
	  else
	    {
	      m = Event->p_mass[idx];
	    }

	  q = m/mbary; //ratio of moon mass to all internal mass
	  double acomb = Event->p_a[idx];
	  double a1 = 1.0/(1.0+q)*acomb;
	  double a2 = q/(1.0+q)*acomb;

	  if(Paramfile->verbosity>1) cout << "Planet idx=" << idx << " acomb=" << acomb << " a1=" << a1 << " abary=" << a2 << " m1=" << m << " mbary=" << mbary << endl;

	  mbary += m;
	  Event->p_period[idx] = sqrt(cube(acomb)/mbary);
	  Event->p_dL[idx] = 360.0/Event->p_period[idx];
		  	    
	  l_elements[idx+1].push_back(orbitalElements(a1, Event->p_e[idx], Event->p_I[idx], Event->p_L0[idx], Event->p_w[idx], Event->p_O[idx], Event->p_dL[idx]));
	  //reflex motion of the main star
	  if(Paramfile->verbosity>1) cout << "Star reflex" << endl;
	  l_elements[0].push_back(orbitalElements(-a2, Event->p_e[idx], Event->p_I[idx], Event->p_L0[idx], Event->p_w[idx], Event->p_O[idx], Event->p_dL[idx]));

	  //If we are dealing with the planet, add the motion to the moons in the system too
	  if(planet_moon==1)
	    {
	      for(auto jdx : orbsize_order)
		{
		  if(Event->p_orbtype[jdx]==3)
		    {
		      if(Paramfile->verbosity>1) cout << "Moon jdx=" << jdx << " with planet idx=" << idx << " +a=" << a1 << endl;
		      l_elements[jdx+1].push_back(orbitalElements(a1, Event->p_e[idx], Event->p_I[idx], Event->p_L0[idx], Event->p_w[idx], Event->p_O[idx], Event->p_dL[idx]));
		    }
		}
	      planet_moon=2; //We've dealt with the system not
	    }
	  
	  //add the reflex motion to any other bodies inside this one's orbit
	  int outside=0;
	  for(auto jdx : orbsize_order)
	    {
	      if(jdx==idx)
		{
		  //we've reached this object, the rest can be skipped unless they are a moon - they are outside this obj
		  break; 
		}

	      if(Event->p_orbtype[jdx]!=3)
		{
		  if(Paramfile->verbosity>1 && Event->p_orbtype[jdx]!=3) cout << "Planet reflex jdx=" << jdx << " with planet idx=" << idx << " +a=" << Event->p_a[idx] << endl;
		  
		  l_elements[jdx+1].push_back(orbitalElements(-a2, Event->p_e[idx], Event->p_I[idx], Event->p_L0[idx], Event->p_w[idx], Event->p_O[idx], Event->p_dL[idx]));
	      
		  if(jdx==1)
		    {
		      //This is the planet, move its moons as well
		      for(auto kdx : orbsize_order)
			{
			  if(Event->p_orbtype[kdx]==3)
			    {
			      if(Paramfile->verbosity>1) cout << "Moon reflex kdx=" << jdx << " with planet jdx=" << idx << " +a=" << Event->p_a[idx] << endl;
		  
			      l_elements[kdx+1].push_back(orbitalElements(-a2, Event->p_e[idx], Event->p_I[idx], Event->p_L0[idx], Event->p_w[idx], Event->p_O[idx], Event->p_dL[idx]));
			    }
			}
		    }
		  
	      //skip if it is a moon whose planet hasn't been reached yet
		}
		
	    }		      
		
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
      int count=0;
      for(auto lorb : l_elements)
	{
	  cout << "Lens " << count << ":" << endl;
	  for(auto lel : lorb) lel.print_elements();
	  count++;
	}
      
      cout << "s_elements.size()=" << s_elements.size();
      for(auto sel : s_elements)
	{
	  cout << " " << sel.size();
	}
      cout << endl;
      count=0;
      for(auto sorb : s_elements)
	{
	  cout << "Source " << count << ":" << endl;
	  for(auto sel : sorb) sel.print_elements();
	  count++;
	}
    }

  double time_elapsed=0;

  //Calculate the lightcurve
  double last_progress=0;

  vector<double> msource(nsrc);


  //Throw away some events which we know to not be realistic or to stretch the
  //capabilities of the MultiLens generator

  int bad_scenario=0;

  if(nlens>=4) bad_scenario+=iPow(2,0);
  if(Event->mixedbinary>0) bad_scenario+=iPow(2,1);
  else if(Event->circumbinary>0)
    {
      double e = Event->p_e[nlens-1];
      double mu_hw = Event->p_q[nlens-1];
      double acrit = 1.60 + 5.10*e - 2.22*sqr(e) + 4.12*mu_hw - 4.27*e*mu_hw
	- 5.09*sqr(mu_hw) + 4.61*sqr(e*mu_hw); //Holman & Wiegert (1999)
      cout << "circum binary acrit = " << acrit << " " << Event->p_a[0] << " " << Event->p_a[nlens-1]*acrit << endl;
      if(Event->p_a[0]<Event->p_a[nlens-1]*acrit) bad_scenario+=iPow(2,2);
      else
	{
	  double period_ratio = Event->p_period[0]/Event->p_period[nlens-1];
	  double prround = round(period_ratio);
	  if(prround<=9 && (fmod(period_ratio,prround)<0.02 || fmod(period_ratio,prround)>0.98))
	    bad_scenario+=iPow(2,3);
	}
    }
  else if(Event->distantbinary>0)
    {
      double e = Event->p_e[nlens-1];
      double mu_hw = Event->p_q[nlens-1];
      double acrit = 0.464 - 0.380*mu_hw - 0.631*e + 0.586*mu_hw*e
	+ 0.150*sqr(e) - 0.198*mu_hw*sqrt(e); //Holman & Wiegert (1999)
      cout << "distant binary acrit = " << acrit << " " << Event->p_a[0] << " " << Event->p_a[nlens-1]*acrit << endl;

      if(Event->p_a[0]>Event->p_a[nlens-1]*acrit) bad_scenario+=iPow(2,4);
    }
  

  if(bad_scenario>0)
    {
      Event->lcerror = 9000 + bad_scenario;
      Event->detected = 0;
      Event->deterror = 0;
      if(Paramfile->verbosity >= 1)
	{
	  cout << "Lightcurve generation skipped due to a bad scenario, event " << Event->id << ", code " << bad_scenario << "" << endl;
	}
      if(logfile_ptr.good())
	{
	  logfile_ptr << "Lightcurve generation skipped due to a bad scenario, event " << Event->id << ", code " << bad_scenario << "" << endl;
	}
      return;
    }


  
  
  for(int idx=0; idx<Event->nepochs; idx++)
    {
      
      double progress = (double(idx)/double(Event->nepochs)*100.0);
      if(Paramfile->verbosity>=1 && floor(progress/10)!=floor(last_progress/10)) cout << "." << flush;
      last_progress=progress;
      if(enforce_timeout)
	{
	  time_t now = time(NULL);
	  time_elapsed = difftime(now, starttime);
	  if (time_elapsed > timeout)
	    {
	      timed_out = true;
	      break;
	    }
	}
      obsidx = Event->obsidx[idx];
      shiftedidx=idx-idxshift[obsidx];
      int filt = World[obsidx].filter;
      double amp = 0.0;

      msource[0] = Sources->mags[sn][filt];
      if(Paramfile->multiple_sources && Event->scompanions.size()>0)
	msource[1] = Sources->mags[sc][filt];

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
	      if(Paramfile->skip_magnification==0)
		mu[is] = Event->vbm->ESPLMag2(u, rho);
	      else mu[is]=1.0;
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

	      if(Paramfile->skip_magnification==0)
		mu[is] = Event->vbm->BinaryMag2(s,q,xsi, ysi, rho);
	      else mu[is] = 1.0;
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
	      if(Paramfile->skip_magnification==0)
		{
		  int check=0;
		  logfile_ptr.precision(16);
		  logfile_ptr << Event->id << " " << Event->epoch[idx] << " ";
		  for(int ilp=0;ilp<nlens*3;ilp++)
		    {
		      logfile_ptr << lens_parameters[ilp] << " ";
		      if(!isfinite(lens_parameters[ilp])) check++;
		    }
		  logfile_ptr << xs[is] << " " << ys[is] << " " << rho << endl;
		  if(!isfinite(xs[is])) check++;
		  if(!isfinite(ys[is])) check++;
		  if(!isfinite(rho)) check++;
		      

		  if(check>0)
		    {
		      Event->lcerror = LCGEN_INPUT_ERR;
		      Event->detected = 0;
		      Event->deterror = 0;
		      if(Paramfile->verbosity >= 1)
			{
			  cout << "Lightcurve generation halted due to bad inputs, event " << Event->id << ", see logfile for parameters." << endl;
			}
		      if(logfile_ptr.good())
			{
			  logfile_ptr << Event->id << "Lightcurve generation halted due to bad inputs" << endl;
			}
		      return;
		    }
		  

		  //thread time_thread(&sleep_thread,1000);

		  //Spool up a new vbm for each calculation
		  VBMicrolensing VBMlocal;

		  VBMlocal.SetLensGeometry(nlens,lens_parameters);
		  VBMlocal.a1 = Event->gamma;
		  VBMlocal.Tol=Paramfile->vbm_tol;
		  VBMlocal.RelTol=Paramfile->vbm_reltol;
		  VBMlocal.SetMethod(VBMicrolensing::Method::Nopoly);

		  double u_min=1e50;
		  for(int i=0;i<nlens;i++)
		    {
		      double u_lens = qAdd(xs[is]-lens_parameters[i*3+0],ys[is]-lens_parameters[i*3+1])/sqrt(lens_parameters[i*3+2]);
		      if(u_lens<u_min) u_min=u_lens;
		    }

		  if(u_min<10 && msource[is]<40)
		    {
		      mu[is] = VBMlocal.MultiMag2(xs[is], ys[is], rho);
		      logfile_ptr << mu[is] << " " << VBMlocal.therr << " " << VBMlocal.NPS << endl;
		    }
		  else
		    {
		      //if the source is far from all lenses, just use the nearest single lens magnification
		      mu[is] = Event->vbm->ESPLMag2(u_min, rho);
		      logfile_ptr << mu[is] << " " << (u_min>=10?"single":"null") << " " << (msource[is]>=40?"faint":"null") << endl;
		    }
		}
	      else mu[is] = 1.0;
	      //handle astrometry
	    }
	}

      for(int is=0;is<nsrc;is++)
	{
	  Event->mu_src[is][idx] = mu[is];
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

  if(Paramfile->verbosity>=2)
    {
      cout << "Lightcurve took " << time_elapsed << " seconds to generate" << endl;
    }
  
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
