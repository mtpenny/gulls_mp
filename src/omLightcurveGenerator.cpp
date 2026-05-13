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
#include<limits>
#include<numeric>
#include<sys/stat.h>
#include<cmath>
//#include<thread>
//#include<chrono>


#define DEBUGVAR 1
// This assumes Event->vbm has already been initialized and configured

//General lightcurve generator that incorporates heirarchical orbital motion for multiple lenses and sources

// Keep this helper local to each generator translation unit.
// The call sites differ by API usage, but error flagging/logging behavior must stay identical.
static bool handle_vbm_api_error(const char* api_name, struct filekeywords* Paramfile, struct event* Event, ofstream& logfile_ptr, VBMicrolensing* vbm)
{
  if(!vbm->HasLastError())
    {
      return false;
    }

  const VBMicrolensing::LastError& err = vbm->GetLastError();
  Event->vbm_error_category = static_cast<int>(err.category);
  Event->vbm_error_source = err.where;
  Event->vbm_error_message = err.message;
  Event->lcerror = (err.category == VBMTimeoutError::TimeoutCategory::Unknown) ? LCGEN_VBM_ERR : LCGEN_TIMEOUT_ERR;
  Event->detected = 0;
  Event->deterror = 0;
  Event->outputthis = 0;

  if(Paramfile->verbosity >= 1)
    {
      cout << "VBM error in " << api_name << ": " << err.message << endl;
      cout << "Timeout category: " << VBMTimeoutError::CategoryName(err.category) << endl;
      if(!err.where.empty()) cout << "Timeout source: " << err.where << endl;
    }
  if(logfile_ptr.good())
    {
      logfile_ptr << "VBM error in " << api_name << ": " << err.message << endl;
      logfile_ptr << "Timeout category: " << VBMTimeoutError::CategoryName(err.category) << endl;
      if(!err.where.empty()) logfile_ptr << "Timeout source: " << err.where << endl;
      logfile_ptr.flush();
    }

  vbm->ClearLastError();
  return true;
}

//void sleep_thread(int n)
//{
//  this_thread::sleep_for(chrono::milliseconds(n));
//}

void lightcurveGenerator(struct filekeywords* Paramfile, struct event *Event, struct obsfilekeywords World[], struct slcat *Sources, struct slcat *Lenses, ofstream& logfile_ptr)
{

  Event->lcerror=0;
  Event->deterror=0;
  Event->vbm_error_category = static_cast<int>(VBMTimeoutError::TimeoutCategory::Unknown);
  Event->vbm_error_source.clear();
  Event->vbm_error_message.clear();
  Event->vbm->ClearLastError();

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
  Event->vbm->astrometry = (Paramfile->astrometry_on != 0);


  Event->astrox1_raw.resize(Event->nsrc);
  Event->astrox2_raw.resize(Event->nsrc);
  // astro*_raw in event frame with per epoch and per source values for debugging and testing. These are the raw outputs of the VBM, before any rotation or translation to event frame coordinates, so they are in the VBM frame and centered on the center of mass of the lens system. We store these for debugging and testing, and then we will apply the appropriate rotations and translations to get the final astrometric centroid in event frame coordinates, which will be stored in Event->astrox1/2. This way we can test the VBM astrometry logic path independently of the coordinate transformations to get to event frame coordinates, and we can also test the coordinate transformations independently by comparing the raw VBM output to the final event frame astrometry.
  Event->astrox_raw.resize(Event->nsrc);
  Event->astroy_raw.resize(Event->nsrc);
  for(int i=0; i<Event->nsrc; i++)
    {
      Event->astrox1_raw[i].assign(Event->nepochs,0.0);
      Event->astrox2_raw[i].assign(Event->nepochs,0.0);
      Event->astrox_raw[i].assign(Event->nepochs,0.0);
      Event->astroy_raw[i].assign(Event->nepochs,0.0);
    }
  Event->xc_srcs_only.assign(Event->nepochs,0.0);  // flux weighted addition of source-image centroids
  Event->yc_srcs_only.assign(Event->nepochs,0.0);
  Event->xc_src_lens.assign(Event->nepochs,0.0);  // flux weighted addition of source-image centroids and luminous lens centroids, for astrometry path only. This is the relevant blended centroid for astrometry, since ambient light does not contribute to astrometric blending by contract.
  Event->yc_src_lens.assign(Event->nepochs,0.0);

  Event->src_flux_total.assign(Event->nepochs,0.0); // total source flux (for calculating blended centroid), in units of the unmagnified source flux, for astrometry path only

  Event->moons = 0;
  Event->circumbinary=0;
  Event->distantbinary=0;
  Event->mixedbinary=0;
  
  //Conventions:
  //Track the apparent motion of the centers of mass of the source and lens, then compute offsets from them
  //for each component

  // Alpha is a scalar angle (not a vector) derived from the lens-source
  // relative proper-motion vector mu_rel = (mu_lambda, mu_beta), where
  // mu_rel = mu_lens - mu_source in ecliptic components.
  // Angle convention here:
  //   phi_pi = atan2(y, x) with y=mu_rel_lambda and x=mu_rel_beta
  // so phi_pi is measured from +mu_beta toward +mu_lambda (CCW in the
  // (mu_beta, mu_lambda) plane). Then:
  //   alpha = phi_pi - pi
  // alpha is used to orient the source trajectory via xs0/ys0.
  double phi_pi = atan2(Event->pllx[0].mulam_r, Event->pllx[0].mubet_r);
  double alpharad = phi_pi - pi;
  Event->alpha = 180.0/pi * alpharad;
    
  double antipode_ra = c.fold(Event->ra + PI,0,twoPi); //Used for computing orbits
  double antipode_dec = c.fold(-Event->dec,-PI,PI);
      
  vector<double> s_delta(3,0.0); //Offset of the chosen source from its center of mass at tref
  vector<double> l_delta(3,0.0); //Offset of the chosen lens from its center of mass at tref

  int sn = Event->source;
  int sc = -1;
  int ln = Event->lens;

  double rEsrc = Event->rE * Sources->data[sn][Sources->DIST]/Lenses->data[ln][Lenses->DIST];

  double cosa = cos(alpharad);
  double sina = sin(alpharad);
  
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
      double a1 = 1.0/(1.0+Event->scomp_q[0])*acomb;
      double a2 = Event->scomp_q[0]/(1.0+Event->scomp_q[0])*acomb;

      //orbitalElements(double a, double e, double I, double L, double w, double O, double dL_, double epoch_=J2000)
      s_elements[1][0] = orbitalElements(a1, Event->scomp_e[0], Event->scomp_I[0], Event->scomp_L0[0], Event->scomp_w[0], Event->scomp_O[0], Event->scomp_dL[0]);
      s_elements[0][0] = orbitalElements(-a2, Event->scomp_e[0], Event->scomp_I[0], Event->scomp_L0[0], Event->scomp_w[0], Event->scomp_O[0], Event->scomp_dL[0]);

      

      //Compute the origin shift relative to the center of mass of the lens - normalization by rEsrc will be done later
      vector<double> xp;      
      for(int j=0;j<int(s_elements[0].size());j++)
	{
	  s_elements[0][j].viewfrom(Event->tref+Paramfile->simulation_zerotime,antipode_ra,antipode_dec,&xp);
	  s_delta[0] += xp[0]; s_delta[1] += xp[1]; s_delta[2] += xp[2];
	}
      
    }
  //Hold the magnifications
  vector<double> mu(nsrc,0.0);
  vector<double> mupeak(nsrc,0.0);

  //Setup orbits for the lens(es)

  vector<vector<orbitalElements> > l_elements; //Orbital elements classes for the lenses
  int nlens = Event->nlens; //1 + Event->lcompanions.size() + Event->p_a.size();
  l_elements.clear();
  l_elements.resize(nlens);
  Event->VBM_function = "none";

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

	  // I don't know, Matt. Is it? -A
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
	      Event->p_x0.push_back(0.0);
	      Event->p_y0.push_back(0.0);
	      Event->p_z0.push_back(0.0);
	      Event->p_dsdt.push_back(0.0);
	      Event->p_dalphadt.push_back(0.0);
	      Event->p_dxdt.push_back(0.0);
	      Event->p_dydt.push_back(0.0);
	      Event->p_dzdt.push_back(0.0);
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

      int barycenters = 0;
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

  //Set up the origin for t0,u0. For all cases this will be relative to the first lens at time t_ref

  //Compute the origin shift relative to the center of mass of the lens - normalization by rEsrc will be done later

  const double toff = 1.0/(24.0*60.0); // one minute, in days
  
  if(nlens>=2)
    {
      vector<double> xp; //orbit contribution to vector
      vector<double> rp(3,0.0); //total position vector at time tref
      vector<double> rp1(3,0.0); //total position vector at time tref+delta
      vector<double> primary_rp(3,0.0); //primary lens position at tref, in AU
      vector<double> primary_rp1(3,0.0); //primary lens position at tref+delta, in AU
      for(int i=0;i<int(l_elements.size());i++)  // loop { star reflex orbits, planet orbits, moon orbits }
	{
	  rp.assign(3,0.0);  // Codex is insisting this needs to be here, but I get it. -A
	  rp1.assign(3,0.0);
	  for(int j=0;j<int(l_elements[i].size());j++)  // loop over orbital contributions for this lens object; each contribution has { a, e, I, L0, w, O, dL }
	    {
	      l_elements[i][j].viewfrom(Event->tref+Paramfile->simulation_zerotime,antipode_ra,antipode_dec,&xp);
	      rp[0] += xp[0]; rp[1] += xp[1]; rp[2] += xp[2];
	      l_elements[i][j].viewfrom(Event->tref+Paramfile->simulation_zerotime+toff,antipode_ra,antipode_dec,&xp);
	      rp1[0] += xp[0]; rp1[1] += xp[1]; rp1[2] += xp[2];
	      
	    }

	  if(i==0)
	    {
	      primary_rp = rp;
	      primary_rp1 = rp1;
	      l_delta[0] = rp[0]/Event->rE;
	      l_delta[1] = rp[1]/Event->rE;
	      l_delta[2] = rp[2]/Event->rE;
	    }
	  
	  // converting from physical units (AU) to Einstein radius units, and applying the origin shift for the lens-system barycenter.
	  rp[0] = (rp[0]-primary_rp[0])/Event->rE;
	  rp[1] = (rp[1]-primary_rp[1])/Event->rE;
	  rp[2] = (rp[2]-primary_rp[2])/Event->rE;
	  rp1[0] = (rp1[0]-primary_rp1[0])/Event->rE;
	  rp1[1] = (rp1[1]-primary_rp1[1])/Event->rE;
	  rp1[2] = (rp1[2]-primary_rp1[2])/Event->rE;

	  if(i>0)
	    {
	      Event->p_s0[i-1] = qAdd(rp[0],rp[1]);  // is p meant to be in AU here?
	      Event->p_x0[i-1] = rp[0];
	      Event->p_y0[i-1] = rp[1];
	      Event->p_z0[i-1] = rp[2];
	      double s1 = qAdd(rp1[0],rp1[1]);
	      Event->p_dsdt[i-1] = (s1-Event->p_s0[i-1])/toff;
	      Event->p_dalphadt[i-1] = (atan2(rp1[1],rp1[0])-atan2(rp[1],rp[0]))/toff;
	      Event->p_dxdt[i-1] = (rp1[0]-rp[0])/toff;
	      Event->p_dydt[i-1] = (rp1[1]-rp[1])/toff;
	      Event->p_dzdt[i-1] = (rp1[2]-rp[2])/toff;
	      
	      //p_dalphadt, p_dxdt, p_dydt;
	    }
	  

	}
    }
  

  

  double time_elapsed=0;
  bool warned_single_lens_nonzero_origin = false;

  //Calculate the lightcurve
  double last_progress=0;

  vector<double> msource(nsrc);


  //Throw away some events which we know to not be realistic or to stretch the
  //capabilities of the MultiLens generator

  int bad_scenario=0;

  if(nlens>Paramfile->num_lens_max)
    {
      bad_scenario+=iPow(2,0);
      cout << "Too many lenses, nlens=" << nlens << endl;
    }
  if(Event->mixedbinary>0)
    {
      bad_scenario+=iPow(2,1);
      cout << "Mixed binary, will skip" << endl;
    }
  else if(Event->circumbinary>0)
    {
      double e = Event->p_e[nlens-2];
      double mu_hw = Event->p_q[nlens-2];
      if(mu_hw>1) mu_hw = 1.0/mu_hw;
      double acrit = 1.60 + 5.10*e - 2.22*sqr(e) + 4.12*mu_hw - 4.27*e*mu_hw
	- 5.09*sqr(mu_hw) + 4.61*sqr(e*mu_hw); //Holman & Wiegert (1999)
      //cout << "circum binary acrit = " << acrit << " " << Event->p_a[0] << " " << Event->p_a[nlens-2]*acrit << endl;
      cout << "circum binary e=" << e << " mu_hw=" << mu_hw << " acrit=" << acrit << " a0=" << Event->p_a[0] << " ab=" << Event->p_a[nlens-2] << " a[n-1]*acrit" << Event->p_a[nlens-2]*acrit << " nlens=" << nlens << " p_e.size=" << Event->p_e.size() << endl;
      if(Event->p_a[0]<Event->p_a[nlens-2]*acrit) bad_scenario+=iPow(2,2);
      else
	{
	  double period_ratio = Event->p_period[0]/Event->p_period[nlens-2];
	  double prround = round(period_ratio);
	  if(prround<=9 && (fmod(period_ratio,prround)<0.02 || fmod(period_ratio,prround)>0.98))
	    bad_scenario+=iPow(2,3);
	}
    }
  else if(Event->distantbinary>0)
    {
      double e = Event->p_e[nlens-2];
      double mu_hw = Event->p_q[nlens-2];
      if(mu_hw>1) mu_hw = 1.0/mu_hw;
      double acrit = 0.464 - 0.380*mu_hw - 0.631*e + 0.586*mu_hw*e
	+ 0.150*sqr(e) - 0.198*mu_hw*sqrt(e); //Holman & Wiegert (1999)
      cout << "distant binary e=" << e << " mu_hw=" << mu_hw << " acrit=" << acrit << " a0=" << Event->p_a[0] << " ab=" << Event->p_a[nlens-2] << " a[n-1]*acrit=" << Event->p_a[nlens-2]*acrit << " nlens=" << nlens << " p_e.size=" << Event->p_e.size() << endl;

      if(Event->p_a[0]>Event->p_a[nlens-2]*acrit) bad_scenario+=iPow(2,4);
    }
  
  if(Paramfile->verbosity>=1)
    {
      for(int i=0;i<nlens-1;i++)
	{
	  cout << "lens " << i;
	  cout << " a=" << Event->p_a[i] << " ";
	  cout << " e=" << Event->p_e[i] << " ";
	  cout << " i=" << Event->p_I[i] << " ";
	  cout << " L0=" << Event->p_L0[i] << " ";
	  cout << " w=" << Event->p_w[i] << " ";
	  cout << " O=" << Event->p_O[i] << " ";
	  cout << " dL=" << Event->p_dL[i] << " ";
	  cout << " period=" << Event->p_period[i] << " ";
	  cout << " mass=" << Event->p_mass[i] << " ";
	  cout << " q=" << Event->p_q[i] << " ";
	  cout << " s0=" << Event->p_s0[i] << " ";
	  cout << " x0=" << Event->p_x0[i] << " ";
	  cout << " y0=" << Event->p_y0[i] << " ";
	  cout << " z0=" << Event->p_z0[i] << " ";
	  cout << " dsdt=" << Event->p_dsdt[i] << " ";
	  cout << " dalphadt=" << Event->p_dalphadt[i] << " ";
	  cout << " dxdt=" << Event->p_dxdt[i] << " ";
	  cout << " dydt=" << Event->p_dydt[i] << " ";
	  cout << " dzdt=" << Event->p_dzdt[i] << " ";
	  cout << " orbtype=" << Event->p_orbtype[i] << " ";
	  cout << endl;
	}
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


  if(Paramfile->verbosity>=1) cout << "Starting lightcurve generation" << endl;
  
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
	      cout << "Timeout reached, time_elapsed=" << time_elapsed << " timeout=" << timeout << endl;
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

      //Compute the source reference point in event-frame coordinates using
      //the requested alpha convention.
      double xs0 = uu * sina + tt * cosa;
      double ys0 = -uu * cosa + tt * sina;

      //Shift from CoM origin to lens1@tref origin
      xs0 += l_delta[0]; //already in units of Einstein ring
      ys0 += l_delta[1];

      vector<double> xs(nsrc,0.0); //source position in the plane of the sky, ecliptic sky coordinates in AU
      vector<double> ys(nsrc,0.0); // AU is a strange unit. Aren't these angles on the sky? -A
      vector<double> ds(nsrc,0.0); //source distance

      if(nsrc>1)
	{
	  for(int i=0;i<nsrc;i++)  // loop through source companions
	    {
	      vector<double> xp;
	      for(int j=0;j<int(s_elements[i].size());j++)
		{
		  s_elements[i][j].viewfrom(Event->jdtimes[obsidx][shiftedidx],antipode_ra,antipode_dec,&xp);
		  xs[i] += xp[0]; ys[i] += xp[1]; ds[i] += xp[2];
		}
	      xs[i] -= s_delta[0]; ys[i] -= s_delta[1]; ds[i] -= s_delta[2];
	      xs[i] /= rEsrc; ys[i] /= rEsrc; ds[i] /= rEsrc;  // now they are angles? -A MP: unitless relative to the Einstein radius
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
		  if(Paramfile->verbosity>=3 || idx==0) l_elements[i][j].print_elements();
		  l_elements[i][j].viewfrom(Event->jdtimes[obsidx][shiftedidx],antipode_ra,antipode_dec,&xp);
		  if(Paramfile->verbosity>=3 || idx==0) cout << setprecision(16) << i << " " << j << " " << Event->jdtimes[obsidx][shiftedidx] << " " << xp[0] << " " << xp[1] << " " << xp[2] << endl;
		  xl[i] += xp[0]; yl[i] += xp[1]; dl[i] += xp[2];
		}
	      //This shift of origin away from the center of mass has been moved to the source position
	      //xl[i] -= l_delta[0]; yl[i] -= l_delta[1]; dl[i] -= l_delta[2];


	      xl[i] /= Event->rE; yl[i] /= Event->rE; dl[i] /= Event->rE;

	      //rotations needed here?
	      lens_parameters[3*i+0] = xl[i];
	      lens_parameters[3*i+1] = yl[i];
	      Event->xlens[i][idx] = xl[i];  // is this lens parallax? -A MP: No, this is in the center of mass frame
	      Event->ylens[i][idx] = yl[i];
	    }

	  Event->vbm->SetLensGeometry(nlens,lens_parameters);
	}
      else
	{
	  xl[0] = 0.0; yl[0] = 0.0;
	  Event->xlens[0][idx] = 0.0; Event->ylens[0][idx] = 0.0; // how come this doesn't need parallax? -A MP: again, center of mass frame
	}

      

      //Finally ready to compute magnifications
      vector<double> astro_x(nsrc,0.0);
      vector<double> astro_y(nsrc,0.0);

      if(nlens==1)
	{
	  if(!warned_single_lens_nonzero_origin)
	    {
	      // whether we have COM or lens 1 at the origin of the event-frame, they would
	      // both mean the same thing for the single-lens case. 
	      // VBM's astrometry is scalar for the single lens case, so it assumes the lens 
	      // is at the origin of the event frame. 
	      const double lens_origin_tol = 1e-12;
	      if(fabs(xl[0]) > lens_origin_tol || fabs(yl[0]) > lens_origin_tol)
		{
		  cerr << "WARNING: single-lens ESPL assumes lens at origin (xl=yl=0), but got "
		       << "xl=" << xl[0] << ", yl=" << yl[0]
		       << " for event " << Event->id
		       << " at epoch " << Event->epoch[idx] << endl;
		  warned_single_lens_nonzero_origin = true;
		}
	    }
	  for(int is=0;is<nsrc;is++)
	    {
	      if(is==0) rho = Event->rs;	
	      else rho = Event->scomp_rs[is-1];
	      u = qAdd(xs[is],ys[is]);  //magnitude of the relative source-lens position vector (per source per epoch)
	      bool used_vbm_astrometry = false;
	      if(Paramfile->skip_magnification==0)
		{
		  mu[is] = Event->vbm->ESPLMag2(u, rho);  //calculate the single-lens magnification (per source per epoch)
		  Event->VBM_function = "ESPLMag2";  
		  if(handle_vbm_api_error("ESPLMag2", Paramfile, Event, logfile_ptr, Event->vbm)) return;
		  used_vbm_astrometry = true;
		}
	      else mu[is]=1.0;
		  
	      // ESPLMag2 provides a radial centroid for single-lens geometry.
	      // Project that radial value onto the source direction in event-frame.
	      if(Paramfile->astrometry_on && used_vbm_astrometry && u>1e-12)
		{
		  astro_x[is] = Event->vbm->astrox1 * xs[is]/u;  // Project the radial centroid onto the source direction in event-frame.
		  astro_y[is] = Event->vbm->astrox1 * ys[is]/u;  // the flux weighted source-image centroid is in the direction of the 
		  // source from the lens, so we can use the source coordinates to get the direction of the centroid shift and apply it 
		  // to the radial value of the centroid shift to get the astrometric centroid in ecliptic coordinates.
		}
	      else
		{
		  astro_x[is] = xs[is];  // if the lensing isn't significant, the "image" centroid is just the source position
		  astro_y[is] = ys[is];  //pr source per epoch, where x and y are in ecliptic coordinates with the lens at the origin
		  // this isn't yet stored in the event, but presumably it gets stored after we blend with any other sources.
		  // check this!!
		  // if we move on later with astrox1_raw and astrox2_raw, we will no longer be in ecliptic
		}
	      if(Paramfile->astrometry_on && used_vbm_astrometry)
		{
		  Event->astrox1_raw[is][idx] = Event->vbm->astrox1;  // x-axis in the ESPL-VBM frame is the lens-source axis, so this is
		  // the raw centroid shift along that axis, which is the only component for a single lens. We can store it incase we
		  // need to debug later.
		  Event->astrox2_raw[is][idx] = 0.0; // ESPL does not have an orthogonal component to the centroid shift, so this is just
		  // set to zero.
		}
	      Event->astrox_raw[is][idx] = astro_x[is];
	      Event->astroy_raw[is][idx] = astro_y[is];
	    } 
		// I feel like we are missing the blending, but lets let her cook.
	}
      else if(nlens==2)  //binary lens
	{
	  double s = qAdd(xl[1]-xl[0],yl[1]-yl[0]);  // scalar angular separation of the two lenses, in units of the Einstein radius
	  double q = 0.0;
	  if (Event->lcompanions.size()>0)
	    {
	      q = Event->lcomp_q[0];  // mass ratio of the two stellar lenses
	      
	    }
	  else
	    {
	      q = Event->p_q[0];  // mass ratio of the two lenses
	    }
	  double rot = atan2(yl[1],xl[1]);  // angle to rotate coordinates into the VBM binary lens frame, which is defined such that 
	  // the two lenses lie on the x-axis. This rotation is needed because the VBM binary lens magnification functions assume the 
	  // binary axis is along the x-axis, which can be at any angle on the sky, but we have the lens positions in ecliptic coordinates 
	  // (xl, yl). 
	  // The rotation (rot) is in the direction to rotate the event-frame lens positions into the VBM frame, which is a counterclockwise 
	  // rotation by the angle of the binary axis. To rotate the source positions into the VBM frame, we need to rotate by -rot.
	  double cr = cos(-rot); double sr = sin(-rot);  // event -> VBM frame coefficients
	  double cr_inv = cos(rot); double sr_inv = sin(rot);   // VBM -> event frame coefficients

	  for(int is=0;is<nsrc;is++)  // loop over sources
	    {
	      if(is==0) rho = Event->rs;  // get primary-source angular radius in units of the Einstein radius from rs
	      else rho = Event->scomp_rs[is-1];  // or companion-source angular radius from scomp_rs (does not include the primary).
	      //rotate coordintates to binary axis

		  // Shift to COM
	      //Binary mag works from the center of mass, so translate source to CoM, then rotate
	      double xs_com_ecl = xs[is] + l_delta[0];  // is this to COM or from?
	      double ys_com_ecl = ys[is] + l_delta[1];  // VBM is COM centered, but what was the event frame centered on?
		  // it must have lens 1 for this to make sense.

	      // rotate from event -> VBM frame
	      double xsi = cr*xs_com_ecl - sr*ys_com_ecl;
	      double ysi = sr*xs_com_ecl + cr*ys_com_ecl;
	      bool used_vbm_astrometry = false;  // this gets set to true if we use the VBM astrometry logic path.
	      if(Paramfile->skip_magnification==0)  // if we aren't skipping the magnification calculation...
		{
		  mu[is] = Event->vbm->BinaryMag2(s,q,xsi, ysi, rho);  //calculate the per source per epoch magnification using VBM
		  Event->VBM_function = "BinaryMag2";  // store the function that was used, for debugging
		  if(handle_vbm_api_error("BinaryMag2", Paramfile, Event, logfile_ptr, Event->vbm)) return;
		  used_vbm_astrometry = true;  // mark this path as executed, for debugging
		}
	      else mu[is] = 1.0;  // if we are skipping the magnification calculation, set the magnification to 1, and we won't use the VBM astrometry
	      if(Paramfile->astrometry_on && used_vbm_astrometry)
		{
		  // BinaryMag2 centroid is in binary-axis coordinates; inverse-rotate
		  // back to the canonical ecliptic axes, preserving barycenter origin.
		  double cx_bin = Event->vbm->astrox1;
		  double cy_bin = Event->vbm->astrox2;
		  astro_x[is] = cr_inv*cx_bin - sr_inv*cy_bin;
		  astro_y[is] = sr_inv*cx_bin + cr_inv*cy_bin;
		  // do we need to shift back to lens 1 as the origin? MP: I don't think so, it is the center of mass that is an inertial frame
		}
	      else
		{
		  astro_x[is] = xs[is];
		  astro_y[is] = ys[is];
		}
	      if(Paramfile->astrometry_on && used_vbm_astrometry) // if we aren't calculating the magnification, we aren't 
		  // calculating the astrometric shift, so they stay zero (I think. Provided they are initialized as such).
		{
		  Event->astrox1_raw[is][idx] = Event->vbm->astrox1;
		  Event->astrox2_raw[is][idx] = Event->vbm->astrox2;
		}
	      Event->astrox_raw[is][idx] = astro_x[is];
	      Event->astroy_raw[is][idx] = astro_y[is];
	    }
	}
      else //We're using multi-body lensing
	{	 
	  for(int is=0;is<nsrc;is++)
	    {
	      if(is==0) rho = Event->rs;
	      else rho = Event->scomp_rs[is-1];
	      bool used_vbm_astrometry = false;
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
		  VBMlocal.SetTimeouts(Event->vbm->GetTimeouts());
		  VBMlocal.SetErrorPolicy(Event->vbm->GetErrorPolicy());
		  VBMlocal.astrometry = Event->vbm->astrometry;

		  double u_min=1e50;
		  int fallback_lens_idx = -1;
		  double fallback_lens_weight = 0.0;
		  for(int i=0;i<nlens;i++)
		    {
		      const double lens_mass_frac = lens_parameters[i*3+2];
		      if(lens_mass_frac <= 0.0) continue;
		      double u_lens = qAdd(xs[is]-lens_parameters[i*3+0],ys[is]-lens_parameters[i*3+1])/sqrt(lens_mass_frac);
		      if(u_lens<u_min)
			{
			  u_min=u_lens;
			  fallback_lens_idx = i;
			  fallback_lens_weight = lens_mass_frac;
			}
		    }

		  VBMicrolensing* astrometry_vbm = nullptr;
		  if(u_min<10 && msource[is]<40)
		    {
		      mu[is] = VBMlocal.MultiMag2(xs[is], ys[is], rho);
		      Event->VBM_function = "MultiMag2";
		      if(handle_vbm_api_error("MultiMag2", Paramfile, Event, logfile_ptr, &VBMlocal)) return;
		      used_vbm_astrometry = true;
		      astrometry_vbm = &VBMlocal;
		      logfile_ptr << mu[is] << " " << VBMlocal.therr << " " << VBMlocal.NPS << endl;
		    }
		  else
		    {
		      // If the source is far from all lenses, approximate the system as the single
		      // component lens chosen by the fallback metric (minimum source-lens separation
		      // in component Einstein-radius units; the most influential lenser). ESPLMag2 works in 
		      // that lens's natural Einstein-radius units, so convert both u and rho before calling VBM.
		      const double rho_lens = (fallback_lens_weight > 0.0 ? rho/sqrt(fallback_lens_weight) : rho);
		      mu[is] = Event->vbm->ESPLMag2(u_min, rho_lens);
		      Event->VBM_function = "ESPLMag2";
		      if(handle_vbm_api_error("ESPLMag2", Paramfile, Event, logfile_ptr, Event->vbm)) return;
		      used_vbm_astrometry = true;
		      astrometry_vbm = Event->vbm;
		      logfile_ptr << mu[is] << " " << (u_min>=10?"single":"null") << " " << (msource[is]>=40?"faint":"null") << endl;
		    }

		  if(Paramfile->astrometry_on && used_vbm_astrometry && astrometry_vbm)
		    {
		      if(Event->VBM_function == "ESPLMag2" && fallback_lens_idx >= 0 && fallback_lens_weight > 0.0)
			{
			  const double dx = xs[is] - lens_parameters[fallback_lens_idx*3+0];
			  const double dy = ys[is] - lens_parameters[fallback_lens_idx*3+1];
			  const double dist = qAdd(dx,dy);
			  if(dist > 1e-12)
			    {
			      // ESPLMag2 returns a scalar centroid shift along the source-lens axis in the
			      // component-lens Einstein units. Project it back onto the 2D event frame and
			      // rescale by sqrt(mass fraction) to recover event-thetaE units.
			      const double ast_shift_event = astrometry_vbm->astrox1 * sqrt(fallback_lens_weight);
			      astro_x[is] = lens_parameters[fallback_lens_idx*3+0] + ast_shift_event * dx/dist;
			      astro_y[is] = lens_parameters[fallback_lens_idx*3+1] + ast_shift_event * dy/dist;
			    }
			  else
			    {
			      astro_x[is] = xs[is];
			      astro_y[is] = ys[is];
			    }
			}
		      else
			{
			  astro_x[is] = astrometry_vbm->astrox1;
			  astro_y[is] = astrometry_vbm->astrox2;
			}
		    }
		  else
		    {
		      astro_x[is] = xs[is];
		      astro_y[is] = ys[is];
		    }
		  if(Paramfile->astrometry_on && used_vbm_astrometry && astrometry_vbm)
		    {
		      Event->astrox1_raw[is][idx] = astrometry_vbm->astrox1;
		      Event->astrox2_raw[is][idx] = astrometry_vbm->astrox2;
		    }
		}

	      else mu[is] = 1.0;
	      Event->astrox_raw[is][idx] = astro_x[is];
	      Event->astroy_raw[is][idx] = astro_y[is];
	    }
	}

      for(int is=0;is<nsrc;is++)
	{
	  Event->mu_src[is][idx] = mu[is];
	  if(mu[is]>mupeak[is])
	    {
	      mupeak[is] = mu[is];
	      Event->tpeak[is] = Event->epoch[idx];
	      Event->upeak[is] = 1.0/mupeak[is];
	    }
	}

      //Here Atrue is magnification, but later it gets converted into fractional flux
      Event->Atrue[idx] = mu[0];
      for(int is=1;is<nsrc;is++)
	{
	  Event->Atrue[idx] += Event->scomp_fsofs1[is-1][filt] * (mu[is]-1.0);
	}

      if(Paramfile->astrometry_on)
	{
	  double src_flux_tot = mu[0];  // just prinary source flux at this point
	  double cx_src = astro_x[0] * mu[0];  // weighted primary source positon
	  double cy_src = astro_y[0] * mu[0];
	  for(int is=1;is<nsrc;is++)
	    {
	      double src_flux = Event->scomp_fsofs1[is-1][filt] * mu[is];  //calculate epoch-wise flux of each companion source
	      src_flux_tot += src_flux;  // iteratively add the companion source fluxes
	      cx_src += astro_x[is] * src_flux;  // weighted centroid addition from each companion source
	      cy_src += astro_y[is] * src_flux;
	    }
	  if(src_flux_tot>0.0)
	    {
	      cx_src /= src_flux_tot;  // normalize the centroid by the total flux to get the flux-weighted centroid position
	      cy_src /= src_flux_tot;
	    }
	  else
	    {
	      cx_src = std::numeric_limits<double>::quiet_NaN();  // if total source flux is negative, something is wrong
	      cy_src = std::numeric_limits<double>::quiet_NaN();
	    }

	  Event->xc_srcs_only[idx] = cx_src;  // blended apparent source centroid
	  Event->yc_srcs_only[idx] = cy_src;
	  Event->src_flux_total[idx] = src_flux_tot;  // so I don't have to recalculate it in photometry.cpp
	}


      // Keep track of highest magnification
      if (Event->Atrue[idx] > Event->Amax)
	{
	  Event->Amax = Event->Atrue[idx];  // this was amp, but I think that was a bug -A
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
