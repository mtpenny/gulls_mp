#include "lightcurveGenerator.h"
#include "astroFns.h"
#include "VBMicrolensingLibrary.h"
#include "constdefs.h"
#include "columnCodes.h"
#include <time.h>
#include <vector>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <cmath>
#include<fstream>

#define DEBUGVAR 0

//Helper functions for VBMicrolensing astrometry parameter conversion
void build_binary_astro_params(struct filekeywords* Paramfile, struct event *Event, 
                              struct slcat *Sources, double t0_abs, double* pr) {
  // Non-orbital Binary Astro parameters for BinaryAstroLightCurve (13 params)
  // pr[0..8]: standard binary parameters
  // pr[9..12]: astrometric parameters (muS_Dec, muS_RA, piS, thetaE)
  
  double s = Event->params[SS];  // separation in Einstein radii
  double q = Event->params[QQ];  // mass ratio
  double u0 = Event->u0;
  double alpha = Event->alpha * TO_RAD;  // convert degrees to radians
  double rho = Event->rs;  // source size
  double tE = Event->tE_r;  // reference frame Einstein time
  
  // Get source catalog data
  int sn = Event->source;
  
  // Convert parallax from ecliptic (piEN,piEE) to equatorial (piN,piE)
  // For Step 1, use simple approximation - will implement full conversion in Step 2
  double piN = Event->piEN;  // North component (simplified)
  double piE = Event->piEE;  // East component (simplified)
  
  // Convert source proper motion from Galactic (MUL/MUB) to equatorial (RA/Dec)
  // For Step 1, use simple approximation - will implement full conversion in Step 2
  double muS_Dec = Sources->data[sn][Sources->MUB];  // mas/yr (simplified)
  double muS_RA = Sources->data[sn][Sources->MUL];   // mas/yr (simplified)
  
  // Source parallax from distance
  double piS = 1000.0 / Sources->data[sn][Sources->DIST]; // mas from kpc
  
  // Einstein angle
  double thetaE = Event->thE;  // mas
  
  // Fill parameter array for BinaryAstroLightCurve (non-orbital)
  pr[0] = log(std::max(1e-12, s));
  pr[1] = log(std::max(1e-12, q));
  pr[2] = u0;
  pr[3] = alpha;
  pr[4] = log(std::max(1e-12, rho));
  pr[5] = log(std::max(1e-12, tE));
  pr[6] = t0_abs;
  pr[7] = piN;
  pr[8] = piE;
  pr[9] = muS_Dec;
  pr[10] = muS_RA;
  pr[11] = piS;
  pr[12] = thetaE;
  
  // Warn about simplified coordinate conversions in Step 1
  static bool warning_shown = false;
  if(!warning_shown && (piS > 0)) {
    std::cout << "Warning: Using simplified coordinate conversions for Step 1. "
              << "Full Galactic<->Equatorial conversion will be implemented in Step 2." << std::endl;
    warning_shown = true;
  }
}

//extern "C"
//{
//  void magfunc_(double *m1, double *a, double *xsCenter,  double *ysCenter, double *rs, double *Gamma, double *amp, double *eps, int *errflag);
//}

//#include "singleLens.h"

void lightcurveGenerator(struct filekeywords* Paramfile, struct event *Event, struct obsfilekeywords World[], struct slcat *Sources, struct slcat *Lenses, ofstream& logfile_ptr)
{
  char str[100];

  double xsCenter, ysCenter, rs, Gamma=0.4;
  double amp, eps=1.0e-3;
  double alpha, cosa, sina;
  
  vector<int> obsoffset(Paramfile->numobservatories,0);

  Event->Amax=-1;
  Event->umin=1e50;
  double lim_gamma=Event->gamma;
  double lcgen=Paramfile->LC_GEN;
  int idx,obsidx;

  int errflag;

  double tt, uu;

  double aorE = Event->params[AA]/Event->rE;
  double cosinc = cos(Event->params[INC]*TO_RAD);
  double phase0 = Event->params[PHASE]*TO_RAD;
  double period = Event->params[TT]*DAYINYR;
  double q=Event->params[QQ];
  double s;
  double a1 = q/(1+q)*aorE; //in rE
  double a2 = aorE - a1;
  double zerotime = Event->t0;
  if(Paramfile->parameterization==1) zerotime = Event->tcroin;

  double x20 = a2 * cos(phase0);
  double y20 = a2 * sin(phase0) * cosinc;
  double s0 = Event->params[SS];
  double xCoM = s0*q/(1+q); //Position of the center of mass in the frame of the primary lens at t0 or tcroin

  
  double phase;
  double x1, y1; //host position in inertial frame in rE
  double x2, y2; //planet position in inertial frame in rE
  double rot, cosrot, sinrot; //rotation angle to subtract to put source in rotating frame
  double xsin, ysin; //source position in the inertial frame
  double xsrot, ysrot; //source position in  the rotating frame

  //work out the event parameters in the fortran parametrization

  rs = Event->rs;	            /* source size */
  alpha = Event->alpha*TO_RAD;	    /* slope of the trajectory */
  Gamma = Event->gamma;	            /* limb-darkening profile */
  Event->vbm->a1 = lim_gamma;             /*  Linear limb-darkening coefficient.*/
 


  cosa = cos(alpha); sina = sin(alpha);
  /*xcom = a*(1-2*m1);*/
  //xcom = -m1*a;  /* Origin is primary lens*/
  //xcom = a*(1-2*m1); /* Origin is the Center of mass */

  Event->lcerror=0;
  errflag=0;

  //if the event is saturated in each band, no need to calculate the lightcurve
  if(Event->nepochs==0 || Event->allsat) 
    {
      return;
    }
  if(Paramfile->verbosity>=4)
    {
            Event->vbm_rootaccuracy.resize(Event->nepochs);
            Event->vbm_squarecheck.resize(Event->nepochs);
            Event->vbm_therr.resize(Event->nepochs);
    }  
    
  // Precompute astrometric centroids for each observatory using VBM Astro API
  vector<vector<double> > centroid_N_obs; // [obsidx][epoch] - North centroids in mas
  vector<vector<double> > centroid_E_obs; // [obsidx][epoch] - East centroids in mas
  
  if(Event->vbm->astrometry) {
    centroid_N_obs.resize(Paramfile->numobservatories);
    centroid_E_obs.resize(Paramfile->numobservatories);
    
    // Build binary astrometry parameters (non-orbital for Step 1)
    double pr[13]; // 13 parameters for BinaryAstroLightCurve
    double t0_abs = Paramfile->simulation_zerotime + Event->t0; // convert to absolute JD
    build_binary_astro_params(Paramfile, Event, Sources, t0_abs, pr);
    
    for(int obs=0;obs<Paramfile->numobservatories;obs++) {
      if(Event->jdtimes[obs].size() == 0) continue;
      
      // Use the per-observatory absolute JD times
      const vector<double>& times = Event->jdtimes[obs];
      int nobs = static_cast<int>(times.size());
      
      // Allocate result arrays (all required by VBM API)
      vector<double> mag(nobs), c1s(nobs), c2s(nobs), c1l(nobs), c2l(nobs), y1(nobs), y2(nobs);
      
      // Call VBM astrometry function with correct signature
      Event->vbm->BinaryAstroLightCurve(
        pr, const_cast<double*>(times.data()),
        mag.data(), c1s.data(), c2s.data(), c1l.data(), c2l.data(),
        y1.data(), y2.data(), nobs);
      
      // Store the sky centroids (already in mas)
      centroid_N_obs[obs] = c1s; // North component
      centroid_E_obs[obs] = c2s; // East component
    }
  }
  
  vector<int> idxshift;
  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    idxshift.push_back(Event->nepochsvec[obsidx]);

  int shiftedidx;

  //Calculate the lightcurve
  for(idx=0;idx<Event->nepochs;idx++)
    {
      obsidx = Event->obsidx[idx];
      shiftedidx = idx-idxshift[obsidx];

      if(Paramfile->identicalSequence && obsidx>0)
	{
	  //lightcurve is identical from observatory to observatory
	  amp = Event->Atrue[shiftedidx];
	  // Copy centroid values from the first observatory (assuming identical observing sequences)
	  if(Event->vbm->astrometry && centroid_N_obs.size() > 0 && 
	     shiftedidx < centroid_N_obs[0].size()) {
	    Event->centroid_N_mas[idx] = centroid_N_obs[0][shiftedidx];
	    Event->centroid_E_mas[idx] = centroid_E_obs[0][shiftedidx];
	  } else {
	    Event->centroid_N_mas[idx] = 0.0;
	    Event->centroid_E_mas[idx] = 0.0;
	  }
	}
      else
	{

	  /* Compute the magnification */

	  //parallax shifts in the fixed reference frame
	  tt = (Event->epoch[idx] - Event->t0) / Event->tE_r;
	  uu = Event->u0;

	  if(Paramfile->pllxMultiplyer)
	    {
	      //if(tt<3)
		//cerr << tt << " " << uu << " " <<Event->pllx[obsidx].tshift(idx-idxshift[obsidx]) << " "  << Event->pllx[obsidx].ushift(idx-idxshift[obsidx]) << endl;
	      //tt += Event->pllx[obsidx].tshift(idx-idxshift[obsidx]);
	      //uu += Event->pllx[obsidx].ushift(idx-idxshift[obsidx]);
	      tt += Event->pllx[obsidx].tshift[shiftedidx];
	      uu += Event->pllx[obsidx].ushift[shiftedidx];
	    }

	  Event->umin=min(Event->umin,qAdd(tt,uu));
	  

	  //orbital calculations
	  phase = phase0 + 2*pi*(Event->epoch[idx] - zerotime) / period;
	  x1 = a1 * cos(phase+pi);
	  y1 = a1 * sin(phase+pi) * cosinc;
	  x2 = a2 * cos(phase);
	  y2 = a2 * sin(phase) * cosinc;
	  Event->xl1[idx] = x1;	  Event->yl1[idx] = y1;
	  Event->xl2[idx] = x2;	  Event->yl2[idx] = y2;
	  s = qAdd(x2-x1,y2-y1);
	  rot = atan2(y2,x2) - atan2(y20,x20);
	  cosrot = cos(-rot); sinrot = sin(-rot);

	  //source position
	  xsin = tt*cosa - uu*sina - xCoM; //as viewed from earth in non-rotating frame
	  ysin = tt*sina + uu*cosa;
	  Event->xs[idx] = xsin; Event->ys[idx] = ysin;
	  xsrot = xsin*cosrot - ysin*sinrot; //in frame rotating with binary
	  ysrot = xsin*sinrot + ysin*cosrot;

	  //VBB uses CoM as the origin
	  //amp = VBBL.BinaryMagDark(s,q,xsrot,ysrot,rs,Gamma,eps);
          //amp = VBBL.BinaryMag2(s, q, xsrot,ysrot,rs);
	  if(Paramfile->verbosity>=3)
	    cout << setprecision(16) << s << " " << q << " " << xsrot
		 << " " << ysrot << " " << rs << setprecision(6)
		 << endl;
	  
	  amp = Event->vbm->BinaryMag2(s,q,xsrot,ysrot,rs);
	  
	  // Retrieve precomputed astrometric centroids
	  if(Event->vbm->astrometry && obsidx < centroid_N_obs.size()) {
	    // Find the index in the observatory-specific time array
	    int obs_epoch_idx = shiftedidx; // This should correspond to jdtimes[obsidx] index
	    if(obs_epoch_idx < centroid_N_obs[obsidx].size()) {
	      Event->centroid_N_mas[idx] = centroid_N_obs[obsidx][obs_epoch_idx];
	      Event->centroid_E_mas[idx] = centroid_E_obs[obsidx][obs_epoch_idx];
	    } else {
	      Event->centroid_N_mas[idx] = 0.0;
	      Event->centroid_E_mas[idx] = 0.0;
	    }
	  } else {
	    // Fallback: set to zero if astrometry not available
	    Event->centroid_N_mas[idx] = 0.0;
	    Event->centroid_E_mas[idx] = 0.0;
	  }
	  
          if(Paramfile->verbosity>=4)
            {
                 Event->vbm_rootaccuracy[idx] = Event->vbm->rootaccuracy;
                 Event->vbm_squarecheck[idx] = Event->vbm->squarecheck;
                 Event->vbm_therr[idx] = Event->vbm->therr;
             }
	  //xsCenter-=xcom; amp = pacAmp(qAdd(xsCenter,ysCenter));//for testing
	}

      Event->Atrue[idx] = amp;
 
      if( errflag != 0) 
	{
	  sprintf(str,"\nerror caught from magfunc_  errval:%d", 
		  Event->lcerror);
	  logfile_ptr << Event->lcerror << endl;
	  logfile_ptr << Event->u0 << " " << Event->tE_r << " " 
		      << Event->t0 << " " << Event->params[QQ] << " " 
		      << Event->params[SS] << " " << Event->rs << " " 
		      << xsrot << " " << ysrot << endl;
	  fmtline(str,WIDTH,"(lightcurveGenerator)");
	  errorHandler(errflag);
	  Event->lcerror=errflag;
	  //break;
	}

      //keep track of highest magnification
      if(amp>Event->Amax) 
	{
	  Event->Amax = amp;
	  Event->peakpoint = idx;
	}
  
    }

}
