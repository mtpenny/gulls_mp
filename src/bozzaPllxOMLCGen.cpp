#include "lightcurveGenerator.h"
#include "astroFns.h"
#include "VBMicrolensingLibrary.h"
#include "constdefs.h"
#include "columnCodes.h"
#include<time.h>
#include<vector>
#include<iomanip>
#include<fstream>

#define DEBUGVAR 0

//Helper functions for VBMicrolensing astrometry
void build_binary_astro_params(struct filekeywords* Paramfile, struct event *Event, 
                              struct slcat *Sources, double t0_abs, double* pr) {
  // Binary Astro parameters for BinaryAstroLightCurveOrbital (20 params)
  // Based on VBM guide and existing GULLS parameters
  
  double s = Event->params[SS];  // separation in Einstein radii
  double q = Event->params[QQ];  // mass ratio
  double u0 = Event->u0;
  double alpha = Event->alpha * TO_RAD;  // convert degrees to radians
  double rho = Event->rs;  // source size
  double tE = Event->tE_r;  // reference frame Einstein time
  
  // For Step 1, we'll use simple placeholders for astrometric parameters
  // These will be properly implemented in Steps 2-3
  double piN = 0.0;  // parallax North component (will implement conversion later)
  double piE = 0.0;  // parallax East component (will implement conversion later)
  double muS_Dec = 0.0;  // source proper motion Dec (will implement conversion later)
  double muS_RA = 0.0;   // source proper motion RA (will implement conversion later)
  double piS = 0.1;      // source parallax in mas (placeholder)
  double thetaE = Event->thE;  // Einstein angle in mas
  
  // Orbital parameters
  double period = Event->params[TT] * DAYINYR;  // period in days
  double a_orb = Event->params[AA];  // semimajor axis
  double inc = Event->params[INC];   // inclination
  double phase0 = Event->params[PHASE]; // initial phase
  double t_per = t0_abs;  // time of periastron (placeholder)
  double eccentricity = 0.0;  // assume circular orbit for now
  double omega = 0.0;  // argument of periastron
  double Omega = 0.0;  // longitude of ascending node
  
  // Fill parameter array for BinaryAstroLightCurveOrbital
  pr[0] = log(s);
  pr[1] = log(q);
  pr[2] = u0;
  pr[3] = alpha;
  pr[4] = log(rho);
  pr[5] = log(tE);
  pr[6] = t0_abs;
  pr[7] = log(period);
  pr[8] = log(a_orb);
  pr[9] = inc;
  pr[10] = phase0;
  pr[11] = t_per;
  pr[12] = eccentricity;
  pr[13] = omega;
  pr[14] = Omega;
  pr[15] = piN;
  pr[16] = piE;
  pr[17] = muS_Dec;
  pr[18] = muS_RA;
  pr[19] = piS;
  pr[20] = thetaE;
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
	  // Copy centroid values from the corresponding epoch
	  Event->centroid_x[idx] = Event->centroid_x[shiftedidx];
	  Event->centroid_y[idx] = Event->centroid_y[shiftedidx];
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
	  
	  // For Step 1: Compute astrometric centroids if astrometry is enabled
	  if(Event->vbm->astrometry) {
	    // Store the current astrometric centroids (VBM computes them during BinaryMag2)
	    // For Step 1, we'll use the source-plane coordinates as placeholders
	    // These are in Einstein units relative to the lens
	    Event->centroid_x[idx] = Event->vbm->astrox1 + xsrot;  // blend source and lens centroid
	    Event->centroid_y[idx] = Event->vbm->astrox2 + ysrot;  // placeholder blending
	  } else {
	    // If astrometry disabled, set to source position
	    Event->centroid_x[idx] = xsrot;
	    Event->centroid_y[idx] = ysrot;
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

