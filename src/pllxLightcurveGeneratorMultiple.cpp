#include "lightcurveGenerator.h"
#include "backupGenerator.h"
#include "astroFns.h"
#include "VBMicrolensingLibrary.h"
#include "singleLens.h"
#include<time.h>
#include<vector>

#include<fstream>
#include <sstream>
#include<iomanip>
#include <sys/stat.h>
#define DEBUGVAR 0

extern "C"
{
  void magfunc_(double *m1, double *a, double *xsCenter,  double *ysCenter, double *rs, double *Gamma, double *amp, double *eps, int *errflag);
}

void lightcurveGenerator(struct filekeywords* Paramfile, struct event *Event, struct obsfilekeywords World[], struct slcat *Sources, struct slcat *Lenses, ofstream& logfile_ptr)
{
  char str[512];

  double m1, a;
  double q;
  double xsCoM,  xsCenter, ysCenter, rs, Gamma;
  double amp, eps=1.0e-3;
  double alpha, cosa, sina,VBM_origin,Mao_origin ;
  int useVBB=1;
  vector<int> obsoffset(Paramfile->numobservatories,0);
  Event->vbm->astrometry = true; // request centroid outputs from VBM
  Event->Amax=-1;
  Event->umin=1e50;
  double ampoldlc, ampvbm, dif_over_amp;
  int idx,obsidx;
  double lim_gamma=Paramfile->LD_GAMMA;
  int errflag;

  double tt, uu;
  double lcgen=Paramfile->LC_GEN;
  if(DEBUGVAR) cout << "LC_GEN: " << lcgen << endl;

  Event->vbm->SetMethod(VBMicrolensing::Method::Multipoly);

  m1 = 1.0/(1 + Event->params[QQ]); /* mass of the first lens m1+m2=1 */
  a = Event->params[SS];            /* separation */
  rs = Event->rs;                   /* source size */
  alpha = Event->alpha*TO_RAD;      /* slope of the trajectory */
  Gamma = Event->gamma;             /* limb-darkening profile */
  Event->vbm->a1 = lim_gamma;              /*  Linear limb-darkening coefficient.*/
  q = Event->params[QQ];
  cosa = cos(alpha); sina = sin(alpha);


  Mao_origin = -m1*a; //Translate from L1 origin to m1z2+m2z1=0 origin
  VBM_origin = (1-m1)*(-a); //Translate from L1 origin to center of mass origin

  fstream fstr;
  fstr.open("test.txt",fstream::out);
  
  
  int nn=4;
  
  double pr[] = {     //parameters
    0.0, 0.0, 1.0,    // First lens: x1_1, x1_2, m1
    1.0, -0.7, 1e-4,  // Second lens: x2_2, x2_2, m2
    2.0, 0.7, 1e-4,   // Third lens: x3_re, x3_im, m3
    0.6, -0.6, 1e-6   // Fourth lens: x4_re, x4_im, m4
  };


  if(Paramfile->verbosity>2) cout << "Set geometry" << endl;
  Event->vbm->SetLensGeometry(nn,pr);
  


  Event->lcerror=0;
  errflag=0;

  //if the event is saturated in each band, no need to calculate the lightcurve
  if(Event->nepochs==0 || Event->allsat)
    {
      return;
    }

  vector<int> idxshift;
  int shiftedidx;
  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    idxshift.push_back(Event->nepochsvec[obsidx]);

  //Calculate the lightcurve
  for(idx=0;idx<Event->nepochs;idx++)
    {
      obsidx = Event->obsidx[idx];
      shiftedidx=idx-idxshift[obsidx];
      double astroX = 0.0;
      double astroY = 0.0;

      if(Paramfile->identicalSequence && obsidx>0)

        {
          //lightcurve is identical from observatory to observatory
          amp = Event->Atrue[shiftedidx];
          astroX = Event->xctrue[shiftedidx];
          astroY = Event->yctrue[shiftedidx];
          Event->Atrue[idx] = Event->Atrue[shiftedidx];
          Event->musrc1[idx] = Event->musrc1[shiftedidx];
          Event->musrc2[idx] = Event->musrc2[shiftedidx];
          Event->vbm_rootaccuracy[idx] = Event->vbm_rootaccuracy[shiftedidx];
          Event->vbm_squarecheck[idx] = Event->vbm_squarecheck[shiftedidx];
          Event->vbm_therr[idx] = Event->vbm_therr[shiftedidx];
          Event->xs[idx] = Event->xs[shiftedidx];
          Event->ys[idx] = Event->ys[shiftedidx];
          Event->xs2[idx] = Event->xs2[shiftedidx];
          Event->ys2[idx] = Event->ys2[shiftedidx];
          Event->xl1[idx] = Event->xl1[shiftedidx];
          Event->yl1[idx] = Event->yl1[shiftedidx];
          Event->xl2[idx] = Event->xl2[shiftedidx];
          Event->yl2[idx] = Event->yl2[shiftedidx];
        }
      else
        {

          /* Compute the magnification */

          tt = (Event->epoch[idx] - Event->t0) / Event->tE_r;
          uu = Event->u0;

          if(Paramfile->pllxMultiplyer)
            {
              //tt += Event->pllx[obsidx].tshift(Event->jdepoch[idx]);
              //uu += Event->pllx[obsidx].ushift(Event->jdepoch[idx]);
              tt += Event->pllx[obsidx].tshift[shiftedidx];
              uu += Event->pllx[obsidx].ushift[shiftedidx];
            }


          Event->umin=min(Event->umin,qAdd(tt,uu));
          
          xsCoM = tt*cosa - uu*sina + VBM_origin; //coordinate shift to center of mass
          xsCenter = tt*cosa - uu*sina + Mao_origin;// coordinate shift to primary lens
          ysCenter = tt*sina + uu*cosa;
          if(Paramfile->verbosity>3) cout << hexfloat << Event->epoch[idx] << " " << Event->t0 << " " << Event->tE_r << " " << tt << " " << xsCoM << " " << ysCenter << " " << rs << endl;
          if(Paramfile->verbosity>3) fstr << hexfloat << Event->epoch[idx] << " " << Event->t0 << " " << Event->tE_r << " " << tt << " " << xsCoM << " " << ysCenter << " " << rs << endl;
          amp = Event->vbm->MultiMag2(xsCoM, ysCenter,rs);
          // In the VBM frame, astrox1 is the X-coordinate and astrox2 is the perpendicular (Y-like) coordinate.
          // They are intentionally mapped to astroX (X) and astroY (Y) in the sky/event frame.
          astroX = Event->vbm->astrox1;
          astroY = Event->vbm->astrox2;
          if(Paramfile->verbosity>3) cout << "Event->vbm->MultiMag2(xsCoM, ysCenter,rs); done" << endl;
          Event->vbm_rootaccuracy[idx] = Event->vbm->rootaccuracy;
          Event->vbm_squarecheck[idx] = Event->vbm->squarecheck;
          Event->vbm_therr[idx] = Event->vbm->therr;
  
          if(Paramfile->verbosity>3) cout << amp << endl;
          Event->musrc1[idx] = amp;
          Event->musrc2[idx] = 0.0;
        }

      // Still needs logic for multiple sources in the future

      Event->xctrue[idx] = astroX;
      Event->yctrue[idx] = astroY;
      Event->xctrueerr[idx] = 0.0;
      Event->yctrueerr[idx] = 0.0;
      Event->xc[idx] = astroX; // alter at the photometry step
      Event->yc[idx] = astroY;
      Event->xcerr[idx] = 0.0;
      Event->ycerr[idx] = 0.0;

      Event->Atrue[idx] = amp;
      if( errflag != 0) 
	{
	  snprintf(str, sizeof(str), "\nerror caught from magfunc_  errval:%d", 
		  Event->lcerror);
	  logfile_ptr << Event->lcerror << endl;
	  logfile_ptr << Event->u0 << " " << Event->tE_r << " " 
		      << Event->t0 << " " << Event->params[QQ] << " " 
		      << Event->params[SS] << " " << Event->rs << " " 
		      << xsCenter << " " << ysCenter << endl;
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

  fstr.close();

}
