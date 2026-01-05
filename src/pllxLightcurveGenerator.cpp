#include "lightcurveGenerator.h"
#include "backupGenerator.h"
#include "astroFns.h"
#include "VBMicrolensingLibrary.h"
#include "singleLens.h"
#include<time.h>
#include<vector>

#include<fstream>
#include<sstream>
#include<iomanip>
#include<sys/stat.h>
#define DEBUGVAR 1
// This assumes Event->vbm has already been initialized and configured

void lightcurveGenerator(struct filekeywords* Paramfile, struct event *Event, struct obsfilekeywords World[], struct slcat *Sources, struct slcat *Lenses, ofstream& logfile_ptr)
{
    double rs = Event->rs;
    double u0 = Event->u0;
    double alpha = Event->alpha * TO_RAD;
    double q = Event->params[QQ];
    double a = Event->params[SS];
    double m1 = 1.0 / (1.0 + q);
    double cosa = cos(alpha);
    double sina = sin(alpha);
    double VBM_origin = (1.0 - m1) * (-a);
    vector<int> obsoffset(Paramfile->numobservatories,0);
    Event->Amax=-1;
    Event->umin=1e50;
    Event->lcerror=0;
    int errflag=0;
    int obsidx;
    double lim_gamma=Event->gamma;
    if(Paramfile->verbosity>=3)
	    cout << "At lightcurveGenerator start, Tol=" 
		 << Event->vbm->Tol 
		 << ", RelTol=" 
		 << Event->vbm->RelTol 
		 << std::endl;


    if(Paramfile->verbosity>=3) cout << "About to resize xsrc" << endl;
    
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

    if(Paramfile->verbosity>=3) cout << "xsrc resized" << endl;

    Event->xlens.clear();
    Event->ylens.clear();
    Event->xlens.resize(Event->nlens);
    Event->ylens.resize(Event->nlens);
    for(int i=0; i<Event->nlens; i++)
      {
	Event->xlens[i].resize(Event->nepochs);
	Event->ylens[i].resize(Event->nepochs);
      }

    if(Paramfile->verbosity>=3) cout << "xlens resized" << endl;
     
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
    for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
       idxshift.push_back(Event->nepochsvec[obsidx]);

    //Calculate the lightcurve
    for (int idx = 0; idx < Event->nepochs; ++idx)
      {
	if (enforce_timeout)
	  {
	    time_t now = time(NULL);
	    if (difftime(now, starttime) > timeout)
	      {
		timed_out = true;
		cout << "Lightcurve generation timed out" << endl;
		break;
	      }
	  }
        obsidx = Event->obsidx[idx];
        shiftedidx=idx-idxshift[obsidx];
        double amp = 0.0;
        if (Paramfile->identicalSequence && obsidx > 0)
	  {
            // Lightcurve is identical from observatory to observatory
            amp = Event->Atrue[shiftedidx];
	  }
	else
	  {
	    double tt = (Event->epoch[idx] - Event->t0) / Event->tE_r;
            double uu = u0;

            if (Paramfile->pllxMultiplyer)
	      {
                tt += Event->pllx[obsidx].tshift[shiftedidx];
                uu += Event->pllx[obsidx].ushift[shiftedidx];
	      }
            Event->umin=min(Event->umin,qAdd(tt,uu));
            
	    double xsCoM = tt * cosa - uu * sina + VBM_origin;
            double ysCenter = tt * sina + uu * cosa;

            Event->xs[idx] = xsCoM;
            Event->ys[idx] = ysCenter;
	    Event->xsrc[0][idx] = xsCoM;
	    Event->ysrc[0][idx] = ysCenter;
            Event->xl1[idx] = Event->xlens[0][idx] = VBM_origin;
            Event->yl1[idx] = Event->ylens[0][idx] = 0.0;
            Event->xl2[idx] = Event->xlens[1][idx] = VBM_origin + a;
            Event->yl2[idx] = Event->ylens[1][idx] = 0.0;
	    Event->vbm->a1 = lim_gamma;
	    amp = Event->vbm->BinaryMag2(a, q, xsCoM, ysCenter, rs);
	    Event->mu_src[0][idx] = amp;

	    cout << 0 << " " << Event->epoch[idx] << " " << amp << endl;

	    Event->vbm_rootaccuracy[idx] = Event->vbm->rootaccuracy;
	    Event->vbm_squarecheck[idx] = Event->vbm->squarecheck;
	    Event->vbm_therr[idx] = Event->vbm->therr;
	    
	    if(Paramfile->multiple_sources && Event->scompanions.size()>0)
	      {
		double xs2CoM, ys2Center;
		cout << "scomp_s.size, scomp_phase.size" << " " << Event->scomp_s[0] << " " << Event->scomp_phase.size() << endl;
		double x2off = Event->scomp_s[0] * cos(Event->scomp_phase[0]*TO_RAD);
		cout << "x2off" << x2off << endl;
		double y2off = Event->scomp_s[0] * sin(Event->scomp_phase[0]*TO_RAD) * cos(Event->scomp_I[0]*TO_RAD);
		cout << "y2off" << y2off << endl;
		xs2CoM = xsCoM + x2off * cosa - y2off * sina;
		cout << "xs2CoM" << xs2CoM << endl;
		ys2Center = ysCenter + x2off * sina + y2off * cosa;
		cout << "ys2Center" << ys2Center << endl;
		double amp2 = Event->vbm->BinaryMag2(a, q, xs2CoM, ys2Center, Event->scomp_rs[0]);
		cout << 1 << " " << Event->epoch[idx] << " " << amp2 << endl;
		Event->xs2[idx] = xs2CoM;
		Event->ys2[idx] = ys2Center;
		Event->xsrc[1][idx] = xs2CoM;
		Event->ysrc[1][idx] = ys2Center;
		Event->mu_src[1][idx] = amp2;
	      
		int filt = World[obsidx].filter;		
		Event->Atrue[idx] = amp + Event->scomp_fsofs1[0][filt] * (amp2-1);

		//Put binary source astrometry here
		//Event->xctrue[idx] = ;
		//Event->yctrue[idx] = ;

	      }
	    else
	      {
		Event->Atrue[idx] = amp;


		//put single source astrometry here
		//Event->xctrue[idx] = ;
		//Event->yctrue[idx] = ;
	      }
	  }

        // Keep track of highest magnification
        if (amp > Event->Amax)
	  {
            Event->Amax = amp;
            Event->peakpoint = idx;
	  }
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
    //if there has been an error - try the backup generator
    if(Event->lcerror)
      {
	if(Event->lcerror == LCGEN_TIMEOUT_ERR)
	  {
	    if(Paramfile->verbosity >= 1)
	      {
		cout << "Skipping backup generator due to lightcurve timeout" << endl;
	      }
	    if(logfile_ptr.good())
	      {
		logfile_ptr << "Skipping backup generator due to lightcurve timeout" << std::endl;
	      }
	    return;
	  }
	Event->lcerror=0;
      backupGenerator(Paramfile, Event, World, Sources, Lenses, logfile_ptr);
    }
}
