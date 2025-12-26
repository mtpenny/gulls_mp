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
#define DEBUGVAR 1
// This assumes Event->vbm has already been initialized and configured

void lightcurveGenerator(struct filekeywords* Paramfile, struct event *Event, struct obsfilekeywords World[], struct slcat *Sources, struct slcat *Lenses, ofstream& logfile_ptr)
{
    char str[512];
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
    Event->vbm->astrometry = true; // ensure centroid information is populated for each call
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
    for (int idx = 0; idx < Event->nepochs; ++idx)  // loop over all epochs
    {
		if (enforce_timeout)
		{
			time_t now = time(NULL);
			if (difftime(now, starttime) > timeout)
			{
				timed_out = true;
				break;
			}
		} // if any epoch exceeds timeout, break loop

        obsidx = Event->obsidx[idx];  // which observatory is this epoch from
        shiftedidx=idx-idxshift[obsidx]; //epochs are stored sequentially for each observatory
		// idxshift is the starting index for each observatory's epochs
        double amp = 0.0;
        double combinedAstroX = 0.0;
        double combinedAstroY = 0.0;
        if (Paramfile->identicalSequence && obsidx > 0)
	  {
            // Lightcurve is identical from observatory to observatory
            amp = Event->Atrue[shiftedidx];  // use previously calculated magnification
			combinedAstroX = Event->xctrue[shiftedidx];  // true centroid x
			combinedAstroY = Event->yctrue[shiftedidx];  // true centroid y
			Event->musrc1[idx] = Event->musrc1[shiftedidx];
			Event->musrc2[idx] = Event->musrc2[shiftedidx];
			Event->Atrue[idx] = Event->Atrue[shiftedidx];
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
	  } else {

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
		Event->xl1[idx] = VBM_origin;
		Event->yl1[idx] = 0.0;
		Event->xl2[idx] = VBM_origin + a;
		Event->yl2[idx] = 0.0;
	    Event->vbm->a1 = lim_gamma;
	    amp = Event->vbm->BinaryMag2(a, q, xsCoM, ysCenter, rs);
	    double src1AstroX = Event->vbm->astrox1;
	    double src1AstroY = Event->vbm->astrox2;
	    combinedAstroX = src1AstroX;  // source 1 only for now
	    combinedAstroY = src1AstroY;

	    Event->vbm_rootaccuracy[idx] = Event->vbm->rootaccuracy;
	    Event->vbm_squarecheck[idx] = Event->vbm->squarecheck;
	    Event->vbm_therr[idx] = Event->vbm->therr;
	    
	    if(Paramfile->multiple_sources && Event->scompanions.size()>0)
	    {
			double xs2CoM, ys2Center;
			double x2off = Event->scomp_s[0] * cos(Event->scomp_phase[0]*TO_RAD);
			double y2off = Event->scomp_s[0] * sin(Event->scomp_phase[0]*TO_RAD) * cos(Event->scomp_inc[0]*TO_RAD);
			xs2CoM = xsCoM + x2off * cos(Event->scomp_alpha[0]*TO_RAD) - y2off * sin(Event->scomp_alpha[0]*TO_RAD);
			ys2Center = ysCenter + x2off * sin(Event->scomp_alpha[0]*TO_RAD) + y2off * cos(Event->scomp_alpha[0]*TO_RAD);

			double amp2 = Event->vbm->BinaryMag2(a, q, xs2CoM, ys2Center, Event->scomp_rs[0]);
			Event->xs2[idx] = xs2CoM;  // actually source 2 position
			Event->ys2[idx] = ys2Center; // actually source 2 position

			// blend centroid shifts from each source using their instantaneous fluxes
			double src2AstroX = Event->vbm->astrox1;
			double src2AstroY = Event->vbm->astrox2;
			double fluxRatio = 0.0;
			int filt = World[obsidx].filter; // which filter is being observed for this epoch
			if(Event->scomp_fsofs1.size()>0 && Event->scomp_fsofs1[0].size()>filt)  // sanity check/data validation
			{
			    fluxRatio = Event->scomp_fsofs1[0][filt]; // flux ratio of source 2/source 1 in this filter
				// scomp_fsofs1 => source companion flux over flux of source 1
			}
			double baseFlux1 = Event->fs[obsidx];
			double baseFlux2 = baseFlux1 * fluxRatio;
			double flux1 = baseFlux1 * amp;
			double flux2 = baseFlux2 * amp2;
			double totalFlux = flux1 + flux2;
			if(totalFlux > 0.0)
			{
			    combinedAstroX = (flux1 * src1AstroX + flux2 * src2AstroX) / totalFlux;
			    combinedAstroY = (flux1 * src1AstroY + flux2 * src2AstroY) / totalFlux;
			}

			// Store individual source magnifications
			Event->musrc1[idx] = amp;
			Event->musrc2[idx] = amp2;
			
			Event->Atrue[idx] = amp + Event->scomp_fsofs1[0][filt] * (amp2-1);

	    } else {
		
			Event->musrc1[idx] = amp;
			Event->musrc2[idx] = 0.0; // No second source
			Event->Atrue[idx] = amp;

		}
	}

	Event->xctrue[idx] = combinedAstroX; //source(s) only, blended centroid
	Event->yctrue[idx] = combinedAstroY;
	Event->xctrueerr[idx] = 0.0;
	Event->yctrueerr[idx] = 0.0;
	Event->xc[idx] = combinedAstroX; // alter at the photometry step
	Event->yc[idx] = combinedAstroY;
	Event->xcerr[idx] = 0.0;
	Event->ycerr[idx] = 0.0;

	// Keep track of highest magnification
	if (amp > Event->Amax) {
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
