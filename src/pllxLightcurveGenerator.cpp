#include "lightcurveGenerator.h"
#include "backupGenerator.h"
#include "astroFns.h"
#include "VBMicrolensingLibrary.h"
#include "coords.h"
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
    char str[100];
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
     
    //if the event is saturated in each band, no need to calculate the lightcurve
    if(Event->nepochs==0 || Event->allsat)
     {
      return;
     }
    vector<int> idxshift;
    int shiftedidx;
    for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
       idxshift.push_back(Event->nepochsvec[obsidx]);

	// Enable astrometry calculation in VBMicrolensing if requested
	Event->vbm->astrometry = (Paramfile->astrometry_on != 0);

	// Pre-compute rotation from lens axes (binary axis) to sky NE frame.
	// Prefer the parallax vector direction; otherwise use μ_rel direction (from Galactic to Equatorial via mulb2ad).
	double dPosAng = 0.0; // static lens here
	double piE_norm = hypot(Event->piEN, Event->piEE);
	bool parallax_on = (Paramfile->pllxMultiplyer && piE_norm > 1e-12);
	double PosAng = 0.0; // default neutral rotation if orientation cannot be inferred
	bool astrom_ok = false;
	if (parallax_on) {
		double phi_pi = atan2(Event->piEE, Event->piEN); // atan2(East, North)
		PosAng = phi_pi - alpha + dPosAng;
		astrom_ok = true;
	} else {
		coords c;
		double mua = 0.0, mud = 0.0; // Equatorial components (East, North)
		c.mulb2ad(Event->l, Event->b, Event->murel_l, Event->murel_b, &mua, &mud);
		double mu_norm = hypot(mua, mud);
		if (mu_norm > 1e-16) {
			double phi_mu = atan2(mua, mud); // atan2(East, North)
			PosAng = phi_mu - alpha + dPosAng;
			astrom_ok = true;
		} else {
			// Unknown sky orientation (no πE and no μ_rel direction); keep lens axes aligned to NE as a neutral fallback.
			if (Paramfile->verbosity >= 2) {
				cout << "[astrometry] Warning: πE off and μ_rel direction unavailable; using neutral rotation (PosAng=0)." << endl;
			}
			if (logfile_ptr.good()) {
				logfile_ptr << "[astrometry] Warning: πE off and μ_rel direction unavailable; using neutral rotation (PosAng=0)." << std::endl;
			}
			PosAng = 0.0;
			astrom_ok = false; // explicitly disable astrometric outputs
		}
	}
	double cosPos = cos(PosAng);
	double sinPos = sin(PosAng);

		//Calculate the lightcurve
    for (int idx = 0; idx < Event->nepochs; ++idx)
      {
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
            Event->xl1[idx] = VBM_origin;
            Event->yl1[idx] = 0.0;
            Event->xl2[idx] = VBM_origin + a;
            Event->yl2[idx] = 0.0;
	    Event->vbm->a1 = lim_gamma;
	    amp = Event->vbm->BinaryMag2(a, q, xsCoM, ysCenter, rs);

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

		// Get centroid for first source (already computed above in BinaryMag2 call)
		double cx1 = Event->vbm->astrox1;
		double cy1 = Event->vbm->astrox2;
		
		double amp2 = Event->vbm->BinaryMag2(a, q, xs2CoM, ys2Center, Event->scomp_rs[0]);
		Event->xs2[idx] = xs2CoM;
		Event->ys2[idx] = ys2Center;
		
		// Get centroid for second source
		double cx2 = Event->vbm->astrox1;
		double cy2 = Event->vbm->astrox2;
	      
		int filt = World[obsidx].filter;		
		Event->Atrue[idx] = amp + Event->scomp_fsofs1[0][filt] * (amp2-1);

		// Flux-weighted centroid for binary source
		// Total flux = amp (source 1) + fsofs1 * (amp2 - 1) (source 2 contribution)
		// Centroid = (amp * c1 + fsofs1 * amp2 * c2) / (amp + fsofs1 * (amp2 - 1))
		double flux1 = amp;
		double flux2 = Event->scomp_fsofs1[0][filt] * amp2;
		double total_flux = flux1 + Event->scomp_fsofs1[0][filt] * (amp2 - 1.0);
		
		// Flux-weighted centroid in lens frame (Einstein radii)
		double cx = (flux1 * cx1 + flux2 * cx2) / (flux1 + flux2);
		double cy = (flux1 * cy1 + flux2 * cy2) / (flux1 + flux2);
		
				// Save lens-frame centroid (Einstein radii) and also compute NE (mas)
						Event->xctrue[idx] = cx; // lens-frame x1
						Event->yctrue[idx] = cy; // lens-frame x2
						double cx_mas = cx * Event->thE;
						double cy_mas = cy * Event->thE;
				
						if (Paramfile->astrometry_on && astrom_ok) {
							double cN = cx_mas * cosPos + cy_mas * sinPos;   // North (mas)
							double cE = -cx_mas * sinPos + cy_mas * cosPos;  // East (mas)
							Event->cNtrue[idx] = cN;
							Event->cEtrue[idx] = cE;
							// initialize observed values here; noise added later in photometry
							Event->cNobs[idx] = cN;
							Event->cEobs[idx] = cE;
						} else {
							Event->cNtrue[idx] = 0.0;
							Event->cEtrue[idx] = 0.0;
							Event->cNobs[idx] = 0.0;
							Event->cEobs[idx] = 0.0;
						}

	      }
	    else
	      {
		Event->Atrue[idx] = amp;

		// Extract centroids from VBM (in Einstein radii in lens frame)
		// and transform to sky frame (mas)
		// VBM centroid is in the lens frame (x1, x2) where:
		//   x1 is along the binary axis
		//   x2 is perpendicular to the binary axis
		// We need to rotate to sky frame (North, East) using the trajectory angle alpha
		
		double cx = Event->vbm->astrox1;  // Centroid x in lens frame (Einstein radii)
		double cy = Event->vbm->astrox2;  // Centroid y in lens frame (Einstein radii)
		
				// Save lens-frame centroid (Einstein radii) and also compute NE (mas)
						Event->xctrue[idx] = cx; // lens-frame x1
						Event->yctrue[idx] = cy; // lens-frame x2
						double cx_mas = cx * Event->thE;
						double cy_mas = cy * Event->thE;
		
		// Rotate from lens frame to sky frame (North, East)
		// The lens frame x-axis is along the trajectory at angle alpha
		// Sky frame: cN = North offset (mas), cE = East offset (mas)
				// Proper rotation: [N, E] = R(PosAng) * [x1, x2]
				if (Paramfile->astrometry_on && astrom_ok) {
					double cN = cx_mas * cosPos + cy_mas * sinPos;   // North (mas)
					double cE = -cx_mas * sinPos + cy_mas * cosPos;  // East (mas)
					Event->cNtrue[idx] = cN;  // sky NE
					Event->cEtrue[idx] = cE;
					Event->cNobs[idx] = cN;   // initialize observed (noise added in photometry)
					Event->cEobs[idx] = cE;
				} else {
					Event->cNtrue[idx] = 0.0;
					Event->cEtrue[idx] = 0.0;
					Event->cNobs[idx] = 0.0;
					Event->cEobs[idx] = 0.0;
				}
	      }
	  }

        // Keep track of highest magnification
        if (amp > Event->Amax) {
            Event->Amax = amp;
            Event->peakpoint = idx;
        }
    }
  //if there has been an error - try the backup generator
  if(Event->lcerror)
    {
      Event->lcerror=0;
      backupGenerator(Paramfile, Event, World, Sources, Lenses, logfile_ptr);
    }
}
