#include "lightcurveGenerator.h"
#include "backupGenerator.h"
#include "astroFns.h"
#include "photometryCodes.h"
#include "random.h"
#include<time.h>

#include<fstream>

#define DEBUGVAR 0

void photometry(struct filekeywords* Paramfile, struct event *Event, struct obsfilekeywords World[], struct slcat *Sources, struct slcat *Lenses, ofstream& logfile_ptr)
{
	(void)Sources;
	(void)Lenses;
	(void)logfile_ptr;
  int idx,obsidx;
  
  int filter;
  int sn = Event->source;

  double baseline;
  double ampmag;
  int satflag;
  double nci, ncs, erri, errs;
  vector<double> phot;

  //if the event is saturated in each band, no need to calculate the lightcurve
  if(Event->nepochs==0 || Event->allsat) 
    {
      return;
    }

  //baseline may be unsaturated, but all photometry may still be
  //will need to test for this
  Event->allsat=1;
  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++) 
    {
      Event->allsatobs[obsidx]=1;
    }

  //Perform the photometry
  for(idx=0;idx<Event->nepochs;idx++)
    {
      obsidx = Event->obsidx[idx];

      filter = World[obsidx].filter;

      if(World[obsidx].photcode==FASTAP)
	{
	  ampmag = Event->Atrue[idx];
	  World[obsidx].im.fast_photometry(ampmag, &nci, &ncs, &erri, &satflag);
	  errs=erri;

	  //store the results
	  baseline = Event->baselineFlux[obsidx] * Event->texp[idx] 
	    * Event->nstack[idx];
	  Event->Atrue[idx] = nci/baseline;
	  Event->Atrueerr[idx] = erri/baseline;
	  Event->Aobs[idx] = ncs/baseline;
	  Event->Aerr[idx] = errs/baseline;

	}
      else
	{
	  //add the background
	  World[obsidx].im.set_background(Event->backmag[idx]);
	  World[obsidx].im.addbg();
	  
	  //add the star
	  ampmag = Sources->mags[sn][filter]-2.5*log10(Event->Atrue[idx]);
	  World[obsidx].im.addstar(Event->xsub[obsidx], Event->ysub[obsidx], 
				   ampmag);

	  //do photometry
	  World[obsidx].im.wis_photometry(Event->xsub[obsidx], 
					  Event->ysub[obsidx],
					  Event->texp[idx], Event->nstack[idx],
					  &phot, &satflag);

	  //store the results
	  baseline = Event->baselineFlux[obsidx] * Event->texp[idx] 
	    * Event->nstack[idx];
	  
	  if(World[obsidx].photcode<2) //aperture photometry
	    {
	      Event->Atrue[idx] = phot[0]/baseline;
	      Event->Atrueerr[idx] = phot[1]/baseline;
	      Event->Aobs[idx] = phot[2]/baseline;
	      Event->Aerr[idx] = phot[3]/baseline;
	    }
	  else //weighted photometry
	    {
	      Event->Atrue[idx] = phot[4]/baseline;
	      Event->Atrueerr[idx] = phot[5]/baseline;
	      Event->Aobs[idx] = phot[6]/baseline;
	      Event->Aerr[idx] = phot[7]/baseline;
	    }

	  
					if(World[obsidx].photcode%2==0) //ideal photometry
	    {
	      Event->Aobs[idx] = Event->Atrue[idx];
	      Event->Aerr[idx] = Event->Atrueerr[idx];
	      Event->xc[idx] = Event->xctrue[idx];
	      Event->yc[idx] = Event->yctrue[idx];
			// For ideal mode, do not add noise to astrometry
			if (Paramfile->astrometry_on) {
				if (Event->cNtrue[idx] != 0.0 || Event->cEtrue[idx] != 0.0) {
					Event->cNobs[idx] = Event->cNtrue[idx];
					Event->cEobs[idx] = Event->cEtrue[idx];
					Event->cNobserr[idx] = Paramfile->astrometry_error_floor_mas;
					Event->cEobserr[idx] = Paramfile->astrometry_error_floor_mas;
				} else {
					Event->cNobs[idx] = 0.0;
					Event->cEobs[idx] = 0.0;
					Event->cNobserr[idx] = 0.0;
					Event->cEobserr[idx] = 0.0;
				}
			}
	    }
	  
	  //subtract the star
	  World[obsidx].im.substar(Event->xsub[obsidx], Event->ysub[obsidx], 
				   ampmag);
	  //subtract the background
	  World[obsidx].im.subbg();
	}

      // Astrometric errors and observed values (sky NE frame in mas)
      // This block runs after Aobs/Aerr are set for both photometry paths
      if (Paramfile->astrometry_on) {
        const double eps = 1e-12;
        double denom = (Event->Aerr[idx] > eps ? Event->Aerr[idx] : eps);
        double snr = fabs(Event->Aobs[idx]) / denom;
        // FWHM in arcsec from PSF; convert to mas
        double fwhm_mas = World[obsidx].im.fwhm * 1000.0;
        // Photon noise term per axis (mas)
        double sigma_photon = (snr > 0.0 ? fwhm_mas / snr : 0.0);
        // Total per-axis sigma
        double floor_mas = max(0.0, Paramfile->astrometry_error_floor_mas);
        double sigma_axis = sqrt(sigma_photon * sigma_photon + floor_mas * floor_mas);

        // Store errors and observed values only if we have a true centroid
        if (Event->cNtrue[idx] != 0.0 || Event->cEtrue[idx] != 0.0) {
          Event->cNobserr[idx] = sigma_axis;
          Event->cEobserr[idx] = sigma_axis;
          // add Gaussian noise
          Event->cNobs[idx] = Event->cNtrue[idx] + sigma_axis * gasdev(Paramfile->seed);
          Event->cEobs[idx] = Event->cEtrue[idx] + sigma_axis * gasdev(Paramfile->seed);
        } else {
          Event->cNobserr[idx] = 0.0;
          Event->cEobserr[idx] = 0.0;
          Event->cNobs[idx] = 0.0;
          Event->cEobs[idx] = 0.0;
        }
      } else {
        // Astrometry disabled
        Event->cNobserr[idx] = 0.0;
        Event->cEobserr[idx] = 0.0;
        Event->cNobs[idx] = 0.0;
        Event->cEobs[idx] = 0.0;
      }

      //Test for saturation
      Event->nosat[idx] = !satflag; //nosat is the oposite of satflag
      if(Event->allsat && !satflag) Event->allsat = 0;
      if(Event->allsatobs[obsidx] && !satflag) Event->allsatobs[obsidx] = 0;
  
    }

}

