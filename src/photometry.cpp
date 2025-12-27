#include "lightcurveGenerator.h"
#include "backupGenerator.h"
#include "astroFns.h"
#include "photometryCodes.h"
#include<time.h>

#include<fstream>

#define DEBUGVAR 0

void photometry(struct filekeywords* Paramfile, struct event *Event, struct obsfilekeywords World[], struct slcat *Sources, struct slcat *Lenses, ofstream& logfile_ptr)
{
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
		//icounts = number of ideal counts
  		//ncounts = poisson realized number of counts
  		//error = error on photometry
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
	    }
	  
	  //subtract the star
	  World[obsidx].im.substar(Event->xsub[obsidx], Event->ysub[obsidx], 
				   ampmag);
	  //subtract the background
	  World[obsidx].im.subbg();
	}
	  // astrometric error from Gould & Yee (2014)
	  //σast = σphot * FWHM / (ln 256)^(1/2) , where σphot is the fractional photometric precision
	  // Astrometric errors and observed values (sky xy frame in mas)
      // This block runs after Aobs/Aerr are set for both photometry paths
      if (Paramfile->astrometry_on) {
        const double eps = 1e-12;
        const double ln256 = log(256.0);
        const double inv_sqrt_ln256 = 1.0 / sqrt(ln256);

        double sigma_phot = 0.0;
        if (Event->Aerr[idx] > eps) {
          sigma_phot = Event->Aerr[idx]/Event->Aobs[idx];
        } else {
          sigma_phot = eps;
        }

        double fwhm_mas = World[obsidx].im.fwhm * 1000.0;
		double fwhm_er = fwhm_mas / Event->thE;  // in einsteins radii
        double sigma_astro = fwhm_er * sigma_phot * inv_sqrt_ln256;
        double floor_mas = max(0.0, Paramfile->astrometry_error_floor_mas);
		double floor_er = floor_mas / Event->thE; // in einsteins radii
		double sigmaAstro = sqrt(sigma_astro * sigma_astro + floor_er * floor_er);
		// blend the source centroid with lens and abient stars
		double fstot = 0.0;
		fstot += Event->fs[obsidx];

		if (Paramfile->multiple_sources && Event->scompanions.size()>0)
		{
			int sc = Event->scompanions[0];
			// flux ratio of source companion to source 1 in this filter
			double fluxRatio = 0.0;
			if(Event->scomp_fsofs1.size()>0 && Event->scomp_fsofs1[0].size()>filter)
			{
				fluxRatio = Event->scomp_fsofs1[0][filter];
			}
			fstot += Event->fs[obsidx] * fluxRatio;
		}

		// blend = baseline - sum(source_fluxes)
		double blend_flux;
		blend_flux = 1 - fstot;
		// blending using flux weighted centroids with the "lens" at (xl1,yl1) and the source(s) at (xctrue,yctrue)
		Event->xctrue[idx] = Event->xctrue[idx]*fstot + Event->xl1[idx]*blend_flux;
		Event->yctrue[idx] = Event->yctrue[idx]*fstot + Event->yl1[idx]*blend_flux;
		Event->xctrueerr[idx] = 0.0;
		Event->yctrueerr[idx] = 0.0;

		// add astrometric noise
		Event->xcerr[idx] = sigmaAstro;
		Event->ycerr[idx] = sigmaAstro;
		Event->xc[idx] = Event->xctrue[idx] + sigmaAstro * gasdev(Paramfile->seed);
		Event->yc[idx] = Event->yctrue[idx] + sigmaAstro * gasdev(Paramfile->seed);

		//Test for saturation
		Event->nosat[idx] = !satflag; //nosat is the opposite of satflag
		if(Event->allsat && !satflag) Event->allsat = 0;
		if(Event->allsatobs[obsidx] && !satflag) Event->allsatobs[obsidx] = 0;
	
    }

}
