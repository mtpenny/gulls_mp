#include "lightcurveGenerator.h"
#include "backupGenerator.h"
#include "astroFns.h"
#include "coords.h"
#include "photometryCodes.h"
#include<time.h>

#include<fstream>

#define DEBUGVAR 0

void photometry(struct filekeywords* Paramfile, struct event *Event, struct obsfilekeywords World[], struct slcat *Sources, struct slcat *Lenses, ofstream& logfile_ptr)
{
  int idx,obsidx;
  
  int filter;
  int sn = Event->source;  // catalog number of the source star
  int ln = Event->lens;  // catalog number of the lens star
  int lc = -1; 
  if(Event->lcompanions.size()>0) lc = Event->lcompanions[0];  // catalog number of the lens companion, if it exists

  double baseline;
  double ampmag;
  int satflag;

  double nci, ncs, erri, errs;
  vector<double> phot;
  coords c;
  const double thE_mas = Event->thE;  // Einstein angles in mas
  double lambda0 = 0.0;
  double beta0 = 0.0;
  c.ad2ecl(Event->ra, Event->dec, &lambda0, &beta0);  // convert the event's ra dec to ecliptic for the astrometry calculations
  const double cos_beta0 = cos(beta0);  // for ecliptic tangent plane to celestial sphere conversions (lambda, beta -> n, e)
  const double safe_cos_beta0 = (fabs(cos_beta0) > 1.0e-12 ? cos_beta0 : (cos_beta0 >= 0 ? 1.0e-12 : -1.0e-12));  // avoid division by zero
  const double cos_dec0 = cos(Event->dec); // for equatorial tangent plane to celestial sphere conversions (ra, dec -> N, E)
  const double safe_cos_dec0 = (fabs(cos_dec0) > 1.0e-12 ? cos_dec0 : (cos_dec0 >= 0 ? 1.0e-12 : -1.0e-12));
  const double mas_to_rad = TO_RAD/(3600.0*1000.0);
  const double mas_to_deg = 1.0/(3600.0*1000.0);
  const double ln256 = log(256.0);  // used in the Gould 2014 astrometric error prescription
  const double inv_sqrt_ln256 = 1.0 / sqrt(ln256);
  const double floor_mas = max(0.0, Paramfile->astrometry_error_floor_mas);  // an additional systematric error for the astrometry
  //double dRAc_from_eE = 0.0;
  //double dDec_from_eE = 0.0;
  //double dRAc_from_eN = 0.0;
  //double dDec_from_eN = 0.0;
  //c.muecl2ad(Event->ra, Event->dec, 1.0, 0.0, &dRAc_from_eE, &dDec_from_eE);
  //c.muecl2ad(Event->ra, Event->dec, 0.0, 1.0, &dRAc_from_eN, &dDec_from_eN);

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

  //Perform the photometry  (epoch loop)
  for(idx=0;idx<Event->nepochs;idx++)
    {
      obsidx = Event->obsidx[idx];
      int shiftedidx = idx-Event->nepochsvec[obsidx];

      filter = World[obsidx].filter;

      // Default astrometry state for per-epoch sky products
      Event->ra_noiseless_deg[idx] =0.0;
      Event->dec_noiseless_deg[idx] = 0.0;
      Event->ra_measured_deg[idx] = 0.0;
      Event->dec_measured_deg[idx] = 0.0;
      Event->sigma_ast_mas[idx] = 0.0;
      Event->ra_err_deg[idx] = 0.0;
      Event->dec_err_deg[idx] = 0.0;

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
		    }
	  
	  //subtract the star
	  World[obsidx].im.substar(Event->xsub[obsidx], Event->ysub[obsidx], 
				   ampmag);
	  //subtract the background
	  World[obsidx].im.subbg();
		}

	      if(Paramfile->astrometry_on)
		{
		  const double murel_ref = (Event->tE_r != 0.0 ? Event->thE/Event->tE_r*DAYINYR : 0.0);
		  const double dt_year = (shiftedidx >= 0 && shiftedidx < int(Event->pllx[obsidx].epochs.size())) // check that the shifted index is within the bounds of the pllx epochs vector for this observatory
		    ? ((Event->pllx[obsidx].epochs[shiftedidx] - Event->pllx[obsidx].tref) / DAYINYR) // if the shifted index is out of bounds, fall back to using the unshifted epoch for this index (which may also be out of bounds, but at least won't be negative)
		    : ((Event->epoch[idx] - Event->tref) / DAYINYR);  // convert the epoch to years relative to the pllx reference epoch (tref) for this observatory, which is used for the proper motion and parallax calculations. Use the shifted index to access the pllx epochs if it's within bounds, otherwise use the unshifted index.
		  // mulam_r and mubet_r are unit-vector components of the reference-frame proper motion
		  // in ecliptic coordinates, so multiply by the scalar magnitude murel_ref to get mas/yr.
		  const double pm_lam_mas = murel_ref * Event->pllx[obsidx].mulam_r * dt_year; // proper motion contribution to the ecliptic eastward centroid shift in mas
		  const double pm_beta_mas = murel_ref * Event->pllx[obsidx].mubet_r * dt_year;

		  const double cx_srcs_thE = Event->xc_srcs_only[idx];  // blended apparent source centroid (without lens light contribution)
		  const double cy_srcs_thE = Event->yc_srcs_only[idx];  // in ecliptic coordinates, in theta E units, from omLightcurveGenerator.cpp
		  const double fstotofs1 = Event->src_flux_total[idx];  // the is the per epoch sum of the magnified source fluxes, divided by fs1

		  // Event->scomp_fsofs1[is-1][filt]
		  // Event->baselineFlux[obsidx]

		  // Calculating the realtive lens fluxes
		  // -2.5log10(Ftot/Fs1) = mag1 - magtot
		  // -2.5log10(Fs1/1) = mag1 - m0
		  // => m0 = mag1 + 2.5log10(Fs1)
		  // where mag1 is the magnitude of the primary source, and m0 is the magnitude zero point for the event,
		  // in a source 1 flux normalized system.
		  const double m0 = Sources->mags[sn][filter] + 2.5*log10(1.0); // zp for a system where the baseline flux of source 1 is 1.0, in the current epochs band
		  // ml - m0 = -2.5log10(Fl/1) => Fl = 10^(-0.4*(ml-m0))
		  const double fl1ofs1 = pow(10.0, -0.4*(Lenses->mags[ln][filter]-m0));  // lens 1 flux relative to source 1 flux in the current band
		  double fl2ofs1 = 0.0;  // lens 2 flux relative to source 1 flux, if there's a luminous companion to the lens
		  //check if there is a lens companion
		  if (Paramfile->multiple_lenses && lc >= 0 && lc < int(Lenses->mags.size()))
		    {
		      fl2ofs1 = pow(10.0, -0.4*(Lenses->mags[lc][filter]-m0));  // lens 2 flux relative to source 1 flux, if there's a luminous companion to the lens
		    }
		  // there's never 3 luminous lenses in the current implementation

		  double cx_blend_thE = cx_srcs_thE; // start with the source centroid, and then add the lens light contribution if there is any luminous lens and if the lens is within the image (i.e. has defined xlens and ylens values) for this epoch
		  double cy_blend_thE = cy_srcs_thE;
		  double flux_sum = fstotofs1;  // start with the total source flux, and then add the lens flux if there is a luminous lens

          // flux weighted addition of the lens 1 centroid
		  if(fl1ofs1 > 0.0 && Event->xlens.size() > 0 && Event->ylens.size() > 0
		     && Event->xlens[0].size() > size_t(idx) && Event->ylens[0].size() > size_t(idx)) // check that the lens has defined positions for this epoch before trying to use them
		    {
		      cx_blend_thE = (cx_blend_thE*flux_sum + Event->xlens[0][idx]*fl1ofs1) / (flux_sum + fl1ofs1);
		      cy_blend_thE = (cy_blend_thE*flux_sum + Event->ylens[0][idx]*fl1ofs1) / (flux_sum + fl1ofs1);
		      flux_sum += fl1ofs1;
		    }

		  // flux weighted addition of the lens 2 centroid, if there's a luminous companion to the lens
		  if(fl2ofs1 > 0.0 && Event->xlens.size() > 1 && Event->ylens.size() > 1
		     && Event->xlens[1].size() > size_t(idx) && Event->ylens[1].size() > size_t(idx))  // check that the second lens has defined positions for this epoch before trying to use them
		    {
		      cx_blend_thE = (cx_blend_thE*flux_sum + Event->xlens[1][idx]*fl2ofs1) / (flux_sum + fl2ofs1);
		      cy_blend_thE = (cy_blend_thE*flux_sum + Event->ylens[1][idx]*fl2ofs1) / (flux_sum + fl2ofs1);
		      flux_sum += fl2ofs1;
		    }

		  Event->xc_src_lens[idx] = cx_blend_thE;  // save the blended centroid with lens contribution in theta E units, in the event structure
		  Event->yc_src_lens[idx] = cy_blend_thE;  // so that we can output it in the lightcurve file, and also use it for diagnostics


		  // Converting into the observable frame

		  // Convert to mas
		  double de_mas = cx_blend_thE * thE_mas;
		  double dn_mas = cy_blend_thE * thE_mas;

		  // move to absolute position in rad
		  double e0 = lambda0;
		  double n0 = beta0;

		  // collect shifts in mas
		  double de_lens_pllx = 0.0;  // TODO: fill from **existing parallax code**
		  double dn_lens_pllx = 0.0;  // mas
		  de_mas += pm_lam_mas + de_lens_pllx;  // adding lens motion shift lens parallax shift and relative centroid shift
		  dn_mas += pm_beta_mas + dn_lens_pllx;
		  double de_rad = de_mas * mas_to_rad;
		  double dn_rad = dn_mas * mas_to_rad;
		  double dlam_rad = de_rad / safe_cos_beta0;
		  double dbet_rad = dn_rad;

		  double lambda_noiseless = e0 + dlam_rad;  // absolute blended apparent source centroid in ecliptic longitude, in radians
		  double beta_noiseless = n0 + dbet_rad;  // absolute blended apparent

		  double ra_noiseless = 0.0;
		  double dec_noiseless = 0.0;
		  c.ecl2ad(lambda_noiseless, beta_noiseless, &ra_noiseless, &dec_noiseless);
		  Event->ra_noiseless_deg[idx] = ra_noiseless * TO_DEG;
		  Event->dec_noiseless_deg[idx] = dec_noiseless * TO_DEG;
		  Event->lambda_noiseless_deg[idx] = lambda_noiseless * TO_DEG;
		  Event->beta_noiseless_deg[idx] = beta_noiseless * TO_DEG;

		  double sigma_phot = 0.0;
		  if(std::isfinite(Event->Aobs[idx]) && std::isfinite(Event->Aerr[idx]) && fabs(Event->Aobs[idx]) > 1.0e-12)
		    sigma_phot = fabs(Event->Aerr[idx]/Event->Aobs[idx]);
		  const double sigma_ast_psf_mas = sigma_phot * (World[obsidx].im.fwhm * 1000.0) * inv_sqrt_ln256;
		  const double sigma_ast_mas = sqrt(sigma_ast_psf_mas*sigma_ast_psf_mas + floor_mas*floor_mas);
		  Event->sigma_ast_mas[idx] = sigma_ast_mas;

		  const double noise_ra_mas = sigma_ast_mas * gasdev(Paramfile->seed);
		  const double noise_dec_mas = sigma_ast_mas * gasdev(Paramfile->seed);
		  Event->ra_measured_deg[idx] = ra_noiseless * TO_DEG + noise_ra_mas * mas_to_deg / safe_cos_dec0;
		  Event->dec_measured_deg[idx] = dec_noiseless * TO_DEG + noise_dec_mas * mas_to_deg;
		  Event->ra_err_deg[idx] = sigma_ast_mas * mas_to_deg / safe_cos_dec0;
		  Event->dec_err_deg[idx] = sigma_ast_mas * mas_to_deg;
		}

	      //Test for saturation
	      Event->nosat[idx] = !satflag; //nosat is the oposite of satflag
	      if(Event->allsat && !satflag) Event->allsat = 0;
      if(Event->allsatobs[obsidx] && !satflag) Event->allsatobs[obsidx] = 0;
  
    }

}
