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
  int sn = Event->source;
  int ln = Event->lens;
  int lc = -1;
  if(Event->lcompanions.size()>0) lc = Event->lcompanions[0];

  double baseline;
  double ampmag;
  int satflag;

  double nci, ncs, erri, errs;
  vector<double> phot;
  coords c;
  double lambda0 = 0.0;
  double beta0 = 0.0;
  c.ad2ecl(Event->ra, Event->dec, &lambda0, &beta0);
  const double cos_beta0 = cos(beta0);
  const double safe_cos_beta0 = (fabs(cos_beta0) > 1.0e-12 ? cos_beta0 : (cos_beta0 >= 0 ? 1.0e-12 : -1.0e-12));
  const double cos_dec0 = cos(Event->dec);
  const double safe_cos_dec0 = (fabs(cos_dec0) > 1.0e-12 ? cos_dec0 : (cos_dec0 >= 0 ? 1.0e-12 : -1.0e-12));
  const double mas_to_rad = TO_RAD/(3600.0*1000.0);
  const double mas_to_deg = 1.0/(3600.0*1000.0);
  const double ln256 = log(256.0);
  const double inv_sqrt_ln256 = 1.0 / sqrt(ln256);
  const double floor_mas = max(0.0, Paramfile->astrometry_error_floor_mas);
  double dRAc_from_eE = 0.0;
  double dDec_from_eE = 0.0;
  double dRAc_from_eN = 0.0;
  double dDec_from_eN = 0.0;
  c.muecl2ad(Event->ra, Event->dec, 1.0, 0.0, &dRAc_from_eE, &dDec_from_eE);
  c.muecl2ad(Event->ra, Event->dec, 0.0, 1.0, &dRAc_from_eN, &dDec_from_eN);

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
      int shiftedidx = idx-Event->nepochsvec[obsidx];

      filter = World[obsidx].filter;

      // Default astrometry state for per-epoch sky products.
      // Keep centroid inputs from omLightcurveGenerator.cpp intact.
      Event->xctrueerr[idx] = 0.0;
      Event->yctrueerr[idx] = 0.0;
      Event->xc[idx] = 0.0;
      Event->yc[idx] = 0.0;
      Event->xcerr[idx] = 0.0;
      Event->ycerr[idx] = 0.0;
      Event->lambda_noiseless_deg[idx] = lambda0 * TO_DEG;
      Event->beta_noiseless_deg[idx] = beta0 * TO_DEG;
      Event->ra_noiseless_deg[idx] = Event->ra * TO_DEG;
      Event->dec_noiseless_deg[idx] = Event->dec * TO_DEG;
      Event->ra_measured_deg[idx] = Event->ra * TO_DEG;
      Event->dec_measured_deg[idx] = Event->dec * TO_DEG;
      Event->sigma_ast_mas[idx] = 0.0;
      Event->ra_err_deg[idx] = 0.0;
      Event->dec_err_deg[idx] = 0.0;
      Event->ra_src_only_deg[idx] = Event->ra * TO_DEG;
      Event->dec_src_only_deg[idx] = Event->dec * TO_DEG;
      Event->ra_src_lens_deg[idx] = Event->ra * TO_DEG;
      Event->dec_src_lens_deg[idx] = Event->dec * TO_DEG;

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
		  const double thE_mas = Event->thE;
		  const double murel_ref = (Event->tE_r != 0.0 ? Event->thE/Event->tE_r*DAYINYR : 0.0);
		  const double dt_year = (shiftedidx >= 0 && shiftedidx < int(Event->pllx[obsidx].epochs.size()))
		    ? ((Event->pllx[obsidx].epochs[shiftedidx] - Event->pllx[obsidx].tref) / DAYINYR)
		    : ((Event->epoch[idx] - Event->tref) / DAYINYR);
		  const double pm_e_mas = murel_ref * Event->pllx[obsidx].mulam_r * dt_year;
		  const double pm_n_mas = murel_ref * Event->pllx[obsidx].mubet_r * dt_year;

		  double cx_src_thE = Event->xctrue[idx];
		  double cy_src_thE = Event->yctrue[idx];
		  if(Event->xc_src_only.size() > size_t(idx))
		    {
		      cx_src_thE = Event->xc_src_only[idx];
		      cy_src_thE = Event->yc_src_only[idx];
		    }

		  double fs_tot = Event->fs[obsidx];
		  if(Paramfile->multiple_sources && Event->scompanions.size()>0)
		    {
		      for(size_t cidx=0; cidx<Event->scompanions.size(); ++cidx)
			{
			  if(Event->scomp_fsofs1.size() > cidx && Event->scomp_fsofs1[cidx].size() > size_t(filter))
			    {
			      fs_tot += Event->fs[obsidx] * Event->scomp_fsofs1[cidx][filter];
			    }
			}
		    }

		  double fl1 = 0.0;
		  double fl2 = 0.0;
		  if(Paramfile->lenslight && ln >= 0 && ln < int(Lenses->mags.size()))
		    {
		      const double f_l1_over_s1 = pow(10.0, -0.4*(Lenses->mags[ln][filter] - Sources->mags[sn][filter]));
		      fl1 = Event->fs[obsidx] * f_l1_over_s1;
		      if(Paramfile->multiple_lenses && lc >= 0 && lc < int(Lenses->mags.size()))
			{
			  const double f_l2_over_s1 = pow(10.0, -0.4*(Lenses->mags[lc][filter] - Sources->mags[sn][filter]));
			  fl2 = Event->fs[obsidx] * f_l2_over_s1;
			}
		    }

		  double cx_blend_thE = cx_src_thE;
		  double cy_blend_thE = cy_src_thE;
		  double flux_sum = fs_tot;
		  if(fl1 > 0.0 && Event->xlens.size() > 0 && Event->ylens.size() > 0
		     && Event->xlens[0].size() > size_t(idx) && Event->ylens[0].size() > size_t(idx))
		    {
		      cx_blend_thE = (cx_blend_thE*flux_sum + Event->xlens[0][idx]*fl1)/(flux_sum + fl1);
		      cy_blend_thE = (cy_blend_thE*flux_sum + Event->ylens[0][idx]*fl1)/(flux_sum + fl1);
		      flux_sum += fl1;
		    }
		  if(fl2 > 0.0 && Event->xlens.size() > 1 && Event->ylens.size() > 1
		     && Event->xlens[1].size() > size_t(idx) && Event->ylens[1].size() > size_t(idx))
		    {
		      cx_blend_thE = (cx_blend_thE*flux_sum + Event->xlens[1][idx]*fl2)/(flux_sum + fl2);
		      cy_blend_thE = (cy_blend_thE*flux_sum + Event->ylens[1][idx]*fl2)/(flux_sum + fl2);
		      flux_sum += fl2;
		    }

		  Event->xc_src_only[idx] = cx_src_thE;
		  Event->yc_src_only[idx] = cy_src_thE;
		  Event->xc_src_lens[idx] = cx_blend_thE;
		  Event->yc_src_lens[idx] = cy_blend_thE;
		  Event->xctrue[idx] = cx_blend_thE;
		  Event->yctrue[idx] = cy_blend_thE;
		  Event->xctrueerr[idx] = 0.0;
		  Event->yctrueerr[idx] = 0.0;

		  const double e_src_only_mas = cx_src_thE * thE_mas + pm_e_mas;
		  const double n_src_only_mas = cy_src_thE * thE_mas + pm_n_mas;
		  const double e_src_lens_mas = cx_blend_thE * thE_mas + pm_e_mas;
		  const double n_src_lens_mas = cy_blend_thE * thE_mas + pm_n_mas;
		  const double lambda_src_only = lambda0 + (e_src_only_mas*mas_to_rad)/safe_cos_beta0;
		  const double beta_src_only = beta0 + n_src_only_mas*mas_to_rad;
		  const double lambda_src_lens = lambda0 + (e_src_lens_mas*mas_to_rad)/safe_cos_beta0;
		  const double beta_src_lens = beta0 + n_src_lens_mas*mas_to_rad;
		  double ra_src_only = Event->ra;
		  double dec_src_only = Event->dec;
		  double ra_src_lens = Event->ra;
		  double dec_src_lens = Event->dec;
		  c.ecl2ad(lambda_src_only, beta_src_only, &ra_src_only, &dec_src_only);
		  c.ecl2ad(lambda_src_lens, beta_src_lens, &ra_src_lens, &dec_src_lens);
		  Event->ra_src_only_deg[idx] = ra_src_only * TO_DEG;
		  Event->dec_src_only_deg[idx] = dec_src_only * TO_DEG;
		  Event->ra_src_lens_deg[idx] = ra_src_lens * TO_DEG;
		  Event->dec_src_lens_deg[idx] = dec_src_lens * TO_DEG;

		  const double e_noiseless_mas = cx_blend_thE * thE_mas + pm_e_mas;
		  const double n_noiseless_mas = cy_blend_thE * thE_mas + pm_n_mas;
		  Event->xc[idx] = e_noiseless_mas;
		  Event->yc[idx] = n_noiseless_mas;

		  const double lambda_noiseless = lambda0 + (e_noiseless_mas*mas_to_rad)/safe_cos_beta0;
		  const double beta_noiseless = beta0 + n_noiseless_mas*mas_to_rad;
		  Event->lambda_noiseless_deg[idx] = lambda_noiseless * TO_DEG;
		  Event->beta_noiseless_deg[idx] = beta_noiseless * TO_DEG;
		  double ra_noiseless = Event->ra;
		  double dec_noiseless = Event->dec;
		  c.ecl2ad(lambda_noiseless, beta_noiseless, &ra_noiseless, &dec_noiseless);
		  Event->ra_noiseless_deg[idx] = ra_noiseless * TO_DEG;
		  Event->dec_noiseless_deg[idx] = dec_noiseless * TO_DEG;

		  double sigma_phot = 0.0;
		  if(std::isfinite(Event->Aobs[idx]) && std::isfinite(Event->Aerr[idx]) && fabs(Event->Aobs[idx]) > 1.0e-12)
		    sigma_phot = fabs(Event->Aerr[idx]/Event->Aobs[idx]);
		  const double sigma_ast_psf_mas = sigma_phot * (World[obsidx].im.fwhm * 1000.0) * inv_sqrt_ln256;
		  const double sigma_ast_mas = sqrt(sigma_ast_psf_mas*sigma_ast_psf_mas + floor_mas*floor_mas);
		  Event->sigma_ast_mas[idx] = sigma_ast_mas;
		  Event->xcerr[idx] = sigma_ast_mas;
		  Event->ycerr[idx] = sigma_ast_mas;

		  const double noise_e_mas = sigma_ast_mas * gasdev(Paramfile->seed);
		  const double noise_n_mas = sigma_ast_mas * gasdev(Paramfile->seed);
		  const double lambda_measured = lambda0 + ((e_noiseless_mas + noise_e_mas)*mas_to_rad)/safe_cos_beta0;
		  const double beta_measured = beta0 + (n_noiseless_mas + noise_n_mas)*mas_to_rad;
		  double ra_measured = Event->ra;
		  double dec_measured = Event->dec;
		  c.ecl2ad(lambda_measured, beta_measured, &ra_measured, &dec_measured);
		  Event->ra_measured_deg[idx] = ra_measured * TO_DEG;
		  Event->dec_measured_deg[idx] = dec_measured * TO_DEG;

		  const double var_dRAc_mas2 = sigma_ast_mas*sigma_ast_mas*(dRAc_from_eE*dRAc_from_eE + dRAc_from_eN*dRAc_from_eN);
		  const double var_dDec_mas2 = sigma_ast_mas*sigma_ast_mas*(dDec_from_eE*dDec_from_eE + dDec_from_eN*dDec_from_eN);
		  Event->ra_err_deg[idx] = sqrt(var_dRAc_mas2) * mas_to_deg / safe_cos_dec0;
		  Event->dec_err_deg[idx] = sqrt(var_dDec_mas2) * mas_to_deg;
		}

	      //Test for saturation
	      Event->nosat[idx] = !satflag; //nosat is the oposite of satflag
	      if(Event->allsat && !satflag) Event->allsat = 0;
      if(Event->allsatobs[obsidx] && !satflag) Event->allsatobs[obsidx] = 0;
  
    }

}
