#include "outputLightcurve.h"
#include "zodiacalLight.h"
#include "astroFns.h"
#include "coords.h"
#include "constants.h"
#include "constdefs.h"
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>

#define DEBUGVAR 0

void outputLightcurve(struct event *Event, struct obsfilekeywords World[], struct filekeywords* Paramfile, struct slcat* Sources, struct slcat* Lenses)
{
  void muVisibility(double *mu, double rs, double z0, double ld1);
  string lcfname;
  //FILE *lcfile_ptr;
  ofstream lcfile;
  FILE *lcdatafile_ptr;
  int fileOpen=0;
  int i, obsidx;
  stringstream data;
  string tmp;
  string extension;
  string lcdatafname;
  double t;
  //double u0_ps,t0_ps,tE_ps,t2_ps,u_ps,psmag,psamp;
  //double u0_fs,t0_fs,tE_fs,t2_fs,u_fs,fsmag,fsamp;

  if(Paramfile->verbosity>=3)
    {
      cout << "Entering output lightcurve" << endl;
    }
  

  if(!Event->outputthis) return;
  else if((Paramfile->outputOnErr || Paramfile->outputOnDet 
	   || Paramfile->outputOnAll)==0) return;
  else
    {
      if((Event->lcerror+Event->deterror))
	{
	  if(!Paramfile->outputOnErr) return;
	}
      else if((Event->detected))
	{
	  if(!Paramfile->outputOnDet) return;
	}
      else if((!Event->detected))
	{
	  if(!Paramfile->outputOnAll) return;
	}

      if(Event->allsat || Event->nepochs==0) return;
    }

  if(Event->detected) extension = "det";
  else if(Event->lcerror||Event->deterror) extension = "err";
  else extension = "all";

  if(Paramfile->choosefield<0)
    {
      lcfname = Paramfile->outputdir + Paramfile->run_name + "_"
	+ to_string(Paramfile->instance) + "_" + to_string(Event->id) + "."
	+ extension + ".lc";
    }
  else
    {
      lcfname = Paramfile->outputdir + Paramfile->run_name + "_"
	+ to_string(Paramfile->instance) + "_" + to_string(Paramfile->choosefield) + "_"
	+ to_string(Event->id) + "." + extension + ".lc";
    }
  if(DEBUGVAR) cout << "lcname: " << lcfname << endl;

  if(Event->nepochs>0)
    {
      //lcfile_ptr = fopen(lcfname.c_str(),"w");
      lcfile.open(lcfname,fstream::out);
      fileOpen=1;
    }
  else return;

  if(Paramfile->verbosity>=4)
    {
      cout << "Writing extra lightcurve file" << endl;
      if(Paramfile->choosefield<0)
        {
	  lcdatafname = Paramfile->outputdir + Paramfile->run_name + "_"
	    + to_string(Paramfile->instance) + "_" + to_string(Event->id) + "."
	    + extension + ".lcdata";
        }
      else
        {
	  lcdatafname = Paramfile->outputdir + Paramfile->run_name + "_"
	    + to_string(Paramfile->instance) + "_" + to_string(Paramfile->choosefield) + "_" +
	    to_string(Event->id) + "." + extension + ".lcdata";
        }
      if(DEBUGVAR) cout << "lcdataname: " << lcdatafname << endl;
      if(Event->nepochs>0)
        {
          lcdatafile_ptr = fopen(lcdatafname.c_str(),"w");
          fileOpen=1;
	  fprintf(lcdatafile_ptr, "time Atrue rootaccuracy squarecheck therr \n");
	  //fprintf(lcdatafile_ptr, "time mag_old_lcgen mag_vbm dif_over_mag\n");
        }
      cout << "Header written" << endl;
     }


  //output header information

  //Is there a source companion
  int sc=-1;
  if(Event->scompanions.size()>0)
    {
      sc = Event->scompanions[0];
    }

  //Is there a lens companion
  int lc=-1;
  if(Event->lcompanions.size()>0)
    {
      lc = Event->lcompanions[0];
    }


  //blending
  //data.str(""); data << "#fs: ";
  lcfile << "#fs: ";
  for(int i=0;i<Paramfile->numobservatories;i++)
    {
      //data << Event->fs[i] << " ";
      lcfile << Event->fs[i] << " ";
    }
  lcfile << endl;
  //fprintf(lcfile_ptr,"%s\n",data.str().c_str());

  if(Paramfile->multiple_sources)
    {
      //data.str(""); data
      lcfile << "#fs2: ";
      for(int i=0;i<Paramfile->numobservatories;i++)
	{
	  if(sc>-1)
	    lcfile << Event->scomp_fsofs1[0][World[i].filter]*Event->fs[i] << " ";
	  else
	    lcfile << 0 << " ";
	}
      lcfile << endl;
      //fprintf(lcfile_ptr,"%s\n",data.str().c_str());
    }  

  //Source magnitudes
  //data.str(""); data
  lcfile << "#Sourcemag: ";
  int sn = Event->source;
  for(int i=0;i<Paramfile->Nfilters;i++)
    {
      lcfile << Sources->mags[sn][i] << " ";
    }
  lcfile << endl;
  //fprintf(lcfile_ptr,"%s\n",data.str().c_str());

  if(Paramfile->multiple_sources)
    {
      //data.str(""); data
      lcfile << "#Source2mag: ";
      for(int i=0;i<Paramfile->Nfilters;i++)
	{
	  if(sc>-1) lcfile << Sources->mags[sc][i] << " ";
	  else lcfile << 99 << " ";
	}
      lcfile << endl;
      //fprintf(lcfile_ptr,"%s\n",data.str().c_str());
    }


  //Source data
  //data.str(""); data
  lcfile << "#Sourcedata: ";
  lcfile << sn << " ";
  for(int i=0;i<Sources->data[sn].size();i++)
    {
      lcfile << Sources->data[sn][i] << " ";
    }
  lcfile << endl;
  //fprintf(lcfile_ptr,"%s\n",data.str().c_str());

  if(Paramfile->multiple_sources)
    {
      //data.str(""); data
      lcfile << "#Source2data: ";
      lcfile << sc << " ";
      for(int i=0;i<Sources->data[sn].size();i++)
	{
	  if(sc>-1) lcfile<< Sources->data[sc][i] << " ";
	  else lcfile << 1e-50 << " ";
	}
      lcfile << endl;
      //fprintf(lcfile_ptr,"%s\n",data.str().c_str());
    }
  

  //Observatory magnitudes
  //data.str(""); data
  lcfile << "#Obssrcmag: ";
  for(int i=0;i<Paramfile->numobservatories;i++)
    {
      lcfile << Sources->mags[sn][World[i].filter] << " ";
    }
  lcfile << endl;
  //fprintf(lcfile_ptr,"%s\n",data.str().c_str());

  if(Paramfile->multiple_sources)
    {
      //data.str(""); data
      lcfile << "#Obssrc2mag: ";
      for(int i=0;i<Paramfile->numobservatories;i++)
	{
	  if(sc>-1) lcfile << Sources->mags[sc][World[i].filter] << " ";
	  else lcfile << 99 << " ";
	}
      lcfile << endl;
      //fprintf(lcfile_ptr,"%s\n",data.str().c_str());
    }


  //Lens magnitudes
  //data.str(""); data
  lcfile << "#Lensmag: ";
  int ln = Event->lens;
  for(int i=0;i<Paramfile->Nfilters;i++)
    {
      lcfile << Lenses->mags[ln][i] << " ";
    }
  lcfile << endl;
  //fprintf(lcfile_ptr,"%s\n",data.str().c_str());

  if(Paramfile->multiple_lenses)
    {
      //data.str(""); data
      lcfile << "#Lens2mag: ";
      int ln = Event->lens;
      for(int i=0;i<Paramfile->Nfilters;i++)
	{
	  if(lc>-1) lcfile << Lenses->mags[lc][i] << " ";
	  else lcfile << 99 << " ";
	}
      lcfile << endl;
      //fprintf(lcfile_ptr,"%s\n",data.str().c_str());
    }

  
  //Lens data
  //data.str(""); data
  lcfile << "#Lensdata: ";
  lcfile << ln << " ";
  for(int i=0;i<Lenses->data[ln].size();i++)
    {
      lcfile << Lenses->data[ln][i] << " ";
    }
  lcfile << endl;
  //fprintf(lcfile_ptr,"%s\n",data.str().c_str());

  if(Paramfile->multiple_lenses)
    {
      //data.str(""); data
      lcfile << "#Lens2data: ";
      lcfile << lc << " ";
      for(int i=0;i<Lenses->data[ln].size();i++) //yes, this is meant to be [ln]
	{
	  if(lc>-1) lcfile << Lenses->data[lc][i] << " ";
	  else lcfile << 1e-50 << " ";
	}
      lcfile << endl;
      //fprintf(lcfile_ptr,"%s\n",data.str().c_str());
    }

  //Observatory magnitudes
  //data.str(""); data
  lcfile << "#Obslensmag: ";
  for(int i=0;i<Paramfile->numobservatories;i++)
    {
      lcfile << Lenses->mags[ln][World[i].filter] << " ";
    }
  lcfile << endl;
  //fprintf(lcfile_ptr,"%s\n",data.str().c_str());

  if(Paramfile->multiple_lenses)
    {
      //data.str(""); data
      lcfile << "#Obslensmag: ";
      for(int i=0;i<Paramfile->numobservatories;i++)
	{
	  if(lc>-1) lcfile << Lenses->mags[lc][World[i].filter] << " ";
	  else lcfile << 99 << " ";
	}
      lcfile << endl;
      //fprintf(lcfile_ptr,"%s\n",data.str().c_str());
    }
  
  //Planet data
  //data.str("");
  lcfile << "#Planet: ";
  //for(int i=0;i<NPLANETINPUT+NPLANETDERIV;i++)
  for(auto param : Event->params)
    {
      lcfile << param << " ";
    }
  lcfile << endl;
  //fprintf(lcfile_ptr,"%s\n",data.str().c_str());

  //Microlensing data
  //data.str(""); data
  lcfile << "#Event: ";
  lcfile << Event->u0 << " " << Event->alpha << " " << setprecision(12)
	 << Event->t0 << " " << Event->tcroin << " " << setprecision(6) << " "
	 << Event->ucroin << " " << Event->rcroin << " " << Event->tE_r << " " << Event->rs;
  //fprintf(lcfile_ptr,"%s\n",data.str().c_str());
  lcfile << endl;

  // Astrometry reference-frame definitions used for output columns.
  // Internal centroid arrays (Event->xc/yc and related diagnostics) are in an
  // event frame whose orientation is set by alpha. We rotate those coordinates
  // into observer-centric ecliptic EN before writing the public centroid columns.
  // RA/Dec columns are generated from ecliptic EN with the local ecliptic->ICRS
  // tangent-plane transform.
  const double ra_base_rad = Event->ra;
  const double dec_base_rad = Event->dec;
  const double ra_base_deg = ra_base_rad * r2d;
  const double dec_base_deg = dec_base_rad * r2d;
  const double alpha_rad = Event->alpha * TO_RAD;
  const double cos_dec_eq = cos(dec_base_rad);
  const double inv_cos_dec_eq = (fabs(cos_dec_eq) > 1.0e-12 ? 1.0 / cos_dec_eq : 0.0);
  const double mas_to_deg = 1.0 / (3600.0 * 1000.0);
  const double days_to_years = 1.0 / DAYINYR;

  coords c;
  double dRAc_from_eE = 0.0;
  double dDec_from_eE = 0.0;
  double dRAc_from_eN = 0.0;
  double dDec_from_eN = 0.0;
  c.muecl2ad(ra_base_rad, dec_base_rad, 1.0, 0.0, &dRAc_from_eE, &dDec_from_eE);
  c.muecl2ad(ra_base_rad, dec_base_rad, 0.0, 1.0, &dRAc_from_eN, &dDec_from_eN);

  // Per-observatory rotation from internal event-frame x/y to ecliptic EN.
  // We anchor the event-frame trajectory direction (set by alpha) to the
  // reference-frame lens-source proper-motion unit vector in ecliptic coords.
  vector<double> evt_to_ecl_cos(Paramfile->numobservatories, 1.0);
  vector<double> evt_to_ecl_sin(Paramfile->numobservatories, 0.0);
  for(int obsi = 0; obsi < Paramfile->numobservatories; ++obsi)
    {
      double mu_ref_e = 0.0;
      double mu_ref_n = 0.0;
      if(obsi < int(Event->pllx.size()))
	{
	  mu_ref_e = Event->pllx[obsi].mulam_r;
	  mu_ref_n = Event->pllx[obsi].mubet_r;
	}

      const double mu_ref_norm = hypot(mu_ref_e, mu_ref_n);
      if(mu_ref_norm > 1.0e-12)
	{
	  const double phi_ref = atan2(mu_ref_n, mu_ref_e);
	  const double gamma = phi_ref - alpha_rad;
	  evt_to_ecl_cos[obsi] = cos(gamma);
	  evt_to_ecl_sin[obsi] = sin(gamma);
	}
    }
  const double gamma0_deg = atan2(evt_to_ecl_sin[0], evt_to_ecl_cos[0]) * r2d;

  // Lens proper motion in RA/Dec (derived from l/b).
  const int lens_idx = Event->lens;
  double lens_mu_ra = NAN;
  double lens_mu_dec = NAN;
  bool lens_mu_ok = false;
  if(lens_idx >= 0 && lens_idx < int(Lenses->data.size()))
    {
      if(int(Lenses->data[lens_idx].size()) > Lenses->BB)
	{
	  const double lens_mul = Lenses->data[lens_idx][Lenses->MUL];
	  const double lens_mub = Lenses->data[lens_idx][Lenses->MUB];
	  const double lens_l_rad = Lenses->data[lens_idx][Lenses->LL] * d2r;
	  const double lens_b_rad = Lenses->data[lens_idx][Lenses->BB] * d2r;
	  if(std::isfinite(lens_mul) && std::isfinite(lens_mub) && std::isfinite(lens_l_rad) && std::isfinite(lens_b_rad))
	    {
	      c.mulb2ad(lens_l_rad, lens_b_rad, lens_mul, lens_mub, &lens_mu_ra, &lens_mu_dec);
	      lens_mu_ok = std::isfinite(lens_mu_ra) && std::isfinite(lens_mu_dec);
	    }
	}
    }

  lcfile << "#Astrometry_Frame: ";
  lcfile << setprecision(12)
	 << "RA_rad=" << ra_base_rad << " Dec_rad=" << dec_base_rad
	 << " RA_deg=" << ra_base_deg << " Dec_deg=" << dec_base_deg
	 << " thE_mas=" << Event->thE << " t0=" << Event->t0
	 << " xy=ecliptic_EN(observer-centric) origin=canonical_pointing@tref";
  lcfile << endl;
  lcfile << "#Astrometry_EventToEcl: "
	 << "E_ecl_mas=cos(gamma_obs)*x_evt_mas-sin(gamma_obs)*y_evt_mas "
	 << "N_ecl_mas=sin(gamma_obs)*x_evt_mas+cos(gamma_obs)*y_evt_mas "
	 << "gamma_obs_deg=atan2(murel_ref_beta,murel_ref_lambda)-alpha_deg "
	 << "gamma_obs0_deg=" << gamma0_deg;
  lcfile << endl;
  lcfile << "#Astrometry_Transform: "
	 << "dRAcosDec_mas=(" << dRAc_from_eE << ")*E_ecl_mas+(" << dRAc_from_eN << ")*N_ecl_mas "
	 << "dDec_mas=(" << dDec_from_eE << ")*E_ecl_mas+(" << dDec_from_eN << ")*N_ecl_mas"
	 << endl;
  lcfile << "#Astrometry_BAGLE: x_E_arcsec=dRAcosDec_mas/1000 y_N_arcsec=dDec_mas/1000 "
	 << "model_frame=lens_relative quantity=centroid_minus_lens "
	 << "blendless_columns=RA_centroid_src_only_deg,Dec_centroid_src_only_deg "
	 << "lens_columns=RA_lens_primary_deg,Dec_lens_primary_deg" << endl;

  //Observatory groups
  for(int obsgroup=0; obsgroup<int(Event->obsgroups.size()); obsgroup++)
    {
      //data.str(""); data
      lcfile << "#Obsgroup: ";
      lcfile << obsgroup << " " << Event->flatchi2[obsgroup] << " "
	   << Event->flag_needFS[obsgroup] << " "
	   << (Event->flag_needFS[obsgroup]?Event->FSPL[obsgroup].chisq:
	       Event->PSPL[obsgroup].chisq); 
      //Members
      for(int grpidx=0; grpidx<int(Event->obsgroups[obsgroup].size()); grpidx++)
	{
	  lcfile << Event->obsgroups[obsgroup][grpidx] << " "; 
	}
      lcfile << endl;
      //fprintf(lcfile_ptr,"%s\n",data.str().c_str());
    }


  int ndF = int(Event->dF.size()) / Event->nepochs;
  // ─────────────── header ───────────────
  
  //static const char* baseCols[] = {
  lcfile << 
    "Simulation_time" << " " << "measured_relative_flux" << " " <<
    "measured_relative_flux_error" << " " << "true_relative_flux" << " " <<
    "true_relative_flux_error" << " " << "observatory_code" << " " <<
    "saturation_flag" << " " << "best_single_lens_fit" << " " <<
    "x_centroid_mas" << " " << "x_centroid_error_mas" << " " <<
    "y_centroid_mas" << " " << "y_centroid_error_mas" << " " <<
    "true_x_centroid_mas" << " " << "true_x_centroid_error_mas" << " " <<
    "true_y_centroid_mas" << " " << "true_y_centroid_error_mas" << " " <<
    "RA_centroid_deg" << " " << "Dec_centroid_deg" << " " <<               // observed ICRS centroid
    "RA_centroid_true_deg" << " " << "Dec_centroid_true_deg" << " " <<     // true ICRS centroid
    "RA_centroid_src_only_deg" << " " << "Dec_centroid_src_only_deg" << " " << // blendless source-only ICRS centroid
    "RA_centroid_src_lens_deg" << " " << "Dec_centroid_src_lens_deg" << " " << // source+lens ICRS centroid
    "RA_lens_primary_deg" << " " << "Dec_lens_primary_deg" << " " <<       // primary-lens ICRS astrometry
    "lens_pm_parallax_dRAcosDec_mas" << " " << "lens_pm_parallax_dDec_mas" << " " << // lens PM+parallax offsets
    "RA_centroid_lpllx_deg" << " " << "Dec_centroid_lpllx_deg" << " " <<   // observed + lens-parallax term
    "RA_true_lpllx_deg" << " " << "Dec_true_lpllx_deg" << " " <<           // true + lens-parallax term
    "lens_dist_kpc" << " " <<                                              // lens distance
    "lens_parallax_x_mas" << " " << "lens_parallax_y_mas" << " " <<       // lens parallax (ecl E/N, mas)
    "parallax_shift_t" << " " << "parallax_shift_u" << " " <<    "BJD" << " " <<
    "parallax_shift_x" << " " << "parallax_shift_y" << " " <<    "parallax_shift_z" << " " <<
    "observer_x_ecl_AU" << " " << "observer_y_ecl_AU" << " " << "observer_z_ecl_AU" << " ";
  // Astrometry diagnostic columns (mas unless noted)
  // Raw VBM output (VBM's internal frame, units: theta_E for debugging)
  lcfile << "vbm_astrox1_raw_thE" << " " << "vbm_astrox2_raw_thE" << " ";
  // Centroid at each blending step (mas)
  lcfile << "centroid_src_x_mas" << " " << "centroid_src_y_mas" << " ";           // sources only
  lcfile << "centroid_src_lens_x_mas" << " " << "centroid_src_lens_y_mas" << " "; // + lenses
  lcfile << "centroid_final_x_mas" << " " << "centroid_final_y_mas" << " ";       // + ambient = true
  
  // Per-source and per-lens positions remain in the internal event frame
  // (units: theta_E). Public centroid columns above are rotated to ecliptic EN.
  for(int i=0;i<Event->nsrc;i++)
    {
      lcfile << "source" << i << "_x_thE" << " " << "source" << i << "_y_thE" << " " << "source" << i << "_mu" << " ";
    }
  // Per-lens positions (event frame, units: theta_E)
  for(int i=0;i<Event->nlens;i++)
    {
      lcfile << "lens" << i << "_x_thE" << " " << "lens" << i << "_y_thE" << " ";
    }

  if (ndF > 0)
    {
      string dF_nopllx = string("dF_t0 dF_tE dF_u0 dF_alpha dF_s dF_q dF_rs");
      string dF_pllx   = string(" dF_piEN dF_piEE");
      vector<string> parstrings;
      int nfixpar = 7 + (Paramfile->pllxMultiplyer ? 2 : 0);
      if (Paramfile->pllxMultiplyer)
	split(dF_nopllx + dF_pllx, parstrings);
      else
	split(dF_nopllx, parstrings);
	
      //stringstream ss;
      for(size_t obsgroup = 0; obsgroup < Event->obsgroups.size(); ++obsgroup)
	{
	  int nobs    = int(Event->obsgroups[obsgroup].size());
	  int nparams = nfixpar + 2*nobs;
	
	  for(int idx = 0; idx < nparams; ++idx)
	    {
	      //ss.str("");  ss.clear();
	      lcfile << "ObsGroup_" << obsgroup << "_";
		
	      if (idx < nfixpar)
		{
		  lcfile << parstrings[idx] << " ";
		}
	      else
		{
		  int grpidx = (idx - nfixpar) / 2;
		  int obsidx = Event->obsgroups[obsgroup][grpidx];
		  
		  if (((idx - nfixpar) % 2) == 0)
		    lcfile << "dF_Fbase" << obsidx << " ";
		  else
		    lcfile << "dF_fs"     << obsidx << " ";
		}  
     
	      //fprintf(lcfile_ptr, "%s ", ss.str().c_str());
	    }
	}
    }
  // finish the header line
  lcfile << endl;
  //fprintf(lcfile_ptr, "\n");

  if(Paramfile->verbosity>=3)
    {
      cerr << "Lightcurve header printed. Will print " << Event->nepochs << " lines of data." << endl;
    }
    
  //output the lightcurve
  int shiftedidx;

  //if(lcfile_ptr!=NULL && fileOpen==1)
  //  {
  for(i=0;i<Event->nepochs;i++)
    {
      t=Event->epoch[i];
      obsidx=Event->obsidx[i];
      shiftedidx = i-Event->nepochsvec[obsidx];
	    
	      double thE = Event->thE;  // Angular Einstein radius (mas) for unit conversion
	      lcfile << setprecision(16) << Event->epoch[i] << " " << Event->Aobs[i] << " " << Event->Aerr[i] << " " << flush;
	      lcfile << Event->Atrue[i] << " " << Event->Atrueerr[i] << " " << obsidx << " " << flush; 
	      lcfile << (Event->nosat[i]?0:1) << " " << Event->Afit[i] << " " << flush;

	      // Rotate internal event-frame centroids into observer-centric ecliptic EN.
	      const double evt2ecl_c = evt_to_ecl_cos[obsidx];
	      const double evt2ecl_s = evt_to_ecl_sin[obsidx];

	      const double xc_obs_evt_mas = Event->xc[i];
	      const double yc_obs_evt_mas = Event->yc[i];
	      const double xctrue_evt_mas = Event->xctrue[i] * thE;
	      const double yctrue_evt_mas = Event->yctrue[i] * thE;
	      const bool have_stage_centroids = (
		Event->xc_src_only.size() > (size_t)i &&
		Event->yc_src_only.size() > (size_t)i &&
		Event->xc_src_lens.size() > (size_t)i &&
		Event->yc_src_lens.size() > (size_t)i &&
		Event->xc_src_lens_amb.size() > (size_t)i &&
		Event->yc_src_lens_amb.size() > (size_t)i
	      );

	      double src_only_e = xctrue_evt_mas;
	      double src_only_n = yctrue_evt_mas;
	      double src_lens_e = xctrue_evt_mas;
	      double src_lens_n = yctrue_evt_mas;
	      double src_lens_amb_e = xctrue_evt_mas;
	      double src_lens_amb_n = yctrue_evt_mas;
	      if(have_stage_centroids)
		{
		  const double src_only_evt_x = Event->xc_src_only[i] * thE;
		  const double src_only_evt_y = Event->yc_src_only[i] * thE;
		  const double src_lens_evt_x = Event->xc_src_lens[i] * thE;
		  const double src_lens_evt_y = Event->yc_src_lens[i] * thE;
		  const double src_lens_amb_evt_x = Event->xc_src_lens_amb[i] * thE;
		  const double src_lens_amb_evt_y = Event->yc_src_lens_amb[i] * thE;

		  src_only_e = evt2ecl_c * src_only_evt_x - evt2ecl_s * src_only_evt_y;
		  src_only_n = evt2ecl_s * src_only_evt_x + evt2ecl_c * src_only_evt_y;
		  src_lens_e = evt2ecl_c * src_lens_evt_x - evt2ecl_s * src_lens_evt_y;
		  src_lens_n = evt2ecl_s * src_lens_evt_x + evt2ecl_c * src_lens_evt_y;
		  src_lens_amb_e = evt2ecl_c * src_lens_amb_evt_x - evt2ecl_s * src_lens_amb_evt_y;
		  src_lens_amb_n = evt2ecl_s * src_lens_amb_evt_x + evt2ecl_c * src_lens_amb_evt_y;
		}

	      const double xc_obs_mas = evt2ecl_c * xc_obs_evt_mas - evt2ecl_s * yc_obs_evt_mas;
	      const double yc_obs_mas = evt2ecl_s * xc_obs_evt_mas + evt2ecl_c * yc_obs_evt_mas;
	      const double xctrue_mas = evt2ecl_c * xctrue_evt_mas - evt2ecl_s * yctrue_evt_mas;
	      const double yctrue_mas = evt2ecl_s * xctrue_evt_mas + evt2ecl_c * yctrue_evt_mas;

	      const double xcerr_evt_mas = Event->xcerr[i];
	      const double ycerr_evt_mas = Event->ycerr[i];
	      const double xctrueerr_evt_mas = Event->xctrueerr[i];
	      const double yctrueerr_evt_mas = Event->yctrueerr[i];
	      const double xcerr_mas = sqrt(evt2ecl_c*evt2ecl_c*xcerr_evt_mas*xcerr_evt_mas
					    + evt2ecl_s*evt2ecl_s*ycerr_evt_mas*ycerr_evt_mas);
	      const double ycerr_mas = sqrt(evt2ecl_s*evt2ecl_s*xcerr_evt_mas*xcerr_evt_mas
					    + evt2ecl_c*evt2ecl_c*ycerr_evt_mas*ycerr_evt_mas);
	      const double xctrueerr_mas = sqrt(evt2ecl_c*evt2ecl_c*xctrueerr_evt_mas*xctrueerr_evt_mas
						+ evt2ecl_s*evt2ecl_s*yctrueerr_evt_mas*yctrueerr_evt_mas);
	      const double yctrueerr_mas = sqrt(evt2ecl_s*evt2ecl_s*xctrueerr_evt_mas*xctrueerr_evt_mas
						+ evt2ecl_c*evt2ecl_c*yctrueerr_evt_mas*yctrueerr_evt_mas);

	      lcfile << xc_obs_mas << " " << xcerr_mas << " " << yc_obs_mas << " " << ycerr_mas << " " << flush; 
	      lcfile << xctrue_mas << " " << xctrueerr_mas << " " << yctrue_mas << " " << yctrueerr_mas << " " << flush; 

	      // Transform ecliptic EN offsets to ICRS tangent-plane offsets.
	      const double dra_cosdec_obs_mas = dRAc_from_eE * xc_obs_mas + dRAc_from_eN * yc_obs_mas;
	      const double ddec_obs_mas = dDec_from_eE * xc_obs_mas + dDec_from_eN * yc_obs_mas;
	      const double dra_cosdec_true_mas = dRAc_from_eE * xctrue_mas + dRAc_from_eN * yctrue_mas;
	      const double ddec_true_mas = dDec_from_eE * xctrue_mas + dDec_from_eN * yctrue_mas;
	      const double dra_cosdec_src_only_mas = dRAc_from_eE * src_only_e + dRAc_from_eN * src_only_n;
	      const double ddec_src_only_mas = dDec_from_eE * src_only_e + dDec_from_eN * src_only_n;
	      const double dra_cosdec_src_lens_mas = dRAc_from_eE * src_lens_e + dRAc_from_eN * src_lens_n;
	      const double ddec_src_lens_mas = dDec_from_eE * src_lens_e + dDec_from_eN * src_lens_n;

	      const double ra_obs_deg = ra_base_deg + dra_cosdec_obs_mas * mas_to_deg * inv_cos_dec_eq;
	      const double dec_obs_deg = dec_base_deg + ddec_obs_mas * mas_to_deg;
	      const double ra_true_deg = ra_base_deg + dra_cosdec_true_mas * mas_to_deg * inv_cos_dec_eq;
	      const double dec_true_deg = dec_base_deg + ddec_true_mas * mas_to_deg;
	      const double ra_src_only_deg = ra_base_deg + dra_cosdec_src_only_mas * mas_to_deg * inv_cos_dec_eq;
	      const double dec_src_only_deg = dec_base_deg + ddec_src_only_mas * mas_to_deg;
	      const double ra_src_lens_deg = ra_base_deg + dra_cosdec_src_lens_mas * mas_to_deg * inv_cos_dec_eq;
	      const double dec_src_lens_deg = dec_base_deg + ddec_src_lens_mas * mas_to_deg;
	      double ra_lens_primary_deg = NAN;
	      double dec_lens_primary_deg = NAN;
	      if(Event->xlens.size() > 0 && Event->ylens.size() > 0
		 && Event->xlens[0].size() > (size_t)i && Event->ylens[0].size() > (size_t)i)
		{
		  const double lens_evt_x = Event->xlens[0][i] * thE;
		  const double lens_evt_y = Event->ylens[0][i] * thE;
		  const double lens_e = evt2ecl_c * lens_evt_x - evt2ecl_s * lens_evt_y;
		  const double lens_n = evt2ecl_s * lens_evt_x + evt2ecl_c * lens_evt_y;
		  const double dra_cosdec_lens_mas = dRAc_from_eE * lens_e + dRAc_from_eN * lens_n;
		  const double ddec_lens_mas = dDec_from_eE * lens_e + dDec_from_eN * lens_n;
		  ra_lens_primary_deg = ra_base_deg + dra_cosdec_lens_mas * mas_to_deg * inv_cos_dec_eq;
		  dec_lens_primary_deg = dec_base_deg + ddec_lens_mas * mas_to_deg;
		}

	      // Lens-parallax term in ecliptic EN (mas), based on observer displacement
	      // relative to the event reference frame.
	      double D_L_kpc = NAN;
	      double lens_pllx_e_mas = 0.0;
	      double lens_pllx_n_mas = 0.0;
	      if(lens_idx >= 0 && lens_idx < int(Lenses->data.size()))
		{
		  D_L_kpc = Lenses->data[lens_idx][Lenses->DIST];
		  if(D_L_kpc > 0.0)
		    {
		      lens_pllx_e_mas = -Event->pllx[obsidx].Eshift[shiftedidx] / D_L_kpc;
		      lens_pllx_n_mas = -Event->pllx[obsidx].Nshift[shiftedidx] / D_L_kpc;
		    }
		}

	      // Expected lens sky position from base RA/Dec + proper motion + parallax.
	      const double bjd = Event->pllx[obsidx].epochs[shiftedidx];
	      double lens_pm_parallax_dra_cosdec_mas = NAN;
	      double lens_pm_parallax_ddec_mas = NAN;
	      if(lens_mu_ok && std::isfinite(bjd) && std::isfinite(Event->pllx[obsidx].tref))
		{
		  const double dt_years = (bjd - Event->pllx[obsidx].tref) * days_to_years;
		  const double dra_cosdec_pm_mas = lens_mu_ra * dt_years;
		  const double ddec_pm_mas = lens_mu_dec * dt_years;
		  const double dra_cosdec_pllx_mas = dRAc_from_eE * lens_pllx_e_mas + dRAc_from_eN * lens_pllx_n_mas;
		  const double ddec_pllx_mas = dDec_from_eE * lens_pllx_e_mas + dDec_from_eN * lens_pllx_n_mas;
		  lens_pm_parallax_dra_cosdec_mas = dra_cosdec_pm_mas + dra_cosdec_pllx_mas;
		  lens_pm_parallax_ddec_mas = ddec_pm_mas + ddec_pllx_mas;
		}

	      lcfile << setprecision(12) << ra_obs_deg << " " << dec_obs_deg << " " << flush;
	      lcfile << ra_true_deg << " " << dec_true_deg << " " << flush;
	      lcfile << ra_src_only_deg << " " << dec_src_only_deg << " " << flush;
	      lcfile << ra_src_lens_deg << " " << dec_src_lens_deg << " " << flush;
	      lcfile << ra_lens_primary_deg << " " << dec_lens_primary_deg << " " << flush;
	      lcfile << lens_pm_parallax_dra_cosdec_mas << " " << lens_pm_parallax_ddec_mas << " " << flush;

	      const double dra_cosdec_obs_lpllx_mas = dRAc_from_eE * (xc_obs_mas + lens_pllx_e_mas)
		+ dRAc_from_eN * (yc_obs_mas + lens_pllx_n_mas);
	      const double ddec_obs_lpllx_mas = dDec_from_eE * (xc_obs_mas + lens_pllx_e_mas)
		+ dDec_from_eN * (yc_obs_mas + lens_pllx_n_mas);
	      const double dra_cosdec_true_lpllx_mas = dRAc_from_eE * (xctrue_mas + lens_pllx_e_mas)
		+ dRAc_from_eN * (yctrue_mas + lens_pllx_n_mas);
	      const double ddec_true_lpllx_mas = dDec_from_eE * (xctrue_mas + lens_pllx_e_mas)
		+ dDec_from_eN * (yctrue_mas + lens_pllx_n_mas);

	      const double ra_obs_lpllx = ra_base_deg + dra_cosdec_obs_lpllx_mas * mas_to_deg * inv_cos_dec_eq;
	      const double dec_obs_lpllx = dec_base_deg + ddec_obs_lpllx_mas * mas_to_deg;
	      const double ra_true_lpllx = ra_base_deg + dra_cosdec_true_lpllx_mas * mas_to_deg * inv_cos_dec_eq;
	      const double dec_true_lpllx = dec_base_deg + ddec_true_lpllx_mas * mas_to_deg;
	      lcfile << ra_obs_lpllx << " " << dec_obs_lpllx << " " << flush;
	      lcfile << ra_true_lpllx << " " << dec_true_lpllx << " " << flush;

	      double obs_x_AU = Event->pllx[obsidx].sslocation[shiftedidx][0];
	      double obs_y_AU = Event->pllx[obsidx].sslocation[shiftedidx][1];
	      double obs_z_AU = Event->pllx[obsidx].sslocation[shiftedidx][2];

	      // Output lens distance, lens parallax term, parallax shifts, epochs,
	      // and observer ecliptic cartesian position (AU).
	      lcfile << setprecision(16) << D_L_kpc << " " << lens_pllx_e_mas << " " << lens_pllx_n_mas << " " << flush;
	      lcfile << Event->pllx[obsidx].tshift[shiftedidx] << " " << flush;
	      lcfile << Event->pllx[obsidx].ushift[shiftedidx] << " " << flush;
	      lcfile << bjd << " " << flush; 
	      lcfile << obs_x_AU << " " << obs_y_AU << " " << obs_z_AU << " " << flush;
	      lcfile << obs_x_AU << " " << obs_y_AU << " " << obs_z_AU << " " << flush;

      // Output astrometry diagnostics
	      if(Event->astrox1_raw.size() > (size_t)i)
		{
		  // Raw VBM output (theta_E - for coordinate system debugging)
		  lcfile << Event->astrox1_raw[i] << " " << Event->astrox2_raw[i] << " " << flush;
		  lcfile << src_only_e << " " << src_only_n << " " << flush;
		  lcfile << src_lens_e << " " << src_lens_n << " " << flush;
		  lcfile << src_lens_amb_e << " " << src_lens_amb_n << " " << flush;
	}
      else
	{
	  lcfile << 0.0 << " " << 0.0 << " " << flush;
	  lcfile << 0.0 << " " << 0.0 << " " << flush;
	  lcfile << 0.0 << " " << 0.0 << " " << flush;
	  lcfile << 0.0 << " " << 0.0 << " " << flush;
	}

      if(Paramfile->verbosity>=4)
	{
	  cout << "xsrc size " << Event->xsrc.size() << " " << Event->xsrc[0].size() << endl;
	  cout << "ysrc size " << Event->ysrc.size() << " " << Event->ysrc[0].size() << endl;
	  cout << "xlens size " << Event->xlens.size() << " " << Event->xlens[0].size() << endl;
	  cout << "ylens size " << Event->ylens.size() << " " << Event->ylens[0].size() << endl;
	}
	  
      for(int s=0;s<Event->nsrc;s++)
	{
	  lcfile << Event->xsrc[s][i] << " " << Event->ysrc[s][i] << " " << Event->mu_src[s][i] << " " << flush;
	}
      lcfile << flush;
      for(int l=0;l<Event->nlens;l++)
	{
	  lcfile << Event->xlens[l][i] << " " << Event->ylens[l][i] << " " << flush;
	}
      lcfile << flush;
	  
      if(ndF>0)
	{
	  for(int j=0;j<ndF;j++)
	    {
	      //fprintf(lcfile_ptr,"%.6g ",Event->dF[i+j*Event->nepochs]);
	      lcfile << Event->dF[i+j*Event->nepochs] << " ";
	    }
	  //for(int j=0;j<ndF;j++)
	  //{
	  //fprintf(lcfile_ptr,"%.6g ",Event->dF_debug[i+j*Event->nepochs]);
	  // }
	  //for(int j=0;j<ndF;j++)
	  // {
	  // fprintf(lcfile_ptr,"%.6g ",Event->dF_diff[i+j*Event->nepochs]);
	  // }
	}
      lcfile << endl;
      //fprintf(lcfile_ptr,"\n");
    } //end for nepochs
  lcfile.close();
      //  } //if lcfileptr
  //fclose(lcfile_ptr);
    
  if(Paramfile->verbosity>=4)
    {
      cout << "Writing extra lightcurve data" << endl;
      //output the lightcurve data
      if(lcdatafile_ptr!=NULL && fileOpen==1)
        {
	  for(i=0;i<Event->nepochs;i++)
            {
	      t=Event->epoch[i];
	      obsidx=Event->obsidx[i];
	      fprintf(lcdatafile_ptr, "%.11g %.11g %.12g %.12g %.12g\n ",
		      Event->epoch[i], Event->Atrue[i], Event->vbm_rootaccuracy[i], Event->vbm_squarecheck[i], Event->vbm_therr[i]);//0, 1, 2, 3,4
	    }
	    fclose(lcdatafile_ptr);
	}

      cout << "Writing extra lightcurve data done" << endl;
    }
}

void outputImages(struct event *Event, struct obsfilekeywords World[], struct slcat* Sources, struct filekeywords* Paramfile)
{
  string tmp1;
  string basefname;
  string extension;
  string oname;
  string imtype;
  int filter;
  double mag;
  double peaktime;
  double background;

  if(!Event->outputthis || !Paramfile->outputImages) return;
  else if((Paramfile->outputOnErr || Paramfile->outputOnDet 
	   || Paramfile->outputOnAll)==0) return;
  else
    {
      if((Event->lcerror+Event->deterror))
	{
	  if(!Paramfile->outputOnErr) return;
	}
      else if((Event->detected))
	{
	  if(!Paramfile->outputOnDet) return;
	}
      else if((!Event->detected))
	{
	  if(!Paramfile->outputOnAll) return;
	}
    }

  if(Event->detected) extension=string(".det");
  else if(Event->lcerror||Event->deterror) extension=string(".err");
  else extension=string(".all");

  if(Paramfile->choosefield<0)
    {
      tmp1 = Paramfile->outputdir +  Paramfile->run_name + "_"
	+ to_string(Paramfile->instance) + "_" + to_string(Event->id);
    }
  else
    {
      // Build the same base name as the other branch — avoid stray format literal and
      // do not use the comma operator. Keep it as a plain concatenation.
      tmp1 = Paramfile->outputdir + Paramfile->run_name + "_"
        + to_string(Paramfile->instance) + "_" + to_string(Paramfile->choosefield) + "_"
        + to_string(Event->id);
    }
  basefname=tmp1;

  //output test images
  for(int obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    {
      if(Event->nepochsvec[obsidx+1]-Event->nepochsvec[obsidx]<=0) continue;
      // do not mutate tmp1/basefname; create a small suffix for this observation
      string obs_suffix = "." + to_string(obsidx) + "_";

      filter = World[obsidx].filter;

      //first the baseline image
      imtype="base";
      oname = basefname + obs_suffix + imtype + extension + ".fits";

      mag = Sources->mags[Event->source][filter];

      //calculate the background at random time - may as well be first point
      background = Event->backmag[Event->nepochsvec[obsidx]];
      World[obsidx].im.set_background(background);
      World[obsidx].im.addbg();

      //add the baseline source
      World[obsidx].im.addstar(Event->xsub[obsidx], Event->ysub[obsidx], mag);
      World[obsidx].im.reset_detector();
      World[obsidx].im.expose(World[obsidx].exptime[0], 
			      World[obsidx].nstack[0]);
      World[obsidx].im.write_fits(oname,true);
      World[obsidx].im.substar(Event->xsub[obsidx], Event->ysub[obsidx], mag);
      World[obsidx].im.subbg();

      //last the peak image
      imtype=string("peak");
      oname = basefname + tmp1 + imtype + extension + string(".fits");
 
      //initialize to constant background specified in the detector file

      if(Event->Amax < 0) 
	{
	  //if the peak wasn't recorded
	  Event->Amax = (sqr(Event->u0)+2.0)/sqrt(sqr(Event->u0)*(sqr(Event->u0)+4.0));
	  peaktime = Event->t0;
	}
      else
	{
	  //cout << "peaktime at epoch " << Event->peakpoint << endl;
	  peaktime = Event->epoch[Event->peakpoint];
	}
 
      //background
      World[obsidx].im.set_background(Event->backmag[Event->peakpoint]);
      World[obsidx].im.addbg();

      //add the peak source
      World[obsidx].im.addstar(Event->xsub[obsidx], Event->ysub[obsidx], 
				   mag-2.5*log10(Event->Amax));
      World[obsidx].im.reset_detector();
      World[obsidx].im.expose(World[obsidx].exptime[0], 
			      World[obsidx].nstack[0]);
      World[obsidx].im.write_fits(oname,true);
      World[obsidx].im.substar(Event->xsub[obsidx], Event->ysub[obsidx], 
			       mag-2.5*log10(Event->Amax));
      World[obsidx].im.subbg();
      //cout << "Amax = " << Event->Amax << endl;
    }

    if(Paramfile->verbosity>=3)
    {
      cout << "Exiting output lightcurve" << endl;
    }
}
