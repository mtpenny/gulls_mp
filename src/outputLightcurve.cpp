#include "outputLightcurve.h"
#include "zodiacalLight.h"
#include "astroFns.h"
#include "coords.h"
#include <iomanip>
#include <sstream>
#include <fstream>
#include <limits>
#include <utility>

#define DEBUGVAR 0

// STRUCTURES.H
//-------------------------------
// Astrometry stage diagnostics (event frame, theta_E)
// omLightcurveGenerator.cpp
/*
vector<double> xc_srcs_only;  // blended apparent source centroid (without lens light contribution), in theta E units
vector<double> yc_srcs_only;
vector<double> xc_src_lens;  // blended apparent source centroid (with lens light contribution), in theta E units
vector<double> yc_src_lens;
vector<vector<double> > astrox1_raw;  // raw vbm frame outputs (per source), in theta E units
vector<vector<double> > astrox2_raw;
vector<vector<double> > astrox_raw;  // vbm frame outputs converted to event_frame astrometric shifts (per source), in theta E units
vector<vector<double> > astroy_raw;
vector<double> src_flux_total; // total source flux (for calculating blended centroid), in units of the unmagnified source flux
// photometry.cpp
vector<double> n_src_lens_mas; // blended ecliptic centroid shift in mas (includes lens proper motion and parallax)
vector<double> e_src_lens_mas;
vector<double> lambda_noiseless_deg; // absolute blended apparent source centroid in ecliptic coordinates, in degrees
vector<double> beta_noiseless_deg;

// Public sky astrometry products (degrees / mas)
// photometry.cpp
vector<double> ra_noiseless_deg;  // absolute blended apparent source centroid in RA/Dec, in degrees
vector<double> dec_noiseless_deg;
vector<double> ra_measured_deg;  // added measurement noise, in the tangent plane, in degrees
vector<double> dec_measured_deg;  // added measurement noise in degrees
vector<double> sigma_ast_mas;  // symmetric tangent plane astrometric error (1-sigma), in mas
vector<double> ra_err_deg;  // RA error component in degrees on the celestial sphere
vector<double> dec_err_deg;  //dec error component in degrees
*/

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
  lcfile << "#VBM_function: name="
	 << (Event->VBM_function.empty() ? "unknown" : Event->VBM_function)
	 << endl;

  // Astrometry metadata
  coords c;
  double lambda0 = 0.0;
  double beta0 = 0.0;
  c.ad2ecl(Event->ra, Event->dec, &lambda0, &beta0);
  const double murel_ref = (Event->tE_r != 0.0 ? Event->thE/Event->tE_r*DAYINYR : 0.0);
  lcfile << "#Astrometry_Frame: "
	 << setprecision(12)
	 << "RA_rad=" << Event->ra << " Dec_rad=" << Event->dec
	 << " RA_deg=" << Event->ra*TO_DEG << " Dec_deg=" << Event->dec*TO_DEG
	 << " lambda0_deg=" << lambda0*TO_DEG << " beta0_deg=" << beta0*TO_DEG
	 << " thE_mas=" << Event->thE << " tref=" << Event->tref
	 << " frame=ecliptic_EN_barycenter";
  lcfile << endl;
  lcfile << "#Astrometry_PM: units=mas_per_yr";
  for(int i=0;i<Paramfile->numobservatories;i++)
    {
      double pm_e = 0.0;
      double pm_n = 0.0;
      if(i < int(Event->pllx.size()))
	{
	  pm_e = murel_ref * Event->pllx[i].mulam_r;
	  pm_n = murel_ref * Event->pllx[i].mubet_r;
	}
      lcfile << " obs" << i << "_pmE=" << pm_e << " obs" << i << "_pmN=" << pm_n;
    }
  lcfile << endl;
  lcfile << "#Astrometry_Contract: version=v2 model_frame=absolute event_true_unit=thetaE noise_frame=ecliptic_tangent" << endl;
  lcfile << "#Astrometry_Columns:"
    << " event_x_thetaE=blended_sources_lenses_x_thetaE"
    << " event_y_thetaE=blended_sources_lenses_y_thetaE"
   << " ecliptic_lambda_noiseless_deg=ecliptic_lambda_noiseless_deg"
   << " ecliptic_beta_noiseless_deg=ecliptic_beta_noiseless_deg"
   << " sky_ra_noiseless_deg=true_RA_deg"
   << " sky_dec_noiseless_deg=true_Dec_deg"
   << " sky_ra_measured_deg=measured_RA_deg"
   << " sky_dec_measured_deg=measured_Dec_deg"
   << " sky_sigma_mas=sigma_astrometric_mas"
   << " sky_ra_err_deg=measured_RA_error_deg"
   << " sky_dec_err_deg=measured_Dec_error_deg"
   << " source_blend_total_flux=fractional_total_flux";
  lcfile << endl;
  lcfile << "#Astrometry_BAGLE: model_frame=absolute blendless_columns=none lens_columns=lens0_x,lens0_y lens_frame=xy_thetaE" << endl;

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
    "blended_sources_only_x_thetaE" << " " << "blended_sources_only_y_thetaE" << " " <<
    "blended_sources_lenses_x_thetaE" << " " << "blended_sources_lenses_y_thetaE" << " " <<
    "ecliptic_lambda_noiseless_deg" << " " << "ecliptic_beta_noiseless_deg" << " " <<
    "true_RA_deg" << " " << "true_Dec_deg" << " " <<
    "measured_RA_deg" << " " << "measured_Dec_deg" << " " <<
    "sigma_astrometric_mas" << " " << "measured_RA_error_deg" << " " << "measured_Dec_error_deg" << " " <<
    "fractional_total_source_flux" << " " <<  // fstot/fs1
    "parallax_shift_t" << " " << "parallax_shift_u" << " " << "BJD" << " " <<
    "parallax_shift_x" << " " << "parallax_shift_y" << " " << "parallax_shift_z" << " ";
  for(int i=0;i<Event->nsrc;i++)
    {
      lcfile << "vbm_astrox1_source" << i << "_thE" << " " << "vbm_astrox2_source" << i << "_thE" << " ";
    }
  for(int i=0;i<Event->nsrc;i++)
    {
      lcfile << "lenses_source" << i << "_x_thE" << " " << "lenses_source" << i << "_y_thE" << " ";
    }
  for(int i=0;i<Event->nsrc;i++)
    {
      lcfile << "source" << i << "_x_thE" << " " << "source" << i << "_y_thE" << " " << "source" << i << "_mu" << " ";  // mu is magnification
    }
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

  auto safe_get = [](const vector<double>& v, int idx) -> double
  {
    return (idx >= 0 && idx < int(v.size()))
      ? v[idx]
      : std::numeric_limits<double>::quiet_NaN();
  };

  auto safe_get2d = [](const vector<vector<double> >& v, const pair<int,int>& idx2) -> double
  {
    const int i0 = idx2.first;
    const int i1 = idx2.second;
    return (i0 >= 0 && i0 < int(v.size()) && i1 >= 0 && i1 < int(v[i0].size()))
      ? v[i0][i1]
      : std::numeric_limits<double>::quiet_NaN();
  };

  //if(lcfile_ptr!=NULL && fileOpen==1)
  //  {
  for(i=0;i<Event->nepochs;i++)
    {
      t=Event->epoch[i];
      obsidx=Event->obsidx[i];
      shiftedidx = i-Event->nepochsvec[obsidx];
	    
	      lcfile << setprecision(16) << Event->epoch[i] << " " << Event->Aobs[i] << " " << Event->Aerr[i] << " " << flush;
	      lcfile << Event->Atrue[i] << " " << Event->Atrueerr[i] << " " << obsidx << " " << flush; 
	      lcfile << (Event->nosat[i]?0:1) << " " << Event->Afit[i] << " " << flush;
        lcfile << safe_get(Event->xc_srcs_only, i) << " " << safe_get(Event->yc_srcs_only, i) << " " << flush;
        lcfile << safe_get(Event->xc_src_lens, i) << " " << safe_get(Event->yc_src_lens, i) << " " << flush;
        lcfile << safe_get(Event->lambda_noiseless_deg, i) << " " << safe_get(Event->beta_noiseless_deg, i) << " " << flush;
        lcfile << safe_get(Event->ra_noiseless_deg, i) << " " << safe_get(Event->dec_noiseless_deg, i) << " " << flush;
        lcfile << safe_get(Event->ra_measured_deg, i) << " " << safe_get(Event->dec_measured_deg, i) << " " << flush;
        lcfile << safe_get(Event->sigma_ast_mas, i) << " " << safe_get(Event->ra_err_deg, i) << " " << safe_get(Event->dec_err_deg, i) << " " << flush;
        lcfile << safe_get(Event->src_flux_total, i) << " " << flush;
	      lcfile << Event->pllx[obsidx].tshift[shiftedidx] << " " << flush;
	      lcfile << Event->pllx[obsidx].ushift[shiftedidx] << " " << flush;
	      lcfile << Event->pllx[obsidx].epochs[shiftedidx] << " " << flush; 
	      lcfile << Event->pllx[obsidx].sslocation[shiftedidx][0] << " " << flush;
	      lcfile << Event->pllx[obsidx].sslocation[shiftedidx][1] << " " << flush; 
	      lcfile << Event->pllx[obsidx].sslocation[shiftedidx][2] << " " << flush;
      for(int s=0;s<Event->nsrc;s++)
  {
    lcfile << safe_get2d(Event->astrox1_raw, make_pair(s,i)) << " " << safe_get2d(Event->astrox2_raw, make_pair(s,i)) << " " << flush;
  }
      for(int s=0;s<Event->nsrc;s++)
  {
    lcfile << safe_get2d(Event->astrox_raw, make_pair(s,i)) << " " << safe_get2d(Event->astroy_raw, make_pair(s,i)) << " " << flush;
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
      tmp1 = Paramfile->outputdir + Paramfile->run_name + "_"
	+ to_string(Paramfile->instance) + "_" + to_string(Paramfile->choosefield) + "_"
	+  to_string(Event->id);
    }
  basefname=tmp1;

  //output test images
  for(int obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    {
      if(Event->nepochsvec[obsidx+1]-Event->nepochsvec[obsidx]<=0) continue;
      const string obs_prefix = basefname + "." + to_string(obsidx) + "_";

      filter = World[obsidx].filter;

      //first the baseline image
      imtype="base";
      oname = obs_prefix + imtype + extension + ".fits";

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
      oname = obs_prefix + imtype + extension + string(".fits");
 
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
