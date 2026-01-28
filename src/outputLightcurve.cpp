#include "outputLightcurve.h"
#include "zodiacalLight.h"
#include "astroFns.h"
#include "constants.h"
#include <iomanip>
#include <sstream>
#include <fstream>

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

  // Astrometry reference frame definition
  // Event frame origin is at catalog (RA, Dec) - does NOT move with lens proper motion
  // Event frame (x,y) orientation relative to sky (E,N) is UNCERTAIN - validate with plots!
  // RA/Dec columns assume x~East, y~North but this may be wrong depending on alpha convention
  // Two RA/Dec versions: without lens parallax (_deg) and with lens parallax attempt (_lpllx_deg)
  lcfile << "#Astrometry_Frame: ";
  lcfile << setprecision(12) << "RA_rad=" << Event->ra << " Dec_rad=" << Event->dec 
         << " RA_deg=" << Event->ra * r2d << " Dec_deg=" << Event->dec * r2d
         << " thE_mas=" << Event->thE << " t0=" << Event->t0;
  lcfile << endl;
  lcfile << "#Astrometry_Frame: origin=catalog_position, xy_orientation=UNCERTAIN(validate!)" << endl;

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
    "RA_centroid_deg" << " " << "Dec_centroid_deg" << " " <<             // observed, no lens parallax
    "RA_centroid_true_deg" << " " << "Dec_centroid_true_deg" << " " <<   // true, no lens parallax
    "RA_centroid_lpllx_deg" << " " << "Dec_centroid_lpllx_deg" << " " << // observed, WITH lens parallax
    "RA_true_lpllx_deg" << " " << "Dec_true_lpllx_deg" << " " <<         // true, WITH lens parallax
    "lens_dist_kpc" << " " <<                                            // lens distance for user validation
    "lens_parallax_x_mas" << " " << "lens_parallax_y_mas" << " " <<       // lens parallax shift (mas) in event-frame x/y
    "parallax_shift_t" << " " << "parallax_shift_u" << " " <<    "BJD" << " " <<
    "parallax_shift_x" << " " << "parallax_shift_y" << " " <<    "parallax_shift_z" << " ";
  // Astrometry diagnostic columns (mas unless noted)
  // Raw VBM output (VBM's internal frame, units: theta_E for debugging)
  lcfile << "vbm_astrox1_raw_thE" << " " << "vbm_astrox2_raw_thE" << " ";
  // Centroid at each blending step (mas)
  lcfile << "centroid_src_x_mas" << " " << "centroid_src_y_mas" << " ";           // sources only
  lcfile << "centroid_src_lens_x_mas" << " " << "centroid_src_lens_y_mas" << " "; // + lenses
  lcfile << "centroid_final_x_mas" << " " << "centroid_final_y_mas" << " ";       // + ambient = true
  
  // Per-source positions (event frame, units: theta_E) and magnifications
  // WARNING: x,y orientation relative to E,N is uncertain - validate with plots!
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
	//fprintf(lcfile_ptr, "%.12g %.8g %g %.12g %g %d %d %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.6g %.6g %16.7f %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g ",
      lcfile << setprecision(16) << Event->epoch[i] << " " << Event->Aobs[i] << " " << Event->Aerr[i] << " " << flush;
      lcfile << Event->Atrue[i] << " " << Event->Atrueerr[i] << " " << obsidx << " " << flush; 
      lcfile << (Event->nosat[i]?0:1) << " " << Event->Afit[i] << " " << flush;
      // Observed centroid (xc, yc) and errors are in mas; true centroid (xctrue, yctrue) in θ_E needs conversion
      lcfile << Event->xc[i] << " " << Event->xcerr[i] << " " << Event->yc[i] << " " << Event->ycerr[i] << " " << flush; 
      lcfile << Event->xctrue[i] * thE << " " << Event->xctrueerr[i] << " " << Event->yctrue[i] * thE << " " << Event->yctrueerr[i] << " " << flush; 
      
      // Convert centroid from event frame (mas) to RA/Dec (degrees)
      // WARNING: Event frame orientation (x,y) relative to (E,N) is uncertain - validate with plots!
      // Assuming x ~ East, y ~ North for now. If wrong, these RA/Dec values will be nonsense.
      // Origin at (Event->ra, Event->dec) which is the catalog position at t_ref (NOT current lens position)
      double ra_base_deg = Event->ra * r2d;  // radians to degrees
      double dec_base_deg = Event->dec * r2d;
      double cos_dec = cos(Event->dec);
      double mas_to_deg = 1.0 / (3600.0 * 1000.0);
      
      // --- VERSION 1: No lens parallax (centroid offset relative to catalog position) ---
      double ra_obs_deg = ra_base_deg + Event->xc[i] * mas_to_deg / cos_dec;
      double dec_obs_deg = dec_base_deg + Event->yc[i] * mas_to_deg;
      lcfile << setprecision(12) << ra_obs_deg << " " << dec_obs_deg << " " << flush;
      
      double xctrue_mas = Event->xctrue[i] * thE;
      double yctrue_mas = Event->yctrue[i] * thE;
      double ra_true_deg = ra_base_deg + xctrue_mas * mas_to_deg / cos_dec;
      double dec_true_deg = dec_base_deg + yctrue_mas * mas_to_deg;
      lcfile << ra_true_deg << " " << dec_true_deg << " " << flush;
      
      // --- VERSION 2: With lens parallax (observer position shifts apparent lens position) ---
      // sslocation = sun-to-observer vector projected on sky (AU), assuming [0]=x, [1]=y in event frame
      // Lens parallax: lens appears shifted by -(observer_pos) / D_lens
      // Sign convention: if observer is at +x from sun, lens at finite distance appears at +x relative to infinity
      int ln = Event->lens;
      double D_L_kpc = Lenses->data[ln][Lenses->DIST];
      double D_L_AU = D_L_kpc * 206265.0;  // 1 kpc ≈ 206265 AU
      
      double obs_x_AU = Event->pllx[obsidx].sslocation[shiftedidx][0];
      double obs_y_AU = Event->pllx[obsidx].sslocation[shiftedidx][1];
      
      // Lens parallax shift (radians): Δθ = -obs_pos / D_lens (standard parallax convention)
      // Actually, closer objects shift WITH observer motion, so Δθ = +obs_pos / D_lens? 
      // Outputting both signs would be silly - let's use standard: parallax = baseline/distance
      // where baseline points FROM observer TO sun, so shift = -obs_pos / D
      double pllx_x_rad = -obs_x_AU / D_L_AU;
      double pllx_y_rad = -obs_y_AU / D_L_AU;
      double pllx_x_mas = pllx_x_rad * r2d * 3600.0 * 1000.0;
      double pllx_y_mas = pllx_y_rad * r2d * 3600.0 * 1000.0;
      
      // Add lens parallax to base position, then add centroid offset
      double ra_obs_lpllx = ra_base_deg + (pllx_x_mas + Event->xc[i]) * mas_to_deg / cos_dec;
      double dec_obs_lpllx = dec_base_deg + (pllx_y_mas + Event->yc[i]) * mas_to_deg;
      lcfile << ra_obs_lpllx << " " << dec_obs_lpllx << " " << flush;
      
      double ra_true_lpllx = ra_base_deg + (pllx_x_mas + xctrue_mas) * mas_to_deg / cos_dec;
      double dec_true_lpllx = dec_base_deg + (pllx_y_mas + yctrue_mas) * mas_to_deg;
      lcfile << ra_true_lpllx << " " << dec_true_lpllx << " " << flush;
      
      // Output lens distance and lens parallax shift for user validation
      lcfile << D_L_kpc << " " << pllx_x_mas << " " << pllx_y_mas << " " << setprecision(6) << flush;
      
      lcfile << Event->pllx[obsidx].tshift[shiftedidx] << " " << flush;
      lcfile << Event->pllx[obsidx].ushift[shiftedidx] << " " << flush;
      lcfile << Event->pllx[obsidx].epochs[shiftedidx] << " " << flush; 
      lcfile << Event->pllx[obsidx].sslocation[shiftedidx][0] << " " << flush;
      lcfile << Event->pllx[obsidx].sslocation[shiftedidx][1] << " " << flush; 
      lcfile << Event->pllx[obsidx].sslocation[shiftedidx][2] << " " << flush;

      // Output astrometry diagnostics
      if(Event->astrox1_raw.size() > (size_t)i)
	{
	  // Raw VBM output (theta_E - for coordinate system debugging)
	  lcfile << Event->astrox1_raw[i] << " " << Event->astrox2_raw[i] << " " << flush;
	  // Centroid at each blending step (converted to mas)
	  lcfile << Event->xc_src_only[i] * thE << " " << Event->yc_src_only[i] * thE << " " << flush;
	  lcfile << Event->xc_src_lens[i] * thE << " " << Event->yc_src_lens[i] * thE << " " << flush;
	  lcfile << Event->xc_src_lens_amb[i] * thE << " " << Event->yc_src_lens_amb[i] * thE << " " << flush;
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
