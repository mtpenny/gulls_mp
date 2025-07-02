#include "outputLightcurve.h"
#include "zodiacalLight.h"
#include "astroFns.h"

#define DEBUGVAR 0

void outputLightcurve(struct event *Event, struct obsfilekeywords World[], struct filekeywords* Paramfile, struct slcat* Sources, struct slcat* Lenses)
{
  void muVisibility(double *mu, double rs, double z0, double ld1);
  char lcfname[1000];
  FILE *lcfile_ptr;
  FILE *lcdatafile_ptr;
  int fileOpen=0;
  int i, obsidx;
  char data[1000];
  char tmp[30];
  char extension[4];
  char lcdatafname[1000];
  double t;
  //double u0_ps,t0_ps,tE_ps,t2_ps,u_ps,psmag,psamp;
  //double u0_fs,t0_fs,tE_fs,t2_fs,u_fs,fsmag,fsamp;

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

  if(Event->detected) strcpy(extension,"det");
  else if(Event->lcerror||Event->deterror) strcpy(extension,"err");
  else strcpy(extension,"all");

  if(Paramfile->choosefield<0)
    {
      sprintf(lcfname, "%s%s_%d_%d.%s.lc", Paramfile->outputdir,
	      Paramfile->run_name, Event->instance, Event->id, extension);
    }
  else
    {
      sprintf(lcfname, "%s%s_%d_%d_%d.%s.lc", Paramfile->outputdir,
	      Paramfile->run_name, Event->instance, Paramfile->choosefield, 
	      Event->id, extension);
    }
  if(DEBUGVAR) printf("lcname: %s\n",lcfname);

  if(Event->nepochs>0)
    {
      lcfile_ptr = fopen(lcfname,"w");
      fileOpen=1;
    }
  else return;

  if(Paramfile->verbosity>=4)
    {
      if(Paramfile->choosefield<0)
        {
          sprintf(lcdatafname, "%s%s_%d_%d.%s.lcdata", Paramfile->outputdir,
              Paramfile->run_name, Event->instance, Event->id, extension);
        }
      else
        {
          sprintf(lcdatafname, "%s%s_%d_%d_%d.%s.lcdata", Paramfile->outputdir,
              Paramfile->run_name, Event->instance, Paramfile->choosefield,
              Event->id, extension);
        }
      if(DEBUGVAR) printf("lcdataname: %s\n",lcdatafname);
      if(Event->nepochs>0)
        {
          lcdatafile_ptr = fopen(lcdatafname,"w");
          fileOpen=1;
	  fprintf(lcdatafile_ptr, "time Atrue rootaccuracy squarecheck therr \n");
	  //fprintf(lcdatafile_ptr, "time mag_old_lcgen mag_vbm dif_over_mag\n");
        }
     }


  //output header information

  //blending
  strcpy(data,"#fs: ");
  for(int i=0;i<Paramfile->numobservatories;i++)
    {
      sprintf(tmp,"%g ",Event->fs[i]);
      strcat(data,tmp);
    }
  fprintf(lcfile_ptr,"%s\n",data);

  //Source magnitudes
  strcpy(data,"#Sourcemag: ");
  int sn = Event->source;
  for(int i=0;i<Paramfile->Nfilters;i++)
    {
      sprintf(tmp,"%g ",Sources->mags[sn][i]);
      strcat(data,tmp);
    }
  fprintf(lcfile_ptr,"%s\n",data);

  //Source data
  strcpy(data,"#Sourcedata: ");
  sprintf(tmp,"%d ",sn); strcat(data,tmp);
  //for(int i=0;i<sOutputCols;i++)
  //  {
  //    sprintf(tmp,"%g ", Sources->data[sn][sOutputColumns[i]]);
  //    strcat(data,tmp);
  //}
  for(int i=0;i<Sources->data[sn].size();i++)
	{
	  sprintf(tmp,"%g ", Sources->data[sn][i]);
      strcat(data,tmp);
	}
  fprintf(lcfile_ptr,"%s\n",data);

  //Observatory magnitudes
  strcpy(data,"#Obssrcmag: ");
  for(int i=0;i<Paramfile->numobservatories;i++)
    {
      sprintf(tmp,"%g ",Sources->mags[sn][World[i].filter]);
      strcat(data,tmp);
    }
  fprintf(lcfile_ptr,"%s\n",data);


  //Lens magnitudes
  strcpy(data,"#Lensmag: ");
  int ln = Event->lens;
  for(int i=0;i<Paramfile->Nfilters;i++)
    {
      sprintf(tmp,"%g ",Lenses->mags[ln][i]);
      strcat(data,tmp);
    }
  fprintf(lcfile_ptr,"%s\n",data);

  //Lens data
  strcpy(data,"#Lensdata: ");
  sprintf(tmp,"%d ",ln); strcat(data,tmp);
  //for(int i=0;i<lOutputCols;i++)
  //  {
  //    sprintf(tmp,"%g ", Lenses->data[ln][lOutputColumns[i]]);
  //    strcat(data,tmp);
  //  }
  for(int i=0;i<Lenses->data[ln].size();i++)
	{
	  sprintf(tmp,"%g ", Lenses->data[ln][i]);
      strcat(data,tmp);
	}
  fprintf(lcfile_ptr,"%s\n",data);

  //Observatory magnitudes
  strcpy(data,"#Obslensmag: ");
  for(int i=0;i<Paramfile->numobservatories;i++)
    {
      sprintf(tmp,"%g ",Lenses->mags[ln][World[i].filter]);
      strcat(data,tmp);
    }
  fprintf(lcfile_ptr,"%s\n",data);


  //Planet data
  strcpy(data,"#Planet: ");
  for(int i=0;i<NPLANETINPUT+NPLANETDERIV;i++)
    {
      sprintf(tmp,"%g ",Event->params[i]);
      strcat(data,tmp);
    }
  fprintf(lcfile_ptr,"%s\n",data);

  //Microlensing data
  strcpy(data,"#Event: ");
  sprintf(tmp,"%g ",Event->u0); strcat(data,tmp);
  sprintf(tmp,"%g ",Event->alpha); strcat(data,tmp);
  sprintf(tmp,"%.12g ",Event->t0); strcat(data,tmp);
  sprintf(tmp,"%.12g ",Event->tcroin); strcat(data,tmp);
  sprintf(tmp,"%g ",Event->ucroin); strcat(data,tmp);
  sprintf(tmp,"%g ",Event->rcroin); strcat(data,tmp);
  sprintf(tmp,"%g ",Event->tE_r); strcat(data,tmp);
  sprintf(tmp,"%g ",Event->rs); strcat(data,tmp);

  fprintf(lcfile_ptr,"%s\n",data);

  //Observatory groups
  for(int obsgroup=0; obsgroup<int(Event->obsgroups.size()); obsgroup++)
    {
      strcpy(data,"#Obsgroup: ");
      sprintf(tmp,"%d ",obsgroup); strcat(data,tmp);
      sprintf(tmp,"%g ",Event->flatchi2[obsgroup]); strcat(data,tmp);
      sprintf(tmp,"%d ",Event->flag_needFS[obsgroup]); strcat(data,tmp);
      sprintf(tmp,"%g ", 
	      (Event->flag_needFS[obsgroup]?Event->FSPL[obsgroup].chisq:Event->PSPL[obsgroup].chisq)); strcat(data,tmp);
      //Members
      for(int grpidx=0; grpidx<int(Event->obsgroups[obsgroup].size()); grpidx++)
	{
	  sprintf(tmp,"%d ",Event->obsgroups[obsgroup][grpidx]); strcat(data,tmp);
	}
      fprintf(lcfile_ptr,"%s\n",data);
    }


  int ndF = int(Event->dF.size()) / Event->nepochs;
  // ─────────────── header ───────────────
  {
    static const char* baseCols[] = {
      "Simulation_time", "measured_relative_flux", "measured_relative_flux_error",
      "true_relative_flux",  "true_relative_flux_error",    "observatory_code",
      "saturation_flag",     "best_single_lens_fit",        "parallax_shift_t",
      "parallax_shift_u",    "BJD",                         "source_x",
      "source_y",            "lens1_x",                     "lens1_y",
      "lens2_x",             "lens2_y",                     "parallax_shift_x",
      "parallax_shift_y",    "parallax_shift_z"
    };
    int nBase = sizeof(baseCols) / sizeof(baseCols[0]);
    for(int i = 0; i < nBase; ++i) {
      fprintf(lcfile_ptr, "%s ", baseCols[i]);
    }
    if (ndF > 0) {
      string dF_nopllx = string("dF_t0 dF_tE dF_u0 dF_alpha dF_s dF_q dF_rs");
      string dF_pllx   = string(" dF_piEN dF_piEE");
      vector<string> parstrings;
      int nfixpar = 7 + (Paramfile->pllxMultiplyer ? 2 : 0);
      if (Paramfile->pllxMultiplyer)
        split(dF_nopllx + dF_pllx, parstrings);
      else
        split(dF_nopllx, parstrings);
      
      stringstream ss;
      for(size_t obsgroup = 0; obsgroup < Event->obsgroups.size(); ++obsgroup) {
        int nobs    = int(Event->obsgroups[obsgroup].size());
        int nparams = nfixpar + 2*nobs;

        for(int idx = 0; idx < nparams; ++idx) {
          ss.str("");  ss.clear();
          ss << "ObsGroup_" << obsgroup << "_";

          if (idx < nfixpar) {
            ss << parstrings[idx];
          } else {
            int grpidx = (idx - nfixpar) / 2;
            int obsidx = Event->obsgroups[obsgroup][grpidx];

            if (((idx - nfixpar) % 2) == 0)
              ss << "dF_Fbase" << obsidx;
            else
              ss << "dF_fs"     << obsidx;
          }  
     
     	  fprintf(lcfile_ptr, "%s ", ss.str().c_str());
        }
      }
    }
    // finish the header line
    fprintf(lcfile_ptr, "\n");
  }
  //output the lightcurve
  int shiftedidx;

  if(lcfile_ptr!=NULL && fileOpen==1)
    {
      for(i=0;i<Event->nepochs;i++)
	{
	  t=Event->epoch[i];
	  obsidx=Event->obsidx[i];
	  shiftedidx = i-Event->nepochsvec[obsidx];
	  
	  fprintf(lcfile_ptr, "%.12g %.8g %g %.12g %g %d %d %.8g %.6g %.6g %16.7f %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g ",
		  Event->epoch[i], Event->Aobs[i], Event->Aerr[i], //0, 1, 2
		  Event->Atrue[i], Event->Atrueerr[i], obsidx, //3, 4, 5
		  (Event->nosat[i]?0:1), Event->Afit[i], //6, 7
		  //Event->pllx[obsidx].tshift(Event->jdepoch[i]), //8
		  //Event->pllx[obsidx].ushift(Event->jdepoch[i]), //9
		  //Event->pllx[obsidx].epochs[Event->jdepoch[i]],  //10
		  Event->pllx[obsidx].tshift[shiftedidx],
		  Event->pllx[obsidx].ushift[shiftedidx],
		  Event->pllx[obsidx].epochs[shiftedidx],
		  Event->xs[i], Event->ys[i], Event->xl1[i],     //11, 12, 13
		  Event->yl1[i], Event->xl2[i], Event->yl2[i], //14, 15, 16
		  //Event->pllx[obsidx].sslocation[Event->jdepoch[i]][0], //17
		  //Event->pllx[obsidx].sslocation[Event->jdepoch[i]][1], //18
		  //Event->pllx[obsidx].sslocation[Event->jdepoch[i]][2]); //19
		  Event->pllx[obsidx].sslocation[shiftedidx][0], //17
		  Event->pllx[obsidx].sslocation[shiftedidx][1], //18
		  Event->pllx[obsidx].sslocation[shiftedidx][2]); //19
		    
	  
	  if(ndF>0)
	    {
	      for(int j=0;j<ndF;j++)
		{
		  fprintf(lcfile_ptr,"%.6g ",Event->dF[i+j*Event->nepochs]);
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
	  fprintf(lcfile_ptr,"\n");
	}
      
      fclose(lcfile_ptr);
    }
  if(Paramfile->verbosity>=4)
    {
      //output the lightcurve data
      if(lcdatafile_ptr!=NULL && fileOpen==1)
        {
	  for(i=0;i<Event->nepochs;i++)
            {
		t=Event->epoch[i];
		obsidx=Event->obsidx[i];
          	fprintf(lcdatafile_ptr, "%.11g %.11g %.12g %.12g %g\n ",
                  Event->epoch[i], Event->Atrue[i], Event->vbm_rootaccuracy[i], Event->vbm_squarecheck[i], Event->vbm_therr[i]);//0, 1, 2, 3,4
	    }
	    fclose(lcdatafile_ptr);
	}}
}

void outputImages(struct event *Event, struct obsfilekeywords World[], struct slcat* Sources, struct filekeywords* Paramfile)
{
  char tmp1[1000];
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
      sprintf(tmp1,"%s%s_%d_%d",Paramfile->outputdir, Paramfile->run_name,
	      Event->instance,  Event->id);
    }
  else
    {
      sprintf(tmp1,"%s%s_%d_%d_%d",Paramfile->outputdir, Paramfile->run_name,
	      Event->instance, Paramfile->choosefield, Event->id);
    }
  basefname=tmp1;

  //output test images
  for(int obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    {
      if(Event->nepochsvec[obsidx+1]-Event->nepochsvec[obsidx]<=0) continue;
      sprintf(tmp1,".%d_",obsidx);

      filter = World[obsidx].filter;

      //first the baseline image
      imtype=string("base");
      oname = basefname + tmp1 + imtype + extension + string(".fits");

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
}
