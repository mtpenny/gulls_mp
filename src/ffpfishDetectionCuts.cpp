#include "lightcurveFitter.h"
#include "fisher.h"
#include "fisherinversion.h"
#include "mderrors.h"

#include<iostream>
#include<sstream>
#include<string>

#define DEBUGVAR 0

using namespace std;

void detectionCuts(struct filekeywords* Paramfile, struct event *Event, struct obsfilekeywords World[], struct slcat *Sources, struct slcat *Lenses)
{

  int obsgroup, grpidx, obsidx;
  double pllx=Paramfile->pllxMultiplyer;

  //dont bother if there was an error generating the lightcurve
  if(Event->lcerror) return;

  double chi2point, chi2model;

  stringstream ss;
  string fmfname;
  string header;
  string values;

  Event->deterror=0;
  Event->detected=0;
  
  //Determine if an event is detected

  for(obsgroup=0; obsgroup<int(Event->obsgroups.size()); obsgroup++)
    {
      double t0_diff =1E50;
      double t0_diff_min = 1E50;
      int t0_diff_min_idx=0;
      Event->currentgroup=obsgroup;
	  int nobs = int(Event->obsgroups[obsgroup].size());
      //Event->obsgroupoutput[obsgroup] = string("999999 999999 999999 999999 999999 999999 999999 999999 999999 999999 999999 999999 999999 999999 999999 ");
      //CHECK

      //find the number of column values we need to fill if the event goes undetected
      //this will be {t0 te u0 rs} or if pllx {t0 te u0 rs pien piee} 
      int nfixpar = 4 + (pllx?2:0);
      //also need to add on the twice the number of observatories in the group for {F0 fs} in each filter
      int nparams = nfixpar + 2*int(Event->obsgroups[obsgroup].size());

      //for each of the nparams, add a column to the obs group data with value 999
      //this string will be replaced the same number of columns but real values if the event is detected
      stringstream tempss;
      for(int i=0;i<nparams;i++)
	{
	  tempss << "999 ";
	}
      //give tempss value to the obsgroupdata
      Event->obsgroupoutput[obsgroup] = tempss.str();

      //Write out the header columns too
      vector<string> fmheaders;
      fmheaders.push_back("t0"); fmheaders.push_back("logtE"); fmheaders.push_back("u0");
      fmheaders.push_back("logrho");
      
      if(pllx)
	{
	  fmheaders.push_back("piEN");
	  fmheaders.push_back("piEE");
	}
      for(int obs=0;obs<nobs;obs++)
	{
	  stringstream tmphdr;
	  tmphdr << "F0" << obs;
	  fmheaders.push_back(tmphdr.str());
	  tmphdr.str(string());
	  tmphdr << "logfs" << obs;
	  fmheaders.push_back(tmphdr.str());
	  tmphdr.str(string());
	}

      stringstream mdhss;
      
      for(int j=0;j<nparams;j++)
	{
	  for(int i=0;i<=j;i++)
	    {
	      mdhss << "ObsGroup_" << obsgroup << "_cov_" << fmheaders[i] << "_" << fmheaders[j] << " ";
	    }
	}
	  
      // mdhss << "ObsGroup_" <<  obsgroup << "_sigma_ ";
      // mdhss << "ObsGroup_" <<  obsgroup << "_sigma_logtE ";
      // mdhss << "ObsGroup_" <<  obsgroup << "_sigma_u0 ";
      // mdhss << "ObsGroup_" <<  obsgroup << "_sigma_logrho ";
      // if(pllx)
      // 	{
      // 	  mdhss << "ObsGroup_" <<  obsgroup << "_sigma_piEN ";
      // 	  mdhss << "ObsGroup_" <<  obsgroup << "_sigma_piEE ";
      // 	}
      // for(int obs=0;obs<nobs;obs++)
      // 	{
      // 	  mdhss << "ObsGroup_" <<  obsgroup << "_sigma_F0_" << obs << " ";
      // 	  mdhss << "ObsGroup_" <<  obsgroup << "_sigma_logfs_" << obs << " ";
      // 	}

      Event->obsgroupoutputheader[obsgroup] = mdhss.str();
      
      //this is the old way, hopefully more adaptable now
      //Event->obsgroupoutput[obsgroup] = string("999 999 999 999 999 999 999 999 999 999 999 999 999 999 999 ");

      Event->flag_needFS[obsgroup]=0;
      Event->PSPL[obsgroup].chisq=0;

      int N3sig=0;
      int nconsec=0;
      double lastt=1e30;
      
      //for each observatory group
      for(grpidx=0;grpidx<nobs;grpidx++)
	{
	  //The observatory index is now the group index?
	  //the nepochsvec is sequential for all groups, so obsidx is the start using the 
	  obsidx=Event->obsgroups[obsgroup][grpidx];
	  //for 
	  for(int idx=Event->nepochsvec[obsidx];idx<Event->nepochsvec[obsidx+1];idx++)
	    {
	      //Are the points detected at > 3 sigma? (must be consecutive here)
	      chi2point = sqr((Event->Aobs[idx]-1.0)/Event->Aerr[idx]);
	      chi2model = sqr((Event->Aobs[idx]-Event->Atrue[idx])/Event->Aerr[idx]);
              t0_diff = (Event->epoch[idx]-Event->t0);
			  
	      if(abs(t0_diff)<abs(t0_diff_min) && grpidx==0)
		{
		  t0_diff_min=t0_diff;
		  t0_diff_min_idx = idx;
		}
			  
			  
	      if(Event->nosat[idx])
		{
		  //store the total chi^2 in the right place
		  Event->PSPL[obsgroup].chisq += chi2point-chi2model;
				  
		  if(chi2point>9)
		    {
		      if(Event->epoch[idx]-lastt>0 && Event->epoch[idx]-lastt<1)
			{
			  nconsec++;
			  if(nconsec>N3sig) N3sig=nconsec;
			}
		      else
			{
			  nconsec=1;
			}
		      //store chi_{3+}^2 in the observatory chi^2 vector
		      Event->PSPL[obsgroup].chisqvec[obsidx] += chi2point-chi2model;
		    }
		  else
		    {
		      nconsec=0;
		    }
		}
			  
	      lastt=Event->epoch[idx];
	    }
	}
	  
      //Doing scaling calculations if in the paramfile
      if(Paramfile->error_scaling)
	{
	  //Find uniform scaling factor that would have \delta\chi^2 > 300
	  Event->scale_factor_300 = sqrt(Event->PSPL[obsgroup].chisq/300);
	  //Find the nth point away from t0 such that we can find the scaling required to reduce n3sigma below threshold
	  int t0_idx;
	  //First we will do n3sig>=3
	  /*In this case
	    t0
	    |          |  |       |
	    t1         t2         t3
	    -->tmin-t0<0
	    -->thrid furthest point is -1 away from tmin_idx
	  */
	  if (t0_diff_min<0)
	    {
	      t0_idx = t0_diff_min_idx -1;
	    }
	  /*In this case
	    t0
	    |       |  |          |
	    t1         t2         t3
	    -->tmin-t0>0
	    -->thrid furthest point is +1 away from tmin_idx
	  */
	  else
	    {
	      t0_idx = t0_diff_min_idx +1;
	    }
	 
		  
          //Then for a SNR of 3
	  Event->scale_factor_n3sig3 = (Event->Atrue[t0_idx]-1)/(3*Event->Atrueerr[t0_idx]);
		  
	  //Next  we will do n3sig>=6
	  /*In this case
	    t0
	    |          |          |  |       |         |          |
	    t1         t2         t3         t4        t5         t6
	    -->tmin-t0<0
	    -->sixth furthest point is +3 away from tmin_idx
	  */
	  if (t0_diff_min<0)
	    {
	      t0_idx = t0_diff_min_idx +3;
	    }
	  /*In this case
	    t0
	    |          |          |       |  |         |          |
	    t1         t2         t3         t4        t5         t6
	    -->tmin-t0>0
	    -->thrid furthest point is +1 away from tmin_idx
	  */
	  else
	    {
	      t0_idx = t0_diff_min_idx -3;
	    }
	  //Then for a SNR of 3
	  Event->scale_factor_n3sig6 = (Event->Atrue[t0_idx]-1)/(3*Event->Atrueerr[t0_idx]);
	}
      
      
      //Store the number of 3 sigma detections in the flatlc flag
      Event->flatlc[obsgroup] = N3sig;
      
      if(Event->PSPL[obsgroup].chisq>Paramfile->minChiSquared) // && Event->flatlc[obsgroup] >= 3)
	{
	  Event->detected=1;
	}
      
      Event->flatchi2[obsgroup] = Event->flatlc[obsgroup];
      if(Event->flatlc[obsgroup]>0) Event->flatlc[obsgroup]=0;
      else Event->flatlc[obsgroup]=1;
      
    }

  if(Event->detected)
    {
      if(DEBUGVAR) cout << "Entering fisherMatrix" << endl;
      fisherMatrix(Paramfile, Event, World, Sources, Lenses);
      if(DEBUGVAR) cout << "Leaving fisherMatrix" << endl;

      for(obsgroup=0; obsgroup<int(Event->obsgroups.size()); obsgroup++)
	{
	  Event->currentgroup=obsgroup;
	  int nobs = int(Event->obsgroups[obsgroup].size());
		  
	  //Generate the header information
	  ss.str(string("")); 
	  ss << Paramfile->outputdir << Paramfile->run_name << "_" 
	     << Paramfile->instance << "_";
	  if(Paramfile->choosefield>=0) ss << Paramfile->choosefield << "_";
	  ss << Event->id << ".det.fm." << obsgroup;
	  fmfname = ss.str();
		  
	  if(!Event->outputthis) fmfname=string("");
		  
	  ss.str(string("")); 
	  ss << "t0\tlogtE\tu0\tlogrs\t";
	  if(pllx) ss << "piEN\tpiEE\t";
	  ss << "{F0\tlogfs}";
	  header = ss.str();
		  
	  ss.str(string(""));
	  //ss.precision(16);
	  ss << Event->t0 << "\t" << log10(Event->tE_r)  << "\t" << Event->u0 << "\t" 
	     << log10(Event->rs) << "\t";
	  if(pllx) ss << Event->piEN << "\t" << Event->piEE << "\t";	  
	  for(grpidx=0;grpidx<nobs;grpidx++)
	    {
	      obsidx = Event->obsgroups[obsgroup][grpidx];
	      ss << 1.0 << "\t" << log10(Event->fs[obsidx]) << "\t";
	    }
	  values = ss.str();
		  
	  //Compile the differential vector
	  vector<double> dF;
	  vector<double> err;
		  
	  int nfixpar = 4 + (pllx?2:0);
	  int nparams = nfixpar + 2*nobs;
		  
		  
	  int nepochs=0;
	  for(grpidx=0;grpidx<nobs;grpidx++)
	    {
	      obsidx=Event->obsgroups[obsgroup][grpidx];
	      int startidx = Event->nepochsvec[obsidx];
	      int endidx = Event->nepochsvec[obsidx+1];
	      nepochs += (endidx - startidx);
	    }
		  
	  if(Paramfile->verbosity>1) cout << "\nFisher matrix, group " << obsgroup << endl;
		  
	  for(int param=0;param<nparams;param++)
	    {
	      int parobsidx = Event->obsgroups[obsgroup][(param-nfixpar)/2];
	      int ofallparams = (param<nfixpar ? 
				 param : 
				 nfixpar + 2*parobsidx + (param-nfixpar)%2);
	      int allshift = Event->nepochs*ofallparams;
			  
	      for(grpidx=0;grpidx<nobs;grpidx++)
		{
		  obsidx=Event->obsgroups[obsgroup][grpidx];
				  
		  int startidx = Event->nepochsvec[obsidx];
		  int endidx = Event->nepochsvec[obsidx+1];
				  
		  if(Paramfile->verbosity>2) cout << "nparams, nepochs, param, nobs, ofallparams, allshift, startidx, endidx, dFsize, Aerrsize = " << nparams << " " << nepochs << " " << param << " " << nobs << " " << ofallparams << " " << allshift << " " << startidx << " " << endidx << " " << Event->dF.size() << " " << Event->Aerr.size() << endl;
				  
		  for (int idx=startidx;idx<endidx;idx++)
		    {
		      //if(idx+allshift>=Event->dF.size()) cout << "Here's the segfault" << endl;
		      dF.push_back(Event->dF[idx + allshift]);
		      if(param==0) err.push_back(Event->Aerr[idx]);
		    } //for idx
		  if(Paramfile->verbosity>2) cout << "dF+Aerr copied for parameter " << param << endl;
				  
		} //for grpidx
	    } //for param
		  
	  if(Paramfile->verbosity>1)
	    {
	      cout << "Computing fisher inversion for group " << obsgroup << " with " << nparams << " params and " << nepochs << " datapoints" << endl;
	    }

	  gsl_matrix* cov = gsl_matrix_alloc(nparams,nparams);
	  fisherInversion(dF, err, nparams, nepochs, fmfname, header, values, 
			  &cov);

	  if(Paramfile->verbosity>1)
	    {
	      cout << "Transforming parameters to mass/pirel"  << endl << endl;
	    }
		  
	  vector<double> bestmdcov(7,999999);
	  vector<int> bestfilters(2,999999);
	  vector<double> piEdir(3,999999);

	  if(Event->obsgroups[obsgroup].size()>=2)
	    {
	      piEdir = piEdirerrors(cov, nparams, Event, Sources, Lenses);
		  
	      for(int obs1=0;obs1<int(Event->obsgroups[obsgroup].size());obs1++)
		{
		  for(int obs2=0;obs2<int(Event->obsgroups[obsgroup].size());obs2++)
		    {
		      int filter1 = World[Event->obsgroups[obsgroup][obs1]].filter;
		      int filter2 = World[Event->obsgroups[obsgroup][obs2]].filter;
		      int worstalpha=0;
		      double k2alpha=10; //only use approved filters!
		      int goodfilter=0;
					  
		      //alphas from fitting ax+b to Dartmouth fehp00afep2 
		      //isochrone for logg<4.5. Assumes that redder filter
		      //sets the zero magnitude angular diameter
		      if(filter1==1 && filter2==4)
			{
			  //g-z
			  k2alpha = 0.217;
			  goodfilter=0;
			}
		      if(filter1==1 && filter2==3)
			{
			  //g-i
			  k2alpha = 0.294;
			  goodfilter=0;
			}
		      if(filter1==1 && filter2==2)
			{
			  //g-r
			  k2alpha = 0.511;
			  goodfilter=0;
			}
		      if(filter1==2 && filter2==3)
			{
			  //r-i
			  k2alpha = 0.962;
			  goodfilter=1;
			}
		      if(filter1==2 && filter2==4)
			{
			  //r-z
			  k2alpha = 0.523;
			  goodfilter=1;
			}
		      if(filter1==3 && filter2==4)
			{
			  //i-z
			  k2alpha = 1.375;
			  goodfilter=1;
			}
		      if(filter1==6 && filter2==8) 
			{
			  //V-I
			  k2alpha = 0.427; 
			  goodfilter=1;
			}
		      

		      //vector<double> mdcov = mderrors(fmfname, Event->obsgroup[obsgroup][obs1], Event->obsgroups[obsgroup][obs2], k2alpha, 1, Event, Sources, Lenses);
		      vector<double> mdcov;
		      if(goodfilter)
			{
			  //what is mdcov?
			  mdcov = mderrors(cov, nparams, Event->obsgroups[obsgroup][obs1], Event->obsgroups[obsgroup][obs2], k2alpha, 1, Event, Sources, Lenses);

			  if(mdcov[0]<bestmdcov[0] && mdcov[0]>0 && mdcov[2]>0)
			    {
			      bestmdcov = mdcov;
			      bestfilters[0] = Event->obsgroups[obsgroup][obs1];
			      bestfilters[1] = Event->obsgroups[obsgroup][obs2];
			    } //if better filter pair

			} //if goodfilter
		    } //for obs2
		} //for obs1

	      //cout << "mderr group " << obsgroup << " " << bestmdcov[0] << " " << bestmdcov[1] << " " << bestmdcov[2] << " " << bestfilters[0] << " " << bestfilters[1] << endl;
	    } //if 2 obs


	  //for(int i=0;i<Paramfile->numobservatories;i++)
	  //  ofile << Event->fs[i] << " ";
	  //ofile << "| ";


	  stringstream mdss;
	  //md - mass distance
	  /*
	    Where I am at: Figuring out how to be compliant if pllx is not requested. The mdss string is going to contain the errors for all nparams. I want to record both piEE and piEN, but 
	    this is probably not the appropriate place to do that. 
	    
	    I need to move this somewhere else.
	  */

	  // 76-85
	  // ALL SIGMAS ARE OUTPUT XXX  sigt0 sigte sigu0 sigrs sigpien sigpiee {sigF0 sigfs}
	  for(int j=0;j<nparams;j++)
	    {
	      for(int i=0;i<=j;i++)
		{
		  mdss << gsl_matrix_get(cov,i,j) << " ";
		  if(Paramfile->verbosity>0) cout << "cov(" << i << "," << j << ")=" << sqrt(gsl_matrix_get(cov,i,j)) << " ";
		  
		}
	    }
	  if(Paramfile->verbosity>0) cout << endl;


	  //Old way, should still work
	  //mdss << bestmdcov[0] << " " << bestmdcov[1] << " " << bestmdcov[2]  << " " << bestfilters[0] << " " << bestfilters[1] << " " << bestmdcov[3] << " " << bestmdcov[4] << " " << bestmdcov[5] << " " << bestmdcov[6] << " " << sqrt(gsl_matrix_get(cov,3,3))/Event->rs/ln10 << " " << sqrt(gsl_matrix_get(cov,4,4))/abs(Event->piEN)/ln10 << " " << sqrt(gsl_matrix_get(cov,5,5))/abs(Event->piEE)/ln10 << " " << piEdir[0] << " " << piEdir[1] << " " << piEdir[2] << " ";
	  
	  Event->obsgroupoutput[obsgroup] = mdss.str();

	  //cout << "before free" << endl;
	  gsl_matrix_free(cov);
	  //cout << "after free" << endl;

	} //for obsgroup

	  

    } //if detected
  
}
