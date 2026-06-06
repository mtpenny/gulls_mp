#include "lightcurveFitter.h"
#include "detectionCuts.h"

void detectionCuts(struct filekeywords* Paramfile, struct event *Event, struct obsfilekeywords World[], struct slcat *Sources, struct slcat *Lenses)
{

  int obsidx;
  
  //dont bother if there was an error generating the lightcurve
  if(Event->lcerror) return;

  Event->deterror=0;
  Event->detected=0;

  for(int obsgroup=0; obsgroup<int(Event->obsgroups.size()); obsgroup++)
    {

      Event->currentgroup=obsgroup;

      //Determine if an event is detected 

      //Fit a PS lightcurve
      lightcurveFitter(Paramfile, World, Event);

      //reset parallax parameters back to nominal values if necessary
      for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
	{
	  Event->pllx[obsidx].provide_murel_h_lb(Event->murel_l,Event->murel_b,
						 Event->piE, Event->thE);
	  Event->pllx[obsidx].compute_tushifts();
	}


      /* If finite source star effects important, fit FS lightcurve */
      if(abs(Event->umin) < 10*Event->rs && Event->PSPL[obsgroup].chisq>Paramfile->minChiSquared)
	{
	  Event->flag_needFS[obsgroup]=1;
	  lightcurveFitter_FS(Paramfile, World, Event);

	  //reset parallax parameters back to nominal values if necessary
	  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
	    {
	      Event->pllx[obsidx].provide_murel_h_lb(Event->murel_l,Event->murel_b,
						     Event->piE, Event->thE);
	      Event->pllx[obsidx].compute_tushifts();
	    }
	}

      //did we detect it?
      double chi2;

      if(Paramfile->outputOnDet==2) chi2 = Event->flatchi2[obsgroup];
      else
	{
	  if(Event->flag_needFS[obsgroup])
	    chi2 = Event->FSPL[obsgroup].chisq;
	  else
	    chi2 = Event->PSPL[obsgroup].chisq;
	}
	  
      if(chi2 > Paramfile->minChiSquared) Event->detected=1;
    
    } //for each obsgroup
}
