#include "lightcurveGenerator.h"
#include "lightcurveFitter.h"
#include "astroFns.h"
#include "integerPowers.h"
#include "wittFSPL.h"
#include "singleLens.h"
#include "fisherinversion.h"

#define DEBUGVAR 1

void fisherMatrix(struct filekeywords* Paramfile, struct event *Event, struct obsfilekeywords World[], struct slcat *Sources, struct slcat *Lenses)
{
  double rs, Gamma;
  double t0, tE, u0;
  double tt,uu;
  double logtE=0;
  double pllx=Paramfile->pllxMultiplyer;
  //vector<parallax> parlx;

  int idx,obsidx;
  
  double piEN=0, piEE=0, piE;
  
  if(DEBUGVAR) cout << "Calculating fisher matrix for event number " << Event->id << endl;

  int Nparams=6;
  if(pllx==0) Nparams=4;
  int Nobsparams = 2*Paramfile->numobservatories;
  int Ntotparams = Nparams + Nobsparams;
  // 0 = t0  1 = tE  2 = u0  3 = rs  4 = piEN 5 = piEE 6 = F0  7 = fs
  Event->dF.resize(Ntotparams*Event->nepochs);
  double* fs = new double[Paramfile->numobservatories];
  double u;

  //vector<int> idxshift;
  //for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
  //  idxshift.push_back(Event->nepochsvec[obsidx]);

  //for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
  //  {
  //    parlx.push_back(Event->pllx[obsidx]);
  //  }

  //first initialize the derivatives to zero so that we can accumulate the sum
  for(int param=0; param<Ntotparams; param++)
    {
      int pshift = param*Event->nepochs;
      for(idx=0; idx<Event->nepochs; idx++) Event->dF[pshift+idx]=0;
    }
	 
  t0 = Event->t0;
  tE = Event->tE_r;
  logtE = log10(tE);
  u0 = Event->u0;
  rs = Event->rs; 
  Gamma = Event->gamma;
  if(pllx)
    {
      piEN = Event->piEN;
      piEE = Event->piEE;
      piE = qAdd(piEN,piEE);
    }

  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    {
      fs[obsidx] = Event->fs[obsidx];
    }

  Event->deterror=0;

  double dmudrho;
  double dmudu;
  double tshift=0, ushift=0;
  double deltaN, deltaE;
  double tau;
  int shiftedidx;

  //Calculate the differential of the lightcurve
  for(idx=0;idx<Event->nepochs;idx++)
    {
      obsidx = Event->obsidx[idx];
      shiftedidx = idx-Event->nepochsvec[obsidx];
      
      tau = (Event->epoch[idx] - t0) / tE;
      tt = tau;
      uu = u0;
      
      if(pllx)
	{
	  //tshift = parlx[obsidx].tshift(Event->jdepoch[idx]);
	  //ushift = parlx[obsidx].ushift(Event->jdepoch[idx]);
	  tshift = Event->pllx[obsidx].tshift[shiftedidx];
	  ushift = Event->pllx[obsidx].ushift[shiftedidx];

	  tt += tshift;
	  uu += ushift;
	}
      
      u = qAdd(tt,uu);

      dmududr_witt(u, rs, &dmudu, &dmudrho);
      //dmudu = dAdu(u);
      //dmudrho=0;
      
      // 0 = t0  1 = tE  2 = u0  3 = rs  4 = piE 5 = phi_pi 6 = F0  7 = fs
      //now  0 = t0  1 = logtE  2 = u0  3 = logrs  4 = piEN 5 = piEE 6 = F0  7 = fs
      Event->dF[idx] = fs[obsidx] * dmudu * tt/u * (-1.0/tE);
      Event->dF[idx+Event->nepochs] = fs[obsidx] * dmudu * tt/u * (-tau/tE) * tE*log(10);
      Event->dF[idx+2*Event->nepochs] = fs[obsidx] * dmudu * uu/u;
      Event->dF[idx+3*Event->nepochs] = fs[obsidx] * dmudrho * rs*log(10);
      if(pllx)
	{
	  //piE,phi_piE
	  //dF[idx+4*Event->nepochs] = fs[obsidx] 
	  //  * dmudu/u/piE * (uu*ushift + tt*tshift);
	  //dF[idx+5*Event->nepochs] = fs[obsidx] 
	  //  * dmudu/u*(uu*tshift + tt*ushift);
	  //piEN,piEE
	  //deltaN = Event->pllx[obsidx].NEshift[shiftedidx][0];
	  //deltaE = Event->pllx[obsidx].NEshift[shiftedidx][1];
	  deltaN = Event->pllx[obsidx].Nshift[shiftedidx];
	  deltaE = Event->pllx[obsidx].Eshift[shiftedidx];
	  Event->dF[idx+4*Event->nepochs] = fs[obsidx] * dmudu/u 
	    * (-tt*deltaN - uu*deltaE);
	  Event->dF[idx+5*Event->nepochs] = fs[obsidx] * dmudu/u 
	    * (-tt*deltaE + uu*deltaN);
	}  
    } //end for idx

  //For the linear flux parameters we can work analytically
  int F0idx, fsidx;

  //Now calculate dF/dF0 and dF/dlogfs
  for(idx=0;idx<Event->nepochs;idx++)
    {
      obsidx = Event->obsidx[idx];
      F0idx = (Nparams+obsidx*2)*Event->nepochs+idx;
      fsidx = F0idx + Event->nepochs;
      Event->dF[F0idx]=Event->Atrue[idx];
      Event->dF[fsidx]=(Event->Atrue[idx]-1)/fs[obsidx] * fs[obsidx]*log(10); 
    }     
  
  delete[] fs;

}

