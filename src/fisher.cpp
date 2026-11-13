#include <sstream>
#include <sys/stat.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <vector>
#include <gsl/gsl_deriv.h>

#include "lightcurveGenerator.h"
#include "backupGenerator.h"
#include "astroFns.h"
#include "VBMicrolensingLibrary.h"
#include "fs.h"
#include "src_cld.h"
#include "lens_binary.h"
#include "cd.h"

#define DEBUGVAR 0

using namespace std;

// Context struct for GSL derivative wrapper
struct DerivContext {
    int param;
    int idx;
    filekeywords* Paramfile;
    event* Event;
    obsfilekeywords* World;
};

// GSL derivative wrapper using context struct
static double deriv_wrapper(double x, void* p)
{
  DerivContext* ctx = static_cast<DerivContext*>(p);
  
  filekeywords* Paramfile = ctx->Paramfile;
  event* Event = ctx->Event;
  obsfilekeywords* World = ctx->World;
  int param = ctx->param;
  int idx = ctx->idx;

  double logq = log10(Event->p_q[0]);
  double logs = log10(Event->p_s0[0]);
  double logtE = log10(Event->tE_r);
  double logrs = log10(Event->rs);
  double piEN = Event->piEN;
  double piEE = Event->piEE;

  double t0 = Event->t0, u0 = Event->u0;

  double rot = atan2(Event->p_y0[0],Event->p_x0[0]);
  
    
  switch (param)
    {
    case 0: t0 += x; break;
    case 1: logtE += x; break;
    case 2: u0 += x; break;
    case 3: rot += x; break;
    case 4: logs += x; break;
    case 5: logq += x; break;
    case 6: logrs += x; break;
    case 7: piEN += x; break;
    case 8: piEE += x; break;
    }

  double tE = pow(10, logtE);
  double q = pow(10, logq);
  double s = pow(10, logs);
  double rs = pow(10, logrs); 
  double m1 = 1.0 / (1.0 + q);
  //double origin = (1.0 - m1) * (-s);
  double a1 = -q*s/(1+q);

  double cr = cos(-rot);
  double sr = sin(-rot);

  int obsidx = Event->obsidx[idx];
  int shiftedidx = idx - Event->nepochsvec[obsidx];

  // Without parallax, use the event's configured trajectory angle directly.
  double alpha = Event->alpha * TO_RAD;
  if (Paramfile->pllxMultiplyer)
    {
      // Each model call may follow a different finite-difference perturbation.
      // Rebuild this observer's direction and this epoch's shift from the
      // current parameters so an earlier call cannot leave parallax state behind.
      Event->pllx[obsidx].provide_observables_NE(piEN, piEE, tE, Event->thE);
      double phi_pi = atan2(Event->pllx[obsidx].mulam_r, Event->pllx[obsidx].mubet_r);
      alpha = phi_pi - pi;
      Event->pllx[obsidx].compute_tushifts(shiftedidx);
    }
  double cosa = cos(alpha), sina = sin(alpha);

  double tt = (Event->epoch[idx] - t0) / tE;
  double uu = u0;
  if (Paramfile->pllxMultiplyer) {
    tt += Event->pllx[obsidx].tshift[shiftedidx];
    uu += Event->pllx[obsidx].ushift[shiftedidx];
  }

  //double xs = tt * cosa - uu * sina + origin;
  //double ys = tt * sina + uu * cosa;

  // The source track starts at the primary, while BinaryMag2 uses the binary
  // center of mass as its origin. Translate the primary to that origin once.
  double xl1 = a1 * cos(rot);
  double yl1 = a1 * sin(rot);

  double xs0 = uu * sina + tt * cosa + xl1;
  double ys0 = -uu * cosa + tt * sina + yl1;


  double xsi = cr*xs0 - sr*ys0;
  double ysi = sr*xs0 + cr*ys0;
  
  double limb_gamma=Paramfile->LD_GAMMA;
  Event->vbm->a1=limb_gamma;

  double amp=Event->vbm->BinaryMag2(s, q, xsi, ysi, rs);

  // Flux-only evaluations use param=-1 and do not need per-call debug output.
  if (param >= 0)
    cout << Event->id << " " << param << " " << x << " " << xl1 << " " << yl1 << " " << xs0 << " " << ys0 << " " << xsi << " " << ysi << " " << rot << " " << alpha << " " << a1 << " " << t0 << " " << tE << " " << u0 << " " << s << " " << q << " " << rs << " " << piEN << " " << piEE << " " << amp << endl;

  return amp;
}

void fisherMatrix(filekeywords* Paramfile, event* Event,
                  obsfilekeywords World[], slcat* Sources, slcat* Lenses)
{
  if (Paramfile->verbosity>0) 
    cout << "Calculating Fisher matrix for event " << Event->id << endl;

  const int Nparams = (Paramfile->pllxMultiplyer == 0 ? 7 : 9);
  const int Nobsparams = 2 * Paramfile->numobservatories;
  const int Ntotparams = Nparams + Nobsparams;

  if (Paramfile->verbosity>0)
    cout << "Nparams=" << Nparams << " Ntotparams=" << Ntotparams << endl;

  //case 0: t0 += x; break;
  //case 1: logtE += x; tE = pow(10, logtE); break;
  //case 2: u0 += x; break;
  //case 3: rot += x; break;
  //case 4: logs += x; break;
  //case 5: logq += x; break;
  //case 6: rs = pow(10, log10(Event->rs) + x); break;
  //case 7: piEN += x; break;
  //case 8: piEE += x; break;
  vector<double> step(Ntotparams);
  for (int i = 0; i < Nparams; ++i)
    {
      //double base = (i < Event->params.size()) ? std::max(fabs(Event->params[i]), 1.0) : 1.0; //MP: I have no idea what this was doing
      
      if (i == 6) step[i] = 1.0e-3;
      else if (i == 7 || i == 8) step[i] = 1.0e-4*Event->piE;
      else if (i == 3) { step[i] = 1.0e-4; } // alpha is stored in degrees; multiply by TO_RAD so δalpha is in radians 
      else step[i] = 1e-4;
    }

  Event->dF.clear();
  Event->dF.resize(Ntotparams * Event->nepochs, 0.0);

  for (int param = 0; param < Nparams; ++param)
    {
      DerivContext ctx{param, 0, Paramfile, Event, World};
      gsl_function F;
      F.function = &deriv_wrapper;
      F.params   = &ctx;

      if (Paramfile->verbosity>0)
	cout << "Calculating derivatives for param " << param << endl;

      vector<bool> suppressDerivative(Event->nepochs, false);
      for (int idx = 0; idx < Event->nepochs; ++idx)
	{
	  ctx.idx = idx;
	  double result = 0.0, abserr = 0.0;
	  gsl_deriv_central(&F, 0.0, step[param], &result, &abserr);
	  int obsidx = Event->obsidx[idx];
	  double fs = Event->fs[obsidx];
	  double snr = fabs(result) / (abserr + 1e-12);
	  bool suppress = std::isnan(result) || std::isnan(abserr) ||
	    ((param == 6) ? (snr < 1.0) : (snr < 3.0));
	  int index = param * Event->nepochs + idx;
	  // Preserve the existing signal-to-noise and anomaly-window decision; defer
	  // replacement until every epoch has its raw derivative.
	  suppressDerivative[idx] = suppress &&
	    fabs(Event->epoch[idx] - Event->t0) >= 0.2 * Event->tE_r;
	  Event->dF[index] = result * fs;
	} //end for idx

      // Fill every epoch before interpolation so future neighbours are real values.
      // Keep neighbours on the same observer: adjacent array entries can cross
      // observer boundaries, where parallax tracks and sampling differ.
      for (int idx = 0; idx < Event->nepochs; ++idx)
	{
	  if (!suppressDerivative[idx]) continue;
	  double interp = 0.0;
	  int count = 0;
	  for (int offset = 1; offset <= 3; ++offset)
	    for (int j = -1; j <= 1; j += 2)
	      {
		int neighbor = idx + j * offset;
		if (neighbor >= 0 && neighbor < Event->nepochs &&
		    Event->obsidx[neighbor] == Event->obsidx[idx])
		  {
		    double val = Event->dF[param * Event->nepochs + neighbor];
		    if (!std::isnan(val))
		      {
			interp += val;
			++count;
		      }
		  }
	      }
	  if (count > 0)
	    Event->dF[param * Event->nepochs + idx] = interp / count;
	}

      if (Paramfile->pllxMultiplyer && (param == 7 || param == 8 || param == 1))
	{
	  // These parameters change the parallax projection during differentiation.
	  // Restore all observers to the event values before the next parameter pass.
	  for(int obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
	    {
	      Event->pllx[obsidx].provide_murel_h_lb(Event->murel_l,Event->murel_b,
						     Event->piE, Event->thE);
	      Event->pllx[obsidx].compute_tushifts();
	    }
	}
      
    } //end for param

  if (Paramfile->verbosity>0)
    cout << "Computing flux derivatives" << endl;
      

  for (int idx = 0; idx < Event->nepochs; ++idx)
    {
      int obsidx = Event->obsidx[idx];
      int base = (Nparams + obsidx * 2) * Event->nepochs + idx;
      DerivContext ctx{-1, idx, Paramfile, Event, World};
      // Fisher flux terms use its fixed projected lens geometry. Event->Atrue
      // may instead contain the general generator's live orbital motion.
      const double amp = deriv_wrapper(0.0, &ctx);
      Event->dF[base] = 1.0 + Event->fs[obsidx] * (amp - 1.0);
      Event->dF[base + Event->nepochs] = amp - 1.0;
    }

  if (Paramfile->verbosity>0)
    cout << "All fisher derivatives calculated" << endl;

}

