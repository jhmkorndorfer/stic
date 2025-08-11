/* ------- file: -------------------------- initscatter.c -----------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Apr  1 09:53:47 2009 --

       --------------------------                      ----------RH-- */

// #include "rh.h"
#include "rhf1d.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "accelerate.h"
#include "error.h"
#include "statistics.h"
#include "inputs.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];


/* ------- begin -------------------------- initScatter.c ----------- */

void initScatter()
{
  const char routineName[] = "initScatter";
  register int n, kr, niter;

  double dJmax;
  Atom *atom;
  AtomicLine *line;

  /* --- Fill the radiative rates from previous solution if 
         PRD lines are present --                      -------------- */

  if (atmos.NPRDactive > 0  && input.startJ == OLD_J) {
    for (n = 0;  n < atmos.Natom;  n++) {
      atom = &atmos.atoms[n];
      if (atom->active) {
	readRadRate(atom);
	for (kr = 0;  kr < atom->Nline;  kr++) {
	  line = atom->line + kr;
	  if (line->PRD) {
	    switch (input.PRD_angle_dep) {
	    case PRD_ANGLE_INDEP:
	      PRDScatter(line, LINEAR);
	      break;
	    case PRD_ANGLE_APPROX:
	      PRDAngleApproxScatter(line, LINEAR);
	      break;
	    case PRD_ANGLE_DEP:
	      PRDAngleScatter(line, LINEAR);
	      break;
	    }
	  }
	}
      }
    }
  }
  /* --- Iterate the scattering in the background. Only iterate
         when we are also planning to do main iterations -- --------- */

  if (spectrum.updateJ) {
    sprintf(messageStr,
	    " %s: Lambda iterating background scattering...\n\n",
	    routineName);
    Error(MESSAGE, routineName, messageStr);

    niter = 0;
    while (input.NmaxIter  &&  niter < input.NmaxScatter) {
      dJmax = solveSpectrum(FALSE, FALSE, 0, TRUE);
      if (dJmax < input.iterLimit) break;
      niter++;
    }
  }
  getCPU(2, TIME_POLL, "Initial solution");
}
/* ------- end ---------------------------- initScatter.c ----------- */


/* ------- begin -------------------------- initScatter.c ----------- */

void initScatter_ctx(RHContext *ctx)
{
  const char routineName[] = "initScatter_ctx";
  register int n, kr, niter;

  Atmosphere *atmosLocal = &ctx->atmos;
  InputData *inputLocal = &ctx->input;
  Spectrum *spectrumLocal = &ctx->spectrum;

  double dJmax;
  Atom *atom;
  AtomicLine *line;

  /* --- Fill the radiative rates from previous solution if 
         PRD lines are present --                      -------------- */

  if (atmosLocal->NPRDactive > 0  && inputLocal->startJ == OLD_J) {
    for (n = 0;  n < atmosLocal->Natom;  n++) {
      atom = &atmosLocal->atoms[n];
      if (atom->active) {
	readRadRate_ctx(atom, ctx);
	for (kr = 0;  kr < atom->Nline;  kr++) {
	  line = atom->line + kr;
	  if (line->PRD) {
	    switch (inputLocal->PRD_angle_dep) {
	    case PRD_ANGLE_INDEP:
	      PRDScatter_ctx(line, LINEAR, ctx);
	      break;
	    case PRD_ANGLE_APPROX:
	      PRDAngleApproxScatter_ctx(line, LINEAR, ctx);
	      break;
	    case PRD_ANGLE_DEP:
	      PRDAngleScatter_ctx(line, LINEAR, ctx);
	      break;
	    }
	  }
	}
      }
    }
  }
  /* --- Iterate the scattering in the background. Only iterate
         when we are also planning to do main iterations -- --------- */

  if (spectrumLocal->updateJ) {
    sprintf(messageStr,
	    " %s: Lambda iterating background scattering...\n\n",
	    routineName);
    Error(MESSAGE, routineName, messageStr);

    niter = 0;
    while (inputLocal->NmaxIter  &&  niter < inputLocal->NmaxScatter) {
      dJmax = solveSpectrum_ctx(FALSE, FALSE, 0, TRUE, ctx);
      if (dJmax < inputLocal->iterLimit) break;
      niter++;
    }
  }
  getCPU(2, TIME_POLL, "Initial solution");
}
/* ------- end ---------------------------- initScatter_ctx.c ----------- */
