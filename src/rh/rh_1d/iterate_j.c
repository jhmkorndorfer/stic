/* ------- file: -------------------------- iterate.c ---------------

       Version:       rh2.0
       Author:        Han Uitenbroek  (huitenbroek@nso.edu)
       Last modified: Tue Nov 16 15:31:48 2010 --

       --------------------------                      ----------RH-- */

/* --- Main iteration routine --                       -------------- */

 
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "background.h"
#include "accelerate.h"
#include "error.h"
#include "statistics.h"
#include "inputs.h"
#include "rhf1d.h"

typedef struct {
  bool_t eval_operator, redistribute;
  int    nspect, iter;
  double dJ;
} threadinfo;

/* --- Function prototypes --                          -------------- */

void *Formal_pthread(void *argument);


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];
extern MPI_t mpi;


/* ------- begin -------------------------- Iterate.c --------------- */

void Iterate_j(int NmaxIter, double iterLimit, double *dpopmax_out)
{
  const char routineName[] = "Iterate";
  register int niter, nact;

  bool_t    eval_operator, write_analyze_output, equilibria_only, old_ne_flag = FALSE;
  int       Ngorder, nsum = 0;
  double    dpopsmax, PRDiterlimit, cswitch;
  Atom     *atom;
  Molecule *molecule;

  if (NmaxIter <= 0) return;
  getCPU(1, TIME_START, NULL);
  
  /* --- Initialize structures for Ng acceleration of population
         convergence --                                  ------------ */

  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    atom = atmos.activeatoms[nact];
    atom->Ng_n = NgInit(atom->Nlevel*atmos.Nspace, input.Ngdelay,
			input.Ngorder, input.Ngperiod, atom->n[0]);
  }
  for (nact = 0;  nact < atmos.Nactivemol;  nact++) {
    molecule = atmos.activemols[nact];
    //Ngorder  = (input.accelerate_mols) ? input.Ngorder : 0;

    molecule->Ng_nv = NgInit(molecule->Nv*atmos.Nspace, input.Ngdelay,
			     input.Ngorder, input.Ngperiod, molecule->nv[0]);
  }
  
  if(input.solve_ne >= ITERATION_EOS)
    atmos.ng_ne = NgInit(atmos.Nspace, input.Ngdelay, input.Ngorder,input.Ngperiod, atmos.ne );

  
  /* --- Start of the main iteration loop --             ------------ */

  niter = 1;

 /* Collisional-radiative switching ? */
  if (input.crsw != 0.0)
    cswitch = input.crsw_ini;
  else
    cswitch = 1.0;
    
  /* PRD switching ? */
  if (input.prdsw > 0.0)
    input.prdswitch = 0.0;
  else
    input.prdswitch = 1.0;

  
  while ((niter <= NmaxIter || niter < 3)) {
    getCPU(2, TIME_START, NULL);
    mpi.iter = niter;
    
    for (nact = 0;  nact < atmos.Nactiveatom;  nact++)
      initGammaAtom(atmos.activeatoms[nact], cswitch);
    for (nact = 0;  nact < atmos.Nactivemol;  nact++)
      initGammaMolecule(atmos.activemols[nact]);

    /* --- Formal solution for all wavelengths --      -------------- */

    solveSpectrum(eval_operator=TRUE, FALSE, niter, FALSE);

    /* --- Solve statistical equilibrium equations --  -------------- */
    
    if(((niter+1) == input.Ngdelay) && (dpopsmax > input.ng_start_limit)){
      nsum = 1;
      if(input.solve_ne >= ITERATION_EOS && atmos.ne_flag) nsum = 4;
      
      for(nact = 0;  nact < atmos.Nactiveatom;  nact++){
	atmos.activeatoms[nact]->Ng_n->Ndelay += nsum;	
      }
      if(input.solve_ne >= ITERATION_EOS) atmos.ng_ne->Ndelay += nsum;
      input.Ngdelay+=nsum;
    }
    
    sprintf(messageStr, "\n -- Iteration %3d\n", niter);
    Error(MESSAGE, routineName, messageStr);
    dpopsmax = updatePopulations(niter);
    if (mpi.stop) return;

    
    old_ne_flag = atmos.ne_flag;
    
    if (atmos.NPRDactive > 0) {
      
      /* --- Redistribute intensity in PRD lines if necessary -- ---- */

      if (input.PRDiterLimit < 0.0)
	PRDiterlimit = MAX(dpopsmax, -input.PRDiterLimit);
      else
	PRDiterlimit = input.PRDiterLimit;

      Redistribute_j(input.PRD_NmaxIter, PRDiterlimit, dpopsmax*0.1);
      if (mpi.stop) return;
	  
    }

    sprintf(messageStr, "Total Iteration %3d", niter);
    getCPU(2, TIME_POLL, messageStr);

    if ((dpopsmax < iterLimit) && (cswitch <= 1.0) && (input.prdswitch >= 1.0)) {
      
      /* --- Make sure that in the last iteration the PRD is converged all the 
	 wav for the response functions --- */
      
      if (atmos.NPRDactive > 0) {
	
	/* --- Redistribute intensity in PRD lines if necessary -- ---- */
	
	if (input.PRDiterLimit < 0.0)
	  PRDiterlimit = MAX(dpopsmax, -input.PRDiterLimit);
	else
	  PRDiterlimit = input.PRDiterLimit;
       
	Redistribute_j(input.PRD_NmaxIter*2+1, PRDiterlimit, dpopsmax);
	if (mpi.stop) return;
      }
      break;
    }
    niter++;
    

    /* Update collisional radiative switching */
    if (input.crsw > 0)
      cswitch = MAX(1.0, cswitch * pow(0.1, 1./input.crsw));
      
    /* Update PRD switching */ 
    if (input.prdsw > 0.0) 
      input.prdswitch = MIN(1.0, input.prdsw * (double) (niter * niter) ); 

    
    if (atmos.hydrostatic) {
      if (!atmos.atoms[0].active) {
	sprintf(messageStr, "Can only perform hydrostatic equilibrium"
                            " for hydrogen active");
	Error(ERROR_LEVEL_2, routineName, messageStr);
      }
      
      // if(atmos.atoms[0].mxchange > 1.e-2)
      Hydrostatic(N_MAX_HSE_ITER, HSE_ITER_LIMIT);
    }
  }
  
  *(dpopmax_out) = dpopsmax; 
  
  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    atom = atmos.activeatoms[nact];
    freeMatrix((void **) atom->Gamma);
    atom->Gamma = NULL;
    NgFree(atom->Ng_n);
  } 
  for (nact = 0;  nact < atmos.Nactivemol;  nact++) {
    molecule = atmos.activemols[nact];
    freeMatrix((void **) molecule->Gamma);
    molecule->Gamma = NULL;
    NgFree(molecule->Ng_nv);
  }
  if(input.solve_ne >= ITERATION_EOS){
    NgFree(atmos.ng_ne);
    atmos.ng_ne = NULL;
  }
  
  getCPU(1, TIME_POLL, "Iteration Total");
}
/* ------- end ---------------------------- Iterate.c --------------- */

/* ------- begin -------------------------- solveSpectrum.c --------- */

double solveSpectrum(bool_t eval_operator, bool_t redistribute, int iter, bool_t synth_all)
{
  register int nspect, n, nt, k;

  int         Nthreads, lambda_max;
  double      dJ, dJmax;
  pthread_t  *thread_id;
  threadinfo *ti;

  /* --- Administers the formal solution for each wavelength. When
         input.Nthreads > 1 the solutions are performed concurrently
         in Nthreads threads. These are POSIX style threads.

    See: - David R. Butenhof, Programming with POSIX threads,
           Addison & Wesley.

         - Multithreaded Programming Guide, http://sun.docs.com
           (search for POSIX threads).

         When solveSpectrum is called with redistribute == TRUE only
         wavelengths that have an active PRD line are solved. The
         redistribute key is passed to the addtoRates routine via
         Formal so that only the radiative rates of PRD lines are
         updated. These are needed for the emission profile ratio \rho.
         --                                            -------------- */

  getCPU(3, TIME_START, NULL);

  /* --- First zero the radiative rates --             -------------- */

  zeroRates(redistribute);
  lambda_max = 0;
  dJmax = 0.0;

  // zero out J in gas parcel's frame
  if (spectrum.updateJ && input.PRD_angle_dep == PRD_ANGLE_APPROX
      && atmos.Nrays > 1  && atmos.NPRDactive > 0){
      memset(&spectrum.Jgas[0][0],0,spectrum.nJlam*atmos.Nspace*sizeof(double));

  }

  
    /* --- Else call the solution for wavelengths sequentially -- --- */
      
    for (nspect = 0;  nspect < spectrum.Nspect;  nspect++) {
      if (!redistribute ||
	  (redistribute && containsPRDline(&spectrum.as[nspect]))) {
	  dJ = Formal(nspect, eval_operator, redistribute, iter);
	if (dJ > dJmax) {
	  dJmax = dJ;
	  lambda_max = nspect;
	}
      }
    }
  

  sprintf(messageStr, " Spectrum max delta J = %6.4E (lambda#: %d)\n",
	  dJmax, lambda_max);
  Error(MESSAGE, NULL, messageStr);

  getCPU(3, TIME_POLL,
	 (eval_operator) ? "Spectrum & Operator" : "Solve Spectrum");

  return dJmax;
}
/* ------- end ---------------------------- solveSpectrum.c --------- */

/* ------- begin -------------------------- Formal_pthread.c -------- */

void *Formal_pthread(void *argument)
{
  threadinfo *ti = (threadinfo *) argument;
  
  /* --- Threads wrapper around Formal --              -------------- */

  ti->dJ = Formal(ti->nspect, ti->eval_operator, ti->redistribute, (int)ti->iter);

  return (NULL);
}
/* ------- end ---------------------------- Formal_pthread.c -------- */
