

/* ------- file: -------------------------- rhf1d.c -----------------

       Version:       rh2.0, 1-D plane-parallel
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Feb 24 16:40:14 2011 --

       --------------------------                      ----------RH-- */

/* --- Main routine of 1D plane-parallel radiative transfer program.
       MALI scheme formulated according to Rybicki & Hummer

  See: G. B. Rybicki and D. G. Hummer 1991, A&A 245, p. 171-181
       G. B. Rybicki and D. G. Hummer 1992, A&A 263, p. 209-215

       Formal solution is performed with Feautrier difference scheme
       in static ctx->atmospheres, and with piecewise quadratic integration
       in moving ctx->atmospheres.

       --                                              -------------- */
#include "rhf1d.h"
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <errno.h>
#include <unistd.h>
#include <sys/stat.h>


// #include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "background.h"
#include "math.h"
#include "statistics.h"
#include "error.h"
#include "inputs.h"
#include "xdr.h"

#include "dummyatmos.h"
#include "sortlambda_j.h"
#include "initial_j.h"
#include "scatter_j.h"
#include "iterate_j.h"

/* --- Function prototypes --                          -------------- */

//extern void hermitian_interpolation(int n, double *x, double *y, int nn, double *xp, double *yp);


/* --- Global variables --                             -------------- */

enum Topology topology = ONE_D_PLANE;

// typedef struct {
//   Atmosphere atmos;
//   Geometry geometry;
//   Spectrum spectrum;
//   InputData input;
//   CommandLine commandline;
//   ProgramStats stats;
//   rhinfo io;
//   rhbgmem *bmem;
//   crhpop *save_popp;
//   MPI_t mpi;

//   // These were static before
//   bool initialized;
//   int save_Nrays;
//   double save_muz, save_mux, save_muy, save_wmu;
//   enum StokesMode oldMode;

// } RHContext;

Atmosphere atmos;
Geometry geometry;
Spectrum spectrum;
InputData input;
CommandLine commandline;

ProgramStats stats;

char messageStr[MAX_MESSAGE_LENGTH];
rhinfo io;
BackgroundData bgdat;
rhbgmem *bmem; // To store background opac in mem
crhpop *save_popp;
MPI_t mpi;



#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

/* ---- Check for directory --- */

int bdir_exists(const char *name){
  DIR* dir = opendir(name);
  if (dir){
    closedir(dir);
    return 1;
  }else return 0;
}


void adjust_vlos(double *vlos, int ndep, double mu)
{
  int register ii;
  for(ii=0; ii<ndep; ++ii)
    vlos[ii] /= mu;
  
}


/* ------- begin -------------------------- rhf1d.c ----------------- */

bool_t rhf1d(RHContext* ctx, float muz, int rhs_ndep, double *rhs_T, double *rhs_rho, 
	     double *rhs_nne, double *rhs_vturb, double *rhs_v, 
	     double *rhs_B, double *rhs_inc, double *rhs_azi,
	     double *rhs_z, double *rhs_nhtot, double *rhs_tau ,
	     double *rhs_cmass, double gravity, bool_t stokes, ospec *sp,
	     crhpop *save_pop, int mynw, double *mylambda, int myrank, int savpop,
	     int iverbose, int *hydrostat, int computing_derivatives)
{
  
  bool_t write_analyze_output, equilibria_only, quiet = ((iverbose <= 1)? TRUE : FALSE);
  int    niter, nact, i, sNgperiod, sNgdelay, sPRDNITER,k;
  static int save_Nrays;
  static double save_muz, save_mux, save_muy, save_wmu;
  static enum StokesMode oldMode;

  double dpopmax, *ne_lte = NULL;
  static Atom *atom;
  static Molecule *molecule;
  
 
  static bool_t firsttime = TRUE;
 
  int argc = 1;
  char *argv[] = {"rhf1d",NULL};


  adjust_vlos(rhs_v, rhs_ndep, muz);
  
  if(firsttime){
    memset(&ctx->atmos, 0, sizeof(Atmosphere));
    memset(&ctx->geometry, 0, sizeof(Geometry));
    memset(&ctx->spectrum, 0, sizeof(Spectrum));
    ctx->atmos.cos_gamma = NULL;
    ctx->atmos.cos_2chi = NULL;
    ctx->atmos.sin_2chi = NULL;
    
    if(bdir_exists("scratch/") == 0) mkdir("scratch", 0700);

  }
  
  ctx->atmos.moving = TRUE;
  ctx->atmos.Stokes = FALSE;
  ctx->atmos.Nspace = rhs_ndep;
  
  /* --- Read ctx->input data and initialize --             -------------- */

  
  setOptions(argc, argv, myrank, quiet);
  if(iverbose == 0) ctx->commandline.quiet = quiet;
  
  
  if(firsttime){
    ctx->atmos.nH = NULL;
    getCPU(0, TIME_START, NULL);
    SetFPEtraps();
    ctx->geometry.Ndep = rhs_ndep;
    ctx->atmos.gravity = gravity;
    readInput_ctx(ctx);
    // readAbundance(&ctx->atmos);
    readAbundance_ctx(ctx);
    DUMMYatmos_ctx(ctx);
    oldMode = ctx->input.StokesMode;
    mpi.rank = myrank;
    ctx->mpi.rank = myrank;
  }
  //POSSIBLE BUGS HERE!!
  mpi.stop = false;
  ctx->mpi.stop = false;


  /* --- Store Ng values --- */

  sNgperiod = ctx->input.Ngperiod;
  sNgdelay = ctx->input.Ngdelay;
  sPRDNITER = ctx->input.PRD_NmaxIter;
  
  /* --- Allocate space for arrays that define structure -- --------- */
  
  ctx->geometry.tau_ref = rhs_tau;//(double *) malloc(Ndep * sizeof(double));
  ctx->geometry.cmass   = rhs_cmass;//(double *) malloc(Ndep * sizeof(double));
  ctx->geometry.height  = rhs_z;//(double *) malloc(Ndep * sizeof(double));
  ctx->atmos.T      = rhs_T;//(double *) malloc(Ndep * sizeof(double));
  ctx->atmos.ne     = rhs_nne;//(double *) malloc(Ndep * sizeof(double));
  ctx->atmos.vturb  = rhs_vturb;//(double *) malloc(Ndep * sizeof(double));
  ctx->geometry.vel = rhs_v;//(double *) malloc(Ndep * sizeof(double));
  ctx->atmos.B = rhs_B;
  ctx->atmos.gamma_B = rhs_inc;
  ctx->atmos.chi_B = rhs_azi;
  ctx->atmos.H_LTE = TRUE;
  ctx->atmos.nHtot = rhs_nhtot;
  ctx->atmos.rho = rhs_rho;
  ctx->save_popp = save_pop;
  ctx->input.StokesMode = oldMode;
  ctx->spectrum.updateJ = TRUE;  
  getCPU(1, TIME_START, NULL);
  if (ctx->input.StokesMode > NO_STOKES)
    ctx->atmos.Stokes = TRUE;
  hydrostat[0] = (int)ctx->atmos.hydrostatic;
  if(ctx->input.solve_ne >= ITERATION_EOS) hydrostat[0] = 1;
  
  ne_lte = (double*)calloc(ctx->atmos.Nspace,sizeof(double));
  memcpy(ne_lte, ctx->atmos.ne, ctx->atmos.Nspace*sizeof(double));

 

  
  /* --- Init atomic/molecular models, extra lambda 
     positions and background opac. Allocated only in the first call --- */
  
  if(firsttime){
    // readAtomicModels();
    readAtomicModels_ctx(ctx);
    // readMolecularModels();
    readMolecularModels_ctx(ctx);
    // SortLambda_j(mynw, mylambda);
    SortLambda_j_ctx(mynw, mylambda, ctx);
    // init_Background_j();
    init_Background_j_ctx(ctx);
    // Initvarious_ctx();
    Initvarious_ctx(ctx);

    /* --- Save ctx->geometry values to change back after --    ------------ */
    save_Nrays = ctx->atmos.Nrays;   save_wmu = ctx->geometry.wmu[0];
    save_muz = ctx->geometry.muz[0]; save_mux = ctx->geometry.mux[0]; save_muy = ctx->geometry.muy[0];
  }
  
  firsttime = FALSE;
  if(ctx->input.solve_ne >= ITERATION_EOS){
    if(ctx->atmos.atoms[0].active) ctx->atmos.ne_flag = TRUE;
    if(ctx->save_popp && ctx->save_popp->ne_dep){
      double *tmp1 = (double*)calloc(ctx->atmos.Nspace,sizeof(double));
      hermitian_interpolation((int)ctx->atmos.Nspace, ctx->save_popp->tau_ref, ctx->save_popp->ne_dep,
			      (int)ctx->atmos.Nspace, ctx->geometry.tau_ref, tmp1, 1);
      
      
      for(k=0;k<ctx->atmos.Nspace;++k) ctx->atmos.ne[k] *= tmp1[k];
      free(tmp1);
    }
  }
  /* --- Reallocate stuff and compute background opac --- */


  UpdateAtmosDep_ctx(ctx);
  Background_j_ctx(write_analyze_output=FALSE, equilibria_only=FALSE, ctx);
  //convertScales(&ctx->atmos, &ctx->geometry);
  //if(ctx->atmos.H->NLTEpops) hydrostat[0] = TRUE;
  printf("I GOT HERE 271!!!!\n");

  //for(i=0;i<ctx->atmos.Nspace;i++) printf("[%3d] %e\n", i, ctx->geometry.tau_ref[i]);
  
  if(!mpi.stop){
    
    /* --- Init profiles, populations and scattering --- */
    getProfiles_ctx(ctx);
    printf("I GOT HERE 277!!!!\n");
    initSolution_j_ctx( myrank, savpop, ctx);
    printf("I GOT HERE 279!!!!\n");

    if(computing_derivatives || (ctx->input.solve_ne < ITERATION_EOS)){
       read_populations_ctx(ctx->save_popp, 0, ctx);
    }
    
    
    if((savpop == 0) && 1){
      ctx->input.Ngdelay = min(15,ctx->input.Ngdelay) ;
      ctx->input.Ngperiod = min(13,ctx->input.Ngperiod) ;
      //ctx->input.PRD_NmaxIter = min(3,ctx->input.PRD_NmaxIter) ;
    }

    // for(niter=0; niter<ctx->spectrum.nPRDlines; niter++)
    // fprintf(stderr,"prdline->frac[0][0]=%e\n", ctx->spectrum.PRDlines[niter]->frac[0][0]);

    initScatter_ctx(ctx);
    printf("I GOT HERE 294!!!!\n");
    
    //getCPU(1, TIME_POLL, "Total Initialize");
    
    
    /* --- Solve radiative transfer for active ingredients -- --------- */

    Iterate_j_ctx(ctx->input.NmaxIter, ctx->input.iterLimit, &dpopmax, ctx);
    fprintf(stderr, "I GOT HERE 302!!!!\n");
    if(isnan(dpopmax) || isinf(dpopmax) || dpopmax < 0){
      mpi.stop = true;
    }

    ctx->input.Ngdelay=sNgdelay, ctx->input.Ngperiod = sNgperiod, ctx->input.PRD_NmaxIter = sPRDNITER;
    
    fprintf(stderr, "I GOT HERE 310!!!!\n");
    /* --- Adjust stokes mode in case we are running POLARIZATION_FREE --- */
    
    if(!mpi.stop){
      fprintf(stderr, "I GOT HERE 311!!!!\n");
      adjustStokesMode_ctx(ctx);
      fprintf(stderr, "I GOT HERE 312!!!!\n");
      niter = 0;
      
      while ((niter < ctx->input.NmaxScatter)) {
        fprintf(stderr, "I GOT HERE 313!!!!\n");
        if (solveSpectrum_ctx(FALSE, FALSE, 0, TRUE, ctx) <= ctx->input.iterLimit)
        {
          fprintf(stderr, "I GOT HERE 314!!!!\n");
          break;
        } 
        fprintf(stderr, "I GOT HERE 315!!!!\n");
        niter++;
      }
    } else dpopmax = 1.0e13;
  } else {
    dpopmax = 1.0e13;
    printf("I GOT HERE 321!!!!\n");
  }
  
  
  bool_t converged = dpopmax < ctx->input.iterLimit;

  /* --- Store populations if needed --- */

  if(savpop > 0 && converged){
    fprintf(stderr, "I GOT HERE 316!!!!\n");
    save_populations_ctx(ctx->save_popp, ne_lte, ctx);
    fprintf(stderr, "I GOT HERE 316!!!!\n");
  }

  fprintf(stderr, "I GOT HERE 317!!!!\n");

  free(ne_lte); ne_lte = NULL;
  
  /* --- Compute output ray --- */
  if(converged){
    ctx->atmos.Nrays     = 1;
    ctx->geometry.Nrays  = 1;
    ctx->geometry.muz[0] = muz;
    ctx->geometry.mux[0] = sqrt(1.0 - SQ(ctx->geometry.muz[0]));
    ctx->geometry.muy[0] = 0.0;
    ctx->geometry.wmu[0] = 1.0;
    ctx->spectrum.updateJ = FALSE;

    fprintf(stderr, "I GOT HERE 360!!!!\n");
    
    calculateRay();
    
    fprintf(stderr, "I GOT HERE 364!!!!\n");
    
    
    /* --- Put back previous values for ctx->geometry  --- */
    
    ctx->atmos.Nrays     = save_Nrays;
    ctx->geometry.Nrays = save_Nrays;
    ctx->geometry.muz[0] = save_muz;
    ctx->geometry.mux[0] = save_mux;
    ctx->geometry.muy[0] = save_muy;
    ctx->geometry.wmu[0] = save_wmu;
    ctx->spectrum.updateJ = TRUE;


  }

  fprintf(stderr, "I GOT HERE 380!!!!\n");

  ctx->input.StokesMode = oldMode;

  
  /* --- Copy desired ray to output arrays---*/
  
  sp->nrays = ctx->atmos.Nrays;
  sp->nlambda = ctx->spectrum.Nspect;
  if(sp->I != NULL) free(sp->I);
  if(sp->Q != NULL) free(sp->Q);
  if(sp->U != NULL) free(sp->U);
  if(sp->V != NULL) free(sp->V);
  if(sp->lambda != NULL) free(sp->V);
  sp->I = calloc(ctx->spectrum.Nspect, sizeof(double));
  sp->Q = calloc(ctx->spectrum.Nspect, sizeof(double));
  sp->U = calloc(ctx->spectrum.Nspect, sizeof(double));
  sp->V = calloc(ctx->spectrum.Nspect, sizeof(double));
  sp->lambda = calloc(ctx->spectrum.Nspect, sizeof(double));

  memcpy(sp->lambda, &ctx->spectrum.lambda[0], ctx->spectrum.Nspect * sizeof(double));

  fprintf(stderr, "I GOT HERE 390!!!!\n");
  if(converged){
    memcpy(sp->I, ctx->spectrum.I[0], ctx->spectrum.Nspect * sizeof(double));
    if(ctx->atmos.Stokes && stokes){
      memcpy(sp->Q, ctx->spectrum.Stokes_Q[0], ctx->spectrum.Nspect * sizeof(double));
      memcpy(sp->U, ctx->spectrum.Stokes_U[0], ctx->spectrum.Nspect * sizeof(double));
      memcpy(sp->V, ctx->spectrum.Stokes_V[0], ctx->spectrum.Nspect * sizeof(double));	
    }
  }

  fprintf(stderr, "I GOT HERE 400!!!!\n");


  /* --- Deallocate n & nstar --- */
  
  for (nact = 0; nact < ctx->atmos.Natom; nact++) {
    atom = &ctx->atmos.atoms[nact];
    if(atom->n == atom->nstar){
      freeMatrix((void **) atom->nstar);
      atom->nstar = NULL;
      atom->n = NULL;
    }else{
      if (atom->nstar != NULL) freeMatrix((void **) atom->nstar);
      if (atom->n != NULL) freeMatrix((void **) atom->n);
      atom->nstar = NULL;
      atom->n = NULL;
    }
  }

  if(!quiet){
    fclose(ctx->commandline.logfile);
    //if(!mpi.stop)
    remove(ctx->commandline.logfileName);
    //else exit(0);
  }

  return converged;
}


void clean_saved_populations(crhpop *save_pop){

  int ii, kr;
  
  if(save_pop->nactive > 0){
    // fprintf(stderr,"clean_saved_populations: cleaning pops\n");
    for(ii = 0;ii < save_pop->nactive;ii++){
      
      free(save_pop->pop[ii].n);
      free( save_pop->pop[ii].ntotal);
      
      save_pop->pop[ii].ntotal = NULL;
      save_pop->pop[ii].n = NULL;
      
      if(save_pop->pop[ii].nprd > 0){
	
	for(kr=0;kr<save_pop->pop[ii].nprd; kr++){
	  // fprintf(stderr,"cleaning rho [%d] -> %p \n", kr, save_pop->pop[ii].line[kr]);
	  
	  free(save_pop->pop[ii].line[kr].rho);
	}

	free(save_pop->pop[ii].line);
	
      } // nprd > 0

      save_pop->pop[ii].nprd = 0;
      
    } // nactive
    
    
    free(save_pop->pop);
    free(save_pop->lambda);
    free(save_pop->J);
    if(input.backgr_pol) free(save_pop->J20);
    free(save_pop->tau_ref);
    if(save_pop->ne_dep){
      free(save_pop->ne_dep);
      save_pop->ne_dep = NULL;
    }
    
    save_pop->pop = NULL;
    save_pop->lambda = NULL;
    save_pop->J = NULL;
    save_pop->J20 = NULL;
    save_pop->tau_ref = NULL;
  }

  save_pop->nactive = 0;
}


void clean_saved_populations_ctx(crhpop *save_pop, RHContext *ctx)
{

  int ii, kr;
  InputData *inputLocal = &ctx->input;
  // crhpop   *save_pop   = ctx->save_popp;
  
  if(save_pop->nactive > 0){
    // fprintf(stderr,"clean_saved_populations: cleaning pops\n");
    for(ii = 0;ii < save_pop->nactive;ii++){
      
      free(save_pop->pop[ii].n);
      free( save_pop->pop[ii].ntotal);
      
      save_pop->pop[ii].ntotal = NULL;
      save_pop->pop[ii].n = NULL;
      
      if(save_pop->pop[ii].nprd > 0){
	
	for(kr=0;kr<save_pop->pop[ii].nprd; kr++){
	  // fprintf(stderr,"cleaning rho [%d] -> %p \n", kr, save_pop->pop[ii].line[kr]);
	  
	  free(save_pop->pop[ii].line[kr].rho);
	}

	free(save_pop->pop[ii].line);
	
      } // nprd > 0

      save_pop->pop[ii].nprd = 0;
      
    } // nactive
    
    
    free(save_pop->pop);
    free(save_pop->lambda);
    free(save_pop->J);
    if(inputLocal->backgr_pol) free(save_pop->J20);
    free(save_pop->tau_ref);
    if(save_pop->ne_dep){
      free(save_pop->ne_dep);
      save_pop->ne_dep = NULL;
    }
    
    save_pop->pop = NULL;
    save_pop->lambda = NULL;
    save_pop->J = NULL;
    save_pop->J20 = NULL;
    save_pop->tau_ref = NULL;
  }

  save_pop->nactive = 0;
}


void read_populations(crhpop *save_pop, int flag){
  Atom *atom;
  int    niter, nact, save_Nrays, nactotal, ii, nprd, kr, kkr, copied = 0;
  AtomicLine *line;
  register int k, j, la;
  double ratio, tmp;
  double *tmp1 = (double*)malloc(atmos.Nspace*sizeof(double));
  
  if(save_pop->pop == NULL || save_pop->nactive != atmos.Nactiveatom){
    // fprintf(stderr,"read_population: atmos->Nactiveatom[%d] != save_pop->nactive[%d]\n",
    //atmos.Nactiveatom, save_pop->nactive);
    return;
  }

  
  for(nact=0;nact < save_pop->nactive;nact++){
    atom = atmos.activeatoms[nact];
    
    /* --- Check dimensions --- */
    if(atom->Nlevel != save_pop->pop[nact].nlevel || atmos.Nspace != save_pop->ndep){
      continue;
    }
    copied += 1;
    
    
    for(j= 0; j < atom->Nlevel; j++ ){
      hermitian_interpolation((int)atmos.Nspace, save_pop->tau_ref, &save_pop->pop[nact].n[j*atmos.Nspace],
      		      (int)atmos.Nspace, geometry.tau_ref, tmp1, 1);      
      
      for(k = 0; k < atmos.Nspace; k++){
	
	atom->n[j][k] = tmp1[k] * atom->nstar[j][k];
      }
    }

    // --- Check that we don't exeed total number of particles --- //
    
    for(k = 0; k < atmos.Nspace; k++){
      tmp = 0.0;
      for(j= 0; j < atom->Nlevel; j++ ) tmp += atom->n[j][k];

      tmp = atom->ntotal[k] / tmp;
      for(j= 0; j < atom->Nlevel; j++ ) atom->n[j][k] *= tmp;      
    }

    
    
    /* --- Copy rho for PRD lines? --- */
    
    if(save_pop->pop[nact].nprd > 0){
      
      /* --- Number of PRD lines in atom --- */
      nprd = 0;
      for (kr = 0;  kr < atom->Nline;  kr++){
	line = &atom->line[kr];
	if (line->PRD) nprd++;
      }
      
      if(nprd == save_pop->pop[nact].nprd){
	for(kr = 0; kr<save_pop->pop[nact].nprd; kr++){
	  
	  line = &atom->line[save_pop->pop[nact].line[kr].idx];
	  
	  if(line->PRD && (save_pop->pop[nact].line[kr].nlambda == line->Nlambda)){
	    // //   fprintf(stderr, "Copying rho array for nact=%d, line=%d\n", nact, kr);
	    //memcpy(&line->rho_prd[0][0], save_pop->pop[nact].line[kr].rho,
	    // 		   line->Nlambda*atmos.Nspace*sizeof(double));
	    for(la=0;la<line->Nlambda;la++){
	      hermitian_interpolation((int)atmos.Nspace, save_pop->tau_ref, &save_pop->pop[nact].line[kr].rho[la*atmos.Nspace],
	      			      (int)atmos.Nspace, geometry.tau_ref, line->rho_prd[la],0);
	    }
	  }else{
	    fprintf(stderr,"read_populations: BAD BOOK-KEEPING, not a PRD line, not copying rho, idx=%d, kkr=%d \n", kr, kkr );
	    
	  }
	  
	}
	
      }else{
	fprintf(commandline.logfile,"read_populations: nprd[%d] != save_pop.pop.nprd[%d], not copying PRD rho!\n", nprd,save_pop->pop[nact].nprd );
      }
      
      //}//else{
      // fprintf(stderr,"read_populations, no PRD lines for ATOM=%d\n",nact);
    }
  } // nact
  
  
  /* --- Copy radiation field --- */
  
  if(spectrum.Nspect == save_pop->nw || atmos.Nspace == save_pop->ndep){
    for(la=0;la<spectrum.Nspect;la++){
      hermitian_interpolation((int)atmos.Nspace, save_pop->tau_ref, &save_pop->J[la*atmos.Nspace],
      			      (int)atmos.Nspace, geometry.tau_ref, spectrum.J[la],0);
      //memcpy(spectrum.J[la], &save_pop->J[la*atmos.Nspace], atmos.Nspace*sizeof(double));
    
    
      if(input.backgr_pol){
	hermitian_interpolation((int)atmos.Nspace, save_pop->tau_ref, &save_pop->J20[la*atmos.Nspace],
				      (int)atmos.Nspace, geometry.tau_ref, spectrum.J20[la],0);
	//memcpy(spectrum.J20[la], &save_pop->J20[la*atmos.Nspace], atmos.Nspace*sizeof(double));
	
	//memcpy(&spectrum.J20[0][0], &save_pop->J20[0],
	//	     spectrum.Nspect*atmos.Nspace*sizeof(double));
	
	//else fprintf(stderr, "NOT READING J20!\n");
      }
    }
  }
  
  free(tmp1);
}

void read_populations_ctx(crhpop *save_pop, int flag, RHContext *ctx)
{
  Atom *atom;
  int    niter, nact, save_Nrays, nactotal, ii, nprd, kr, kkr, copied = 0;
  AtomicLine *line;
  register int k, j, la;
  double ratio, tmp;
  Geometry *geometryLocal = &ctx->geometry;
  Atmosphere *atmosLocal = &ctx->atmos;
  Spectrum *spectrumLocal = &ctx->spectrum;
  InputData *inputLocal = &ctx->input;
  // crhpop   *save_pop   = ctx->save_popp;

  double *tmp1 = (double*)malloc(atmosLocal->Nspace*sizeof(double));
  
  if(save_pop->pop == NULL || save_pop->nactive != atmosLocal->Nactiveatom){
    // fprintf(stderr,"read_population: atmos->Nactiveatom[%d] != save_pop->nactive[%d]\n",
    //atmosLocal->Nactiveatom, save_pop->nactive);
    return;
  }

  
  for(nact=0;nact < save_pop->nactive;nact++){
    atom = atmosLocal->activeatoms[nact];
    
    /* --- Check dimensions --- */
    if(atom->Nlevel != save_pop->pop[nact].nlevel || atmosLocal->Nspace != save_pop->ndep){
      continue;
    }
    copied += 1;
    
    
    for(j= 0; j < atom->Nlevel; j++ ){
      hermitian_interpolation((int)atmosLocal->Nspace, save_pop->tau_ref, &save_pop->pop[nact].n[j*atmosLocal->Nspace],
      		      (int)atmosLocal->Nspace, geometryLocal->tau_ref, tmp1, 1);      
      
      for(k = 0; k < atmosLocal->Nspace; k++){
	
	atom->n[j][k] = tmp1[k] * atom->nstar[j][k];
      }
    }

    // --- Check that we don't exeed total number of particles --- //
    
    for(k = 0; k < atmosLocal->Nspace; k++){
      tmp = 0.0;
      for(j= 0; j < atom->Nlevel; j++ ) tmp += atom->n[j][k];

      tmp = atom->ntotal[k] / tmp;
      for(j= 0; j < atom->Nlevel; j++ ) atom->n[j][k] *= tmp;      
    }

    
    
    /* --- Copy rho for PRD lines? --- */
    
    if(save_pop->pop[nact].nprd > 0){
      
      /* --- Number of PRD lines in atom --- */
      nprd = 0;
      for (kr = 0;  kr < atom->Nline;  kr++){
	line = &atom->line[kr];
	if (line->PRD) nprd++;
      }
      
      if(nprd == save_pop->pop[nact].nprd){
	for(kr = 0; kr<save_pop->pop[nact].nprd; kr++){
	  
	  line = &atom->line[save_pop->pop[nact].line[kr].idx];
	  
	  if(line->PRD && (save_pop->pop[nact].line[kr].nlambda == line->Nlambda)){
	    // //   fprintf(stderr, "Copying rho array for nact=%d, line=%d\n", nact, kr);
	    //memcpy(&line->rho_prd[0][0], save_pop->pop[nact].line[kr].rho,
	    // 		   line->Nlambda*atmosLocal->Nspace*sizeof(double));
	    for(la=0;la<line->Nlambda;la++){
	      hermitian_interpolation((int)atmosLocal->Nspace, save_pop->tau_ref, &save_pop->pop[nact].line[kr].rho[la*atmosLocal->Nspace],
	      			      (int)atmosLocal->Nspace, geometryLocal->tau_ref, line->rho_prd[la],0);
	    }
	  }else{
	    fprintf(stderr,"read_populations: BAD BOOK-KEEPING, not a PRD line, not copying rho, idx=%d, kkr=%d \n", kr, kkr );
	    
	  }
	  
	}
	
      }else{
	fprintf(commandline.logfile,"read_populations: nprd[%d] != save_pop.pop.nprd[%d], not copying PRD rho!\n", nprd,save_pop->pop[nact].nprd );
      }
      
      //}//else{
      // fprintf(stderr,"read_populations, no PRD lines for ATOM=%d\n",nact);
    }
  } // nact
  
  
  /* --- Copy radiation field --- */
  
  if(spectrumLocal->Nspect == save_pop->nw || atmosLocal->Nspace == save_pop->ndep){
    for(la=0;la<spectrumLocal->Nspect;la++){
      hermitian_interpolation((int)atmosLocal->Nspace, save_pop->tau_ref, &save_pop->J[la*atmosLocal->Nspace],
      			      (int)atmosLocal->Nspace, geometryLocal->tau_ref, spectrumLocal->J[la],0);
      //memcpy(spectrumLocal->J[la], &save_pop->J[la*atmosLocal->Nspace], atmosLocal->Nspace*sizeof(double));
    
    
      if(inputLocal->backgr_pol){
	hermitian_interpolation((int)atmosLocal->Nspace, save_pop->tau_ref, &save_pop->J20[la*atmosLocal->Nspace],
				      (int)atmosLocal->Nspace, geometryLocal->tau_ref, spectrumLocal->J20[la],0);
	//memcpy(spectrumLocal->J20[la], &save_pop->J20[la*atmosLocal->Nspace], atmosLocal->Nspace*sizeof(double));
	
	//memcpy(&spectrumLocal->J20[0][0], &save_pop->J20[0],
	//	     spectrumLocal->Nspect*atmosLocal->Nspace*sizeof(double));
	
	//else fprintf(stderr, "NOT READING J20!\n");
      }
    }
  }
  
  free(tmp1);
}


void save_populations(crhpop *save_pop, double *ne_lte){

  Atom *atom;
  int    niter, nact, save_Nrays, nactotal, ii, kr, nprd;
  AtomicLine *line;
  register int k,j;

  
  
  // fprintf(stderr,"save_pop: nactive=%d\n", save_pop->nactive);
  
  if(save_pop->nactive > 0) clean_saved_populations(save_pop);

  save_pop->nactive = atmos.Nactiveatom;
  save_pop->ndep = atmos.Nspace;
  save_pop->nw = spectrum.Nspect;
  //
  save_pop->pop =     (crhatom*) malloc(atmos.Nactiveatom * sizeof(crhatom));
  save_pop->lambda =  (double*) malloc(spectrum.Nspect*sizeof(double));
  //
  save_pop->J =   (double*) malloc(spectrum.Nspect* atmos.Nspace*sizeof(double));
  save_pop->tau_ref = (double*) malloc(atmos.Nspace*sizeof(double));
  save_pop->ne_dep = (double*) calloc(atmos.Nspace, sizeof(double));

  
  if(input.backgr_pol)
    save_pop->J20 = (double*) malloc(spectrum.Nspect* atmos.Nspace*sizeof(double));
  else
    save_pop->J20 = NULL;

  memcpy(save_pop->tau_ref, geometry.tau_ref, sizeof(double)*atmos.Nspace);
  for(k=0;k<atmos.Nspace; ++k) save_pop->ne_dep[k] = atmos.ne[k]/ne_lte[k];
  
  
  /* --- copy J, J20 and lambda ---*/
  
  memcpy(&save_pop->J[0], &spectrum.J[0][0],
	 spectrum.Nspect*atmos.Nspace*sizeof(double));
  
  if(input.backgr_pol){
    memcpy(&save_pop->J20[0], &spectrum.J20[0][0],
	   spectrum.Nspect*atmos.Nspace*sizeof(double));
  }// else  fprintf(stderr, "NOT STORING J20!\n");
    
  memcpy(&save_pop->lambda[0], spectrum.lambda, spectrum.Nspect * sizeof(double));

  /* --- Copy populations for each active atom --- */
  
  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    
    atom = atmos.activeatoms[nact];
    save_pop->pop[nact].nlevel = atom->Nlevel;
    save_pop->pop[nact].converged  = true;

    /* --- Allocate arrays ---*/
    save_pop->pop[nact].n = (double*) malloc( atom->Nlevel * atmos.Nspace*sizeof(double));
    save_pop->pop[nact].ntotal = (double*) malloc(atmos.Nspace*sizeof(double));


    /* --- copy populations and ntotal---*/
    // memcpy(&save_pop->pop[nact].n[0], &atom->n[0][0],
    //atom->Nlevel*atmos.Nspace*sizeof(double));
    memcpy(&save_pop->pop[nact].ntotal[0], &atom->ntotal[0], atmos.Nspace*sizeof(double));

    for(j=0;j<atom->Nlevel;j++)
      for(k=0;k<atmos.Nspace;k++)
	save_pop->pop[nact].n[j*atmos.Nspace+k] = atom->n[j][k] / atom->nstar[j][k]; //atom->ntotal[k];
    

    
    /* --- Check number of PRD lines in atom --- */
    
    nprd = 0;
    for (kr = 0;  kr < atom->Nline;  kr++)
      if (atom->line[kr].PRD) nprd++;
    

    /* --- Allocate arrays to store rho for each PRD line --- */
    
    save_pop->pop[nact].nprd = nprd;
    
    if(nprd >0){
      
      save_pop->pop[nact].line = (crhprd*) malloc(nprd * sizeof(crhprd));
      ii = 0;
      
      for (kr = 0;  kr < atom->Nline;  kr++){
	
	line = &atom->line[kr];
	
	if (line->PRD){
	  save_pop->pop[nact].line[ii].idx = kr;
	  save_pop->pop[nact].line[ii].nlambda = line->Nlambda;
	  save_pop->pop[nact].line[ii].rho = (double*)
	    malloc(line->Nlambda * atmos.Nspace * sizeof(double));
	  //
	  //  fprintf(stderr,"ii=%d, srho=%p, lrho=%p\n",
	  //	  ii, &save_pop->pop[nact].line[ii].rho[0], &line->rho_prd[0][0]);
	  //
	  memcpy(save_pop->pop[nact].line[ii].rho, &line->rho_prd[0][0],
		 line->Nlambda*atmos.Nspace*sizeof(double));
	  ii++;
	} // PRD line
      } // kr
    }// nprd > 0
    
  } // nact

  
}

void save_populations_ctx(crhpop *save_pop, double *ne_lte, RHContext *ctx)
{

  Atom *atom;
  int    niter, nact, save_Nrays, nactotal, ii, kr, nprd;
  AtomicLine *line;
  register int k,j;
  InputData *inputLocal = &ctx->input;
  Atmosphere *atmosLocal = &ctx->atmos;
  Geometry *geometryLocal = &ctx->geometry;
  Spectrum *spectrumLocal = &ctx->spectrum;
  // crhpop   *save_pop   = ctx->save_popp;

  
  
  // fprintf(stderr,"save_pop: nactive=%d\n", save_pop->nactive);

  if(save_pop->nactive > 0) clean_saved_populations_ctx(save_pop, ctx);

  save_pop->nactive = atmosLocal->Nactiveatom;
  save_pop->ndep = atmosLocal->Nspace;
  save_pop->nw = spectrumLocal->Nspect;
  //
  save_pop->pop =     (crhatom*) malloc(atmosLocal->Nactiveatom * sizeof(crhatom));
  save_pop->lambda =  (double*) malloc(spectrumLocal->Nspect*sizeof(double));
  //
  save_pop->J =   (double*) malloc(spectrumLocal->Nspect* atmosLocal->Nspace*sizeof(double));
  save_pop->tau_ref = (double*) malloc(atmosLocal->Nspace*sizeof(double));
  save_pop->ne_dep = (double*) calloc(atmosLocal->Nspace, sizeof(double));

  
  if(inputLocal->backgr_pol)
    save_pop->J20 = (double*) malloc(spectrumLocal->Nspect* atmosLocal->Nspace*sizeof(double));
  else
    save_pop->J20 = NULL;

  memcpy(save_pop->tau_ref, geometryLocal->tau_ref, sizeof(double)*atmosLocal->Nspace);
  for(k=0;k<atmosLocal->Nspace; ++k) save_pop->ne_dep[k] = atmosLocal->ne[k]/ne_lte[k];
  
  
  /* --- copy J, J20 and lambda ---*/
  
  memcpy(&save_pop->J[0], &spectrumLocal->J[0][0],
	 spectrumLocal->Nspect*atmosLocal->Nspace*sizeof(double));
  
  if(inputLocal->backgr_pol){
    memcpy(&save_pop->J20[0], &spectrumLocal->J20[0][0],
	   spectrumLocal->Nspect*atmosLocal->Nspace*sizeof(double));
  }// else  fprintf(stderr, "NOT STORING J20!\n");
    
  memcpy(&save_pop->lambda[0], spectrumLocal->lambda, spectrumLocal->Nspect * sizeof(double));

  /* --- Copy populations for each active atom --- */
  
  for (nact = 0;  nact < atmosLocal->Nactiveatom;  nact++) {
    
    atom = atmosLocal->activeatoms[nact];
    save_pop->pop[nact].nlevel = atom->Nlevel;
    save_pop->pop[nact].converged  = true;

    /* --- Allocate arrays ---*/
    save_pop->pop[nact].n = (double*) malloc( atom->Nlevel * atmosLocal->Nspace*sizeof(double));
    save_pop->pop[nact].ntotal = (double*) malloc(atmosLocal->Nspace*sizeof(double));


    /* --- copy populations and ntotal---*/
    // memcpy(&save_pop->pop[nact].n[0], &atom->n[0][0],
    //atom->Nlevel*atmosLocal->Nspace*sizeof(double));
    memcpy(&save_pop->pop[nact].ntotal[0], &atom->ntotal[0], atmosLocal->Nspace*sizeof(double));

    for(j=0;j<atom->Nlevel;j++)
      for(k=0;k<atmosLocal->Nspace;k++)
	save_pop->pop[nact].n[j*atmosLocal->Nspace+k] = atom->n[j][k] / atom->nstar[j][k]; //atom->ntotal[k];
    

    
    /* --- Check number of PRD lines in atom --- */
    
    nprd = 0;
    for (kr = 0;  kr < atom->Nline;  kr++)
      if (atom->line[kr].PRD) nprd++;
    

    /* --- Allocate arrays to store rho for each PRD line --- */
    
    save_pop->pop[nact].nprd = nprd;
    
    if(nprd >0){
      
      save_pop->pop[nact].line = (crhprd*) malloc(nprd * sizeof(crhprd));
      ii = 0;
      
      for (kr = 0;  kr < atom->Nline;  kr++){
	
	line = &atom->line[kr];
	
	if (line->PRD){
	  save_pop->pop[nact].line[ii].idx = kr;
	  save_pop->pop[nact].line[ii].nlambda = line->Nlambda;
	  save_pop->pop[nact].line[ii].rho = (double*)
	    malloc(line->Nlambda * atmosLocal->Nspace * sizeof(double));
	  //
	  //  fprintf(stderr,"ii=%d, srho=%p, lrho=%p\n",
	  //	  ii, &save_pop->pop[nact].line[ii].rho[0], &line->rho_prd[0][0]);
	  //
	  memcpy(save_pop->pop[nact].line[ii].rho, &line->rho_prd[0][0],
		 line->Nlambda*atmosLocal->Nspace*sizeof(double));
	  ii++;
	} // PRD line
      } // kr
    }// nprd > 0
    
  } // nact

  
}






/* ------- end ---------------------------- rhf1d.c ----------------- */
