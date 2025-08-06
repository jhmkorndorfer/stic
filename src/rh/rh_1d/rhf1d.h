#ifndef RHF1DF_H
#define RHF1DF_H
#ifdef __cplusplus
extern "C" {
#endif
#include "rh.h"
#include <rpc/types.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <errno.h>
#include <unistd.h>
#include <sys/stat.h>

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
// #include "rhf1d.h"
#include "dummyatmos.h"
#include "sortlambda_j.h"
#include "initial_j.h"
#include "scatter_j.h"
#include "iterate_j.h"

  
  typedef struct{
    long *atom_file_pos;
  } rhinfo;
  
  typedef struct{
    bool_t allocated;
    double **chi_b, **eta_b, **sca_b, **chip_b;
  } rhbgmem;
  
  typedef struct{
    int nlambda, nrays;
    double *lambda, *I, *Q, *U, *V;
  } ospec;

  typedef struct{
    int nlambda, idx;
    double *rho;
  } crhprd;
  
  typedef struct{
    int nlevel, converged, nprd;
    double *n, *ntotal;
    crhprd *line;
  } crhatom;
  
  typedef struct{
    int nactive, ndep, nw;
    double *J, *J20, *lambda, *tau_ref, *ne_dep;
    crhatom *pop;
  } crhpop;
  
  typedef struct{
    int rank, verb, iter;
    bool_t stop;
    FILE *logfile;
    char filename[300];
  } MPI_t;

  typedef struct {
  Atmosphere atmos;
  Geometry geometry;
  Spectrum spectrum;
  InputData input;
  CommandLine commandline;
  ProgramStats stats;
  rhinfo io;
  rhbgmem *bmem;
  crhpop *save_popp;
  MPI_t mpi;

  // These were static before
  bool initialized;
  int save_Nrays;
  double save_muz, save_mux, save_muy, save_wmu;
  enum StokesMode oldMode;
  BackgroundData bgdat;

} RHContext;

  
  
  void save_populations(RHContext* ctx, crhpop *save_pop, double *ne_lte);
  void read_populations(RHContext* ctx, crhpop *save_pop, int flag);
  void readInput_ctx(RHContext *ctx);
  void clean_saved_populations(RHContext* ctx, crhpop *save_pop_ref);
  void UpdateAtmosDep(void);
  void Initvarious();
  void calculateRay(void);

  //  For now here I will add the tentative functions for RHF1D thread safe version...
  void DUMMYatmos_ctx(RHContext *ctx);
  void getAngleQuad_ctx(RHContext *ctx);
  void readAbundance_ctx(RHContext *ctx);
  void readAtomicModels_ctx(RHContext *ctx);
  void distribute_nH_ctx(RHContext *ctx);
  void readAtom_ctx(Atom *atom, char *atomFileName, bool_t active, RHContext *ctx);
  bool_t getBarklemactivecross_ctx(AtomicLine *line, RHContext *ctx);
  void readPopulations_ctx(Atom *atom, RHContext *ctx);
  void readMolecularModels_ctx(RHContext *ctx);
  void readMolecule_ctx(Molecule *molecule, char *fileName, bool_t active, RHContext *ctx);
  void LTEmolecule_ctx(Molecule *molecule, RHContext *ctx);
  void readMolecularLines_ctx(struct Molecule *molecule, char *line_data, RHContext *ctx);
  void SortLambda_j_ctx(int mynw, double *mylambda, RHContext *ctx);
  void init_Background_j_ctx(RHContext *ctx);
  void readKuruczLines_ctx(char *inputFile, RHContext *ctx);
  bool_t RLKdet_level_ctx(char* label, RLK_level *level, RHContext *ctx);
  void getUnsoldcross_ctx(RLK_Line *rlk, RHContext *ctx);
  bool_t getBarklemcross_ac_ctx(Barklemstruct *bs, RLK_Line *rlk, RHContext *ctx);






  bool_t rhf1d(RHContext* ctx, float muz, int rhs_ndep, double *rhs_T, double *rhs_rho, 
	       double *rhs_nne, double *rhs_vturb, double *rhs_v, 
	       double *rhs_B, double *rhs_inc, double *rhs_azi,
	       double *rhs_z, double *rhs_nhtot, double *rhs_ltau,
	       double *rhs_cmass, double gravity, bool_t stokes, ospec *sp,
	       crhpop *save_pop, int mynw, double *mylambda, int myrank, int savpop,
	       int iverbose, int *hydrostat, int computing_derivatives);

  void   Redistribute_j(int NmaxIter, double iterLimit, double iprec);
  void hermitian_interpolation(int n, double *x, double *y, int nn, double *xp, double *yp, int lo);
  

  
#ifdef __cplusplus
}
# endif

#endif

