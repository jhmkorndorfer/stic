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

  // For now here I will add the tentative functions for RHF1D thread safe version...
  // Functions below are in order of appearnce in the rhf1d function from rhf1d.c. This also includes subfunctions, for instance: readMolecularModels_ctx calls readMolecule_ctx and readMolecularLines_ctx etc.
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
  void Initvarious_ctx(RHContext *ctx);
  void UpdateAtmosDep_ctx(RHContext *ctx);
  void Background_j_ctx(bool_t write_analyze_output, bool_t equilibria_only, RHContext *ctx);
  void Solve_ne_ctx(double *ne, bool_t fromscratch, RHContext *ctx);
  double getKuruczpf_ctx(Element *element, int stage, int k, RHContext *ctx);
  void SetLTEQuantities_ctx(RHContext *ctx);
  void LTEpops_ctx(Atom *atom, bool_t Debeye, RHContext *ctx);
  void CollisionRate_ctx(struct Atom *atom, char **fp_atom, RHContext *ctx);
  double ar85cea_ctx(int i, int j, int k, struct Atom *atom, RHContext *ctx);
  double summers_ctx(int i, int j, double nne, struct Atom *atom, RHContext *ctx);
  int atomnr_ctx(char ID[ATOM_ID_WIDTH+1], RHContext *ctx);
  void atomnm_ctx(int anr,char *cseq, RHContext *ctx);
  void FixedRate_ctx(Atom *atom, RHContext *ctx);
  void FixedRateOne_ctx(Atom *atom, int k, RHContext *ctx);
  void readMolecules_ctx(char *fileName, RHContext *ctx);
  void ChemicalEquilibrium_ctx(int NmaxIter, double iterLimit, RHContext *ctx);
  void getfjk_ctx(Element *element, double ne, int k, double *fjk, double *dfjk, RHContext *ctx);
  void Thomson_ctx(double *chi, RHContext *ctx);
  bool_t Hminus_bf_ctx(double lambda, double *chi, double *eta, RHContext *ctx);
  bool_t Hminus_ff_ctx(double lambda, double *chi, RHContext *ctx);
  bool_t Hminus_ff_long_ctx(double lambda, double *chi, RHContext *ctx);
  bool_t OH_bf_opac_ctx(double lambda, double *chi, double *eta, RHContext *ctx);
  bool_t CH_bf_opac_ctx(double lambda, double *chi, double *eta, RHContext *ctx);
  bool_t Hydrogen_bf_ctx(double lambda, double *chi, double *eta, RHContext *ctx);
  void Hydrogen_ff_ctx(double lambda, double *chi, RHContext *ctx);
  bool_t Rayleigh_ctx(double lambda, Atom *atom, double *scatt, RHContext *ctx);
  bool_t H2plus_ff_ctx(double lambda, double *chi, RHContext *ctx);
  bool_t Rayleigh_H2_ctx(double lambda, double *scatt, RHContext *ctx);
  bool_t H2minus_ff_ctx(double lambda, double *chi, RHContext *ctx);
  bool_t Metal_bf_ctx(double lambda, int Nmetal, struct Atom *metals, double *chi, double *eta, RHContext *ctx);
  flags passive_bb_ctx(double lambda, int nspect, int mu, bool_t to_obs, double *chi, double *eta, double *chip, RHContext *ctx);
  void Damping_ctx(AtomicLine *line, double *adamp, RHContext *ctx);
  double vproject_ctx(int k, int mu, RHContext *ctx);
  void VanderWaals_ctx(AtomicLine *line, double *GvdW, RHContext *ctx);
  void Stark_ctx(AtomicLine *line, double *GStark, RHContext *ctx);
  void StarkLinear_ctx(AtomicLine *line, double *GStark, RHContext *ctx);
  flags rlk_opacity_ctx(double lambda, int nspect, int mu, bool_t to_obs, double *chi, double *eta, double *scatt, double *chip, RHContext *ctx);
  void LTEpops_elem_ctx(Element *element, RHContext *ctx);
  double RLKProfile_ctx(RLK_Line *rlk, int k, int mu, bool_t to_obs, double lambda, double *phi_Q, double *phi_U, double *phi_V, double *psi_Q, double *psi_U, double *psi_V, RHContext *ctx);
  flags MolecularOpacity_ctx(double lambda, int nspect, int mu, bool_t to_obs, double *chi, double *eta, double *chip, RHContext *ctx);
  double MolProfile_ctx(MolecularLine *mrt, int k, int mu, bool_t to_obs, double lambda, double *phi_Q, double *phi_U, double *phi_V, double *psi_Q, double *psi_U, double *psi_V, RHContext *ctx);
  int writeBackground_j_ctx(int la, int mu, bool_t to_obs, double *chi_c, double *eta_c, double *sca_c, double *chip_c, RHContext *ctx);
  void allocateBack_ctx(int nspect, RHContext *ctx);
  void getProfiles_ctx(RHContext *ctx);
  void Profile_ctx(AtomicLine *line, RHContext *ctx);
  void writeProfile_ctx(AtomicLine *line, int lamu, double *phi, RHContext *ctx);
  void MolecularProfile_ctx(MolecularLine *mrt, RHContext *ctx);
  void MolecularDamping_ctx(MolecularLine *mrt, double *adamp, RHContext *ctx);
  void initSolution_j_ctx( int myrank, int savepop, RHContext *ctx);
  void initSolution_alloc2_ctx(int myrank, RHContext *ctx);
  void zeroRadiation_ctx(Atom *atom, int nact, RHContext *ctx);
  void initGammaAtom_ctx(Atom *atom, double cswitch, RHContext *ctx);
  void statEquil_ctx(Atom *atom, int isum, RHContext *ctx);
  void COcollisions_ctx(struct Molecule *molecule, RHContext *ctx);
  void H2collisions_ctx(struct Molecule *molecule, RHContext *ctx);
  void initScatter_ctx(RHContext *ctx);
  bool_t readRadRate_ctx(Atom *atom, RHContext *ctx);
  bool_t xdr_radrate_ctx(XDR *xdrs, Atom *atom, RHContext *ctx);
  void PRDScatter_ctx(AtomicLine *PRDline, enum Interpolation representation, RHContext *ctx);
  void PRDAngleApproxScatter_ctx(AtomicLine *PRDline, enum Interpolation representation, RHContext *ctx);
  void PRDAngleScatter_ctx(AtomicLine *PRDline, enum Interpolation representation, RHContext *ctx);
  void readImu_ctx(int nspect, int mu, bool_t to_obs, double *I, RHContext *ctx);
  double solveSpectrum_ctx(bool_t eval_operator, bool_t redistribute, int iter, bool_t synth_all, RHContext *ctx);
  void zeroRates_ctx(bool_t redistribute, RHContext *ctx);
  double Formal_ctx(int nspect, bool_t eval_operator, bool_t redistribute, int iter, RHContext *ctx);
  void alloc_as_ctx(int nspect, bool_t crosscoupling, RHContext *ctx);
  bool_t containsPolarized_ctx(ActiveSet *as, RHContext *ctx);
  bool_t containsBoundBound_ctx(ActiveSet *as, RHContext *ctx);
  bool_t containsPRDline_ctx(ActiveSet *as, RHContext *ctx);
  void readJlambda_ctx(int nspect, double *J, RHContext *ctx);
  void readJ20lambda_ctx(int nspect, double *J20, RHContext *ctx);
  int readBackground_j_ctx(int la, int mu, bool_t to_obs, RHContext *ctx);
  void Opacity_ctx(int nspect, int mu, bool_t to_obs, bool_t initialize, RHContext *ctx);
  void readProfile_ctx(AtomicLine *line, int lamu, double *phi, RHContext *ctx);
  void addtoCoupling_ctx(int nspect, RHContext *ctx);
  void PiecewiseStokesBezier3_ctx(int nspect, int mu, bool_t to_obs, double *chi, double **S, double **I, double *Psi, RHContext *ctx);
  void StokesK_ctx(int nspect, int k, double chi_I, double K[4][4], RHContext *ctx);
  void PiecewiseStokes_ctx(int nspect, int mu, bool_t to_obs, double *chi_I, double **S, double **I, double *Psi, RHContext *ctx);
  void Piecewise_Hermite_1D_ctx(int nspect, int mu, bool_t to_obs, double *chi, double *S, double *I, double *Psi, RHContext *ctx);
  void Piecewise_1D_ctx(int nspect, int mu, bool_t to_obs, double *chi, double *S, double *I, double *Psi, RHContext *ctx);
  void addtoGamma_ctx(int nspect, double wmu, double *I, double *Psi, RHContext *ctx);
  bool_t containsActive_ctx(ActiveSet *as, RHContext *ctx);
  void addtoRates_ctx(int nspect, int mu, bool_t to_obs, double wmu, double *I, bool_t redistribute, RHContext *ctx);
  void writeImu_ctx(int nspect, int mu, bool_t to_obs, double *I, RHContext *ctx);
  double Feautrier_ctx(int nspect, int mu, double *chi, double *S, enum FeautrierOrder F_order, double *P, double *Psi, RHContext *ctx);
  void writeJlambda_ctx(int nspect, double *J, RHContext *ctx);
  void writeJ20lambda_ctx(int nspect, double *J20, RHContext *ctx);


















































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

