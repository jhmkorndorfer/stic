/* ------- file: -------------------------- rh.h --------------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Fri Feb 24 08:24:45 2012 --

       --------------------------                      ----------RH-- */

/* --- General include file for Rybicki & Hummer family of radiative
       transfer codes. --                              -------------- */
 

/* --- Include stdio.h here because it defines NULL
     - Include rpc/types.h here because it defines bool_t
     - Include pthread.h for POSIX threads 
       --                                              -------------- */

#include <stdio.h> 
#include <rpc/types.h>
#include <pthread.h>

enum Topology       {ONE_D_PLANE, TWO_D_PLANE, SPHERICAL_SYMMETRIC,
                     THREE_D_PLANE};
enum FeautrierOrder {STANDARD, HERMITE};
enum Interpolation  {LINEAR, SPLINE, EXP_SPLINE};
enum solution       {UNKNOWN=-1, LTE_POPULATIONS, ZERO_RADIATION,
                     OLD_POPULATIONS, PESC, NEW_J, OLD_J};
enum StokesMode     {NO_STOKES, FIELD_FREE, POLARIZATION_FREE, FULL_STOKES};
enum PRDangle       {PRD_ANGLE_INDEP, PRD_ANGLE_APPROX, PRD_ANGLE_DEP};
enum VoigtAlgorithm {ARMSTRONG, RYBICKI, HUI_ETAL, HUMLICEK, LOOKUP};


#define  MAX_LINE_SIZE      512
#define  MAX_MESSAGE_LENGTH 512
#define  MAX_KEYWORD_SIZE   32
#define PRD_FILE_TEMPLATE1  "scratch/PRD_%.1s_%d-%d_p%d.dat"
#define PRD_FILE_TEMPLATE   "scratch/PRD_%s_%d-%d_p%d.dat"

#define  LG10  2.30258509299404568402

#ifndef MAX 
#define  MAX(x, y)    ((x) > (y) ? (x):(y))
#endif
#ifndef MIN
#define  MIN(x, y)    ((x) < (y) ? (x):(y))
#endif

#define SWAPPOINTER(a, b)  {void *tmp;  tmp = a;  a = b;  b = tmp;}
#define SWAPDOUBLE(a, b)   {double tmp;  tmp = a;  a = b;  b = tmp;}

#define  SQ(x)        ((x)*(x))
#define  CUBE(x)      ((x)*(x)*(x))
#define  POW10(x)     exp(LG10 * (x))

#define MODULO(n, N) (((n) >= 0) ? ((n)%(N)) : ((((n)%(N))+(N)) % (N)))

#define FLUX_DOT_OUT  "flux.out"
#define J20_DOT_OUT   "J20.out"

#define PERMISSIONS  0666


/* --- Function prototypes --                          -------------- */

void   SetFPEtraps(void);
void   checkNread(int Nread, int Nrequired, const char *routineName,
		  int checkPoint);
void   setOptions(int argc, char *argv[], int iproc, int quiet);
double vproject(int k, int mu);
void   Bproject(void);
bool_t StopRequested(void);
char **getWords(char *label, char *separator, int *count);


/* --- Matrix manipulation --                          -------------- */

char   **matrix_char(int Nrow, int Ncol);
int    **matrix_int(int Nrow, int Ncol);
double **matrix_double(int Nrow, int Ncol);
short  **matrix_short(int Nrow, int Ncol);
unsigned char **matrix_uchar(int Nrow, int Ncol);
float **matrix_float(int Nrow, int Ncol);
double **d2dim(int x1l,int x1h,int x2l,int x2h);
float **f2dim(int x1l,int x1h,int x2l,int x2h);
void del_d2dim(double **p,int x1l, int x2l);
void del_f2dim(float **p,int x1l, int x2l);
float *f1dim(int x1l,int x1h);
double *d1dim(int x1l,int x1h);
int *i1dim(int x1l,int x1h);
void del_f1dim(float *p,int x1l);
void del_d1dim(double *p,int x1l);
void del_i1dim(int *p,int x1l);




void   freeMatrix(void **Matrix);
void   SolveLinearEq(int N, double **A, double *b, bool_t improve);
void   SolveLinearSvd(int n, double **A, double *b);


/* --- Interpolation routines --                       -------------- */

int   qsascend(const void *v1, const void *v2);
int   qsdescend(const void *v1, const void *v2);
     
void  Hunt(int n, double *array, double value, int *iLower);
void  Locate(int n, double *array, double value, int *theIndex);
     
void   Linear(int Ntable, double *xtable, double *ytable,
	      int N, double *x, double *y, bool_t hunt);
double BiLinear(int Na, double *a_table, double a,
		int Nb, double *b_table, double b,
		double **f, bool_t hunt);

void  splineCoef(int Ntable, double *xtable, double *ytable);
void  splineEval(int N, double *x, double *y, bool_t hunt);
void  exp_splineCoef(int Ntable, double *xtable, double *ytable,
		     double tension);
void  exp_splineEval(int N, double *x, double *y, bool_t hunt);

void  cc_kernel(double s, double *u);
double cubeconvol(int Nx, int Ny, double *f, double x, double y);

/* --- Special functions --                            -------------- */

void   GaussLeg(double x1, double x2, double *x, double *w, int n);
double Voigt(double a, double v, double *F, enum VoigtAlgorithm algorithm);
double gammln(double xx);
     
void   w2(double dtau, double *w);
void   w3(double dtau, double *w);

void   U3(double dtau, double *U);

/* ------- end ---------------------------- rh.h -------------------- */
