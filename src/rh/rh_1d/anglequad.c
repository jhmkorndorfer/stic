/* ------- file: -------------------------- anglequad.c -------------

       Version:       rh2.0, 1-D plane-parallel
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Jan 15 15:47:11 2004 --

       --------------------------                      ----------RH-- */

/* --- Get angle quadrature. --                        -------------- */

 
#include <stdlib.h>
#include <math.h>

// #include "rh.h"
#include "rhf1d.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;


/* ------- begin -------------------------- getAngleQuad.c ---------- */

void getAngleQuad(Geometry *geometry)
{
  register int mu;

  /* --- Copy the number of rays to the geometry structure.
   Note: Nrays is read into the atmos structure in readinput.c -- --- */

  geometry->Nrays = atmos.Nrays;

  geometry->mux = (double *) malloc(geometry->Nrays * sizeof(double));
  geometry->muy = (double *) malloc(geometry->Nrays * sizeof(double));
  geometry->muz = (double *) malloc(geometry->Nrays * sizeof(double));
  geometry->wmu = (double *) malloc(geometry->Nrays * sizeof(double));

  GaussLeg(0.0, 1.0, geometry->muz, geometry->wmu, geometry->Nrays);

  /* --- In the 1-D plane case assume that all the rays lie in the z-x
         plane and are pointed in the direction of positive x. -- --- */

  for (mu = 0;  mu < geometry->Nrays;  mu++) {
    geometry->muy[mu] = 0.0;
    geometry->mux[mu] = sqrt(1.0 - SQ(geometry->muz[mu]));
  }
}
/* ------- end ---------------------------- getAngleQuad.c ---------- */


/* ------- begin -------------------------- getAngleQuad_ctx.c ---------- */

void getAngleQuad_ctx(RHContext *ctx)
{
  Atmosphere *atmosLocal = &ctx->atmos;
  Geometry *geometryLocal = &ctx->geometry;
  register int mu;

  /* --- Copy the number of rays to the geometry structure.
   Note: Nrays is read into the atmos structure in readinput.c -- --- */

  geometryLocal->Nrays = atmosLocal->Nrays;

  geometryLocal->mux = (double *) malloc(geometryLocal->Nrays * sizeof(double));
  geometryLocal->muy = (double *) malloc(geometryLocal->Nrays * sizeof(double));
  geometryLocal->muz = (double *) malloc(geometryLocal->Nrays * sizeof(double));
  geometryLocal->wmu = (double *) malloc(geometryLocal->Nrays * sizeof(double));

  GaussLeg(0.0, 1.0, geometryLocal->muz, geometryLocal->wmu, geometryLocal->Nrays);

  /* --- In the 1-D plane case assume that all the rays lie in the z-x
         plane and are pointed in the direction of positive x. -- --- */

  for (mu = 0;  mu < geometryLocal->Nrays;  mu++) {
    geometryLocal->muy[mu] = 0.0;
    geometryLocal->mux[mu] = sqrt(1.0 - SQ(geometryLocal->muz[mu]));
  }
}
/* ------- end ---------------------------- getAngleQuad_ctx.c ---------- */
