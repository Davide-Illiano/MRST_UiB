/* Copyright 2011 (c) Jostein R. Natvig <Jostein.R.Natvig at sintef.no> */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#ifdef MATLAB_MEX_FILE
#include "grid.h"
#include "reordersequence.h"
#include "tarjan.h"
#include "nlsolvers.h"
#include "twophase.h"
#include "twophasetransport.h"

#include <mex.h>
extern int interrupt_signal;
#define print   mexPrintf
#define malloc  mxMalloc
#define calloc  mxCalloc
#define realloc mxRealloc
#define free    mxFree

#else
#include <opm/core/grid.h>
#include <opm/core/transport/reorder/reordersequence.h>
#include <opm/core/transport/reorder/tarjan.h>
#include <opm/core/transport/reorder/nlsolvers.h>
#include <opm/core/transport/reorder/twophase.h>
#include <opm/core/transport/reorder/twophasetransport.h>

#define print printf

#endif




void twophasetransport(
    const double            *porevolume,
    const double            *source,
    double                   dt,
    struct UnstructuredGrid *grid,
    const double            *darcyflux,
    const int               *satnum,
    double                  *saturation)
{
    int    i;

    int    ncomponents;
    int    *sequence;
    int    *components;

    struct SolverData *data = init_solverdata(grid, darcyflux,
                                              porevolume, source, satnum, dt, saturation);


    struct NonlinearSolverCtrl ctrl;


    /* Compute sequence of single-cell problems */
    sequence   = malloc(  grid->number_of_cells    * sizeof *sequence);
    components = malloc(( grid->number_of_cells+1) * sizeof *components);

    compute_sequence(grid, darcyflux, sequence, components, &ncomponents);
    assert(ncomponents == grid->number_of_cells);



    ctrl.maxiterations = 20;
    ctrl.nltolerance   = 1e-9;
    ctrl.method        = REGULAFALSI;
    ctrl.initialguess  = 0.5;

    /* Assume all strong components are single-cell domains. */
    for(i=0; i<grid->number_of_cells; ++i)
    {
#ifdef MATLAB_MEX_FILE
        if (interrupt_signal)
        {
            mexPrintf("Reorder loop interrupted by user: %d of %d "
                      "cells finished.\n", i, grid->number_of_cells);
            break;
        }
#endif
        solvecell(data, &ctrl, sequence[i]);
        /*print("iterations:%d residual:%20.16f\n",
         *  ctrl.iterations, ctrl.residual);*/
    }

    destroy_solverdata(data);

    free(sequence);
    free(components);
}

/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
