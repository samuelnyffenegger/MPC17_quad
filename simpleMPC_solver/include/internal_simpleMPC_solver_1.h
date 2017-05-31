/*
internal_simpleMPC_solver_1 : A fast customized optimization solver.

Copyright (C) 2013-2016 EMBOTECH GMBH [info@embotech.com]. All rights reserved.


This software is intended for simulation and testing purposes only. 
Use of this software for any commercial purpose is prohibited.

This program is distributed in the hope that it will be useful.
EMBOTECH makes NO WARRANTIES with respect to the use of the software 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
PARTICULAR PURPOSE. 

EMBOTECH shall not have any liability for any damage arising from the use
of the software.

This Agreement shall exclusively be governed by and interpreted in 
accordance with the laws of Switzerland, excluding its principles
of conflict of laws. The Courts of Zurich-City shall have exclusive 
jurisdiction in case of any dispute.

*/

#include <stdio.h>

#ifndef __internal_simpleMPC_solver_1_H__
#define __internal_simpleMPC_solver_1_H__

/* DATA TYPE ------------------------------------------------------------*/
typedef double internal_simpleMPC_solver_1_FLOAT;

typedef double internal_simpleMPC_solver_1INTERFACE_FLOAT;

/* SOLVER SETTINGS ------------------------------------------------------*/
/* print level */
#ifndef internal_simpleMPC_solver_1_SET_PRINTLEVEL
#define internal_simpleMPC_solver_1_SET_PRINTLEVEL    (2)
#endif

/* timing */
#ifndef internal_simpleMPC_solver_1_SET_TIMING
#define internal_simpleMPC_solver_1_SET_TIMING    (1)
#endif

/* Numeric Warnings */
/* #define PRINTNUMERICALWARNINGS */

/* maximum number of iterations  */
#define internal_simpleMPC_solver_1_SET_MAXIT         (200)	

/* scaling factor of line search (affine direction) */
#define internal_simpleMPC_solver_1_SET_LS_SCALE_AFF  (internal_simpleMPC_solver_1_FLOAT)(0.9)      

/* scaling factor of line search (combined direction) */
#define internal_simpleMPC_solver_1_SET_LS_SCALE      (internal_simpleMPC_solver_1_FLOAT)(0.95)  

/* minimum required step size in each iteration */
#define internal_simpleMPC_solver_1_SET_LS_MINSTEP    (internal_simpleMPC_solver_1_FLOAT)(1E-08)

/* maximum step size (combined direction) */
#define internal_simpleMPC_solver_1_SET_LS_MAXSTEP    (internal_simpleMPC_solver_1_FLOAT)(0.995)

/* desired relative duality gap */
#define internal_simpleMPC_solver_1_SET_ACC_RDGAP     (internal_simpleMPC_solver_1_FLOAT)(0.0001)

/* desired maximum residual on equality constraints */
#define internal_simpleMPC_solver_1_SET_ACC_RESEQ     (internal_simpleMPC_solver_1_FLOAT)(1E-06)

/* desired maximum residual on inequality constraints */
#define internal_simpleMPC_solver_1_SET_ACC_RESINEQ   (internal_simpleMPC_solver_1_FLOAT)(1E-06)

/* desired maximum violation of complementarity */
#define internal_simpleMPC_solver_1_SET_ACC_KKTCOMPL  (internal_simpleMPC_solver_1_FLOAT)(1E-06)


/* RETURN CODES----------------------------------------------------------*/
/* solver has converged within desired accuracy */
#define internal_simpleMPC_solver_1_OPTIMAL      (1)

/* maximum number of iterations has been reached */
#define internal_simpleMPC_solver_1_MAXITREACHED (0)

/* no progress in line search possible */
#define internal_simpleMPC_solver_1_NOPROGRESS   (-7)

/* fatal internal error - nans occurring */
#define internal_simpleMPC_solver_1_NAN  (-10)


/* PARAMETERS -----------------------------------------------------------*/
/* fill this with data before calling the solver! */
typedef struct internal_simpleMPC_solver_1_params
{
    /* vector of size 1 */
    internal_simpleMPC_solver_1_FLOAT p_1[1];

    /* vector of size 14 */
    internal_simpleMPC_solver_1_FLOAT p_2[14];

    /* vector of size 24 */
    internal_simpleMPC_solver_1_FLOAT p_3[24];

    /* vector of size 18 */
    internal_simpleMPC_solver_1_FLOAT p_4[18];

} internal_simpleMPC_solver_1_params;


/* OUTPUTS --------------------------------------------------------------*/
/* the desired variables are put here by the solver */
typedef struct internal_simpleMPC_solver_1_output
{
    /* vector of size 1 */
    internal_simpleMPC_solver_1_FLOAT o_1[1];

    /* vector of size 1 */
    internal_simpleMPC_solver_1_FLOAT o_2[1];

    /* vector of size 1 */
    internal_simpleMPC_solver_1_FLOAT o_3[1];

    /* vector of size 1 */
    internal_simpleMPC_solver_1_FLOAT o_4[1];

} internal_simpleMPC_solver_1_output;


/* SOLVER INFO ----------------------------------------------------------*/
/* diagnostic data from last interior point step */
typedef struct internal_simpleMPC_solver_1_info
{
    /* iteration number */
    int it;

	/* number of iterations needed to optimality (branch-and-bound) */
	int it2opt;
	
    /* inf-norm of equality constraint residuals */
    internal_simpleMPC_solver_1_FLOAT res_eq;
	
    /* inf-norm of inequality constraint residuals */
    internal_simpleMPC_solver_1_FLOAT res_ineq;

    /* primal objective */
    internal_simpleMPC_solver_1_FLOAT pobj;	
	
    /* dual objective */
    internal_simpleMPC_solver_1_FLOAT dobj;	

    /* duality gap := pobj - dobj */
    internal_simpleMPC_solver_1_FLOAT dgap;		
	
    /* relative duality gap := |dgap / pobj | */
    internal_simpleMPC_solver_1_FLOAT rdgap;		

    /* duality measure */
    internal_simpleMPC_solver_1_FLOAT mu;

	/* duality measure (after affine step) */
    internal_simpleMPC_solver_1_FLOAT mu_aff;
	
    /* centering parameter */
    internal_simpleMPC_solver_1_FLOAT sigma;
	
    /* number of backtracking line search steps (affine direction) */
    int lsit_aff;
    
    /* number of backtracking line search steps (combined direction) */
    int lsit_cc;
    
    /* step size (affine direction) */
    internal_simpleMPC_solver_1_FLOAT step_aff;
    
    /* step size (combined direction) */
    internal_simpleMPC_solver_1_FLOAT step_cc;    

	/* solvertime */
	internal_simpleMPC_solver_1_FLOAT solvetime;   

} internal_simpleMPC_solver_1_info;









/* SOLVER FUNCTION DEFINITION -------------------------------------------*/
/* examine exitflag before using the result! */
#ifdef _cplusplus
extern "C" {
#endif
int internal_simpleMPC_solver_1_solve(internal_simpleMPC_solver_1_params* params, internal_simpleMPC_solver_1_output* output, internal_simpleMPC_solver_1_info* info, FILE* fs);

#ifdef _cplusplus
}
#endif

#endif