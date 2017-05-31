/*
Header file containing definitions for C interface of simpleMPC_solver,
 a fast costumized optimization solver.
*/

#include <stdio.h>

#ifndef __simpleMPC_solver_H__
#define __simpleMPC_solver_H__

/* DATA TYPE ------------------------------------------------------------*/
typedef double simpleMPC_solver_FLOAT;

/* SOLVER SETTINGS ------------------------------------------------------*/
/* print level */
#ifndef simpleMPC_solver_SET_PRINTLEVEL
#define simpleMPC_solver_SET_PRINTLEVEL    (2)
#endif

/* PARAMETERS -----------------------------------------------------------*/
/* fill this with data before calling the solver! */
typedef struct simpleMPC_solver_params
{
	/* column vector of length 11 */
	simpleMPC_solver_FLOAT inputvariables[11];

} simpleMPC_solver_params;


/* OUTPUTS --------------------------------------------------------------*/
/* the desired variables are put here by the solver */
typedef struct simpleMPC_solver_output
{
	/* column vector of length 4 */
	simpleMPC_solver_FLOAT u0_variable[4];

} simpleMPC_solver_output;


/* SOLVER INFO ----------------------------------------------------------*/
/* diagnostic data from last interior point step */
typedef struct simpleMPC_solver_info
{
	/* iteration number */
	int it;

	/* number of iterations needed to optimality (branch-and-bound) */
	int it2opt;

	/* inf-norm of equality constraint residuals */
	simpleMPC_solver_FLOAT res_eq;

	/* inf-norm of inequality constraint residuals */
	simpleMPC_solver_FLOAT res_ineq;

	/* primal objective */
	simpleMPC_solver_FLOAT pobj;

	/* dual objective */
	simpleMPC_solver_FLOAT dobj;

	/* duality gap := pobj - dobj */
	simpleMPC_solver_FLOAT dgap;

	/* relative duality gap := |dgap / pobj | */
	simpleMPC_solver_FLOAT rdgap;

	/* duality measure */
	simpleMPC_solver_FLOAT mu;

	/* duality measure (after affine step) */
	simpleMPC_solver_FLOAT mu_aff;

	/* centering parameter */
	simpleMPC_solver_FLOAT sigma;

	/* number of backtracking line search steps (affine direction) */
	int lsit_aff;

	/* number of backtracking line search steps (combined direction) */
	int lsit_cc;

	/* step size (affine direction) */
	simpleMPC_solver_FLOAT step_aff;

	/* step size (combined direction) */
	simpleMPC_solver_FLOAT step_cc;

	/* solvertime */
	simpleMPC_solver_FLOAT solvetime;

} simpleMPC_solver_info;


/* SOLVER FUNCTION DEFINITION -------------------------------------------*/
/* examine exitflag before using the result! */
int simpleMPC_solver_solve(simpleMPC_solver_params* params, simpleMPC_solver_output* output, simpleMPC_solver_info* info, FILE* fs);

#endif
