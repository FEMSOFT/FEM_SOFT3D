#include <phg.h>

#if (PHG_VERSION_MAJOR <= 0 && PHG_VERSION_MINOR < 9)
# undef ELEMENT
typedef SIMPLEX ELEMENT;
#endif

#include <string.h>
#include <math.h>

void
build_linear_system(SOLVER *solver, DOF *u_h, DOF *f_h, DISCRETEFORM *discreteform, RHST  *rhs, GRID *g)
{
    int i, j;
    ELEMENT *e;
    assert(u_h->dim == 1);
    ForAllElements(g, e) {
	int N = DofGetNBas(u_h, e);	/* number of bases in the element */
	FLOAT A[N][N], rhs[N], buffer[N];
	INT I[N];
	/* compute \int \grad\phi_j \cdot \grad\phi_i making use of symmetry */
	for (i = 0; i < N; i++) {
	    I[i] = phgSolverMapE2L(solver, 0, e, i);//local number of basis
	    for (j = 0; j <= i; j++) {
		A[j][i] = A[i][j] =
		    /* stiffness */
		    phgQuadGradBasDotGradBas(e, u_h, j, u_h, i, QUAD_DEFAULT) +
		    /* mass */
		    a * phgQuadBasDotBas(e, u_h, j, u_h, i, QUAD_DEFAULT);
	    }
	}

	/* loop on basis functions */
	for (i = 0; i < N; i++) {
	    if (phgDofDirichletBC(u_h, e, i, func_u, buffer, rhs+i,
					DOF_PROJ_NONE)) {
		phgSolverAddMatrixEntries(solver, 1, I + i, N, I, buffer); 
	    }
	    else {	/* interior node */
		/* right hand side = \int f * phi_i */
		phgQuadDofTimesBas(e, f_h, u_h, i, QUAD_DEFAULT, rhs + i);
		phgSolverAddMatrixEntries(solver, 1, I + i, N, I, A[i]); 
	    }
	}
	phgSolverAddRHSEntries(solver, N, I, rhs);
    }//for all element
}



int
PSolver(DISCRETEFORM *discreteform, RHST *rhs, GRID *g,int nev, double *evals, double **evecs)
{
    INT periodicity = 0 /* X_MASK | Y_MASK | Z_MASK */;
    INT i;
    ELEMENT *e;
    DOF *u_h, *f_h, *grad_u, *error, *u;
    SOLVER *solver;

    /* The discrete solution */
    if (FALSE) {
	/* u_h is h-p type */
	HP_TYPE *hp = phgHPNew(g, HP_HB);
	ForAllElements(g, e)
	e->hp_order = DOF_DEFAULT->order + GlobalElement(g, e->index) % 3;
	phgHPSetup(hp, FALSE);
	u_h = phgHPDofNew(g, hp, 1, "u_h", DofInterpolation);
	phgHPFree(&hp);
    }
    else {
	/* u_h is non h-p type */
	u_h = phgDofNew(g, DOF_DEFAULT, 1, "u_h", DofInterpolation);
    }
    phgDofSetDataByValue(u_h, 0.0);

    /* RHS function */
    f_h = phgDofNew(g, DOF_DEFAULT, 1, "f_h",  func_f);

    /* The analytic solution */
    u = phgDofNew(g, DOF_ANALYTIC, 1, "u", func_u);

    while (TRUE) {
	double t0 = phgGetTime(NULL);
	if (phgBalanceGrid(g, 1.2, 1, NULL, 0.))
	    phgPrintf("Repartition mesh, load imbalance: %lg\n",
			(double)g->lif);
	phgPrintf("%"dFMT" DOF, %"dFMT" elements, %"dFMT
		  " submeshes, load imbalance: %lg\n",
			DofGetDataCountGlobal(u_h), g->nleaf_global, g->nprocs,
			(double)g->lif);

	solver = phgSolverCreate(SOLVER_DEFAULT, u_h, NULL);
	phgPrintf("  DOF: %"dFMT", unknowns: %"dFMT
		  ", Dirichlet bdry: %"dFMT"\n",
		DofGetDataCountGlobal(u_h), solver->rhs->map->nglobal,
		solver->rhs->map->bdry_nglobal);
	build_linear_system(solver, u_h, f_h,discreteform,rhs,g);
	phgSolverSolve(solver, TRUE, u_h, NULL);
	phgSolverDestroy(&solver);
   }

#if 0
    phgPrintf("Final mesh written to \"%s\".\n",
	phgExportDX(g, "simplest.dx", u_h, error, NULL));
#elif 0
    phgPrintf("Final mesh written to \"%s\".\n",
	phgExportVTK(g, "simplest.vtk", u_h, error, NULL));
#endif

    phgDofFree(&u);
    phgDofFree(&u_h);
    phgDofFree(&f_h);
    phgDofFree(&error);

    return 0;
}

