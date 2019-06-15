/*
 * =====================================================================================
 *
 *       Filename:  Mesh.h
 *
 *    Description:  Save the 2D mesh for finite element discretization
 *
 *        Version:  1.0
 *        Created:  2012/04/01 03时09分20秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  hhxie@lsec.cc.ac.cn
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef __EIGENPAIRSOLVER__
#define __EIGENPAIRSOLVER__


#include <Matrix.h>
#include <Enumerations.h>

/** Solve the Eigenvalue problem for the Laplace equation: double symmetric matrix */
void DSEigenSolver(MATRIX *A,  MATRIX *M, int nev, int indictor, double tau, 
		 double *Evals, double **Evecs, bool Change);

/** Solve the Eigenvalue problem for the Laplace equation: double non-symmetric matrix */
void DNEigenSolver(MATRIX *A,  MATRIX *M, int nev, int indictor, double tau, 
		double *Evals, double **Evecs, bool Change);	

/** Solve the Eigenvalue problem for the Laplace equation: double non-symmetric matrix */
/** which comes from the Eigenvalue Multigrid method */
// void DNEigenSolver(MATRIX *D, MATRIX *M, double **vector_D, double **vector_M,
// 	      double** D22,double** M22,int nev, int indictor, double tau,
// 	      double *Evals, double **Evecs, bool Change);

/** Solve the Eigenvalue problem for the Laplace equation: complex non-symmetric matrix */
// void ZNEigenSolver(MATRIX *A,  MATRIX *M, int nev, int indictor, double tau, 
// 		 std::complex<double> *Evals, std::complex<double> **Evecs);

/** Solve the Eigenvalue problem for the General equation: double non-symmetric matrix */
// void DNEigenSolver(int N_matrices_rows, int N_matrices_column, int n_sqmatrices, 
// 			 MATRIX **sqmatrices, int n_matrices, MATRIX **matrices,
// 			 int *ordering, int *diag_ordering, int n_sqmassmatrices, 
// 			 MATRIX **sqmassmatrices,int n_massmatrices, 
// 			 MATRIX **massmatrices, int *massordering, int *massdiag_ordering,
// 			 bool Destroy, bool Singular, int nev, int indictor, double tau, 
// 		         double *Evals, double **Evecs);
			
/** Solve the Eigenvalue problem for the Stokes equation with LPS methd: 
    double non-symmetric matrix */
/*int StokesDNEigenSolver(MATRIX *A, MATRIX *B1, MATRIX *B2, MATRIX *C,
			MATRIX *M, int nev, int indictor,
			double tau, double *Evals, double **Evecs);*/		      
#endif
