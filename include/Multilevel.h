#ifndef __LINEARALG__
#define __LINEARALG__
#include "Matrix.h"
#include "DiscreteForm3D.h"
#include "Mesh3D.h"
#include "Fespace3D.h"
#include "MG.h"
#include "Rhst.h"
#include "Matrix.h"

void AssembleMatrixMix(MATRIX *Matrix,MATRIX *MatrixB);
void ScfIteration(DISCRETEFORM3D *discreteformA, DISCRETEFORM3D *discreteformB,int nev,double **vecs,
                  double *evals,double tol,int maxinum);
void BuildProlongRestrictFa(Fespace3D *coarse_spsace,Fespace3D *finer_space, MATRIX *Prolong, MATRIX *Restrict);
void BuildProlongRestrictAn(Fespace3D *coarse_spsace,Fespace3D *finer_space, MATRIX *Prolong, MATRIX *Restrict);
void SetSubMatrix(MATRIX *A,double *AA, int n_row,int *row,int n_col,int *col);
void GetLocalProlong(ELEMENT3D *coarse_elem3D,BASEFUNCTION3D *coarse_base,int n_coarse_base,
		     ELEMENT3D *finer_elem3D, BASEFUNCTION3D *finer_base,int n_finer_base,
		     ELEMENT3D *elem3D,double *Local_Prolong);
void ProduceRelativeElem3D(ELEMENT3D *coarse_elem3D, ELEMENT3D *finer_elem3D, ELEMENT3D *elem3D);
void ProduceRelativeCoord(ELEMENT3D *coarse_elem3D, double dx, double dy, double dz, double *elem3D);
void AssembleRHSTAllMix(RHST *Rhs);
void BuildProlongFa(Fespace3D *coarse_space,Fespace3D *finer_space, MATRIX *Prolong);
void BuildRestrictFa(Fespace3D *coarse_spsace,Fespace3D *finer_space, MATRIX *Restrict);
void BuildProlongAn(Fespace3D *coarse_spsace,Fespace3D *finer_space, MATRIX *Prolong);
void BuildRestrictAn(Fespace3D *coarse_spsace,Fespace3D *finer_space, MATRIX *Restrict);
MATRIX* ExtensionMatrix(MATRIX *matrix,int num);
MATRIX* ExtensionMatrixSim(MATRIX *matrix,int num);

double* ProlongH2h(MULTILEVEL *multilevel,int num_H,int num_h,double *vector);
void ProlongH2h2(MULTILEVEL *multilevel,int num_H,int num_h,double *vector,double *tmp_vector2);
double* Restricth2H(MULTILEVEL *multilevel,int num_h,int num_H,double *vector);
void Restricth2H2(MULTILEVEL *multilevel,int num_h,int num_H,double *vector,double *vector_H);
MATRIX* BuildMatrixMul(DISCRETEFORM3D *discreteform);
void AssembleMatrixMul(MATRIX *Matrix,Fefunction3D *Fefunct);//matrix是要合成的矩阵，matrixA是细网格刚度矩阵
void AssembleMatrixMulSimplify(MATRIX *Matrix,Fefunction3D *Fefunct);//matrix是要合成的矩阵，matrixA是细网格刚度矩阵
void AssembleMatrixHinh(MATRIX *Matrix,Fefunction3D *Fefunct);//当利用到细网格单元上的信息时，在细网格上操作

void ScfIterationCoupleBEC(DISCRETEFORM3D *discreteformA1, DISCRETEFORM3D *discreteformA2, DISCRETEFORM3D*discreteformB,
		                  int nev,double *eigenvectors1,double *evals1,double *eigenvectors2,double *evals2,double tol, int maxinum);
void ScfIterationCoupleBECGeneral(DISCRETEFORM3D **discreteformA, DISCRETEFORM3D **discreteformB, int nev, double ***vecs, 
		                          double **evals, int number, double tol, int maxinum, double alpha);
void ScfIterationCoupleBECCorrection(DISCRETEFORM3D **discreteformA, DISCRETEFORM3D **discreteformB, RHST **Rhs, int nev, double ***vecs,
                            double **evals, int number, Fefunction3D **Fefunct, double tol, int maxinum, double alpha);
void ScfIterationMix( DISCRETEFORM3D *discreteformA, DISCRETEFORM3D*discreteformB, RHST *Rhs, int nev, double **vecs, double *evals,double tol, int maxinum);
void ScfIterationMul(DISCRETEFORM3D *discreteformA, DISCRETEFORM3D *discreteformB,int nev,double *vecs,
                     double *evals,Fefunction3D *Fefunct,double tol,int maxinum);
void ScfIterationCorrection(DISCRETEFORM3D *discreteformA, DISCRETEFORM3D *discreteformB, RHST *Rhs, int nev, double *vecs,
                            double *evals, Fefunction3D *Fefunct, double tol, int maxinum);

void NonScfIteration(DISCRETEFORM3D *discreteform,RHST *rhs,double *sol,double tol,int maxinum);
void NonScfIterationWithMatrix(MATRIX *matrix,RHST *rhs,double *sol,double tol,int maxinum);
void NonScfCorrection(MATRIX *stiff_matrix, MATRIX *matrix,DISCRETEFORM3D *discreteformH, RHST *rhs, double *result, Fefunction3D *fefuncbase, double tol, int maxinum);

#endif
