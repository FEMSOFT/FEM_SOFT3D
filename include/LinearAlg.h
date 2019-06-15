#ifndef __LINEARALG__
#define __LINEARALG__
#include "Matrix.h"
void MaxVec(double *a, int dim,double *value,int *index);
void MaxAbsVec(double *a,int dim,double *value,int *index);
void ScalVec(double alpha, double *p, int size);
void SumVecSelf(double *a, double alpha, double *b, double beta, int size);
void SumVec(double *a, double alpha, double *b, double beta,double *c,int size);
void ScalVec(double alpha, double *p, int size);
void NumVec(double alpha, double *p, double *q, int size);
double NormVec(double *b, int size);
void MatrixDotVec(MATRIX *Matrix,double *x, double *r);
void MatMatVec(MATRIX *MatrixA, MATRIX *MatrixB,double *x, double *r);
void TripleMatrix(MATRIX *matrixa,MATRIX *matrixb,MATRIX *matrixc,int num, double *r);
int power(int a,int b);
double VecDotVec(double *a, double *b, int size);
void VecDotVecs(double *a, double *b, int size,double *c);
void ShowVec(double *p, int size);
void AssignVec(double *p, double *q, int size);
void TransposeMatrix(MATRIX *A, MATRIX *B);
void MatMatVec(MATRIX *MatrixA, MATRIX *MatrixB,double *x, double *r);
void VecMatVec(double *r,MATRIX *MatrixA,double *c);
void  VecMatVecs(double *r,MATRIX *MatrixA,double *s,double *c);
void MatMatMat(double *a,double *b,double *d);
void VecArrayArray(double localval[4], double v1,double v2,double v3,double v4,double A[4][4],double B[4][4]);
void VecArrayVec(double lastvalstiff[1], double v1,double v2,double v3,double v4,double A[4][4]);
void VecArrayVecs(double lastvalstiff[1], double v1,double v2,double v3,double v4,double A[4][4], double w1,double w2,double w3,double w4);
#endif
