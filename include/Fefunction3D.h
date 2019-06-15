/*
 * =====================================================================================
 *
 *       Filename:  Fefunction3D.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2014年10月20日 14时22分45秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Xu Fei (Fly), xufei@lsec.cc.ac.cn
 *        Company:  
 *
 * =====================================================================================
 */
#ifndef __FEFUNCTION3D__
#define __FEFUNCTION3D__
#include "Fespace3D.h"
#include "Mesh3D.h"
#include "Enumerations.h"
#include "Base3D.h"


typedef struct Fefunction3D_ {  
  //自由度的个数
  int DOF_ALL;
  //有限元函数在每个自由度的值
  double *Values;
  //相对应的有限元空间
  Fespace3D *Fespace;  
} Fefunction3D;


Fefunction3D* BuildFefunction3D(Fespace3D *fesp,double *values);
// compute the value at a point on an element
double GetFefunctValue(Fefunction3D *Fefunct,ELEMENT3D *elem3D,
                     double Quad_x,double Quad_y,double Quad_z, 
                     double Quad_xi,double Quad_eta,double Quad_zeta,
                     MultiIndex3D MultiIndex);
// compute the values for multi derivatives at a point on an element
void GetFefunctValues(Fefunction3D *Fefunct,ELEMENT3D *elem3D,
                     double Quad_x,double Quad_y,double Quad_z, 
                     double Quad_xi,double Quad_eta,double Quad_zeta,
		     int N_MultiIndex,MultiIndex3D *MultiIndex,
		     double *Values);

// Interpolation of the Fefunction and Function2D
void Interpolation3D(Fefunction3D *fefun3D,Functions3D *funct3D);
// Error estimate of the exact solution 
// N_Funct2d = N_Derives, one function corresponding to one derives
// L^2 error: N_Funct2d = N_Derives =1, multiindex2D = D00
// H^1 error: N_Funct2d = N_Derives =3, multiindex2D = [D00,D10,D01]
double ErrorFefunction(Fefunction3D *fefunction1,Fefunction3D *fefunction2,int N_MultiIndex, MultiIndex3D *multiindex3D);
/*
double ErrorEstimate3D(Fefunction3D *fefun3D,Functions3D *Exact3D, 
		       int N_Derives,MultiIndex3D *multiindex3D,int Quad_Points);
double ErrorEstimateEigen3D(Fefunction3D *fefun3D,Functions3D *Exact3D, 
		       int N_Derives,MultiIndex3D *multiindex3D,int Quad_Points);
double ErrorEstimateEigen3DExactSolution(Fefunction3D *fefun3D,Functions3D *Exact3D, 
		       int N_Derives,MultiIndex3D *multiindex3D, int Quad_Points);
double ErrorEstimateEigen3DVecExactSolution(double *values,Fespace3D *fesp,Functions3D *Exact3D, 
		       int N_Derives,MultiIndex3D *multiindex3D,double *val);
double DifferenceEstimate3D(Fefunction3D *fefun1_3D,Fefunction3D *fefun2_3D,
			    int N_Derives,MultiIndex3D *multiindex1_3D, 
			    MultiIndex3D *multiindex2_3D, int Quad_Points);
double ErrorFefunction(Fefunction3D *fefunction1,Fefunction3D *fefunction2,int N_MultiIndex, MultiIndex3D *multiindex3D);
double ErrorFefunction1(Fefunction3D *fefunction1,Fefunction3D *fefunction2,int N_MultiIndex, MultiIndex3D *multiindex3D);
void ErrorEstimateMul(MULTILEVEL *Multilevel,Fefunction3D **fefunction, double *error,int N_MultiIndex, MultiIndex3D *multiindex3D);//计算所有的误差
void ErrorEstimate4NumerSolution(MATRIX **matrix, double **fefunction, Fespace3D *fesp,int num,int N_MultiIndex,
                                 MultiIndex3D *multiindex3D,double *error);
void ErrorEstimateMultiEigen(MATRIX **matrix, double **fefunction,Fespace3D *fesp,int num,int N_MultiIndex,
                             MultiIndex3D *multiindex3D,double *error);//计算所有的误差
void ErrorEstimateMultiEigen1(MATRIX **matrix, double **fefunction,Fespace3D *fesp,int num,int N_MultiIndex,
                             MultiIndex3D *multiindex3D,double *error);//计算所有的误差
							 void OutPutFefunction3D(Fefunction3D * Fefunct,char *filename);
*/
void FreeFefunction3D(Fefunction3D *funct);
void Normalization(Fefunction3D *fefunction);
void Normalization1(Fefunction3D *fefunction);
void NormalizationUserfunc(Fefunction3D *fefunction);
void NormalizationVector1(int dof, double *values);
void WriteSolution(double *solution, int num, char *file);
double*  QuadFefunction(Fefunction3D *fefunction,double *val);
void GetMaxValue(double *vector,int length, double *norm);
void Quad_Vector(Fespace3D *fespace,double *vector, double *val);
void Normalization_Vector(Fespace3D *fespace,double *vector);


#endif
