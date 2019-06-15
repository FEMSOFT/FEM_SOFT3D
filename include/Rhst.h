/*
 * =====================================================================================
 *
 *       Filename:  Rhst.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2014年10月20日 14时28分15秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Xu Fei (Fly), xufei@lsec.cc.ac.cn
 *        Company:  
 *
 * =====================================================================================
 */
#ifndef __RHST__
#define __RHST__
#include "Fespace3D.h"
#include "Fefunction3D.h"
#include "Mesh3D.h"
#include "Base3D.h"
#include "Enumerations.h"
#include "Constants.h"
#include "Quadrature3D.h"
#include "Quadrature2D.h"
#include "Quadrature1D.h"


typedef struct RHST_ {
  /** test finite element space */
  Fespace3D *Test_Space;  
  /** the Multiindex for testspace */
  int N_Test_MultiIndex;  
  MultiIndex3D *Test_MultiIndex;  
  //辅助有限元函数的个数
  int N_AuxFeFun;   
  //存储有限元函数的向执政指针  
  Fefunction3D **AuxFeFun;    
  //每个有限元函数需要求解的导数个数
  int *N_AuxFeFun_MultiIndex;
  //具体的有限元函数的导数编号
  MultiIndex3D *AuxFeFun_MultiIndex;
  //所有的有限元函数导数的向量的维数
  int N_AuxFeFun_Values;
  //具体用来存储有限元函数导数的向量
  double *AuxFeFun_Values;   
  /** The information of the quadrature scheme */
  /**用户函数*/
  int N_UserFunction;
  Functions3D  **UserFunction;
  //变化的常数
  int N_FeConst;
  double *FeConst;
  //三维的积分信息
  Quadrature3D *Quad3D;  
  //二维的积分信息
  Quadrature2D *Quad2D;  
  //一维的积分信息
  Quadrature1D *Quad1D;
  /** Discrete form for RHST */
  //在体上的离散变分形式函数
  DiscreteFormRHSVolu *DiscreteFormVolu;  
  //在面上的离散变分形式函数
  DiscreteFormRHSFace *DiscreteFormFace;  
  //在线上的离散变分形式函数
  DiscreteFormRHSLine *DiscreteFormLine;
  //相对于边界条件的变分形式函数
  //DiscreteFormRHSLine *DiscreteFormBD;
  /** the dimension of this rhs vector */
  int DOF_ALL;  
  /** save the entries of this vector */
  double *Entries;
  
} RHST;


RHST *BuildRHST(Fespace3D *test_fesp,int n_test_MultiIndex,MultiIndex3D *test_MultiIndex,
			    DiscreteFormRHSVolu *discreteformvolu,int Quad_Points3D);
RHST *BuildRHSTBD(Fespace3D *test_fesp,int n_test_MultiIndex,MultiIndex3D *test_MultiIndex,
			    DiscreteFormRHSVolu *discreteformvolu,DiscreteFormRHSFace *discreteformface,
			    int Quad_Points3D ,int Quad_Points2D, int Quad_Points1D);
RHST *BuildRHSTAuxFeFunct(Fespace3D *test_fesp,int n_test_MultiIndex,MultiIndex3D *test_MultiIndex,
                          int n_auxFefunct,Fefunction3D **auxFefunct,int *n_auxFefun_MultiIndex,
                          MultiIndex3D *auxFefun_MultiIndex, DiscreteFormRHSVolu *discreteformvolu,
			  int Quad_Points3D);
RHST *BuildRHSTAll(Fespace3D *test_fesp,int n_test_MultiIndex,
                          MultiIndex3D *test_MultiIndex,                        
                        int n_auxFefunct,Fefunction3D **auxFefunct,int *n_auxFefun_MultiIndex,
                        MultiIndex3D *auxFefun_MultiIndex,DiscreteFormRHSVolu *discreteformvolu, 
		        DiscreteFormRHSFace *discreteformface,int Quad_Points3D,int Quad_Points2D,
		        int Quad_Points1D);
RHST *BuildRHSTFeConst(  Fespace3D *test_fesp,int n_test_MultiIndex,
                        MultiIndex3D *test_MultiIndex,                        
                        int n_auxFefunct,Fefunction3D **auxFefunct,int *n_auxFefun_MultiIndex,
                        MultiIndex3D *auxFefun_MultiIndex,
                        int N_FeConst,double *FeConst,
		        DiscreteFormRHSVolu *discreteformvolu,int Quad_Points3D);
RHST *BuildRHSTAllMix(  Fespace3D *test_fesp,int n_test_MultiIndex,
                        MultiIndex3D *test_MultiIndex,                        
                        int n_auxFefunct,Fefunction3D **auxFefunct,int *n_auxFefun_MultiIndex,
                        MultiIndex3D *auxFefun_MultiIndex,
                        int N_UserFunction,Functions3D **UserFunction,int N_FeConst,double *FeConst,
		        DiscreteFormRHSVolu *discreteformvolu,int Quad_Points3D,int Quad_Points2D,
		        int Quad_Points1D);

void AssembleRHST(RHST *Rhs);
//Sub functions for the building 
void SubAssembleRhs3D(ELEMENT3D *elem3D,BASEFUNCTION3D *test_base, double *Test_Values,
		      int n_test_multiindex,MultiIndex3D *test_multiindex, double *Test_MultiIndex_Values, int *Dim_Values_AuxFeFun, 
		      int N_AuxFeFun, Fefunction3D **AuxFeFuns, int *N_AuxFeFun_MultiIndex, MultiIndex3D *AuxFeFun_MultiIndex,
		      int N_AuxFeFun_Values, double *AuxFeFun_Values,int N_FeConst,double *FeConst,double  Quad_xi,double Quad_eta,
		       double Quad_zeta,double  Quad_W, 
		      DiscreteFormRHSVolu *discreteformvolu, double *Rhs_Entries);

void SubAssembleRhs3DSimple(RHST *RHS, ELEMENT3D *elem3D,double *Test_Values,
		      double *Test_MultiIndex_Values, int *Dim_Values_AuxFeFun, 
		      double  Quad_xi,double Quad_eta,double Quad_zeta,double Quad_W,  double *Rhs_Entries);

void SubAssembleRhs3DLine(int ID_bd,double *normal, double *tangent, double local_direction,
			                RHST *RHS,ELEMENT3D *elem3D,double *Test_Values,
			                double *Test_MultiIndex_Values,  double  Quad_xi, 
			                double Quad_eta,double Quad_W, double length, 
				        DiscreteFormRHSLine *discreteformline, 
					double *Rhs_Entries);
void SubAssembleRhs3DAll(ELEMENT3D *elem3D,BASEFUNCTION3D *test_base, double *Test_Values,
		      int n_test_multiindex,MultiIndex3D *test_multiindex, double *Test_MultiIndex_Values,
		      int *Dim_Values_AuxFeFun, 
		      int N_AuxFeFun, Fefunction3D **AuxFeFuns, int *N_AuxFeFun_MultiIndex,
		      MultiIndex3D *AuxFeFun_MultiIndex,
		      int N_AuxFeFun_Values, double *AuxFeFun_Values,
                      int N_UserFunction,Functions3D **UserFunction, int N_FeConst,double *FeConst,
                      double  Quad_xi,double Quad_eta,
		      double Quad_zeta,double Quad_W, 
		      DiscreteFormRHSVolu *discreteformvolu, double *Rhs_Entries);

void AddSubRhs(RHST *Rhs, double *subRhs,int n_row, int *row);
void OutPutRHST(RHST *Rhs);
void FreeRHST(RHST *Rhs);
#endif

