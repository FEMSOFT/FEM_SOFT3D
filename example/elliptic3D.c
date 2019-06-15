/*
 * =====================================================================================
 *
 *       Filename:  elliptic3D.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2014年10月26日 20时17分41秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Xu Fei (Fly), xufei@lsec.cc.ac.cn
 *        Company:  
 *
 * =====================================================================================
 */
#include<stdio.h>	
#include<stdlib.h>
#include<string.h>
#include <math.h>
#include <time.h>
#include "Mesh3D.h"
#include "Fespace3D.h"
#include "Quadrature3D.h"
#include "Matrix.h"
#include "Quadrature2D.h"
#include "Quadrature1D.h"
#include "Enumerations.h"
#include "Constants.h"
#include "Rhst.h"
#include "DIRECTSOLVER.h"
#include "Fefunction3D.h"
#include "DiscreteForm3D.h"
//#include "MultiLevel.h"

typedef void Fun_Test(int num);
double StiffMatrix(double Quad_x, double Quad_y,double Quad_z,int N_Test,double *Test_Value,int N_Anasatz,
		           double *Anasatz_Value,int N_AuxFeFun,double *AuxFeFun);
void Diffusion_Matrix(double x,double y,double z, int N, double *values);           
void ExactU(double x,double y, double z,int n, double *values);
double RHSTerm(double x, double y,double z, int N_Test,double *Test_Value,int N_AuxFeFun,double *AuxFeFun);
// kind of boundary condition
void BoundCondition(int i, double t, BoundType *);
//void BoundValue(int BdComp, double Param, double x, double y, double *value);
//void BoundFun(int ID_Bound,Functions3D *boundfun);
Functions3D * BoundFun(int ID_Bound);
//void BoundValue(double x, double y, int n, double *value);
void  BoundValue(double x, double y,double z, int dim, double *values);
void GradExactU(double x,double y,double z,  int n, double *values);
double  DiscreteFormBD(int ID_bd,double  x,double y,double z,double *normal, double *tangent, 
                       double local_direction, int n_test, double *test_values,int N_AuxFeFun,
		       double*AuxFeFun);


int 
main(int argc, char **argv)
{
  //网格部分
  MESH *mesh;
  mesh=malloc(sizeof(MESH));
  int i,j;
  printf ( "please input the iteration number\n" );
  scanf("%d",&j);
  /*
  ReadMesh(mesh,"../dat/new.dat");
  for(i=0;i<j;++i)
  {
     AdjustMesh(mesh);
     MeshRefineConsistent(mesh);
  }
 */
// ReadMeshBis(mesh,"../dat/adaptive-one.dat");
  double time1,time2;
  ReadMeshBis(mesh,"../dat/cubeadaptive.dat");
  for(i=0;i<j;++i)
  {
     MeshRefineUniformBisectionNew1(mesh);
     MeshFullInfoL(mesh);
     MeshRenew(mesh);
  }
  MeshRefineUniformBisectionNew1(mesh);
  AdjustMeshForRHF(mesh);
  
    double t1, t2;

  //WriteMesh(mesh,"../dat/new.vtk");
  WriteMesh(mesh,"../dat/ellipticbisection.vtk");
  /** 建立有限元空间 */  
  Fespace3D *fesp3D;
  //include文件中用typedef声明的函数
  BoundCondFunct3D *BoundCond;
  BoundCond = BoundCondition;
  fesp3D = BUILDFEMSPACE3D(mesh,C_T_P1_3D,BoundCond);
  int Quad_Points3D; 
  Quad_Points3D = 15;
   int Quad_Points2D; 
  Quad_Points2D = 4;
   int Quad_Points1D; 
  Quad_Points1D = 4;
/** 建立刚度矩阵*/
  MATRIX *matrix;
  int n_test_MultiIndex,n_anasatz_MultiIndex;
  n_test_MultiIndex = 3;
  n_anasatz_MultiIndex = 3;
  MultiIndex3D test_MultiIndex[3] =  {D100,D010,D001};
  MultiIndex3D anasatz_MultiIndex[3]={D100,D010,D001};
  DiscreteFormMatrix *DiscreteA;
  DiscreteA = StiffMatrix;

  DISCRETEFORM3D *discreteform;  
  discreteform = BuildDiscreteForm3D(fesp3D,fesp3D,n_test_MultiIndex,test_MultiIndex,n_anasatz_MultiIndex,
	                             anasatz_MultiIndex,DiscreteA,
	                             Quad_Points3D, Quad_Points2D,Quad_Points1D);
  /** Assemble the matrix */  
  matrix = BuildMatrix(discreteform);
  printf("Assembling matrix\n");
  t1=GetTime();
  AssembleMatrix(matrix);
  t2=GetTime();
  printf("Assembling time:%f\n",t2-t1);
  //OutPutMatrix(matrix);
  printf ( "matirx information:%d,%d,%d\n",matrix->N_Rows,matrix->N_Columns,matrix->N_Entries); 
  WriteMatrixPS( matrix, "./matAA.ps");
  /**建立右端项 */
  RHST *rhs;
  int n_rhs_MultiIndex =1;
  MultiIndex3D rhs_MultiIndex[1] = {D000};
  DiscreteFormRHSVolu *discreteformrhs;
  discreteformrhs = RHSTerm;
  DiscreteFormRHSFace *discreteformrhsface;
  discreteformrhsface = DiscreteFormBD;
  
  printf("Assemblingm rhs\n");
  t1 = GetTime();
  rhs = BuildRHSTBD(fesp3D,n_rhs_MultiIndex,rhs_MultiIndex,discreteformrhs,discreteformrhsface,
		    Quad_Points3D,Quad_Points2D,Quad_Points1D);
  AssembleRHST(rhs);
  t2=GetTime();
  printf("Assembling rhs time: %f\n",t2-t1);
  /**求解线性方程组 */
  int dof_all;
  dof_all = rhs->DOF_ALL;
  double *sol;
  sol = malloc(dof_all*sizeof(double));
  memset(sol,0.0,dof_all*sizeof(double));
  //treat the Dirichlet boundary condition
  BoundaryFunction3D *boundaryfunction;
  boundaryfunction = BoundFun;
  //处理边界条件
  BoundaryTreatment(matrix, rhs->Entries, boundaryfunction);
  //求解

  t1 = GetTime();
  DIRECTSOLVER(matrix,rhs->Entries,sol);
  t2 = GetTime();
  printf("Solving time: %f\n",t2-t1);
  //CG(matrix,rhs->Entries,sol,1e-6);
  //Error estimate 
  /**建立有限元函数空间 */
  Fefunction3D *Fefunct;
  Fefunct = BuildFefunction3D(fesp3D,sol);

  Functions3D *Exact;
  Exact = ExactU;
  MultiIndex3D ErrIndex[1] = { D000 };
  double Error;
  //compute the error in L2 norm
  Error = ErrorEstimate3D(Fefunct,Exact, 1,ErrIndex,Quad_Points3D);
  //Error = ErrorEstimate3D(Fefunct,Exact, 1,ErrIndex,15);
  Functions3D *GradExact;
  GradExact = GradExactU;
  MultiIndex3D GradErrIndex[3] = { D100,D010,D001 };
  double GradError;
  //Compute the error in the seminorm H1
  GradError = ErrorEstimate3D(Fefunct,GradExactU, 3,GradErrIndex,Quad_Points3D);
  int Num_Elems;
  Num_Elems = Fefunct->Fespace->Mesh->Num_Volus_Global;
  printf("Solving time: %f\n",t2-t1);
  printf("Num_Elems = %d, L2-Norm Error = %18.15e\n",Num_Elems,Error);
  printf("Num_Elems = %d, H1-Norm Error = %18.15e\n",Num_Elems,GradError);

  double *sol_int;
  sol_int = malloc(dof_all*sizeof(double));
  memset(sol_int,0.0,dof_all*sizeof(double)); 
  Fefunction3D *Fefunct_int;
  Fefunct_int = BuildFefunction3D(fesp3D,sol_int);
  //do the interpolation process
  Interpolation3D(Fefunct_int,Exact); 
  Error = ErrorEstimate3D(Fefunct_int,Exact, 1,ErrIndex,Quad_Points3D);
  //Compute the error in the seminorm H1
  GradError = ErrorEstimate3D(Fefunct_int,GradExactU, 3,GradErrIndex,Quad_Points3D);
  printf("Interpolation: Num_Elems = %d, L2-Norm Error = %18.15e\n",Num_Elems,Error);
  printf("Interpolation: Num_Elems = %d, H1-Norm Error = %18.15e\n",Num_Elems,GradError);
 
  exit(0);
//Free the memory
  FreeMatrix(matrix);
  FreeRHST(rhs);
  FreeFespace3D(fesp3D);
  FreeDiscreteForm3D(discreteform);
  free(sol);
  free(sol_int);
  FreeFefunction3D(Fefunct);
  FreeFefunction3D(Fefunct_int);
  FreeMesh(mesh);
  return 0;
}

/**============================================================================ */
/**                                    SubPrograms 
 * */
/**=============================================================================*/
//离散变分函数
//typedef double DiscreteFormMatrix(double,double,double,int,double *,int,double*,int, double *);
double StiffMatrix(double Quad_x,double Quad_y,double Quad_z,int N_Test,double *Test_Value,int N_Anasatz,
                   double *Anasatz_Value, int N_AuxFeFun,double *AuxFeFun)
{
  //return  Test_Value[0]*Anasatz_Value[0]+Test_Value[1]*Anasatz_Value[1]
  //       +Test_Value[2]*Anasatz_Value[2]+Test_Value[3]*Anasatz_Value[3];
  return   Test_Value[0]*Anasatz_Value[0]+Test_Value[1]*Anasatz_Value[1]
          +Test_Value[2]*Anasatz_Value[2];
}
//有端项离散变分形式
double RHSTerm(double x, double y, double z, int N_Test,double *Test_Value,int N_AuxFeFun,double *AuxFeFun)
{
  //return Test_Value[0]*(3.0*PI*PI+1.0)*cos(PI*x)*cos(PI*y)*cos(PI*z);
  return Test_Value[0]*3.0*PI*PI*sin(PI*x)*sin(PI*y)*sin(PI*z);
}
//真实解 
void ExactU(double x,double y,double z, int n, double *values)
{
  /**Example 1*/
  values[0] = sin(PI*x)*sin(PI*y)*sin(PI*z);
  /** Example 2*/
  //values[0] = cos(PI*x)*cos(PI*y)*cos(PI*z);
}
//真实解的导数
void GradExactU(double x,double y, double z, int n, double *values)
{
    /**Example 1*/
   values[0] = PI*cos(PI*x)*sin(PI*y)*sin(PI*z);
   values[1] = PI*sin(PI*x)*cos(PI*y)*sin(PI*z);
   values[2] = PI*sin(PI*x)*sin(PI*y)*cos(PI*z);
  /** Example 2*/
  //printf("Exact value\n");
  //values[0] = -PI*sin(PI*x)*cos(PI*y)*cos(PI*z);
  //values[1] = -PI*cos(PI*x)*sin(PI*y)*cos(PI*z);
  //values[2] = -PI*cos(PI*x)*cos(PI*y)*sin(PI*z);
}
//边界条件类型函数
void BoundCondition(int i, double t, BoundType* cond)
{
  if( i>0)
  {
    cond[0] = DIRICHLETT; 
  }
  else
  {
    cond[0] = NOTBOUNDARY;
  }
}

//定义NEUMANNN 和ROBIN条件的法向导函数
double  DiscreteFormBD(int ID_bd,double  x,double y,double z,double *normal, double *tangent, 
                       double local_direction, int N_Test, double *Test_Values, int N_AuxFeFun,
		       double*AuxFeFun)
{
  switch(ID_bd)
  {    
//     case 1:
//       return -PI*sin(PI*x)*Test_Values[0];
//       break;
//     case 2:
//       return -PI*sin(PI*y)*Test_Values[0];
//       break;
//     case 3:
//       return -PI*sin(PI*x)*Test_Values[0];
//       break;
//     case 4:
//       return -PI*sin(PI*y)*Test_Values[0];
//       break;
    default:
      return 0.0;  
      break;
  }
}
//边界条件：返回的是边界条件的函数
// value of boundary condition
Functions3D * BoundFun(int ID_Bound) //,Functions3D *boundfun)
{
  return  BoundValue;
}
//边界条件函数
void BoundValue(double x, double y,double z, int n, double *value)
{
  /** Example 1 */
  value[0] = 0;
  /** Example 2*/
  //value[0] = cos(PI*x)*cos(PI*y)*cos(PI*z);
}
