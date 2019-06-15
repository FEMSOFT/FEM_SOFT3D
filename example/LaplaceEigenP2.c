/*
 * =====================================================================================
 *       Filename:  LinearEigen.c
 *    Description:    
 *
 *        Version:  1.0
 *        Created:  2014年1月27日 09时50分30秒
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
#include<math.h>
#include<time.h>
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
#include "EigenSolver.h"
#include "Multilevel.h"
#include "MG.h"
typedef void Fun_Test(int num);
double StiffMatrix(double Quad_x, double Quad_y,double Quad_z,int N_Test,double *Test_Value,int N_Anasatz,
		           double *Anasatz_Value,int N_AuxFeFun,double *AuxFeFun,int N_UserFunction,
		   double *UserFunction,int N_FeConst,double *FeConst);
double MassMatrix(double Quad_x, double Quad_y,double Quad_z,int N_Test,double *Test_Value,int N_Anasatz,
		           double *Anasatz_Value,int N_AuxFeFun,double *AuxFeFun,int N_UserFunction,
		   double *UserFunction,int N_FeConst,double *FeConst);
void ExactU(double x,double y, double z,int n, double *values);
void ExactU1(double x,double y, double z,int n, double *values);
void Userfunction(double x,double y,double z,int n,double *values);
void Coeffi(double x,double y,double z,int n,double *values);
//double RHSTerm(double x, double y,double z, int N_Test,double *Test_Value,int N_AuxFeFun,double *AuxFeFun);
void BoundCondition(int i, double t, BoundType *);
Functions3D * BoundFun(int ID_Bound);
void  BoundValue(double x, double y,double z, int dim, double *values);
void GradExactU(double x,double y,double z,  int n, double *values);
void GradExactU1(double x,double y,double z,  int n, double *values);
double  DiscreteFormBD(int ID_bd,double  x,double y,double z,double *normal, double *tangent, 
                       double local_direction, int n_test, double *test_values,int N_AuxFeFun,
		       double*AuxFeFun);
#define min(a, b) ((a>b)? b: a)
int 
main(int argc, char **argv)
{
  //网格部分
  MESH *mesh,*mesh_H;
  mesh=malloc(sizeof(MESH));
  mesh_H=malloc(sizeof(MESH));
  int i,j,iter;
 /* 
  ReadMesh(mesh,"../dat/new.dat");
  MeshConsist(mesh,3);
  for(i=0;i<mesh->Num_Volus_Global;++i)
    mesh->Volus[i]->Ancestor=i;  
 */
  ReadMeshBis(mesh,"../dat/cubeadaptive.dat");  
  MeshBisection(mesh,2);
  //MeshConsist(mesh,2);
  for(i=0;i<mesh->Num_Volus_Global;++i)
    mesh->Volus[i]->Ancestor=i;  
//------------------
  CopyMesh(mesh_H,mesh);
  /** 建立有限元空间 */  
  Fespace3D *fesp3D,*fesp3D_H;
  BoundCondFunct3D *BoundCond;
  BoundCond = BoundCondition;
  fesp3D = BUILDFEMSPACE3D(mesh,C_T_P2_3D,BoundCond);
  fesp3D_H = BUILDFEMSPACE3D(mesh_H,C_T_P2_3D,BoundCond);
  int dof_all = fesp3D->DOF_ALL;
  int Quad_Points3D,Quad_Points2D,Quad_Points1D; 
  Quad_Points3D = 4;
  Quad_Points2D = 4;
  Quad_Points1D = 4;
  /** 建立刚度矩阵*/
  int n_test_MultiIndexA,n_anasatz_MultiIndexA,n_test_MultiIndexB,n_anasatz_MultiIndexB;
  n_test_MultiIndexA = 3;
  n_test_MultiIndexB = 1;
  n_anasatz_MultiIndexA = 3;
  n_anasatz_MultiIndexB = 1;
  MultiIndex3D test_MultiIndexA[3] =  {D100,D010,D001};
  MultiIndex3D test_MultiIndexB[1] =  {D000};
  MultiIndex3D anasatz_MultiIndexA[3]={D100,D010,D001};
  MultiIndex3D anasatz_MultiIndexB[1]={D000};
  DiscreteFormMatrix *DiscreteA,*DiscreteB;
  DiscreteA = StiffMatrix;
  DiscreteB = MassMatrix;
  /**用到一个有限元函数 */
  Fefunction3D  *fefuncbase;
  double time1,time2;
  time1 = GetTime();
  int Num_Outloop ;//网格加密次数
  printf ( "please input the number of refine times\n" );
  scanf("%d",&Num_Outloop);
  int outloop;
  int nev =1;//求解的特征值个数
  double *eigenvalues;
  eigenvalues = malloc(sizeof(double)*(Num_Outloop+1)*nev);
  memset(eigenvalues,0.0,(Num_Outloop+1)*nev*sizeof(double)); 
  //InitialPointer(eigenvalues,Num_Outloop+1,0.0);
  double **eigenvectors;
  eigenvectors = malloc(nev*(Num_Outloop+1)*sizeof(double*)); 
  for(i=0; i<nev; i++)
  {
    eigenvectors[i] = calloc(dof_all, sizeof(double));   
  }
  DISCRETEFORM3D *discreteformA, *discreteformB, *discreteformA_H, *discreteformB_H;
  discreteformA = BuildDiscreteForm3D(fesp3D,fesp3D,n_test_MultiIndexA,test_MultiIndexA,
	                              n_anasatz_MultiIndexA,anasatz_MultiIndexA,DiscreteA,Quad_Points3D,
				      Quad_Points2D,Quad_Points1D);
  discreteformB = BuildDiscreteForm3D(fesp3D,fesp3D,n_test_MultiIndexB,test_MultiIndexB,
	                              n_anasatz_MultiIndexB,anasatz_MultiIndexB,DiscreteB,Quad_Points3D,
				      Quad_Points2D,Quad_Points1D);
  discreteformA_H = BuildDiscreteForm3D(fesp3D_H,fesp3D_H,n_test_MultiIndexA,test_MultiIndexA,
	                              n_anasatz_MultiIndexA,anasatz_MultiIndexA,DiscreteA,Quad_Points3D,
				      Quad_Points2D,Quad_Points1D);
  discreteformB_H = BuildDiscreteForm3D(fesp3D_H,fesp3D_H,n_test_MultiIndexB,test_MultiIndexB,
	                              n_anasatz_MultiIndexB,anasatz_MultiIndexB,DiscreteB,Quad_Points3D,
				      Quad_Points2D,Quad_Points1D);
  MULTILEVEL *multilevel;
  multilevel = BuildMultiLevel(mesh,discreteformA,discreteformB);
  MATRIX *stiff_matrix,*mass_matrix;
  stiff_matrix = multilevel->Stiff_Matrix[0];
  mass_matrix = multilevel->Mass_Matrix[0];
  BoundaryFunction3D *boundaryfunction;
  boundaryfunction = BoundFun;
  BoundaryTreatmentEigen(stiff_matrix, mass_matrix, boundaryfunction);
  //printf("@@@@ the first eigen solve -- stiff->nr = %d\n", stiff_matrix->N_Rows);
  //printf("write matrix->nr = %d\n", stiff_matrix->N_Rows);
  //WriteMatrix(stiff_matrix,"../dat/stiff_matrix-consistent1.dat");
  //WriteMatrix(mass_matrix, "../dat/mass_matrix-consistent1.dat");
  DNEigenSolver(stiff_matrix,mass_matrix,nev,-1,0.0,eigenvalues,eigenvectors+0,1);
  //for(i=0; i<nev; i++) printf("@@@@ nev %d = %18.15f\n", i, eigenvalues[i]);
  //exit(-1);

 //第二步，开始多水平校正迭代--------------------------------------------------------------------------------------------------------------------
  double *val, *solution,*rhs;
  MATRIX *Stiff_Matrix_H,*Mass_Matrix_H,*Prolong;
  double **eigenvector_H;
  eigenvector_H = malloc(sizeof(double*)*nev);
  for(i=0;i<nev;++i)
    eigenvector_H[i] = malloc(sizeof(double)*(multilevel->Stiff_Matrix[0]->N_Rows+1));
  Stiff_Matrix_H = ExtensionMatrix(multilevel->Stiff_Matrix[0],1);
  Mass_Matrix_H = ExtensionMatrix(multilevel->Mass_Matrix[0],1);
  Stiff_Matrix_H->DiscreteForm3D = discreteformA_H;
  Mass_Matrix_H->DiscreteForm3D = discreteformB_H;
  double *line_h,*line_H,*number;  
  number = malloc(sizeof(double));
  val=malloc(sizeof(double));
  rhs=malloc(sizeof(double));
  solution=malloc(sizeof(double));
  line_h=malloc(sizeof(double));
  int correctionstep;
  for(outloop=0;outloop<Num_Outloop;++outloop)
  {
    printf("=================================\n");
    printf ( "begin the %d-th outer loop!\n",outloop );
    AddLevelFEM(multilevel); 
    printf("@@@@ nrow = %d\n", multilevel->Stiff_Matrix[outloop+1]->N_Rows);
    stiff_matrix = multilevel->Stiff_Matrix[outloop+1];
    mass_matrix = multilevel->Mass_Matrix[outloop+1];
    /**建立延拓算子，把上一层网格解投影到当前层*/
    Prolong = multilevel->Prolongs[outloop];
    val=realloc(val,sizeof(double)*stiff_matrix->N_Rows);
    MatrixDotVec(Prolong, eigenvectors[outloop],val);
    rhs = realloc(rhs,sizeof(double)*stiff_matrix->N_Rows);
    solution = realloc(solution,stiff_matrix->N_Rows*sizeof(double));
    line_h = realloc(line_h,sizeof(double)*stiff_matrix->N_Rows);
    //for(correctionstep=0;correctionstep<1;++correctionstep)
    //{
      MatrixDotVec(mass_matrix,val,rhs);
      for(i=0;i<stiff_matrix->N_Rows;++i)
        rhs[i]*=eigenvalues[outloop];
      for(i=0;i<stiff_matrix->N_Rows;++i)
        solution[i] = val[i];
      BoundaryFunction3D *boundaryfunction;
      boundaryfunction = BoundFun;
      BoundaryTreatment(stiff_matrix,rhs,boundaryfunction);
      //求解
      printf("begin to solve");   
      //DIRECTSOLVER(stiff_matrix,rhs->Entries,solution);i
      //tol  = (1/(pow(2,outloop+4)))* (1/(pow(2,outloop+4)));
      CG(stiff_matrix,rhs,solution,1e-15,(int)(pow(2.0,1+2.01*(Num_Outloop-outloop))));
      //GMGSolver(multilevel,rhs,solution,outloop+1,3,1);
      printf ( "end solve end solve \n" );
      //----------2.1：多水平校正里面，开始自洽场迭代求粗空间上的解----------------------------------------------
       //Stiff_Matrix
      MatrixDotVec(stiff_matrix,solution,line_h);
      line_H = Restricth2H(multilevel,outloop+1,0,line_h);
      for(i=0;i<Stiff_Matrix_H->N_Rows-1;++i)
      {
        Stiff_Matrix_H->Entries[Stiff_Matrix_H->RowPtr[i+1]-1]=line_H[i];
      }
      for(i=0;i<Stiff_Matrix_H->N_Rows-1;++i)
        Stiff_Matrix_H->Entries[Stiff_Matrix_H->RowPtr[Stiff_Matrix_H->N_Rows-1]+i]=line_H[i];
      VecMatVec(solution,stiff_matrix,number);
      Stiff_Matrix_H->Entries[Stiff_Matrix_H->RowPtr[Stiff_Matrix_H->N_Rows]-1]=number[0];
      //Mass_Matrix
      MatrixDotVec(mass_matrix,solution,line_h);
      line_H = Restricth2H(multilevel,outloop+1,0,line_h);
      for(i=0;i<Mass_Matrix_H->N_Rows-1;++i)
      Mass_Matrix_H->Entries[Mass_Matrix_H->RowPtr[i+1]-1]=line_H[i];
      for(i=0;i<Mass_Matrix_H->N_Rows-1;++i)
      Mass_Matrix_H->Entries[Stiff_Matrix_H->RowPtr[Mass_Matrix_H->N_Rows-1]+i]=line_H[i];
      VecMatVec(solution,mass_matrix,number);
      Mass_Matrix_H->Entries[Mass_Matrix_H->RowPtr[Mass_Matrix_H->N_Rows]-1]=number[0];
      //特征值求解
      BoundaryTreatmentEigen(Stiff_Matrix_H, Mass_Matrix_H, boundaryfunction);
      // OutPutMatrix(Stiff_Matrix_H);
      // exit(0);
      // memset(eigenvector_H[0],0.0,(multilevel->Stiff_Matrix[0]->N_Rows+1)*sizeof(double));
      //OutPutMatrix(Mass_Matrix_H);
      DNEigenSolver(Stiff_Matrix_H,Mass_Matrix_H, nev, -1, 1.0, eigenvalues+outloop+1, eigenvector_H,1);
      printf("eigenvalues=======%f\n",eigenvalues[outloop+1]); 
      //把求得的解延拓到细空间
      eigenvectors[outloop+1]= ProlongH2h(multilevel,0,outloop+1,eigenvector_H[0]);
      for(i=0;i<stiff_matrix->N_Rows;++i)
      eigenvectors[outloop+1][i]+=eigenvector_H[0][Stiff_Matrix_H->N_Rows-1]*solution[i];
      for(i=0;i<stiff_matrix->N_Rows;++i)
      val[i]=eigenvectors[outloop+1][i];
    //}//correctionstep结束
  }//多水平校正结束
  printf("write matrix->nr = %d\n", multilevel->Stiff_Matrix[Num_Outloop]->N_Rows);
  //WriteMatrix(multilevel->Mass_Matrix[Num_Outloop], "../dat/mass_matrix-consistent.dat");
  //WriteMatrix(multilevel->Stiff_Matrix[Num_Outloop],"../dat/stiff_matrix-consistent.dat");

  time2 = GetTime();
  printf ( "time=%f\n",(time2-time1));
  printf("!!!!!!!!!!!!!!!!!!!!!!\n");
  for(i=0;i<Num_Outloop+1;++i)
    printf ( "the final eigenvalue=%f, %f\n",eigenvalues[i],eigenvalues[i]-3*PI*PI);
  printf("!!!!!!!!!!!!!!!!!!!!!!\n");
  //特征函数估计
  Fespace3D *fesp3D_new; 
  fesp3D_new = BUILDFEMSPACE3D(multilevel->mesh,C_T_P2_3D,BoundCond);

  FreeMultiLevel(multilevel);
  printf("!!!!!!!!!!!!!!!!!!!!!!\n");


  double *error0,*error1,*error01,*error11;
  error0 = malloc(sizeof(double));
  error01 = malloc(sizeof(double));
  error1 = malloc(sizeof(double));
  error11 = malloc(sizeof(double));
  MultiIndex3D ErrIndex[1] = { D000 };
  ErrorEstimateEigen3DVecExactSolution(eigenvectors[Num_Outloop],fesp3D_new,ExactU,1,ErrIndex,error0); 
  ErrorEstimateEigen3DVecExactSolution(eigenvectors[Num_Outloop],fesp3D_new,ExactU1,1,ErrIndex,error01); 
  MultiIndex3D GradErrIndex[3] = { D100,D010,D001 };
  ErrorEstimateEigen3DVecExactSolution(eigenvectors[Num_Outloop],fesp3D_new,GradExactU,3,GradErrIndex,error1); 
  ErrorEstimateEigen3DVecExactSolution(eigenvectors[Num_Outloop],fesp3D_new,GradExactU1,3,GradErrIndex,error11); 
  //for(i=0;i<Num_Outloop+1;++i)
  //printf ( "the final eigenvalue=%f\n",eigenvalues[i]);
  printf ( "-------------------------------------------------------------------------------------\n" );
  printf ( "num_elem=%d,%f,%f\n",mesh->Num_Volus_Global,min(error0[0],error01[0]),min(error1[0],error11[0])); 
  free(solution);
  free(line_h);
  free(rhs);


  /* 
  double *error0,*error1,*error01,*error11;
  error0 = malloc(sizeof(double)*Num_Outloop);
  error01 = malloc(sizeof(double)*Num_Outloop);
  error1 = malloc(sizeof(double)*Num_Outloop);
  error11 = malloc(sizeof(double)*Num_Outloop);
  MultiIndex3D ErrIndex[1] = { D000 };
  ErrorEstimateMultiEigen(multilevel->Prolongs,eigenvectors,fesp3D_new,Num_Outloop,1,ErrIndex,error0); 
  ErrorEstimateMultiEigen1(multilevel->Prolongs,eigenvectors,fesp3D_new,Num_Outloop,1,ErrIndex,error01); 
  MultiIndex3D GradErrIndex[3] = { D100,D010,D001 };
  ErrorEstimateMultiEigen(multilevel->Prolongs,eigenvectors,fesp3D_new,Num_Outloop,3,GradErrIndex,error1); 
  ErrorEstimateMultiEigen1(multilevel->Prolongs,eigenvectors,fesp3D_new,Num_Outloop,3,GradErrIndex,error11); 
  //for(i=0;i<Num_Outloop+1;++i)
  //printf ( "the final eigenvalue=%f\n",eigenvalues[i]);
  printf ( "-------------------------------------------------------------------------------------\n" );
  for(i=0;i<Num_Outloop;++i)
  printf ( "num_elem=%d,%f,%f\n",mesh->Num_Volus_Global,error0[i],error1[i]); 
  printf ( "-------------------------------------------------------------------------------------\n" );
   for(i=0;i<Num_Outloop;++i)
  printf ( "num_elem=%d,%f,%f\n",mesh->Num_Volus_Global,error01[i],error11[i]); 
  */
 


  //printf ( "time=ESH=%f,zhuzhuang%f,qiujie%f,zuihou%f\n",tt2-tt1,tt3-tt2,tt4-tt3,tt5-tt4);


  

  //内存释放
  //FreeFespace3D(fesp3D_new);//会和下一个重复
  for(i=0;i<Num_Outloop+1;++i)
     free(eigenvectors[i]);
  free(eigenvectors);
  free(eigenvalues);
  free(val);
  return 0;
}

/**============================================================================ */
/**                                    SubPrograms 
 * */
/**=============================================================================*/
//离散变分函数
//typedef double DiscreteFormMatrix(double,double,double,int,double *,int,double*,int, double *);
double StiffMatrix(double Quad_x,double Quad_y,double Quad_z,int N_Test,double *Test_Value,int N_Anasatz,
                   double *Anasatz_Value, int N_AuxFeFun,double *AuxFeFun,int N_UserFunction,
		   double *UserFunction,int N_FeConst,double *FeConst)
{
    return Test_Value[0]*Anasatz_Value[0]+Test_Value[1]*Anasatz_Value[1]+Test_Value[2]*Anasatz_Value[2];

}
double MassMatrix(double Quad_x,double Quad_y,double Quad_z,int N_Test,double *Test_Value,int N_Anasatz,
                   double *Anasatz_Value, int N_AuxFeFun,double *AuxFeFun,int N_UserFunction,
		   double *UserFunction,int N_FeConst,double *FeConst)
{
   return  Test_Value[0]*Anasatz_Value[0];
}

//有端项离散变分形式
//double RHSTerm(double x, double y, double z, int N_Test,double *Test_Value,int N_AuxFeFun,double *AuxFeFun)
double RHSTerm(double Quad_x,double Quad_y,double Quad_z,int N_Test,double *Test_Value,
                   int N_AuxFeFun,double *AuxFeFun,int N_FeConst,double *FeConst)
{ 
   return Test_Value[0]*FeConst[0]*AuxFeFun[0];
}
//真实解 
void ExactU(double x,double y,double z, int n, double *values)
{
  /**Example 1*/
  values[0] = sin(PI*x)*sin(PI*y)*sin(PI*z);
  /** Example 2*/
  //values[0] =8*cos(PI*x)*cos(PI*y)*cos(PI*z)/(PI*PI*PI);
}
void ExactU1(double x,double y,double z, int n, double *values)
{
  /**Example 1*/
  values[0] =- sin(PI*x)*sin(PI*y)*sin(PI*z);
  /** Example 2*/
  //values[0] =8*cos(PI*x)*cos(PI*y)*cos(PI*z)/(PI*PI*PI);
}

//真实解的导数
void GradExactU(double x,double y, double z, int n, double *values)
{
    /**Example 1*/
   values[0] = -PI*cos(PI*x)*sin(PI*y)*sin(PI*z);
   values[1] = -PI*sin(PI*x)*cos(PI*y)*sin(PI*z);
   values[2] = -PI*sin(PI*x)*sin(PI*y)*cos(PI*z);
  /** Example 2*/
  //printf("Exact value\n");
  //values[0] = -8*PI*sin(PI*x)*cos(PI*y)*cos(PI*z)/(PI*PI*PI);
  //values[1] = -8*PI*cos(PI*x)*sin(PI*y)*cos(PI*z)/(PI*PI*PI);
  //values[2] = -8*PI*cos(PI*x)*cos(PI*y)*sin(PI*z)/(PI*PI*PI);
}
void GradExactU1(double x,double y, double z, int n, double *values)
{
    /**Example 1*/
   values[0] = PI*cos(PI*x)*sin(PI*y)*sin(PI*z);
   values[1] = PI*sin(PI*x)*cos(PI*y)*sin(PI*z);
   values[2] = PI*sin(PI*x)*sin(PI*y)*cos(PI*z);
  /** Example 2*/
  //printf("Exact value\n");
  //values[0] = -8*PI*sin(PI*x)*cos(PI*y)*cos(PI*z)/(PI*PI*PI);
  //values[1] = -8*PI*cos(PI*x)*sin(PI*y)*cos(PI*z)/(PI*PI*PI);
  //values[2] = -8*PI*cos(PI*x)*cos(PI*y)*sin(PI*z)/(PI*PI*PI);
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
  /** Example 2*/
  value[0] = 0;
}


