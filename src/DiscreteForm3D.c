/*
 * =====================================================================================
 *
 *       Filename:  DiscreteForm3D.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2014年10月26日 22时20分48秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Xu Fei (Fly), xufei@lsec.cc.ac.cn
 *        Company:  
 *
 * =====================================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "DiscreteForm3D.h"

//下面是一些操作函数：最重要的是建立离散变分形式: 本质上这里什么也没有做，只是做了一些登记
/** 建立离散变分形式 （不需要解析函数和有限元函数）*/
DISCRETEFORM3D *BuildDiscreteForm3D(Fespace3D *anasatz_space,Fespace3D *test_space,int n_test_MultiIndex,
			MultiIndex3D *test_MultiIndex,int n_anasatz_MultiIndex,MultiIndex3D *anasatz_MultiIndex,
			DiscreteFormMatrix *discreteformvolu,int Quad_Points3D,int Quad_Points2D,
			int Quad_Points1D)
 {
  DISCRETEFORM3D *discreteform;
  discreteform = malloc(sizeof(DISCRETEFORM3D));
  discreteform->Anasatz_Space = anasatz_space;
  discreteform->Test_Space = test_space;
  discreteform->N_Test_MultiIndex = n_test_MultiIndex;
  discreteform->Test_MultiIndex = test_MultiIndex;
  discreteform->N_Anasatz_MultiIndex = n_anasatz_MultiIndex; 
  discreteform->Anasatz_MultiIndex = anasatz_MultiIndex; 

  discreteform->N_AuxFeFun = 0;
  discreteform->AuxFeFun = NULL; 
  discreteform->N_AuxFeFun_MultiIndex = NULL;
  discreteform->AuxFeFun_MultiIndex = NULL;
  discreteform->AuxFeFun_Values = NULL;
  discreteform->N_AuxFeFun_Values = 0;
 
  discreteform->N_UserFunction=0;
  discreteform->UserFunction=NULL;

  discreteform->N_FeConst=0;
  discreteform->FeConst=NULL;

  discreteform->DiscreteFormVolu = discreteformvolu;
  discreteform->DiscreteFormFace = NULL;
  discreteform->DiscreteFormLine = NULL;
  //discreteform->DiscreteFormBD = NULL;    

  discreteform->Quad3D = BuildQuad2Mesh3D(anasatz_space->Mesh, Quad_Points3D);
  discreteform->Quad2D = BuildQuad2Mesh2D(anasatz_space->Mesh, Quad_Points2D);
  discreteform->Quad1D = BuildQuad1D(Quad_Points1D);

  return discreteform;
}

/** 建立矩阵对象：（需要解析函数和不需要有限元函数）*/
/*DISCRETEFORM3D *Builddiscreteform3DAuxFeFun(Fespace3D *anasatz_space,Fespace3D *test_space,
					  int n_test_MultiIndex,MultiIndex3D *test_MultiIndex,
					  int n_anasatz_MultiIndex,MultiIndex3D *anasatz_MultiIndex,
					  int n_auxfefun, Fefunction3D **auxfefun,
					  int *n_auxFefun_MultiIndex,MultiIndex3D *auxFefun_MultiIndex,
					  DiscreteFormMatrixFace *discreteformface,int Quad_Points3D,
					  int Quad_Points2D, int Quad_Points1D)
{
  DISCRETEFORM3D *discreteform;
  discreteform = malloc(sizeof(DISCRETEFORM3D));
  discreteform->Anasatz_Space = anasatz_space;
  discreteform->Test_Space = test_space;
  discreteform->N_Test_MultiIndex = n_test_MultiIndex;
  discreteform->Test_MultiIndex = test_MultiIndex;
  discreteform->N_Anasatz_MultiIndex = n_anasatz_MultiIndex; 
  discreteform->Anasatz_MultiIndex = anasatz_MultiIndex; 
  
  discreteform->N_AuxFeFun = n_auxfefun;
  discreteform->AuxFeFun = auxfefun; 
  discreteform->N_AuxFeFun_MultiIndex = n_auxFefun_MultiIndex;
  discreteform->AuxFeFun_MultiIndex = auxFefun_MultiIndex;
  int i, N_AuxFeFun_Values;
  N_AuxFeFun_Values = 0; 
  for(i=0;i<n_auxfefun;i++)
    N_AuxFeFun_Values += n_auxFefun_MultiIndex[i]*auxfefun[i]->Fespace->Base->Value_Dim;        
  discreteform->AuxFeFun_Values = malloc(N_AuxFeFun_Values*sizeof(double));
  discreteform->N_AuxFeFun_Values = N_AuxFeFun_Values;
  
  discreteform->DiscreteFormFace = discreteformface; 
  discreteform->DiscreteFormLine = NULL;
  //discreteform->DiscreteFormBD = NULL;
  discreteform->DiscreteFormVert = NULL;
  
  discreteform->Quad3D = BuildQuad2Mesh3D(anasatz_space->Mesh, Quad_Points3D);
  discreteform->Quad2D = BuildQuad2Mesh2D(anasatz_space->Mesh, Quad_Points2D);
  discreteform->Quad1D = BuildQuad1D(Quad_Points1D);
  return discreteform; 

 }
 */


/** 建立矩阵对象 （需要有限元函数和解析函数）*/
DISCRETEFORM3D *BuildDiscreteFormAll3D(Fespace3D *anasatz_space,Fespace3D *test_space,int n_test_MultiIndex,
			  MultiIndex3D *test_MultiIndex,int n_anasatz_MultiIndex,
			  MultiIndex3D *anasatz_MultiIndex,int n_auxfefun,Fefunction3D **auxfefun,
			  int *N_AuxFefun_MultiIndex,MultiIndex3D *AuxFefun_MultiIndex,
			  int N_UserFunction,Functions3D **userfunction,DiscreteFormMatrix *discreteformvolu,
			  int Quad_Points3D,int Quad_Points2D,int Quad_Points1D)
{
  DISCRETEFORM3D *discreteform; 
  discreteform = malloc(sizeof(DISCRETEFORM3D));
  discreteform->Anasatz_Space = anasatz_space;
  discreteform->Test_Space = test_space;
  discreteform->N_Test_MultiIndex = n_test_MultiIndex;
  discreteform->Test_MultiIndex = test_MultiIndex;
  discreteform->N_Anasatz_MultiIndex = n_anasatz_MultiIndex; 
  discreteform->Anasatz_MultiIndex = anasatz_MultiIndex; 
   
  discreteform->N_AuxFeFun = n_auxfefun;//fefunction的个数
  discreteform->AuxFeFun = auxfefun;//存fefunction 
  discreteform->N_AuxFeFun_MultiIndex = N_AuxFefun_MultiIndex;//??????????????
  discreteform->AuxFeFun_MultiIndex = AuxFefun_MultiIndex;//fefunction的求导方式
  int i, N_AuxFeFun_Values;
  N_AuxFeFun_Values = 0; 
  for(i=0;i<n_auxfefun;i++)
    N_AuxFeFun_Values += N_AuxFefun_MultiIndex[i];        
  discreteform->AuxFeFun_Values = malloc(N_AuxFeFun_Values*sizeof(double));
  discreteform->N_AuxFeFun_Values = N_AuxFeFun_Values;
  
  discreteform->N_UserFunction=N_UserFunction;
  discreteform->UserFunction=userfunction;

  discreteform->N_FeConst=0;
  discreteform->FeConst=NULL;



  discreteform->DiscreteFormVolu = discreteformvolu; 
  //discreteform->DiscreteFormLine = discreteformline;
  //discreteform->DiscreteFormBD = discreteformbd;
  //discreteform->DiscreteFormVert = discreteformvert;
  
  discreteform->Quad3D = BuildQuad2Mesh3D(anasatz_space->Mesh, Quad_Points3D);
  discreteform->Quad2D = BuildQuad2Mesh2D(anasatz_space->Mesh, Quad_Points2D);
  discreteform->Quad1D = BuildQuad1D(Quad_Points1D);
  return discreteform;
}



DISCRETEFORM3D *BuildDiscreteFormFeConst3D(Fespace3D *anasatz_space,Fespace3D *test_space,int n_test_MultiIndex,
			  MultiIndex3D *test_MultiIndex,int n_anasatz_MultiIndex,
			  MultiIndex3D *anasatz_MultiIndex,int n_auxfefun,Fefunction3D **auxfefun,
			  int *N_AuxFefun_MultiIndex,MultiIndex3D *AuxFefun_MultiIndex,
			  int N_UserFunction,Functions3D **userfunction,DiscreteFormMatrix *discreteformvolu,
			  int N_FeConst,double *FeConst,int Quad_Points3D,int Quad_Points2D,int Quad_Points1D)
{
  DISCRETEFORM3D *discreteform; 
  discreteform = malloc(sizeof(DISCRETEFORM3D));
  discreteform->Anasatz_Space = anasatz_space;
  discreteform->Test_Space = test_space;
  discreteform->N_Test_MultiIndex = n_test_MultiIndex;
  discreteform->Test_MultiIndex = test_MultiIndex;
  discreteform->N_Anasatz_MultiIndex = n_anasatz_MultiIndex; 
  discreteform->Anasatz_MultiIndex = anasatz_MultiIndex; 
   
  discreteform->N_AuxFeFun = n_auxfefun;//fefunction的个数
  discreteform->AuxFeFun = auxfefun;//存fefunction 
  discreteform->N_AuxFeFun_MultiIndex = N_AuxFefun_MultiIndex;//??????????????
  discreteform->AuxFeFun_MultiIndex = AuxFefun_MultiIndex;//fefunction的求导方式
  int i, N_AuxFeFun_Values;
  N_AuxFeFun_Values = 0; 
  for(i=0;i<n_auxfefun;i++)
    N_AuxFeFun_Values += N_AuxFefun_MultiIndex[i];        
  discreteform->AuxFeFun_Values = malloc(N_AuxFeFun_Values*sizeof(double));
  discreteform->N_AuxFeFun_Values = N_AuxFeFun_Values;
  
  discreteform->N_UserFunction=N_UserFunction;
  discreteform->UserFunction=userfunction;
  
  discreteform->N_FeConst=N_FeConst;
  discreteform->FeConst=FeConst;


  discreteform->DiscreteFormVolu = discreteformvolu; 
  //discreteform->DiscreteFormLine = discreteformline;
  //discreteform->DiscreteFormBD = discreteformbd;
  //discreteform->DiscreteFormVert = discreteformvert;
  
  discreteform->Quad3D = BuildQuad2Mesh3D(anasatz_space->Mesh, Quad_Points3D);
  discreteform->Quad2D = BuildQuad2Mesh2D(anasatz_space->Mesh, Quad_Points2D);
  discreteform->Quad1D = BuildQuad1D(Quad_Points1D);
  return discreteform;
}





void FreeDiscreteForm3D(DISCRETEFORM3D *discreteform)
{
  //printf("delete the DiscreteForm3D!\n");
  if(discreteform->Quad3D)
    FreeQuad3D(discreteform->Quad3D);

  if(discreteform->Quad2D)
    FreeQuad2D(discreteform->Quad2D);

  if(discreteform->Quad1D)
    FreeQuad1D(discreteform->Quad1D);

  if(discreteform->AuxFeFun_Values)
    free(discreteform->AuxFeFun_Values);
   
 free(discreteform);

 
}
