/*
 * =====================================================================================
 *
 *       Filename:  Rhst.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2014年10月31日 12时34分24秒
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
 #include "Rhst.h"
 #include "Base3D.h"

RHST *BuildRHST(Fespace3D *test_fesp,int n_test_MultiIndex,
                 MultiIndex3D *test_MultiIndex,
	      	 DiscreteFormRHSVolu *discreteformvolu,
	         int Quad_Points3D)
{
    RHST *rhst;
    rhst = malloc(1*sizeof(RHST));
    rhst->Test_Space = test_fesp;
    rhst->N_Test_MultiIndex = n_test_MultiIndex;
    rhst->Test_MultiIndex = test_MultiIndex; 

    rhst->N_AuxFeFun = 0;
    rhst->AuxFeFun = NULL;
    rhst->N_AuxFeFun_MultiIndex = NULL;
    rhst->AuxFeFun_MultiIndex = NULL;
    rhst->N_AuxFeFun_Values = 0;
    rhst->AuxFeFun_Values = NULL;
    
    rhst->N_UserFunction=0;
    rhst->UserFunction=NULL;
    rhst->N_FeConst=0;
    rhst->FeConst=NULL;

    rhst->DiscreteFormVolu = discreteformvolu;
    rhst->DiscreteFormFace = NULL;
    rhst->DiscreteFormLine = NULL;
    //rhst->DiscreteFormBD = NULL;

    rhst->Quad3D = BuildQuad2Mesh3D(test_fesp->Mesh, Quad_Points3D); 
    rhst->Quad2D = NULL;
    rhst->Quad1D = NULL;
    rhst->DOF_ALL = test_fesp->DOF_ALL;
    rhst->Entries = malloc(rhst->DOF_ALL*sizeof(double));
    memset(rhst->Entries,0,rhst->DOF_ALL*sizeof(double));
    
    return rhst;
}


RHST *BuildRHSTBD(Fespace3D *test_fesp,int n_test_MultiIndex,
                   MultiIndex3D *test_MultiIndex,                        
	           DiscreteFormRHSVolu *discreteformvolu,
		   DiscreteFormRHSFace *discreteformface,
		   int Quad_Points3D,int Quad_Points2D, int Quad_Points1D)
{
    RHST *rhst;
    rhst = malloc(1*sizeof(RHST));
    rhst->Test_Space = test_fesp;
    rhst->N_Test_MultiIndex = n_test_MultiIndex;
    rhst->Test_MultiIndex = test_MultiIndex; 

    rhst->N_AuxFeFun = 0;
    rhst->AuxFeFun = NULL;
    rhst->N_AuxFeFun_MultiIndex = NULL;
    rhst->AuxFeFun_MultiIndex = NULL;
    rhst->N_AuxFeFun_Values = 0;
    rhst->AuxFeFun_Values = NULL;
        
    rhst->N_UserFunction=0;
    rhst->UserFunction=NULL;
    rhst->N_FeConst=0;
    rhst->FeConst=NULL;
 
    
    rhst->DiscreteFormVolu = discreteformvolu;
    rhst->DiscreteFormFace = discreteformface;
    rhst->DiscreteFormLine = NULL;

    rhst->Quad3D = BuildQuad2Mesh3D(test_fesp->Mesh, Quad_Points3D); 
    rhst->Quad2D = BuildQuad2Mesh2D(test_fesp->Mesh, Quad_Points2D); 
    rhst->Quad1D = BuildQuad1D(Quad_Points1D);
    rhst->DOF_ALL = test_fesp->DOF_ALL;
    rhst->Entries = malloc(rhst->DOF_ALL*sizeof(double));
    memset(rhst->Entries,0.0,rhst->DOF_ALL*sizeof(double));
    return rhst;
}

RHST *BuildRHSTAuxFeFunct(Fespace3D *test_fesp,int n_test_MultiIndex,
                          MultiIndex3D *test_MultiIndex,
                          int n_auxFefunct,Fefunction3D **auxFefunct,
			  int *n_auxFefun_MultiIndex,
                          MultiIndex3D *auxFefun_MultiIndex, 
			  DiscreteFormRHSVolu *discreteformvolu, int Quad_Points3D)
{
    RHST *rhst;
    rhst = malloc(1*sizeof(RHST));
    rhst->Test_Space = test_fesp;
    rhst->N_Test_MultiIndex = n_test_MultiIndex;
    rhst->Test_MultiIndex = test_MultiIndex;   
   
    rhst->N_AuxFeFun = n_auxFefunct;
    rhst->AuxFeFun = auxFefunct;
    rhst->N_AuxFeFun_MultiIndex = n_auxFefun_MultiIndex;
    rhst->AuxFeFun_MultiIndex = auxFefun_MultiIndex;
    int i, N_AuxFeFun_Values;
    N_AuxFeFun_Values = 0; 
    for(i=0;i<n_auxFefunct;i++)
      N_AuxFeFun_Values += n_auxFefun_MultiIndex[i]*test_fesp->Base->Value_Dim;        
    rhst->AuxFeFun_Values = malloc(N_AuxFeFun_Values*sizeof(double));
    memset(rhst->AuxFeFun_Values,0,N_AuxFeFun_Values*sizeof(double));
    
    rhst->N_AuxFeFun_Values = N_AuxFeFun_Values;
    
    rhst->N_UserFunction=0;
    rhst->UserFunction=NULL;
    rhst->N_FeConst=0;
    rhst->FeConst=NULL;

    rhst->DiscreteFormVolu = discreteformvolu;
    rhst->DiscreteFormFace = NULL;
    rhst->DiscreteFormLine = NULL;
    rhst->DOF_ALL = test_fesp->DOF_ALL;
    rhst->Quad3D = BuildQuad2Mesh3D(test_fesp->Mesh, Quad_Points3D); 
    rhst->Quad2D = NULL; 
    rhst->Quad1D = NULL;
    rhst->Entries = malloc(rhst->DOF_ALL*sizeof(double));
    memset(rhst->Entries,0,rhst->DOF_ALL*sizeof(double));
    
    return rhst;
   
}
RHST *BuildRHSTAll(Fespace3D *test_fesp,int n_test_MultiIndex,
                   MultiIndex3D *test_MultiIndex,                        
                   int n_auxFefunct,Fefunction3D **auxFefunct,
		   int *n_auxFefun_MultiIndex,
                   MultiIndex3D *auxFefun_MultiIndex,DiscreteFormRHSVolu *discreteformvolu, 
		   DiscreteFormRHSFace *discreteformface, int Quad_Points3D,int Quad_Points2D, 
		   int Quad_Points1D)
{
    RHST *rhst;
    rhst = malloc(1*sizeof(RHST));
    rhst->Test_Space = test_fesp;
    rhst->N_Test_MultiIndex = n_test_MultiIndex;
    rhst->Test_MultiIndex = test_MultiIndex;
     
    rhst->N_AuxFeFun = n_auxFefunct;
    rhst->AuxFeFun = auxFefunct;
    rhst->N_AuxFeFun_MultiIndex = n_auxFefun_MultiIndex;
    rhst->AuxFeFun_MultiIndex = auxFefun_MultiIndex;
    int i, N_AuxFeFun_Values;
    N_AuxFeFun_Values = 0; 
    for(i=0;i<n_auxFefunct;i++)
      N_AuxFeFun_Values += n_auxFefun_MultiIndex[i]*test_fesp->Base->Value_Dim;        
    rhst->AuxFeFun_Values = malloc(N_AuxFeFun_Values*sizeof(double));
    memset(rhst->AuxFeFun_Values,0,N_AuxFeFun_Values*sizeof(double));
    rhst->N_AuxFeFun_Values = N_AuxFeFun_Values;
     
    rhst->N_UserFunction=0;
    rhst->UserFunction=NULL;
    rhst->N_FeConst=0;
    rhst->FeConst=NULL;

    rhst->DiscreteFormVolu = discreteformvolu;
    rhst->DiscreteFormFace = discreteformface;
    //rhst->DiscreteFormBD = discreteformbd;
    //rhst->DiscreteFormLine = discreteformline;
    rhst->DOF_ALL = test_fesp->DOF_ALL;
    
    rhst->Quad3D = BuildQuad2Mesh3D(test_fesp->Mesh, Quad_Points3D); 
    rhst->Quad2D = BuildQuad2Mesh2D(test_fesp->Mesh, Quad_Points2D); 
    rhst->Quad1D = BuildQuad1D(Quad_Points1D);
    rhst->Entries = malloc(rhst->DOF_ALL*sizeof(double));
    memset(rhst->Entries,0.0,rhst->DOF_ALL*sizeof(double));
    
    return rhst;
}

RHST *BuildRHSTFeConst(Fespace3D *test_fesp,int n_test_MultiIndex,
                      MultiIndex3D *test_MultiIndex,                        
                      int n_auxFefunct,Fefunction3D **auxFefunct,
		      int *n_auxFefun_MultiIndex,
                      MultiIndex3D *auxFefun_MultiIndex, 
                      int N_Feconst, double *Feconst,
		      DiscreteFormRHSVolu *discreteformvolu, int Quad_Points3D)
{
    //printf("1111111111111111\n");
    RHST *rhst;
    rhst = malloc(1*sizeof(RHST));
    rhst->Test_Space = test_fesp;
    rhst->N_Test_MultiIndex = n_test_MultiIndex;
    rhst->Test_MultiIndex = test_MultiIndex;
     
    rhst->N_AuxFeFun = n_auxFefunct;
    rhst->AuxFeFun = auxFefunct;
    rhst->N_AuxFeFun_MultiIndex = n_auxFefun_MultiIndex;
    rhst->AuxFeFun_MultiIndex = auxFefun_MultiIndex;
    int i, N_AuxFeFun_Values;
    N_AuxFeFun_Values = 0; 
    for(i=0;i<n_auxFefunct;i++)
      N_AuxFeFun_Values += n_auxFefun_MultiIndex[i]*test_fesp->Base->Value_Dim;        
    rhst->AuxFeFun_Values = malloc(N_AuxFeFun_Values*sizeof(double));
    memset(rhst->AuxFeFun_Values,0,N_AuxFeFun_Values*sizeof(double));
    rhst->N_AuxFeFun_Values = N_AuxFeFun_Values;
    //printf("2222222222222222222\n");
    rhst->N_UserFunction=0;
    rhst->UserFunction=NULL;
    rhst->N_FeConst=N_Feconst;
    rhst->FeConst=Feconst;
    //printf ( "inner==%d,%d\n", rhst->N_UserFunction,N_UserFunction);
    rhst->DiscreteFormVolu = discreteformvolu;
    //rhst->DiscreteFormFace = discreteformface;
    //rhst->DiscreteFormBD = discreteformbd;
    //rhst->DiscreteFormLine = discreteformline;
    rhst->DOF_ALL = test_fesp->DOF_ALL;
    //printf("33333333333333333333\n");
    rhst->Quad3D = BuildQuad2Mesh3D(test_fesp->Mesh, Quad_Points3D); 
    rhst->Quad2D = NULL; 
    rhst->Quad1D = NULL;
    //printf("333333333334444444444\n");
    rhst->Entries = malloc((rhst->DOF_ALL)*sizeof(double));
    //printf("55555555555555555\n");
    memset(rhst->Entries,0.0,rhst->DOF_ALL*sizeof(double));
    //printf("4444444444444444444\n");
    return rhst;
}



//为了牛顿迭代设计的
RHST *BuildRHSTAllMix(Fespace3D *test_fesp,int n_test_MultiIndex,
                      MultiIndex3D *test_MultiIndex,                        
                      int n_auxFefunct,Fefunction3D **auxFefunct,
		      int *n_auxFefun_MultiIndex,
                      MultiIndex3D *auxFefun_MultiIndex, 
                      int N_Userfunction, Functions3D  **Userfunction, int N_Feconst, double *Feconst,
		      DiscreteFormRHSVolu *discreteformvolu, int Quad_Points3D, int Quad_Points2D, 
		      int Quad_Points1D)
{
    //printf("1111111111111111\n");
    RHST *rhst;
    rhst = malloc(1*sizeof(RHST));
    rhst->Test_Space = test_fesp;
    rhst->N_Test_MultiIndex = n_test_MultiIndex;
    rhst->Test_MultiIndex = test_MultiIndex;
     
    rhst->N_AuxFeFun = n_auxFefunct;
    rhst->AuxFeFun = auxFefunct;
    rhst->N_AuxFeFun_MultiIndex = n_auxFefun_MultiIndex;
    rhst->AuxFeFun_MultiIndex = auxFefun_MultiIndex;
    int i, N_AuxFeFun_Values;
    N_AuxFeFun_Values = 0; 
    for(i=0;i<n_auxFefunct;i++)
      N_AuxFeFun_Values += n_auxFefun_MultiIndex[i]*test_fesp->Base->Value_Dim;        
    rhst->AuxFeFun_Values = malloc(N_AuxFeFun_Values*sizeof(double));
    memset(rhst->AuxFeFun_Values,0,N_AuxFeFun_Values*sizeof(double));
    rhst->N_AuxFeFun_Values = N_AuxFeFun_Values;
    //printf("2222222222222222222\n");
    rhst->N_UserFunction=N_Userfunction;
    rhst->UserFunction=Userfunction;
    rhst->N_FeConst=N_Feconst;
    rhst->FeConst=Feconst;
    //printf ( "inner==%d,%d\n", rhst->N_UserFunction,N_UserFunction);
    rhst->DiscreteFormVolu = discreteformvolu;
    //rhst->DiscreteFormFace = discreteformface;
    //rhst->DiscreteFormBD = discreteformbd;
    //rhst->DiscreteFormLine = discreteformline;
    rhst->DOF_ALL = test_fesp->DOF_ALL+1;
    //printf("33333333333333333333\n");
    rhst->Quad3D = BuildQuad2Mesh3D(test_fesp->Mesh, Quad_Points3D); 
    rhst->Quad2D = BuildQuad2Mesh2D(test_fesp->Mesh, Quad_Points2D); 
    rhst->Quad1D = BuildQuad1D(Quad_Points1D);
    //printf("333333333334444444444\n");
    rhst->Entries = malloc((rhst->DOF_ALL)*sizeof(double));
    //printf("55555555555555555\n");
    memset(rhst->Entries,0.0,rhst->DOF_ALL*sizeof(double));
    //printf("4444444444444444444\n");
    return rhst;
}




// Assemble the RHST (right hand side term)
void AssembleRHST(RHST *Rhs)
{
  Fespace3D *test_space;
  test_space = Rhs->Test_Space;
  MESH *mesh;
  mesh =test_space->Mesh;
  BASEFUNCTION3D *test_base;
  test_base = test_space->Base;
  int n_bas_test, N_Rhs_Entries;
  double *Rhs_Entries;
  //printf ( "rhs11111\n" );
  n_bas_test = test_base->Num_Bas;  
  N_Rhs_Entries = n_bas_test;
  Rhs_Entries = malloc(N_Rhs_Entries*sizeof(double));
  int dim_test_value;
  dim_test_value = test_base->Value_Dim;
  
  int *Test_GlobalNumbers, *Test_BeginIndex;
  Test_GlobalNumbers = test_space->GlobalNumbers;
  Test_BeginIndex = test_space->BeginIndex;
  
  //printf ( "rhs22222\n" );
  int num_volus,i,j;
  num_volus = mesh->Num_Volus_Global;
  Quadrature3D *Quad3D;
  Quadrature2D *Quad2D;
  Quadrature1D *Quad1D;

  Quad3D = Rhs->Quad3D;  
  Quad2D = NULL;  
  Quad1D = NULL;
  if(Rhs->Quad1D)
    Quad1D = Rhs->Quad1D;
  if(Rhs->Quad2D)
    Quad2D = Rhs->Quad2D;
  //printf ( "rhs33333\n" );
  ELEMENT3D *elem3D;
  elem3D = malloc(1*sizeof(ELEMENT3D));
  InitialElem3D(mesh,elem3D); 

  //取出右端项所需要计算的导数值
  int n_multiindex_test;
  n_multiindex_test = Rhs->N_Test_MultiIndex;  
  MultiIndex3D *Test_MultiIndex;
  Test_MultiIndex = Rhs->Test_MultiIndex;
  //printf ( "rhs4444444444444\n" );
  double *Test_Values, *Test_MultiIndex_Values;  
  Test_Values = malloc(n_bas_test*dim_test_value*n_multiindex_test*sizeof(double));
  Test_MultiIndex_Values = malloc(n_multiindex_test*dim_test_value*sizeof(double));
  memset(Test_MultiIndex_Values,0.0,n_multiindex_test*dim_test_value*sizeof(double));
  //在体上的离散变分形式函数
  DiscreteFormRHSVolu *discreteformvolu;  
  //在面上的离散变分形式函数
  DiscreteFormRHSFace *discreteformface;  
  //printf ( "rhs5555555555\n" );
  discreteformvolu = Rhs->DiscreteFormVolu;
  discreteformface = Rhs->DiscreteFormFace;
  
  double Quad_xi, Quad_eta, Quad_zeta, Quad_W;
  int curr_volu, Quad_Points3D, Quad_Points2D, Quad_Points1D;
  Quad_Points3D = Quad3D->N_Points;
  if(Rhs->Quad1D)
    Quad_Points1D = Quad1D->N_Points;
  if(Rhs->Quad2D)
    Quad_Points2D = Quad2D->N_Points;
  //printf ( "rhs666666666666\n" );
  int N_AuxFeFun, *N_AuxFeFun_MultiIndex, N_AuxFeFun_Values;
  MultiIndex3D *AuxFeFun_MultiIndex;
  N_AuxFeFun = Rhs->N_AuxFeFun;
  AuxFeFun_MultiIndex = Rhs->AuxFeFun_MultiIndex;
  N_AuxFeFun_MultiIndex = Rhs->N_AuxFeFun_MultiIndex;
  //printf ( "rhs7777\n" );
  Fefunction3D **AuxFeFuns;
  N_AuxFeFun_Values = Rhs->N_AuxFeFun_Values;
  AuxFeFuns = Rhs->AuxFeFun;
  double *AuxFeFun_Values;
  AuxFeFun_Values = Rhs->AuxFeFun_Values;
  int *Dim_Values_AuxFeFun;
  //printf ( "rhs777778888888888\n" );
  if(N_AuxFeFun)
  {
    Dim_Values_AuxFeFun = malloc(N_AuxFeFun*sizeof(int));
    for(i=0;i<N_AuxFeFun;i++)
      Dim_Values_AuxFeFun[i] = AuxFeFuns[i]->Fespace->Base->Value_Dim;
  }
  //printf ( "rrhs888888888888888\n" );

  int num_feconst = Rhs->N_FeConst;
  double *feconst;
  feconst = Rhs->FeConst;

  double XI,ETA;
  double local_direction;
  int elem_faces, ID_bd;
 
  BoundType bdtype; 
  BoundCondFunct3D *BoundCond;
  BoundCond = test_space->BoundCond;
  
  //N_AuxFeFun_Values = 0;
  //AuxFeFun_Values = NULL;
  printf ( "beign to loop\n" );
  for(curr_volu=0;curr_volu<num_volus;curr_volu++)
  {
   //printf("curr_volu==%d\n",curr_volu);
    //volu = Volus[curr_volu];
    GetElement3D(mesh,curr_volu,elem3D);
    memset(Rhs_Entries,0.0,N_Rhs_Entries*sizeof(double));
    if(discreteformvolu!=NULL)
    {
      for(i=0;i<Quad_Points3D;i++)
      {

	Quad_xi = Quad3D->Quad_X[i];
	Quad_eta = Quad3D->Quad_Y[i];
	Quad_zeta = Quad3D->Quad_Z[i];
	Quad_W = Quad3D->Quad_W[i];
	//先获得所有矩阵和右端项在本积分点的值
	//assemble 
	//每个矩阵包含自己的取值空间和试探空间，后面的微分算子也是针对
	//针对这两个空间的。另外可能还有其他的辅助函数空间和相应的微分 
	//算子，最后也会包含辅助函数来定义矩阵
	//
	SubAssembleRhs3D(elem3D,test_base,Test_Values,n_multiindex_test,Test_MultiIndex,
			  Test_MultiIndex_Values, Dim_Values_AuxFeFun,N_AuxFeFun,AuxFeFuns,
			  N_AuxFeFun_MultiIndex,
			  AuxFeFun_MultiIndex, N_AuxFeFun_Values,AuxFeFun_Values,num_feconst,feconst,
			  Quad_xi,Quad_eta,Quad_zeta,Quad_W,discreteformvolu,Rhs_Entries); 	
         
      } 
    }//end if
    //  printf ( "-----------------------------------\n" );
    //printf ( "-----------------------------------------------\n" );
    //printf ( "%d-th,%lf,%lf,%lf,%lf\n",curr_volu,Rhs_Entries[0], Rhs_Entries[1],Rhs_Entries[2],Rhs_Entries[3]);
    AddSubRhs(Rhs,Rhs_Entries,n_bas_test,Test_GlobalNumbers+Test_BeginIndex[curr_volu]);
 }//单元循环结束    
    FreeElem3D(elem3D);
    free(Rhs_Entries);
    free(Test_Values);
    free(Test_MultiIndex_Values);
    //free(normal);
    //free(tangent);
    if(N_AuxFeFun)
    {
      free(Dim_Values_AuxFeFun);
    }
}

void SubAssembleRhs3D(ELEMENT3D *elem3D,BASEFUNCTION3D *test_base, double *Test_Values,
		      int n_test_multiindex,MultiIndex3D *test_multiindex, double *Test_MultiIndex_Values,
		      int *Dim_Values_AuxFeFun, 
		      int N_AuxFeFun, Fefunction3D **AuxFeFuns, int *N_AuxFeFun_MultiIndex,
		      MultiIndex3D *AuxFeFun_MultiIndex,
		      int N_AuxFeFun_Values, double *AuxFeFun_Values,int N_FeConst,double *FeConst,double  Quad_xi,double Quad_eta,
		      double Quad_zeta,double Quad_W, 
		      DiscreteFormRHSVolu *discreteformvolu, double *Rhs_Entries)
{
    double Quad_x, Quad_y,Quad_z;
    GetElementCoord(elem3D,Quad_xi,Quad_eta, Quad_zeta,&Quad_x,&Quad_y,&Quad_z); 
    int n_bas_test,i,j,k;
    n_bas_test = test_base->Num_Bas;

    //discreteformface = RHS->DiscreteFormFace;
    int n_userfunction;
    double *userfunctionvalue;
    n_userfunction=0;
      userfunctionvalue=NULL;
    int n_feconst;
    double *feconst;
    n_feconst=N_FeConst;
    feconst= FeConst;
    //MultiIndex2D *test_multiindex;
    //BASEFUNCTION2D *test_base;
    //test_base = RHS->Test_Space->Base;
    //test_multiindex = RHS->Test_MultiIndex;
    int dim_test_value;
    dim_test_value = test_base->Value_Dim;
    //compute the values of the test bases
    for(i=0;i<n_test_multiindex;i++)
    {
      GetBaseValues(test_base,test_multiindex[i],elem3D,Quad_x,Quad_y,Quad_z,Quad_xi,Quad_eta,Quad_zeta,
		                                   Test_Values+i*dim_test_value*n_bas_test);	
    }
    int curr_pos, m; 
    Fefunction3D *AuxFeFun;
    if(N_AuxFeFun>0)
    {
        curr_pos = 0;
        for(i=0;i<N_AuxFeFun;i++)
        {
	  AuxFeFun = AuxFeFuns[i];
	  GetFefunctValues(AuxFeFuns[i],elem3D,Quad_x,Quad_y,Quad_z,
                             Quad_xi,Quad_eta,Quad_zeta,N_AuxFeFun_MultiIndex[i],
                             AuxFeFun_MultiIndex+curr_pos,AuxFeFun_Values+curr_pos);
            curr_pos += N_AuxFeFun_MultiIndex[i]*Dim_Values_AuxFeFun[i];            
        }//end for i
    }
    else
    {
	N_AuxFeFun_Values = 0;
    }
    //for the assemble 
    for(i=0;i<n_bas_test;i++)
    {
      // get the values at each element of the test and anasatz with all the multiindices
      //printf("i=%d\n",i);
      for(k=0;k<n_test_multiindex;k++)
      {
	for(m=0;m<dim_test_value;m++)
	{
	  Test_MultiIndex_Values[k*dim_test_value+m] = Test_Values[k*n_bas_test*dim_test_value+i*dim_test_value+m];
	}//end for m
      }//end for k
      
      Rhs_Entries[i] += discreteformvolu(Quad_x,Quad_y,Quad_z,n_test_multiindex,Test_MultiIndex_Values,
				         N_AuxFeFun_Values,AuxFeFun_Values,n_userfunction,userfunctionvalue,
					 n_feconst,feconst)
	                *elem3D->Volu*Quad_W;
       //printf("Rhs_Entries[%d]=%f\n",i,Rhs_Entries[i]);
    }// end for i
}
//for check, which is faster, A simple parameter version
                                                   
void SubAssembleRhs3DSimple(RHST *RHS, ELEMENT3D *elem3D,double *Test_Values,
		      double *Test_MultiIndex_Values, int *Dim_Values_AuxFeFun, 
		      double  Quad_xi,double Quad_eta,double Quad_zeta,double Quad_W,  double *Rhs_Entries)
{
    double Quad_x, Quad_y, Quad_z;    
    GetElementCoord(elem3D,Quad_xi,Quad_eta, Quad_zeta, &Quad_x,&Quad_y,&Quad_z); 
    Fespace3D *test_space; 
    test_space = RHS->Test_Space;
    
    MultiIndex3D *test_multiindex;
    BASEFUNCTION3D *test_base;
    test_base = RHS->Test_Space->Base;
    test_multiindex = RHS->Test_MultiIndex;
 
    int n_bas_test,i,j,k;
    n_bas_test = test_base->Num_Bas;

    int n_test_multiindex;
    n_test_multiindex = RHS->N_Test_MultiIndex; 
    DiscreteFormRHSVolu *discreteformvolu;
    discreteformvolu = RHS->DiscreteFormVolu;

    int dim_test_value;
    dim_test_value = test_base->Value_Dim;
    //compute the values of the test bases 
    for(i=0;i<n_test_multiindex;i++)
    {
      //printf("i=%d\n",i);
        //Base->Base[multiindex[i]](elem2D,Quad_x,Quad_y,Quad_xi,Quad_eta,Test_Values+i*n_bas_test);    
	GetBaseValues(test_base,test_multiindex[i],elem3D,Quad_x,Quad_y,Quad_z,Quad_xi,Quad_eta,Quad_zeta,
		                                   Test_Values+i*dim_test_value*n_bas_test);	
    }
    int curr_pos, m; 
    //N_AuxFeFun,*AuxFeFunPtr, *N_AuxFeFun_MultiIndex;//N_AuxFefun_Values
    //MultiIndex2D *AuxFeFun_MultiIndex;
    //double *AuxFefun_Values;
    Fefunction3D *AuxFeFun, **AuxFeFuns;
    AuxFeFuns = RHS->AuxFeFun;
    int N_AuxFeFun;
    N_AuxFeFun = RHS->N_AuxFeFun;
    
    MultiIndex3D *AuxFeFun_MultiIndex;
    AuxFeFun_MultiIndex = RHS->AuxFeFun_MultiIndex;
    
     int *N_AuxFeFun_MultiIndex;
     N_AuxFeFun_MultiIndex = RHS->N_AuxFeFun_MultiIndex;
     int N_AuxFeFun_Values;
     N_AuxFeFun_Values = RHS->N_AuxFeFun_Values;
     double *AuxFeFun_Values;
     AuxFeFun_Values = RHS->AuxFeFun_Values;	      
		      
    if(N_AuxFeFun>0)
    {
        //printf("ttttttttttttttt\n");
        //N_AuxFeFun = RHS->N_AuxFeFun;
        //AuxFeFun_MultiIndex = RHS->AuxFeFun_MultiIndex;
        //N_AuxFeFun_MultiIndex = RHS->N_AuxFeFun_MultiIndex;
	//AuxFefunct = RHS->AuxFeFunct;
        //N_AuxFefun_Values = 0;   
	//printf("ttttttttttttttt\n");
        //for(i=0;i<N_AuxFefunct;i++)
        //    N_AuxFefun_Values += N_AuxFefun_MultiIndex[i];
        //printf("ttttttttttttttt\n");
        //AuxFefun_Values = malloc(N_AuxFefun_Values*sizeof(double));
        curr_pos = 0;
        for(i=0;i<N_AuxFeFun;i++)
        {
	  AuxFeFun = AuxFeFuns[i];
	  GetFefunctValues(AuxFeFuns[i],elem3D,Quad_x,Quad_y, Quad_z, 
                             Quad_xi,Quad_eta,Quad_zeta,N_AuxFeFun_MultiIndex[i],
                             AuxFeFun_MultiIndex+curr_pos,AuxFeFun_Values+curr_pos);
            curr_pos += N_AuxFeFun_MultiIndex[i]*Dim_Values_AuxFeFun[i];            
        }//end for i
    }
    else
    {
      //printf("No AuxFefunct\n");
      //RHS->N_AuxFeFun_Values = 0;
       //RHS->AuxFeFun_Values = NULL;
	N_AuxFeFun_Values = 0;
        //AuxFeFun_Values = NULL;
    }
    //printf("111111111111\n");
    //for the assemble 
    for(i=0;i<n_bas_test;i++)
    {
      // get the values at each element of the test and anasatz with all the multiindices
      //printf("i=%d\n",i);
      for(k=0;k<n_test_multiindex;k++)
      {
	for(m=0;m<dim_test_value;m++)
	{
	  Test_MultiIndex_Values[k*dim_test_value+m] = Test_Values[k*n_bas_test*dim_test_value+i*dim_test_value+m];
	}//end for m
      }//end for k
      Rhs_Entries[i] += discreteformvolu(Quad_x,Quad_y,Quad_z,n_test_multiindex,Test_MultiIndex_Values,
					 N_AuxFeFun_Values,AuxFeFun_Values,0,NULL,0,NULL)*elem3D->Volu*Quad_W;
       //printf("Rhs_Entries[%d]=%f\n",i,Rhs_Entries[i]);
    }// end for i
}



void SubAssembleRhs3DAll(ELEMENT3D *elem3D,BASEFUNCTION3D *test_base, double *Test_Values,
		      int n_test_multiindex,MultiIndex3D *test_multiindex, double *Test_MultiIndex_Values,
		      int *Dim_Values_AuxFeFun, 
		      int N_AuxFeFun, Fefunction3D **AuxFeFuns, int *N_AuxFeFun_MultiIndex,
		      MultiIndex3D *AuxFeFun_MultiIndex,
		      int N_AuxFeFun_Values, double *AuxFeFun_Values,
                      int N_UserFunction,Functions3D **UserFunction, int N_FeConst,double *FeConst,
                      double  Quad_xi,double Quad_eta,
		      double Quad_zeta,double Quad_W, 
		      DiscreteFormRHSVolu *discreteformvolu, double *Rhs_Entries)
{  
    //printf ( "subsubsubsubsub\n");
    double Quad_x, Quad_y,Quad_z;
    GetElementCoord(elem3D,Quad_xi,Quad_eta,Quad_zeta,&Quad_x,&Quad_y,&Quad_z); 
    int n_bas_test,i,j,k;
    n_bas_test = test_base->Num_Bas;
    int n_userfunction;
    double *userfunctionvalue;
    //printf ( "cccccccccccccccccccccc\n");
    n_userfunction=N_UserFunction;
    userfunctionvalue=malloc(sizeof(double)*n_userfunction);
    //printf ( "vvvvvvvvvvvvvvvvvvv\n");
    //printf ( "%d,%d\n", N_UserFunction,n_userfunction);
    if(n_userfunction)
    {
       for(i=0;i<n_userfunction;++i)
    	  UserFunction[i](Quad_x,Quad_y,Quad_z,1,userfunctionvalue+i);
    }
    else
      userfunctionvalue=NULL;
    //printf ( "zzzzzzzzzzzzzzzzzzzzzzzzzzzz\n" );
    int n_feconst;
    double *feconst;
    n_feconst=N_FeConst;
    feconst=FeConst;
    //printf ( "xxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n" );
    //discreteformface = RHS->DiscreteFormFace;
    int dim_test_value;
    dim_test_value = test_base->Value_Dim;
    //compute the values of the test bases
    for(i=0;i<n_test_multiindex;i++)
    {
      GetBaseValues(test_base,test_multiindex[i],elem3D,Quad_x,Quad_y,Quad_z,Quad_xi,Quad_eta,Quad_zeta,
    		                                   Test_Values+i*dim_test_value*n_bas_test);	
    }
    int curr_pos, m; 
    Fefunction3D *AuxFeFun;
    if(N_AuxFeFun>0)
    {
        curr_pos = 0;
        for(i=0;i<N_AuxFeFun;i++)
        {
	  AuxFeFun = AuxFeFuns[i];
	  GetFefunctValues(AuxFeFuns[i],elem3D,Quad_x,Quad_y,Quad_z,
                             Quad_xi,Quad_eta,Quad_zeta,N_AuxFeFun_MultiIndex[i],
                             AuxFeFun_MultiIndex+curr_pos,AuxFeFun_Values+curr_pos);
            curr_pos += N_AuxFeFun_MultiIndex[i]*Dim_Values_AuxFeFun[i];            
        }//end for i
    }
    else
    {
	N_AuxFeFun_Values = 0;
    }
    //for the assemble 
    for(i=0;i<n_bas_test;i++)
    {
      // get the values at each element of the test and anasatz with all the multiindices
      //printf("i=%d\n",i);
      for(k=0;k<n_test_multiindex;k++)
      {
	for(m=0;m<dim_test_value;m++)
	{
	  Test_MultiIndex_Values[k*dim_test_value+m] = Test_Values[k*n_bas_test*dim_test_value+i*dim_test_value+m];
	}//end for m
      }//end for k
      
      Rhs_Entries[i] += discreteformvolu(Quad_x,Quad_y,Quad_z,n_test_multiindex,Test_MultiIndex_Values,
					 N_AuxFeFun_Values,AuxFeFun_Values,n_userfunction,userfunctionvalue, 
                                         n_feconst,feconst)*elem3D->Volu*Quad_W;
       //printf("Rhs_Entries[%d]=%f\n",i,Rhs_Entries[i]);
    }// end for i
}





/*  
void SubAssembleRhs3DFace(int ID_bd,double *normal, double *tangent, double local_direction,
			                RHST *RHS,ELEMENT3D *elem3D,double *Test_Values,
			                double *Test_MultiIndex_Values,  double  Quad_xi, 
			                double Quad_eta,double Quad_zeta,double Quad_W, 
					double area, DiscreteFormRHSFace *discreteformface, 
					double *Rhs_Entries)
{
    double Quad_x,Quad_y,Quad_z;
    GetElementCoord(elem3D,Quad_xi,Quad_eta,Quad_zeta, &Quad_x,&Quad_y,&Quad_z);
    Fespace3D *test_space; 
    test_space = RHS->Test_Space;
 
    int n_bas_test, n_test_multiindex,i,j,k;
    MultiIndex3D *test_MultiIndex;
    n_bas_test = test_space->Base->Num_Bas;

    n_test_multiindex = RHS->N_Test_MultiIndex;
    
    MultiIndex3D *test_multiindex;
    BASEFUNCTION3D *test_base;
    test_base = RHS->Test_Space->Base;
    test_multiindex = RHS->Test_MultiIndex;
    int dim_test_value;
    dim_test_value = test_base->Value_Dim;
    //compute the values of the test bases
    for(i=0;i<n_test_multiindex;i++)
    {
	GetBaseValues(test_base,test_multiindex[i],elem3D,Quad_x,Quad_y,Quad_z,Quad_xi,Quad_eta,
		      Quad_zeta,Test_Values+i*n_bas_test*dim_test_value);	
    }
    int curr_pos,N_AuxFeFun,*AuxFeFunPtr, *N_AuxFeFun_MultiIndex;
    MultiIndex3D *AuxFeFun_MultiIndex;
    Fefunction3D *AuxFeFun;
    double *AuxFeFun_Values;
    int m;
    if(RHS->N_AuxFeFun)
    {
        //printf("ttttttttttttttt\n");
        N_AuxFeFun = RHS->N_AuxFeFun;
        AuxFeFun_MultiIndex = RHS->AuxFeFun_MultiIndex;
        N_AuxFeFun_MultiIndex = RHS->N_AuxFeFun_MultiIndex;
	AuxFeFun_Values = RHS->AuxFeFun_Values;
        curr_pos = 0;
        for(i=0;i<N_AuxFeFun;i++)
        {
            AuxFeFun = RHS->AuxFeFun[i];
            GetFefunctValues(AuxFeFun,elem3D,Quad_x,Quad_y,Quad_z, 
                             Quad_xi,Quad_eta,Quad_zeta,N_AuxFeFun_MultiIndex[i],
                             AuxFeFun_MultiIndex+curr_pos,AuxFeFun_Values+curr_pos);
            curr_pos += N_AuxFeFun_MultiIndex[i]*AuxFeFun->Fespace->Base->Value_Dim;            
        }//end for i
    }
    else
    {
        RHS->N_AuxFeFun_Values = 0;
    }
    //for the assemble discreteformline
    for(i=0;i<n_bas_test;i++)
    {
      for(k=0;k<n_test_multiindex;k++)
      {
	for(m=0;m<dim_test_value;m++)
	{
	  Test_MultiIndex_Values[k*dim_test_value+m] = Test_Values[k*n_bas_test*dim_test_value+i*dim_test_value+m];
	}//end for m
      }//end for k
      
      Rhs_Entries[i] += discreteformface(ID_bd,Quad_x,Quad_y,Quad_z,normal,tangent,local_direction,
	                                 n_test_multiindex,Test_MultiIndex_Values,RHS->N_AuxFeFun_Values,
	                                 RHS->AuxFeFun_Values)*area*Quad_W;
    }// end for i
}


*/


// add a sub rhs to the RHS
void AddSubRhs(RHST *Rhs, double *subRhs,int n_row, int *row)
{
  int i;
  for(i=0;i<n_row;i++)
  {
    Rhs->Entries[row[i]] += subRhs[i];
  }
}

//输出RHS的内容
void OutPutRHST(RHST *Rhs)
{
  int i,dof_all;
  dof_all = Rhs->DOF_ALL;
  
  printf("\n");
  printf("Output the Rhs: \n");
  for(i=0;i<dof_all;i++)
    printf("Rhs(%d)=%18.15f\n",i+1,Rhs->Entries[i]);
  
  printf("\n");
}
//释放出RHST的内存
void FreeRHST(RHST *Rhs)
{
  printf("delete the RHST!\n");
  if(Rhs->AuxFeFun_Values)
    free(Rhs->AuxFeFun_Values);
   if(Rhs->Quad3D)
    FreeQuad3D(Rhs->Quad3D);
 if(Rhs->Quad2D)
    FreeQuad2D(Rhs->Quad2D);
  if(Rhs->Quad1D)
    FreeQuad1D(Rhs->Quad1D);
  if(Rhs->Entries)
    free(Rhs->Entries);
  free(Rhs);
}







//////////////////////////////////////////////
//
   /* 
    if(discreteformface!=NULL)
     {
      elem_faces = 4;
      for(j=0;j<elem_faces;j++)
      {
	face = Faces[volu->Faces[j]];
	ID_bd = face->ID_Boundary;
	BoundCond(ID_bd,0.0,&bdtype);
	if(bdtype!=DIRICHLETT)//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	{
	  area = FACEE_AREA(mesh,face);//!!!!!!!!!!!!!!!!!!!!!!!!!
	  //首先获得在边界上的积分点和积分的权重
	  for(i=0;i<Quad_Points2D;i++)
	  {
	    XI = Quad2D->Quad_X[i];
	    ETA = Quad2D->Quad_Y[i];
	    Quad_xi = xis[j]*(1.0-XI)  + xis[j+1]*XI;
	    Quad_eta = ets[j]*(1.0-XI) + ets[j+1]*XI;
	    Quad_zeta = ets[j]*(1.0-XI) + ets[j+1]*XI;
	    Quad_W = Quad1D->Quad_W[i];	  
	    //如果是在边界上，并且不是Dirichlet强制边界条件，就需要处理边界条件所带来的对刚度矩阵的处理
	    SubAssembleRhs3DFace(ID_bd,normal,tangent,local_direction,Rhs,elem3D,Test_Values,
		                 Test_MultiIndex_Values,Quad_xi,Quad_eta,Quad_zeta,Quad_W,area,
				 discreteformface,Rhs_Entries);  	    
	  }//end for i	  
	}//end for bytype!=DIRICHLETT
      }//end for j      
    }//end if	
*/

