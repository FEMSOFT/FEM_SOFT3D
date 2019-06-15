/*
 * =====================================================================================
 *
 *       Filename:  Fefunction3D.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2014年10月29日 15时23分09秒
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
#include "Fefunction3D.h"
#include "Quadrature3D.h"
#include "MG.h"
#include "Base3D.h"
#include "Mesh3D.h"

Fefunction3D * BuildFefunction3D(Fespace3D *fesp,double *values)
{
  Fefunction3D *funct3D;
  funct3D = malloc(sizeof(Fefunction3D));
  funct3D->DOF_ALL = fesp->DOF_ALL;
  funct3D->Fespace = fesp;
  funct3D->Values = values;
  
  return funct3D;
}
// compute the value on one point of one element
double GetFefunctValue(Fefunction3D *Fefunct,ELEMENT3D *elem3D,
                     double Quad_x,double Quad_y,double Quad_z, 
                     double Quad_xi,double Quad_eta,double Quad_zeta,
                     MultiIndex3D MultiIndex) 
{
  double value,*Bas_Values;
  int ID_elem,i,j, num_bas, *GlobalNumbers, *BeginIndex;
  int start;
  ID_elem = elem3D->ID_Num;
  num_bas = Fefunct->Fespace->Base->Num_Bas;
  Bas_Values = malloc(num_bas*sizeof(double));
  GetBaseValues(Fefunct->Fespace->Base,MultiIndex,elem3D,Quad_x,Quad_y,Quad_z,
	        Quad_xi,Quad_eta,Quad_zeta,Bas_Values);
  BeginIndex = Fefunct->Fespace->BeginIndex;
  GlobalNumbers = Fefunct->Fespace->GlobalNumbers;
  start = BeginIndex[ID_elem];
  value = 0.0;
  for(i=0;i<num_bas;i++)
    value += Bas_Values[i]*Fefunct->Values[GlobalNumbers[start+i]];
  free(Bas_Values);  
  return value;
}

// compute the values at a point on an element
void GetFefunctValues(Fefunction3D *Fefunct,ELEMENT3D *elem3D,
                     double Quad_x,double Quad_y,double Quad_z, 
                     double Quad_xi,double Quad_eta,double Quad_zeta,
		     int N_MultiIndex,MultiIndex3D *MultiIndex,
		     double *Values)
{
  double value,*Bas_Values;
  int ID_elem,i,j,num_bas,*GlobalNumbers,*BeginIndex;
  int start;
  ID_elem = elem3D->ID_Num;
  num_bas = Fefunct->Fespace->Base->Num_Bas;
  Bas_Values = malloc(num_bas*sizeof(double));
  BeginIndex = Fefunct->Fespace->BeginIndex;;
  GlobalNumbers = Fefunct->Fespace->GlobalNumbers;
  start = BeginIndex[ID_elem];
  for(i=0;i<N_MultiIndex;i++)
  {
    GetBaseValues(Fefunct->Fespace->Base,MultiIndex[i],elem3D,Quad_x,Quad_y,Quad_z,Quad_xi,
	          Quad_eta,Quad_zeta,Bas_Values);
    Values[i] = 0.0;
    for(j=0;j<num_bas;j++)
    {
      Values[i] += Bas_Values[j]*Fefunct->Values[GlobalNumbers[start+j]];
    }
  }
  free(Bas_Values);
}


void Interpolation3D(Fefunction3D *fefun3D,Functions3D *funct3D)
{
  int Num_Volus,curr_volu, i, j, start;
  MESH *Mesh;
  Mesh =  fefun3D->Fespace->Mesh;
  Num_Volus = Mesh->Num_Volus_Global;
  int *BeginIndex, *GlobalNumbers;
  double *Values;
  BeginIndex = fefun3D->Fespace->BeginIndex;
  GlobalNumbers = fefun3D->Fespace->GlobalNumbers;
  Values = fefun3D->Values;
 
  int num_bas;
  double *elem_values;
  num_bas = fefun3D->Fespace->Base->Num_Bas;
  //dim_value = fefun2D->Fespace->Base->Value_Dim;
  elem_values = malloc(num_bas*sizeof(double));
  memset(elem_values,0,num_bas*sizeof(double));
  ELEMENT3D *elem3D;
  elem3D = malloc(1*sizeof(ELEMENT3D));
  // initializing the element
  InitialElem3D(Mesh,elem3D);
  for(curr_volu=0;curr_volu<Num_Volus;curr_volu++)
  {
    GetElement3D(Mesh,curr_volu,elem3D);
    fefun3D->Fespace->Base->NodalF3D(elem3D,funct3D,elem_values);
    start = BeginIndex[curr_volu];
    for(i=0;i<num_bas;i++)
      Values[GlobalNumbers[start+i]] = elem_values[i];
  }//end for curr_face
  free(elem_values);
  FreeElem3D(elem3D);
}



// Interpolation of the Fefunction and Function2D
/*
void Interpolation3D(Fefunction3D *fefun3D,Functions3D *funct3D)
{
  int Num_Faces,curr_volu, i, j, start;
  MESH *Mesh;
  Mesh =  fefun3D->Fespace->Mesh;
  Num_Faces = Mesh->Num_Faces;
  int *BeginIndex, *GlobalNumbers;
  double *Values; 
  BeginIndex = fefun3D->Fespace->BeginIndex;
  GlobalNumbers = fefun3D->Fespace->GlobalNumbers;
  Values = fefun3D->Values;
  
  int num_bas; 
  double *elem_values;
  num_bas = fefun3D->Fespace->Base->Num_Bas;
  //dim_value = fefun2D->Fespace->Base->Value_Dim;
  elem_values = malloc(num_bas*sizeof(double));
  memset(elem_values,0,num_bas*sizeof(double));
  ELEMENT3D *elem3D;
  elem3D = malloc(1*sizeof(ELEMENT3D));
  // initializing the element
  InitialElem3D(Mesh,elem3D);
  
  for(curr_volu=0;curr_volu<Num_Faces;curr_volu++)
  {
    GetElement3D(Mesh,curr_volu,elem3D);
    fefun3D->Fespace->Base->NodalF3D(elem3D,funct3D,elem_values);
    start = BeginIndex[curr_volu];
    for(i=0;i<num_bas;i++)
      Values[GlobalNumbers[start+i]] = elem_values[i];
  }//end for curr_face
  free(elem_values);
  FreeElem2D(elem2D);
}*/
// Error estimate of the exact solution 
// N_Funct2d = N_Derives, one function corresponding to one derives
// L^2 error: N_Funct2d = N_Derives =1, multiindex2D = D000
// H^1 error: N_Funct2d = N_Derives =3, multiindex2D = [D00,D100,D010]
/*
double ErrorEstimate3D(Fefunction3D *fefun3D,Functions3D *Exact3D, 
		       int N_Derives,MultiIndex3D *multiindex3D, int Quad_Points)
{
  printf("Compute the error estimate\n");
  int i,j,num_volu;
  MESH *mesh;
  mesh = fefun3D->Fespace->Mesh;
  num_volu = mesh->Num_Volus_Global;
  ELEMENT3D *elem3D;
  elem3D = malloc(1*sizeof(ELEMENT3D));
  InitialElem3D(mesh,elem3D);
  
  double Quad_xi, Quad_eta,Quad_zeta, Quad_W;
  int curr_volu, kk,ii;

  Quadrature3D *Quad3D;
  Quad3D = BuildQuad2Mesh3D(fefun3D->Fespace->Mesh, Quad_Points);
  
  Quad_Points = Quad3D->N_Points;
  
  double Error, Err_tmp, tmp, Err_K;
  Error = 0.0;
  //Err_tmp = 0.0;
  int dim_value;
  dim_value = fefun3D->Fespace->Base->Value_Dim;
  double *Values_FE, *Values_Exact;
  Values_FE = malloc(N_Derives*dim_value*sizeof(double));
  Values_Exact = malloc(N_Derives*dim_value*sizeof(double));
  double Quad_x,Quad_y,Quad_z;
  
  for(curr_volu=0;curr_volu<num_volu;curr_volu++)
  {
    GetElement3D(mesh,curr_volu,elem3D);
    Err_K = 0.0;
    for(i=0;i<Quad_Points;i++)
    {
      Quad_xi = Quad3D->Quad_X[i];
      Quad_eta = Quad3D->Quad_Y[i];
      Quad_zeta = Quad3D->Quad_Z[i];
      Quad_W = Quad3D->Quad_W[i];
      GetElementCoord(elem3D,Quad_xi,Quad_eta, Quad_zeta,&Quad_x,&Quad_y,&Quad_z);
      //printf ( "%f,%f,%f,%f,%f,%f\n", Quad_xi, Quad_eta, Quad_zeta, Quad_x,Quad_y,Quad_z );
      //求真解的值
      Exact3D(Quad_x,Quad_y,Quad_z,N_Derives*dim_value,Values_Exact);
      //求有限元的值
      GetFefunctValues(fefun3D,elem3D,Quad_x,Quad_y,Quad_z,
		       Quad_xi,Quad_eta,Quad_zeta,N_Derives,multiindex3D,Values_FE);
      Err_tmp = 0.0;
      for(j=0;j<N_Derives*dim_value;j++)
      {
	tmp = Values_Exact[j]-Values_FE[j];
//	printf("i=%d, j=%d, curr_volu=%d, tmp=%f,exact=%f,fevalue=%f\n",i,j,curr_volu,tmp,
//	       Values_Exact[j],Values_FE[j]);	
        Err_tmp += tmp*tmp;
      }//end for j
      Err_K += Err_tmp*Quad_W;
    }//end for i
    Error += Err_K*elem3D->Volu;   
    //printf ( "error=%d,%lf\n",curr_volu,Error );
  }//end for curr_volu  
  
  FreeElem3D(elem3D);  
  free(Values_Exact);
  free(Values_FE);
  FreeQuad3D(Quad3D);
  return sqrt(Error);
}


double ErrorEstimateEigen3DExactSolution(Fefunction3D *fefun3D,Functions3D *Exact3D, 
		       int N_Derives,MultiIndex3D *multiindex3D, int Quad_Points)
{
  printf("Compute the error estimate\n");
  int i,j,num_volu;
  MESH *mesh;
  mesh = fefun3D->Fespace->Mesh;
  num_volu = mesh->Num_Volus_Global;
  ELEMENT3D *elem3D;
  elem3D = malloc(1*sizeof(ELEMENT3D));
  InitialElem3D(mesh,elem3D);
  
  double Quad_xi, Quad_eta,Quad_zeta, Quad_W;
  int curr_volu, kk,ii;

  Quadrature3D *Quad3D;
  Quad3D = BuildQuad2Mesh3D(fefun3D->Fespace->Mesh, Quad_Points);
  
  Quad_Points = Quad3D->N_Points;
  
  double Error, Err_tmp, tmp, Err_K;
  Error = 0.0;
  //Err_tmp = 0.0;
  int dim_value;
  dim_value = fefun3D->Fespace->Base->Value_Dim;
  double *Values_FE, *Values_Exact;
  Values_FE = malloc(N_Derives*dim_value*sizeof(double));
  Values_Exact = malloc(N_Derives*dim_value*sizeof(double));
  double Quad_x,Quad_y,Quad_z;
 
  double *tmper;
  tmper = malloc(sizeof(double));
  GetMaxValue(fefun3D->Values,fefun3D->DOF_ALL, tmper);
  for(i=0;i<fefun3D->DOF_ALL;++i)
     fefun3D->Values[i]/=tmper[0];
  free(tmper);

  for(curr_volu=0;curr_volu<num_volu;curr_volu++)
  {
    GetElement3D(mesh,curr_volu,elem3D);
    Err_K = 0.0;
    for(i=0;i<Quad_Points;i++)
    {
      Quad_xi = Quad3D->Quad_X[i];
      Quad_eta = Quad3D->Quad_Y[i];
      Quad_zeta = Quad3D->Quad_Z[i];
      Quad_W = Quad3D->Quad_W[i];
      GetElementCoord(elem3D,Quad_xi,Quad_eta, Quad_zeta,&Quad_x,&Quad_y,&Quad_z);
      //求真解的值
      Exact3D(Quad_x,Quad_y,Quad_z,N_Derives*dim_value,Values_Exact);
      //求有限元的值
      GetFefunctValues(fefun3D,elem3D,Quad_x,Quad_y,Quad_z,
		       Quad_xi,Quad_eta,Quad_zeta,N_Derives,multiindex3D,Values_FE);
      Err_tmp = 0.0;
      for(j=0;j<N_Derives*dim_value;j++)
      {
	tmp = Values_Exact[j]-Values_FE[j];
//	printf("i=%d, j=%d, curr_volu=%d, tmp=%f,exact=%f,fevalue=%f\n",i,j,curr_volu,tmp,
//	       Values_Exact[j],Values_FE[j]);	
        Err_tmp += tmp*tmp;
      }//end for j
      Err_K += Err_tmp*Quad_W;
    }//end for i
    Error += Err_K*elem3D->Volu;   
    //printf ( "error=%d,%lf\n",curr_volu,Error );
  }//end for curr_volu  
  
  FreeElem3D(elem3D);  
  free(Values_Exact);
  free(Values_FE);
  FreeQuad3D(Quad3D);
  return sqrt(Error);
}


double ErrorEstimateEigen3DVecExactSolution(double *values,Fespace3D *fesp,Functions3D *Exact3D, 
		       int N_Derives,MultiIndex3D *multiindex3D,double *val)
{
   Fefunction3D *fefun3D;
   fefun3D = BuildFefunction3D(fesp,values);
   val[0] = ErrorEstimateEigen3DExactSolution(fefun3D,Exact3D,N_Derives,multiindex3D, 4);
}







double ErrorEstimateEigen3D(Fefunction3D *fefun3D,Functions3D *Exact3D, 
		       int N_Derives,MultiIndex3D *multiindex3D, int Quad_Points)
{
  printf("Compute the error estimate\n");
  int i,j,num_volu;
  MESH *mesh;
  mesh = fefun3D->Fespace->Mesh;
  num_volu = mesh->Num_Volus_Global;
  ELEMENT3D *elem3D;
  elem3D = malloc(1*sizeof(ELEMENT3D));
  InitialElem3D(mesh,elem3D);
  
  double Quad_xi, Quad_eta,Quad_zeta, Quad_W;
  int curr_volu, kk,ii;

  Quadrature3D *Quad3D;
  Quad3D = BuildQuad2Mesh3D(fefun3D->Fespace->Mesh, Quad_Points);
  
  Quad_Points = Quad3D->N_Points;
  
  double Error, Err_tmp, tmp, Err_K;
  Error = 0.0;
  //Err_tmp = 0.0;
  int dim_value;
  dim_value = fefun3D->Fespace->Base->Value_Dim;
  double *Values_FE, *Values_Exact;
  Values_FE = malloc(N_Derives*dim_value*sizeof(double));
  Values_Exact = malloc(N_Derives*dim_value*sizeof(double));
  double Quad_x,Quad_y,Quad_z;
 
  double average[10],coeffi;
  coeffi=0.0;
  int temp=0;
  for(i=0;i<10;++i)
    {
      GetElement3D(mesh,20*i+9,elem3D);
      Quad_xi = elem3D->Vert_X[0];
      Quad_eta = elem3D->Vert_Y[0];
      Quad_zeta = elem3D->Vert_Z[0];
      GetElementCoord(elem3D,Quad_xi,Quad_eta, Quad_zeta,&Quad_x,&Quad_y,&Quad_z);
      //求真解的值
      Exact3D(Quad_x,Quad_y,Quad_z,N_Derives*dim_value,Values_Exact);
      //求有限元的值
      GetFefunctValues(fefun3D,elem3D,Quad_x,Quad_y,Quad_z,
		       Quad_xi,Quad_eta,Quad_zeta,N_Derives,multiindex3D,Values_FE);
     if(Values_Exact[0]!=0)
     {
     average[i]=Values_FE[0]/Values_Exact[0];
     coeffi+=average[i];
     temp++;
     }

   }
  coeffi=coeffi/temp;
  printf ( "coeffi=%f\n",coeffi );

  for(curr_volu=0;curr_volu<num_volu;curr_volu++)
  {
    GetElement3D(mesh,curr_volu,elem3D);
    Err_K = 0.0;
    for(i=0;i<Quad_Points;i++)
    {
      Quad_xi = Quad3D->Quad_X[i];
      Quad_eta = Quad3D->Quad_Y[i];
      Quad_zeta = Quad3D->Quad_Z[i];
      Quad_W = Quad3D->Quad_W[i];
      GetElementCoord(elem3D,Quad_xi,Quad_eta, Quad_zeta,&Quad_x,&Quad_y,&Quad_z);
      //printf ( "%f,%f,%f,%f,%f,%f\n", Quad_xi, Quad_eta, Quad_zeta, Quad_x,Quad_y,Quad_z );
      //求真解的值
      Exact3D(Quad_x,Quad_y,Quad_z,N_Derives*dim_value,Values_Exact);
      //求有限元的值
      GetFefunctValues(fefun3D,elem3D,Quad_x,Quad_y,Quad_z,
		       Quad_xi,Quad_eta,Quad_zeta,N_Derives,multiindex3D,Values_FE);
      Err_tmp = 0.0;
      for(j=0;j<N_Derives*dim_value;j++)
      {
	tmp = Values_Exact[j]*coeffi-Values_FE[j];
//	printf("i=%d, j=%d, curr_volu=%d, tmp=%f,exact=%f,fevalue=%f\n",i,j,curr_volu,tmp,
//	       Values_Exact[j],Values_FE[j]);	
        Err_tmp += tmp*tmp;
      }//end for j
      Err_K += Err_tmp*Quad_W;
    }//end for i
    Error += Err_K*elem3D->Volu;   
    //printf ( "error=%d,%lf\n",curr_volu,Error );
  }//end for curr_volu  
  
  FreeElem3D(elem3D);  
  free(Values_Exact);
  free(Values_FE);
  FreeQuad3D(Quad3D);
  return sqrt(Error);
}
*/

double ErrorFefunction(Fefunction3D *fefunction1,Fefunction3D *fefunction2,int N_MultiIndex, MultiIndex3D *multiindex3D)
{
  //printf("Compute the error estimate\n");
  int i,j,num_volu;
  MESH *mesh;
  mesh = fefunction2->Fespace->Mesh;
  num_volu = mesh->Num_Volus_Global;
  ELEMENT3D *elem3D;
  elem3D = malloc(1*sizeof(ELEMENT3D));
  InitialElem3D(mesh,elem3D);
  
  double Quad_xi, Quad_eta,Quad_zeta, Quad_W;
  int curr_volu, kk,ii;

  Quadrature3D *Quad3D;
  Quad3D = BuildQuad2Mesh3D(fefunction2->Fespace->Mesh, 4);
  
  int Quad_Points = Quad3D->N_Points;
  double Error, Err_tmp, tmp, Err_K;
  Error = 0.0;
  //Err_tmp = 0.0;
  int dim_value;
  dim_value = fefunction2->Fespace->Base->Value_Dim;
  double *Values_FE1, *Values_FE2;
  Values_FE1 = malloc(N_MultiIndex*dim_value*sizeof(double));
  Values_FE2 = malloc(N_MultiIndex*dim_value*sizeof(double));
  double Quad_x,Quad_y,Quad_z;
  
  for(curr_volu=0;curr_volu<num_volu;curr_volu++)
  {
    GetElement3D(mesh,curr_volu,elem3D);
    Err_K = 0.0;
    for(i=0;i<Quad_Points;i++)
    {
      Quad_xi = Quad3D->Quad_X[i];
      Quad_eta = Quad3D->Quad_Y[i];
      Quad_zeta = Quad3D->Quad_Z[i];
      Quad_W = Quad3D->Quad_W[i];
      GetElementCoord(elem3D,Quad_xi,Quad_eta, Quad_zeta,&Quad_x,&Quad_y,&Quad_z);
      //求有限元的值
      GetFefunctValues(fefunction1,elem3D,Quad_x,Quad_y,Quad_z,
		       Quad_xi,Quad_eta,Quad_zeta,N_MultiIndex,multiindex3D,Values_FE1);
      GetFefunctValues(fefunction2,elem3D,Quad_x,Quad_y,Quad_z,
		       Quad_xi,Quad_eta,Quad_zeta,N_MultiIndex,multiindex3D,Values_FE2);
      Err_tmp = 0.0;
      for(j=0;j<N_MultiIndex;j++)
      {
	tmp = Values_FE1[j]-Values_FE2[j];
        Err_tmp += tmp*tmp;
      }//end for j
      Err_K += Err_tmp*Quad_W;
    }//end for i
    Error += Err_K*elem3D->Volu;   
  }//end for curr_volu  
  
  FreeElem3D(elem3D);  
  free(Values_FE1);
  free(Values_FE2);
  FreeQuad3D(Quad3D);
  return sqrt(Error);

}

/*
double ErrorFefunction1(Fefunction3D *fefunction1,Fefunction3D *fefunction2,int N_MultiIndex, MultiIndex3D *multiindex3D)
{
  //printf("Compute the error estimate\n");
  int i,j,num_volu;
  MESH *mesh;
  mesh = fefunction2->Fespace->Mesh;
  num_volu = mesh->Num_Volus_Global;
  ELEMENT3D *elem3D;
  elem3D = malloc(1*sizeof(ELEMENT3D));
  InitialElem3D(mesh,elem3D);
  
  double Quad_xi, Quad_eta,Quad_zeta, Quad_W;
  int curr_volu, kk,ii;

  Quadrature3D *Quad3D;
  Quad3D = BuildQuad2Mesh3D(fefunction2->Fespace->Mesh, 4);
  
  int Quad_Points = Quad3D->N_Points;
  double Error, Err_tmp, tmp, Err_K;
  Error = 0.0;
  //Err_tmp = 0.0;
  int dim_value;
  dim_value = fefunction2->Fespace->Base->Value_Dim;
  double *Values_FE1, *Values_FE2;
  Values_FE1 = malloc(N_MultiIndex*dim_value*sizeof(double));
  Values_FE2 = malloc(N_MultiIndex*dim_value*sizeof(double));
  double Quad_x,Quad_y,Quad_z;
  for(curr_volu=0;curr_volu<num_volu;curr_volu++)
  {
    GetElement3D(mesh,curr_volu,elem3D);
    Err_K = 0.0;
    for(i=0;i<Quad_Points;i++)
    {
      Quad_xi = Quad3D->Quad_X[i];
      Quad_eta = Quad3D->Quad_Y[i];
      Quad_zeta = Quad3D->Quad_Z[i];
      Quad_W = Quad3D->Quad_W[i];
      GetElementCoord(elem3D,Quad_xi,Quad_eta, Quad_zeta,&Quad_x,&Quad_y,&Quad_z);
      //求有限元的值
      GetFefunctValues(fefunction1,elem3D,Quad_x,Quad_y,Quad_z,
		       Quad_xi,Quad_eta,Quad_zeta,N_MultiIndex,multiindex3D,Values_FE1);
      GetFefunctValues(fefunction2,elem3D,Quad_x,Quad_y,Quad_z,
		       Quad_xi,Quad_eta,Quad_zeta,N_MultiIndex,multiindex3D,Values_FE2);
      Err_tmp = 0.0;
      for(j=0;j<N_MultiIndex;j++)
      {
	tmp = Values_FE1[j]+Values_FE2[j];
        Err_tmp += tmp*tmp;
      }//end for j
      Err_K += Err_tmp*Quad_W;
    }//end for i
    Error += Err_K*elem3D->Volu;   
  }//end for curr_volu  
  
  FreeElem3D(elem3D);  
  free(Values_FE1);
  free(Values_FE2);
  FreeQuad3D(Quad3D);
  return sqrt(Error);

}




void ErrorEstimateMul(MULTILEVEL *Multilevel,Fefunction3D **fefunction, double *error,int N_MultiIndex, MultiIndex3D *multiindex3D)//计算所有的误差
 {
   int i,j,k;
   double *tmp_Values;
   tmp_Values = malloc(sizeof(double)*fefunction[Multilevel->Num_Levels-1]->DOF_ALL);
   for(i=0;i<Multilevel->Num_Levels;++i)
    {  
       for(j=i;j<Multilevel->Num_Levels;++j)
        {
          MatrixDotVec(Multilevel->Stiff_Matrix[j],fefunction[i]->Values,tmp_Values);
	  for(k=0;k<fefunction[i]->DOF_ALL;++k)
           fefunction[i]->Values[k]=tmp_Values[k];
	}	   
    error[i]=ErrorFefunction(fefunction[i],fefunction[Multilevel->Num_Levels],N_MultiIndex,multiindex3D);
    }   
   free(tmp_Values);
 }


void ErrorEstimate4NumerSolution(MATRIX **matrix, double **fefunction, Fespace3D *fesp,int num,int N_MultiIndex,
                                 MultiIndex3D *multiindex3D,double *error)
{
  int i,j,k;
   double *tmp_values1,*tmp_values2;
   tmp_values1 = malloc(sizeof(double)*fesp->DOF_ALL);
   tmp_values2 = malloc(sizeof(double)*fesp->DOF_ALL);
   Fefunction3D *fefunc1,*fefunc2;
   fefunc1 = BuildFefunction3D(fesp,fefunction[num]);
   for(i=0;i<num;++i)
    { 
       for(k=0;k<matrix[i]->N_Columns;++k)
         tmp_values1[k]=fefunction[i][k];
       for(j=i;j<num;++j)
        {
          memset(tmp_values2,0.0,sizeof(double)*fesp->DOF_ALL);
          MatrixDotVec(matrix[j],tmp_values1,tmp_values2);
	  for(k=0;k<fesp->DOF_ALL;++k)
           tmp_values1[k]=tmp_values2[k];
	}
       fefunc2 =  BuildFefunction3D(fesp,tmp_values2);
    error[i]=ErrorFefunction(fefunc1,fefunc2,N_MultiIndex,multiindex3D);
    }  

   free(tmp_values1);
   free(tmp_values2);

 }



void ErrorEstimateMultiEigen(MATRIX **matrix, double **fefunction,Fespace3D *fesp,int num,int N_MultiIndex,
                             MultiIndex3D *multiindex3D,double *error)//计算所有的误差
 {
   printf ( "Compute the error!\n" );
   int i,j,k;
   double *tmp_values1,*tmp_values2;
   tmp_values1 = malloc(sizeof(double)*fesp->DOF_ALL);
   tmp_values2 = malloc(sizeof(double)*fesp->DOF_ALL);
   Fefunction3D *fefunc1,*fefunc2;
   fefunc1 = BuildFefunction3D(fesp,fefunction[num]);
   NormalizationUserfunc(fefunc1);
   for(i=0;i<num;++i)
    { 
       for(k=0;k<matrix[i]->N_Columns;++k)
         tmp_values1[k]=fefunction[i][k];
       for(j=i;j<num;++j)
        {
          memset(tmp_values2,0.0,sizeof(double)*fesp->DOF_ALL);
          MatrixDotVec(matrix[j],tmp_values1,tmp_values2);
	  for(k=0;k<fesp->DOF_ALL;++k)
           tmp_values1[k]=tmp_values2[k];
	}
       fefunc2 =  BuildFefunction3D(fesp,tmp_values2);
       NormalizationUserfunc(fefunc2);
    error[i]=ErrorFefunction(fefunc1,fefunc2,N_MultiIndex,multiindex3D);
    }  
   //fefunc2 =  BuildFefunction3D(fesp,fefunction[num]);
   //  Normalization(fefunc2);
   //error[num] = ErrorFefunction(fefunc1,fefunc2,N_MultiIndex,multiindex3D);

   free(tmp_values1);
   free(tmp_values2);
 }

void ErrorEstimateMultiEigen1(MATRIX **matrix, double **fefunction,Fespace3D *fesp,int num,int N_MultiIndex,
                             MultiIndex3D *multiindex3D,double *error)//计算所有的误差
 {
   printf ( "Compute the error!\n" );
   int i,j,k;
   double *tmp_values1,*tmp_values2;
   tmp_values1 = malloc(sizeof(double)*fesp->DOF_ALL);
   tmp_values2 = malloc(sizeof(double)*fesp->DOF_ALL);
   Fefunction3D *fefunc1,*fefunc2;
   fefunc1 = BuildFefunction3D(fesp,fefunction[num]);
   NormalizationUserfunc(fefunc1);
   for(i=0;i<num;++i)
    { 
       for(k=0;k<matrix[i]->N_Columns;++k)
         tmp_values1[k]=fefunction[i][k];
       for(j=i;j<num;++j)
        {
          memset(tmp_values2,0.0,sizeof(double)*fesp->DOF_ALL);
          MatrixDotVec(matrix[j],tmp_values1,tmp_values2);
	  for(k=0;k<fesp->DOF_ALL;++k)
           tmp_values1[k]=tmp_values2[k];
	}
       fefunc2 =  BuildFefunction3D(fesp,tmp_values2);
       NormalizationUserfunc(fefunc2);
    error[i]=ErrorFefunction1(fefunc1,fefunc2,N_MultiIndex,multiindex3D);
    }  

   //fefunc2 =  BuildFefunction3D(fesp,fefunction[num]);
   //  Normalization(fefunc2);
   //error[num] = ErrorFefunction1(fefunc1,fefunc2,N_MultiIndex,multiindex3D);

   free(tmp_values1);
   free(tmp_values2);
 }
*/
double DifferenceEstimate3D(Fefunction3D *fefun1_3D,Fefunction3D *fefun2_D3,
			    int N_Derives,MultiIndex3D *multiindex1_3D, 
			    MultiIndex3D *multiindex2_3D, int Quad_Points)
{
  //compute the error 
  int i,j,num_volus;
  MESH *mesh;
  mesh = fefun1_3D->Fespace->Mesh;
  
  num_volus = mesh->Num_Volus_Global;
  
  ELEMENT3D *elem3D;
  elem3D = malloc(1*sizeof(ELEMENT3D));
  // initializing the element
  InitialElem3D(mesh,elem3D);
  
  double Quad_xi, Quad_eta, Quad_zeta,Quad_W;
  int curr_volu, kk,ii;

  Quadrature3D *Quad3D;
  Quad3D = BuildQuad2Mesh3D(fefun1_3D->Fespace->Mesh, Quad_Points);
  
  Quad_Points = Quad3D->N_Points;
  
  double Error, Err_tmp, tmp, Err_K;
  Error = 0.0;
  Err_tmp = 0.0;
  
  double *Values1_FE, *Values2_FE;
  Values1_FE = malloc(N_Derives*sizeof(double));
  Values2_FE = malloc(N_Derives*sizeof(double));
  
  double Quad_x,Quad_y,Quad_z;
    
  
  for(curr_volu=0;curr_volu<num_volus;curr_volu++)
  {
    GetElement3D(mesh,curr_volu,elem3D);
    
    Err_K = 0.0;
    for(i=0;i<Quad_Points;i++)
    {
      //get the quadrature information 
      Quad_xi = Quad3D->Quad_X[i];
      Quad_eta = Quad3D->Quad_Y[i];
      Quad_zeta = Quad3D->Quad_Z[i];
      Quad_W = Quad3D->Quad_W[i];
      
      GetElementCoord(elem3D,Quad_xi,Quad_eta, Quad_zeta,&Quad_x,&Quad_y,&Quad_z);
      
      //求有限元的值
      GetFefunctValues(fefun1_3D,elem3D,Quad_x,Quad_y,Quad_z,
			      Quad_xi,Quad_eta,Quad_zeta,N_Derives,multiindex1_3D,Values1_FE);
      
      GetFefunctValues(fefun2_D3,elem3D,Quad_x,Quad_y,Quad_z,

			      Quad_xi,Quad_eta,Quad_zeta,N_Derives,multiindex2_3D,Values2_FE);
      Err_tmp = 0.0;
      for(j=0;j<N_Derives;j++)
      {
	tmp = Values1_FE[j]-Values2_FE[j];
	Err_tmp += tmp*tmp;
      }//end for j
      //printf("Quad_W=%18.15f, Err_tmp=%18.15f\n",Quad_W,Err_tmp);
      Err_K += Err_tmp*Quad_W;
      
    }//end for (i=0;i<Quad_Points;i++)
    //printf("Area=%f, Err_K=%f\n",elem2D->Area,Err_K *elem2D->Area);
    
    Error += Err_K *elem3D->Volu;   
    //printf("curr_face=%d, Error=%")
  }//end for curr_face  
  
  FreeElem3D(elem3D);  
  free(Values1_FE);
  free(Values2_FE);
  FreeQuad3D(Quad3D);
  
  return sqrt(Error);
}

void FreeFefunction3D(Fefunction3D *funct)//?????????????????????????????????????????????????????????????????
{
   printf ( "begin free\n" );
   if(funct->Values)
     free(funct->Values);
   if(funct->Fespace)
     free(funct->Fespace);
  free(funct);
   printf ( "end free\n" );
}




/* 有限元函数归一化 */
void Normalization(Fefunction3D *fefunction)
{
  int i,j,num_volu,dof_all;
  MESH *mesh;
  mesh = fefunction->Fespace->Mesh;
  num_volu = mesh->Num_Volus_Global;
  dof_all = fefunction->Fespace->DOF_ALL;
  ELEMENT3D *elem3D;
  elem3D = malloc(1*sizeof(ELEMENT3D));
  InitialElem3D(mesh,elem3D);
  double Quad_xi, Quad_eta,Quad_zeta, Quad_W;
  int curr_volu;

  Quadrature3D *Quad3D;
  Quad3D = BuildQuad2Mesh3D(fefunction->Fespace->Mesh,4);
  int Quad_Points = Quad3D->N_Points;
  
  double Norm,Norm_tmp,tmp,Norm_K;
  Norm = 0.0;
  int dim_value;
  dim_value = fefunction->Fespace->Base->Value_Dim;
  double *Values_FE;
  Values_FE = malloc(dim_value*sizeof(double));
  double Quad_x,Quad_y,Quad_z;
  for(curr_volu=0;curr_volu<num_volu;curr_volu++)
  {
    GetElement3D(mesh,curr_volu,elem3D);
    Norm_K = 0.0;
    for(i=0;i<Quad_Points;i++)
    {
      Quad_xi = Quad3D->Quad_X[i];
      Quad_eta = Quad3D->Quad_Y[i];
      Quad_zeta = Quad3D->Quad_Z[i];
      Quad_W = Quad3D->Quad_W[i];
      GetElementCoord(elem3D,Quad_xi,Quad_eta, Quad_zeta,&Quad_x,&Quad_y,&Quad_z);
      //求有限元的值
      Values_FE[0]=GetFefunctValue(fefunction,elem3D,Quad_x,Quad_y,Quad_z,
		       Quad_xi,Quad_eta,Quad_zeta,D000);
      Norm_tmp = 0.0;
      for(j=0;j<dim_value;j++)
      {
        tmp =Values_FE[j];
        Norm_tmp += tmp*tmp;
      }//end for j
      Norm_K += Norm_tmp*Quad_W;
    }//end for i
    Norm += Norm_K*elem3D->Volu;   
  }//end for curr_volu  
  tmp=sqrt(Norm);
  Norm=tmp;
  for(i=0;i<dof_all;++i)
  {
   fefunction->Values[i]=fefunction->Values[i]/Norm;
  }
  FreeElem3D(elem3D);  
  free(Values_FE);
  FreeQuad3D(Quad3D);
 }


/* 有限元函数归一化 */
void Normalization1(Fefunction3D *fefunction)
{
  int i,j,num_volu,dof_all;
  MESH *mesh;
  mesh = fefunction->Fespace->Mesh;
  num_volu = mesh->Num_Volus_Global;
  dof_all = fefunction->Fespace->DOF_ALL;
  ELEMENT3D *elem3D;
  elem3D = malloc(1*sizeof(ELEMENT3D));
  InitialElem3D(mesh,elem3D);
  double Quad_xi, Quad_eta,Quad_zeta, Quad_W;
  int curr_volu;

  Quadrature3D *Quad3D;
  Quad3D = BuildQuad2Mesh3D(fefunction->Fespace->Mesh,4);
  int Quad_Points = Quad3D->N_Points;
  
  double Norm;
  Norm = 0.0;
  int dim_value;
  dim_value = fefunction->Fespace->Base->Value_Dim;
  double *Values_FE;
  Values_FE = malloc(dim_value*sizeof(double));
  double Quad_x,Quad_y,Quad_z;
  for(curr_volu=0;curr_volu<num_volu;curr_volu++)
  {
    GetElement3D(mesh,curr_volu,elem3D);
    for(i=0;i<Quad_Points;i++)
    {
      Quad_xi = Quad3D->Quad_X[i];
      Quad_eta = Quad3D->Quad_Y[i];
      Quad_zeta = Quad3D->Quad_Z[i];
      Quad_W = Quad3D->Quad_W[i];
      GetElementCoord(elem3D,Quad_xi,Quad_eta, Quad_zeta,&Quad_x,&Quad_y,&Quad_z);
      //求有限元的值
      Values_FE[0]=GetFefunctValue(fefunction,elem3D,Quad_x,Quad_y,Quad_z,
		       Quad_xi,Quad_eta,Quad_zeta,D000);
      if(fabs(Norm)<fabs(Values_FE[0]))
            Norm = fabs(Values_FE[0]);
    }//end for i
  }//end for curr_volu  
  for(i=0;i<dof_all;++i)
  {
   fefunction->Values[i]=fefunction->Values[i]/Norm;
  }
  FreeElem3D(elem3D);  
  free(Values_FE);
  FreeQuad3D(Quad3D);
 }

/* 有限元函数归一化 */
void NormalizationVector1(int dof, double *values)
{
  int i;
  double tmp;
  tmp = 0.0;
  for(i=0;i<dof;i++)
  {
      if(fabs(tmp)<fabs(values[i]))
            tmp = fabs(values[i]);
  }//end for i
  for(i=0;i<dof;++i)
  {
   values[i]=values[i]/tmp;
  }
}






void NormalizationUserfunc(Fefunction3D *fefunction)
{
  int i,j,num_volu,dof_all;
  MESH *mesh;
  mesh = fefunction->Fespace->Mesh;
  num_volu = mesh->Num_Volus_Global;
  dof_all = fefunction->Fespace->DOF_ALL;
  ELEMENT3D *elem3D;
  elem3D = malloc(1*sizeof(ELEMENT3D));
  InitialElem3D(mesh,elem3D);
  double Quad_xi, Quad_eta,Quad_zeta, Quad_W;
  int curr_volu;

  Quadrature3D *Quad3D;
  Quad3D = BuildQuad2Mesh3D(fefunction->Fespace->Mesh,4);
  int Quad_Points = Quad3D->N_Points;
  
  double Norm,Norm_tmp,tmp,Norm_K;
  Norm = 0.0;
  int dim_value;
  dim_value = fefunction->Fespace->Base->Value_Dim;
  double *Values_FE;
  Values_FE = malloc(dim_value*sizeof(double));
  double Quad_x,Quad_y,Quad_z;
  for(curr_volu=0;curr_volu<num_volu;curr_volu++)
  {
    GetElement3D(mesh,curr_volu,elem3D);
    Norm_K = 0.0;
    for(i=0;i<Quad_Points;i++)
    {
      Quad_xi = Quad3D->Quad_X[i];
      Quad_eta = Quad3D->Quad_Y[i];
      Quad_zeta = Quad3D->Quad_Z[i];
      Quad_W = Quad3D->Quad_W[i];
      GetElementCoord(elem3D,Quad_xi,Quad_eta, Quad_zeta,&Quad_x,&Quad_y,&Quad_z);
      //求有限元的值
      Values_FE[0]=GetFefunctValue(fefunction,elem3D,Quad_x,Quad_y,Quad_z,
		       Quad_xi,Quad_eta,Quad_zeta,D000);
      Norm_tmp = 0.0;
      for(j=0;j<dim_value;j++)
      {
        tmp =Values_FE[j];
        Norm_tmp += tmp*tmp;
      }//end for j
      Norm_K += Norm_tmp*Quad_W;
    }//end for i
    Norm += Norm_K*elem3D->Volu;   
  }//end for curr_volu  
  tmp=sqrt(Norm);
  Norm=tmp;
  for(i=0;i<dof_all;++i)
  {
   fefunction->Values[i]=fefunction->Values[i]/Norm;
  }
  FreeElem3D(elem3D);  
  free(Values_FE);
  FreeQuad3D(Quad3D);
 }


//求函数的积分
double*  QuadFefunction(Fefunction3D *fefunction,double *val)
{
  int i,j,num_volu,dof_all;
  MESH *mesh;
  mesh = fefunction->Fespace->Mesh;
  num_volu = mesh->Num_Volus_Global;
  dof_all = fefunction->Fespace->DOF_ALL;
  ELEMENT3D *elem3D;
  elem3D = malloc(1*sizeof(ELEMENT3D));
  InitialElem3D(mesh,elem3D);
  double Quad_xi, Quad_eta,Quad_zeta, Quad_W;
  int curr_volu;

  Quadrature3D *Quad3D;
  Quad3D = BuildQuad2Mesh3D(fefunction->Fespace->Mesh,4);
  int Quad_Points = Quad3D->N_Points;
  
  double Norm,Norm_tmp,tmp,Norm_K;
  Norm = 0.0;
  int dim_value;
  dim_value = fefunction->Fespace->Base->Value_Dim;
  double *Values_FE;
  Values_FE = malloc(dim_value*sizeof(double));
  double Quad_x,Quad_y,Quad_z;
  for(curr_volu=0;curr_volu<num_volu;curr_volu++)
  {
    GetElement3D(mesh,curr_volu,elem3D);
    Norm_K = 0.0;
    for(i=0;i<Quad_Points;i++)
    {
      Quad_xi = Quad3D->Quad_X[i];
      Quad_eta = Quad3D->Quad_Y[i];
      Quad_zeta = Quad3D->Quad_Z[i];
      Quad_W = Quad3D->Quad_W[i];
      GetElementCoord(elem3D,Quad_xi,Quad_eta, Quad_zeta,&Quad_x,&Quad_y,&Quad_z);
      //求有限元的值
      Values_FE[0]=GetFefunctValue(fefunction,elem3D,Quad_x,Quad_y,Quad_z,
		       Quad_xi,Quad_eta,Quad_zeta,D000);
      Norm_tmp = 0.0;
      for(j=0;j<dim_value;j++)
      {
        tmp =Values_FE[j];
        Norm_tmp += tmp*tmp;
      }//end for j
      Norm_K += Norm_tmp*Quad_W;
    }//end for i
    Norm += Norm_K*elem3D->Volu;   
  }//end for curr_volu  
  tmp=sqrt(Norm);
  Norm=tmp;
  FreeElem3D(elem3D);  
  free(Values_FE);
  FreeQuad3D(Quad3D);
 //return Norm;
 val[0]=Norm;
 }



 void WriteSolution(double *solution, int num, char *file)
 {
  FILE *fp = fopen(file, "w");
  int i;
  for(i=0;i<num;++i) 
   fprintf(fp,"%f",solution[i]);
 fclose(fp);
 } 




void GetMaxValue(double *vector,int length, double *norm)
{
  int i;
  double tmp = 0.0;
  for(i=0;i<length;++i)
     {
       if(fabs(tmp)<fabs(vector[i]))
            tmp=vector[i];    
     }
   norm[0] = tmp;
}











//求一个向量的积分值
void Quad_Vector(Fespace3D *fespace,double *vector, double *val)
{
  Fefunction3D *fefunction;
  fefunction = BuildFefunction3D(fespace,vector);
  QuadFefunction(fefunction,val);
}

//对一个向量规范化
void Normalization_Vector(Fespace3D *fespace,double *vector)
{
 double *quad;
 int i;
 quad = malloc(sizeof(double));
 Quad_Vector(fespace,vector,quad);
 for(i=0;i<fespace->DOF_ALL;++i)
   vector[i] = vector[i]/quad[0];

 free(quad);
}









/*
void OutPutFefunction3D(Fefunction3D * Fefunct,char *filename)
{

  
  //OutPutMesh(Fefunct->Fespace->Mesh,filename);
  
  FILE *file; 
  file = fopen(filename,"w");
  if(!file)
  {
    printf("\ncannot open the solution.dat!\n");
    exit(0);
    
  }
  printf("open the output file successfully!\n");
  
  int num_verts,num_faces,n_face_vert;
    int i,j,k;
    char res[100];
    VERT *vert;
    FACEE *face;
    LINE *line;
    MESH *mesh;
    mesh=Fefunct->Fespace->Mesh;
    num_verts = mesh->Num_Verts;
    num_faces = mesh->Num_Faces;
    sprintf(res, "# The mesh information #\n");
    fputs(res, file);
    sprintf(res, "# the vertices information #\n");
    fputs(res, file);
    sprintf(res, "%6d\n", num_verts);
    fputs(res, file);
    for(i=0;i<num_verts;i++)
    {
        vert = mesh->Verts[i];
        sprintf(res,"%6d %18.15e %18.15e %6d\n",vert->Num,vert->Coord[0],
        vert->Coord[1],vert->ID_Boundary);
	fputs(res, file);
    }
    sprintf(res,"\n");
    fputs(res, file);
    sprintf(res,"# elements informations #\n");
    fputs(res, file);
    sprintf(res,"%6d\n",num_faces);
    fputs(res, file);
    for(i=0;i<num_faces;i++)
    {
        face = mesh->Faces[i];
        sprintf(res,"%6d %6d %6d %6d\n",face->Num,face->Verts[0],
                 face->Verts[1],face->Verts[2]);
	fputs(res, file);
    }
    
  int dof_all;
  //char res[100];
  
  dof_all = Fefunct->DOF_ALL;
  
   sprintf(res,"\n");
   fputs(res, file);
    
  sprintf(res, "# the solution information #\n");
  fputs(res, file);  
  sprintf(res,"%6d\n", dof_all);
  fputs(res, file);
  
  for(i=0;i<dof_all;i++)
  {
    sprintf(res,"%18.15e\n",Fefunct->Values[i]);
    fputs(res, file);    
  }
  fclose(file);
  
}
*/
