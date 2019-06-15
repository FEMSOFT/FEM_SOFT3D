/*
 * =====================================================================================
 *
 *       Filename:  C_T_P1_3D.h
 *
 *    Description:  define the base functions in 3D case
 *
 *        Version:  1.0
 *        Created:  2014/10/26 10时44分41秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  xufei@lsec.cc.ac.cn
 *        Company:  
 *
 * =====================================================================================
 */
#ifndef __CTP13D__
#define __CTP13D__

#include <stdbool.h>
#include "Mesh3D.h"
#include "Enumerations.h"
#include "Constants.h"
// ***********************************************************************
// P1 element, conforming, 3D
// ***********************************************************************

int C_T_P1_3D_dof[4] = {1,0,0,0};
double C_T_P1_3D_nodal_points[12] = {0.0,0.0,0.0, 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0}; 
int C_T_P1_3D_Num_Bas = 4;
int C_T_P1_3D_Value_Dim =1;
int C_T_P1_3D_Polydeg =1;
int C_T_P1_3D_Accuracy = 1;
bool C_T_P1_3D_Maptype = Affine;
// base function values
static void C_T_P1_3D_D000(ELEMENT3D *elem, double x,double y,double z,double xi, double eta,double zeta,
                          double *values)
{
  values[0]=1-xi-eta-zeta;
  values[1]=xi;
  values[2]=eta;
  values[3]=zeta;
  
}

// values of the derivatives in xi direction
static void C_T_P1_3D_D100(ELEMENT3D *elem, double x,double y,double z,double xi, double eta,double zeta,
                            double *values)
{
  values[0]=-1;
  values[1]= 1;
  values[2]= 0;  
  values[3]= 0;
}

// values of the derivatives in eta direction
static void C_T_P1_3D_D010(ELEMENT3D *elem, double x,double y, double z,double xi, double eta,double zeta,
                           double *values)
{
  values[0]=-1;
  values[1]= 0;
  values[2]= 1;
  values[3]= 0;  
}
// values of the derivatives in zeta direction
static void C_T_P1_3D_D001(ELEMENT3D *elem, double x,double y, double z,double xi, double eta, double zeta,
                            double *values)
{
  values[0]=-1;
  values[1]= 0;
  values[2]= 0;
  values[3]= 1;  
}

// values of the derivatives in xi-xi  direction
static void C_T_P1_3D_D200(ELEMENT3D *elem, double x,double y, double z,double xi, double eta,double zeta,
                            double *values)
{
  values[0]=0;
  values[1]=0;
  values[2]=0;
  values[3]=0;
}
// values of the derivatives in xi-xi  direction
static void C_T_P1_3D_D020(ELEMENT3D *elem, double x,double y, double z,double xi, double eta,double zeta,
                            double *values)
{
  values[0]=0;
  values[1]=0;
  values[2]=0;
  values[3]=0;
}
// values of the derivatives in xi-xi  direction
static void C_T_P1_3D_D002(ELEMENT3D *elem, double x,double y, double z,double xi, double eta,double zeta,
                            double *values)
{
  values[0]=0;
  values[1]=0;
  values[2]=0;
  values[3]=0;
}
// values of the derivatives in xi-xi  direction
static void C_T_P1_3D_D110(ELEMENT3D *elem, double x,double y, double z,double xi, double eta,double zeta,
                            double *values)
{
  values[0]=0;
  values[1]=0;
  values[2]=0;
  values[3]=0;
}
// values of the derivatives in xi-xi  direction
static void C_T_P1_3D_D101(ELEMENT3D *elem, double x,double y, double z,double xi, double eta,double zeta,
                            double *values)
{
  values[0]=0;
  values[1]=0;
  values[2]=0;
  values[3]=0;
}
// values of the derivatives in xi-xi  direction
static void C_T_P1_3D_D011(ELEMENT3D *elem, double x,double y, double z,double xi, double eta,double zeta,
                            double *values)
{
  values[0]=0;
  values[1]=0;
  values[2]=0;
  values[3]=0;
}
static void C_T_P1_3D_Nodal(ELEMENT3D *elem,Functions3D *fun, double *values)
{
  int dim = 1;  
  fun(elem->Vert_X[0],elem->Vert_Y[0], elem->Vert_Z[0], dim,values);
  fun(elem->Vert_X[1],elem->Vert_Y[1], elem->Vert_Z[1], dim,values+1);
  fun(elem->Vert_X[2],elem->Vert_Y[2], elem->Vert_Z[2], dim,values+2);
  fun(elem->Vert_X[3],elem->Vert_Y[3], elem->Vert_Z[3], dim,values+3);

}
#endif
