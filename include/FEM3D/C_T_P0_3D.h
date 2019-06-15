/*
 * =====================================================================================
 *
 *       Filename:  C_T_P0_3D.h
 *
 *    Description:  define the base functions in 3D case
 *
 *        Version:  1.0
 *        Created:  2017/09/16 10时44分41秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  xufei@lsec.cc.ac.cn
 *        Company:  
 *
 * =====================================================================================
 */
#ifndef __CTP03D__
#define __CTP03D__

#include <stdbool.h>
#include "Mesh3D.h"
#include "Enumerations.h"
#include "Constants.h"
// ***********************************************************************
// P0 element, conforming, 3D
// ***********************************************************************

int C_T_P0_3D_dof[4] = {0,0,0,1};
double C_T_P0_3D_nodal_points[3] = {0.167,0.167,0.167}; 
int C_T_P0_3D_Num_Bas = 1;
int C_T_P0_3D_Value_Dim =1;
int C_T_P0_3D_Polydeg =1;
int C_T_P0_3D_Accuracy = 1;
bool C_T_P0_3D_Maptype = Affine;
// base function values
static void C_T_P0_3D_D000(ELEMENT3D *elem, double x,double y,double z,double xi, double eta,double zeta,
                          double *values)
{
  values[0]=1;
}

// values of the derivatives in xi direction
static void C_T_P0_3D_D100(ELEMENT3D *elem, double x,double y,double z,double xi, double eta,double zeta,
                            double *values)
{
  values[0]=0;
}

// values of the derivatives in eta direction
static void C_T_P0_3D_D010(ELEMENT3D *elem, double x,double y, double z,double xi, double eta,double zeta,
                           double *values)
{
  values[0]=0; 
}
// values of the derivatives in zeta direction
static void C_T_P0_3D_D001(ELEMENT3D *elem, double x,double y, double z,double xi, double eta, double zeta,
                            double *values)
{
  values[0]=0;  
}

// values of the derivatives in xi-xi  direction
static void C_T_P0_3D_D200(ELEMENT3D *elem, double x,double y, double z,double xi, double eta,double zeta,
                            double *values)
{
  values[0]=0;
}
// values of the derivatives in xi-xi  direction
static void C_T_P0_3D_D020(ELEMENT3D *elem, double x,double y, double z,double xi, double eta,double zeta,
                            double *values)
{
  values[0]=0;
}
// values of the derivatives in xi-xi  direction
static void C_T_P0_3D_D002(ELEMENT3D *elem, double x,double y, double z,double xi, double eta,double zeta,
                            double *values)
{
  values[0]=0;
}
// values of the derivatives in xi-xi  direction
static void C_T_P0_3D_D110(ELEMENT3D *elem, double x,double y, double z,double xi, double eta,double zeta,
                            double *values)
{
  values[0]=0;
}
// values of the derivatives in xi-xi  direction
static void C_T_P0_3D_D101(ELEMENT3D *elem, double x,double y, double z,double xi, double eta,double zeta,
                            double *values)
{
  values[0]=0;
}
// values of the derivatives in xi-xi  direction
static void C_T_P0_3D_D011(ELEMENT3D *elem, double x,double y, double z,double xi, double eta,double zeta,
                            double *values)
{
  values[0]=0;
}
static void C_T_P0_3D_Nodal(ELEMENT3D *elem,Functions3D *fun, double *values)
{
  int dim = 1;  
  fun(elem->Vert_X[0],elem->Vert_Y[0], elem->Vert_Z[0], dim,values);
}
#endif
