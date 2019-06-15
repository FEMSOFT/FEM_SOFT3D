/*
 * =====================================================================================
 *
 *       Filename:  C_T_P2_3D.h
 *
 *    Description:  define the base functions in 2D case
 *
 *        Version:  1.0
 *        Created:  2014/10/26 14时39分41秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:   xufei@lsec.cc.ac.cn
 *        Company:  
 *
 * =====================================================================================
 */
#ifndef __CTP23D__
#define __CTP23D__

#include <stdbool.h>
#include "Mesh3D.h"
#include "Enumerations.h"
#include "Constants.h"
// ***********************************************************************
// P2 element, conforming, 3D
// ***********************************************************************

int C_T_P2_3D_dof[4] = {1,1,0,0};
double C_T_P2_3D_nodal_points[30] = {0.0,0.0,0.0,  1.0,0.0,0.0,  0.0,1.0,0.0, 0.0,0.0,1.0,
                                     0.5,0.0,0.0,  0.0,0.5,0.0,  0.0,0.0,0.5,
				     0.5,0.5,0.0,  0.5,0.0,0.5,  0.0,0.5,0.5}; 
int C_T_P2_3D_Num_Bas = 10;
int C_T_P2_3D_Value_Dim =1;
int C_T_P2_3D_Polydeg =2;
int C_T_P2_3D_Accuracy = 2;
bool C_T_P2_3D_Maptype = Affine;

// base function values
static void C_T_P2_3D_D000(ELEMENT3D *elem, double x,double y,double z,double xi, double eta,double zeta,
                           double *values)
{ 
  values[0] = 1+2*xi*xi+4*xi*eta+4*xi*zeta-3*xi+2*eta*eta+4*eta*zeta-3*eta+2*zeta*zeta-3*zeta; 
  values[1] = -xi+2*xi*xi;
  values[2] = -eta+2*eta*eta;
  values[3] = -zeta+2*zeta*zeta;
  values[4] = 4*xi-4*xi*zeta-4*xi*eta-4*xi*xi;
  values[5] = 4*eta-4*eta*zeta-4*xi*eta-4*eta*eta;
  values[6] = 4*zeta-4*zeta*zeta-4*xi*zeta-4*eta*zeta;
  values[7] = 4*xi*eta;
  values[8] = 4*xi*zeta;
  values[9] = 4*eta*zeta;
}

// values of the derivatives in xi direction
static void C_T_P2_3D_D100(ELEMENT3D *elem, double x,double y,double z,double xi, double eta,double zeta,
                           double *values)
{ 
  values[0] = -3+4*xi+4*eta+4*zeta; 
  values[1] = -1+4*xi;
  values[2] =  0;
  values[3] =  0;
  values[4] =  4-4*zeta-4*eta-8*xi;
  values[5] = -4*eta;
  values[6] = -4*zeta;
  values[7] =  4*eta;
  values[8] =  4*zeta;
  values[9] =  0;
}
// values of the derivatives in eta direction
static void C_T_P2_3D_D010(ELEMENT3D *elem, double x,double y,double z,double xi, double eta,double zeta,
                           double *values)
{ 
  values[0] = -3+4*xi+4*eta+4*zeta; 
  values[1] =  0;
  values[2] = -1+4*eta;
  values[3] =  0;
  values[4] = -4*xi;
  values[5] =  4-4*zeta-4*xi-8*eta;
  values[6] = -4*zeta;
  values[7] =  4*xi;
  values[8] =  0;
  values[9] =  4*zeta;
}
// values of the derivatives in zeta direction
static void C_T_P2_3D_D001(ELEMENT3D *elem, double x,double y,double z,double xi, double eta,double zeta,
                           double *values)
{ 
  values[0] = -3+4*xi+4*eta+4*zeta; 
  values[1] =  0;
  values[2] =  0;
  values[3] = -1+4*zeta;
  values[4] = -4*xi;
  values[5] = -4*eta;
  values[6] =  4-4*xi-4*eta-8*zeta;
  values[7] =  0;
  values[8] =  4*xi;
  values[9] =  4*eta;
}
// values of the derivatives in xi direction
static void C_T_P2_3D_D200(ELEMENT3D *elem, double x,double y,double z,double xi, double eta,double zeta,
                           double *values)
{ 
  values[0] =  4; 
  values[1] =  4;
  values[2] =  0;
  values[3] =  0;
  values[4] = -8;
  values[5] =  0;
  values[6] =  0;
  values[7] =  0;
  values[8] =  0;
  values[9] =  0;
}


// values of the derivatives in xi direction
static void C_T_P2_3D_D020(ELEMENT3D *elem, double x,double y,double z,double xi, double eta,double zeta,
                           double *values)
{ 
  values[0] =  4; 
  values[1] =  0;
  values[2] =  4;
  values[3] =  0;
  values[4] =  0;
  values[5] = -8;
  values[6] =  0;
  values[7] =  0;
  values[8] =  0;
  values[9] =  0;
}

// values of the derivatives in xi direction
static void C_T_P2_3D_D002(ELEMENT3D *elem, double x,double y,double z,double xi, double eta,double zeta,
                           double *values)
{ 
  values[0] =  4; 
  values[1] =  0;
  values[2] =  0;
  values[3] =  4;
  values[4] =  0;
  values[5] =  0;
  values[6] = -8;
  values[7] =  0;
  values[8] =  0;
  values[9] =  0;
}

// values of the derivatives in xi direction
static void C_T_P2_3D_D110(ELEMENT3D *elem, double x,double y,double z,double xi, double eta,double zeta,
                           double *values)
{ 
  values[0] =  4; 
  values[1] =  0;
  values[2] =  0;
  values[3] =  0;
  values[4] = -4;
  values[5] = -4;
  values[6] =  0;
  values[7] =  4;
  values[8] =  0;
  values[9] =  0;
}

// values of the derivatives in xi direction
static void C_T_P2_3D_D101(ELEMENT3D *elem, double x,double y,double z,double xi, double eta,double zeta,
                           double *values)
{ 
  values[0] =  4; 
  values[1] =  0;
  values[2] =  0;
  values[3] =  0;
  values[4] = -4;
  values[5] =  0;
  values[6] = -4;
  values[7] =  0;
  values[8] =  4;
  values[9] =  0;
}
// values of the derivatives in eta direction
static void C_T_P2_3D_D011(ELEMENT3D *elem, double x,double y,double z,double xi, double eta,double zeta,
                           double *values)
{ 
  values[0] =  4; 
  values[1] =  0;
  values[2] =  0;
  values[3] =  0;
  values[4] =  0;
  values[5] = -4;
  values[6] = -4;
  values[7] =  0;
  values[8] =  0;
  values[9] =  4;
}



// values of the derivatives in eta-eta direction
// static void C_T_P2_2D_Nodal(ELEMENT2D *elem, Function2D *fun, double *values)
// {
//   // dof on the verts   
//   values[0] = fun(elem->Vert_X[0],elem->Vert_Y[0]);
//   values[1] = fun(elem->Vert_X[1],elem->Vert_Y[1]);
//   values[2] = fun(elem->Vert_X[2],elem->Vert_Y[2]);
//   
//   //dof on the midpoints of the edges
//   values[3] = fun(0.5*(elem->Vert_X[1]+elem->Vert_X[2]),0.5*(elem->Vert_Y[1]+elem->Vert_Y[2]));
//   values[4] = fun(0.5*(elem->Vert_X[2]+elem->Vert_X[0]),0.5*(elem->Vert_Y[2]+elem->Vert_Y[0]));
//   values[5] = fun(0.5*(elem->Vert_X[0]+elem->Vert_X[1]),0.5*(elem->Vert_Y[0]+elem->Vert_Y[1]));
// }

// values of the derivatives in eta-eta direction
static void C_T_P2_3D_Nodal(ELEMENT3D *elem, Functions3D *fun, double *values)
{
  // dof on the verts   
  int dim =1;
 fun(elem->Vert_X[0],elem->Vert_Y[0],elem->Vert_Z[0],dim, values);
 fun(elem->Vert_X[1],elem->Vert_Y[1],elem->Vert_Z[1],dim, values+1);
 fun(elem->Vert_X[2],elem->Vert_Y[2],elem->Vert_Z[2],dim, values+2);
 fun(elem->Vert_X[3],elem->Vert_Y[3],elem->Vert_Z[3],dim, values+3);
  //dof on the midpoints of the edges
 fun(0.5*(elem->Vert_X[0]+elem->Vert_X[1]),0.5*(elem->Vert_Y[0]+elem->Vert_Y[1]),0.5*(elem->Vert_Z[0]+elem->Vert_Z[1]),dim,  values+4);
 fun(0.5*(elem->Vert_X[0]+elem->Vert_X[2]),0.5*(elem->Vert_Y[0]+elem->Vert_Y[2]),0.5*(elem->Vert_Z[0]+elem->Vert_Z[2]),dim,values+5);
 fun(0.5*(elem->Vert_X[0]+elem->Vert_X[3]),0.5*(elem->Vert_Y[0]+elem->Vert_Y[3]),0.5*(elem->Vert_Z[0]+elem->Vert_Z[3]),dim,values+6);
 fun(0.5*(elem->Vert_X[1]+elem->Vert_X[2]),0.5*(elem->Vert_Y[1]+elem->Vert_Y[2]),0.5*(elem->Vert_Z[1]+elem->Vert_Z[2]),dim,  values+7);
 fun(0.5*(elem->Vert_X[1]+elem->Vert_X[3]),0.5*(elem->Vert_Y[1]+elem->Vert_Y[3]),0.5*(elem->Vert_Z[1]+elem->Vert_Z[3]),dim,values+8);
 fun(0.5*(elem->Vert_X[2]+elem->Vert_X[3]),0.5*(elem->Vert_Y[2]+elem->Vert_Y[3]),0.5*(elem->Vert_Z[2]+elem->Vert_Z[3]),dim,values+9);
}
#endif
