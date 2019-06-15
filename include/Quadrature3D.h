/*
 * =====================================================================================
 *
 *       Filename:  Quadrature3D.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2014年10月26日 16时04分48秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Xu Fei (Fly), xufei@lsec.cc.ac.cn
 *        Company:  
 *
 * =====================================================================================
 */
#ifndef __QUADRATURE3D__
#define __QUADRATURE3D__
#include "Mesh3D.h"
#include "Enumerations.h"

typedef struct Quadrature3D_ {
  QuadType3D QuadType;
  int Order;
  int N_Points;
  double *Quad_X;
  double *Quad_Y;
  double *Quad_Z;
  double *Quad_W;
} Quadrature3D;

void BuildQuad2Volu3D(VOLU *volu,int N_Points, Quadrature3D *quad3d);
Quadrature3D *BuildQuad2Mesh3D(MESH *mesh, int N_Points);
void Produce_Quad3D(int N_Lines, int N_Points, Quadrature3D *quad3d);
void OutPutQuad(Quadrature3D *quad);
void SetValueQuad3D(int N_Points,double *Quad_X,double *Quad_Y,
                    double *Quad_Z,double *Quad_W, Quadrature3D *quad3d);
void FreeQuad3D(Quadrature3D *quad);
#endif
