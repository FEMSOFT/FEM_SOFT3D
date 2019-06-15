/*
 * =====================================================================================
 *
 *       Filename:  Quadrature1D.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2014年10月26日 16时03分21秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Xu Fei (Fly), xufei@lsec.cc.ac.cn
 *        Company:  
 *
 * =====================================================================================
 */
#ifndef __QUADRATURE1D__
#define __QUADRATURE1D__
#include "Enumerations.h"

typedef struct Quadrature1D_ {
  QuadType1D QuadType;
  int Order;
  int N_Points;
  double *Quad_X;
  double *Quad_W;
} Quadrature1D;

//void BuildQuad2Face2D(FACEE *face,int N_Points, Quadrature2D *quad2d);
Quadrature1D *BuildQuad1D(int N_Points);
void Produce_Quad1D(int N_Points, Quadrature1D *quad1d);
void OutPutQuad1D(Quadrature1D *quad);
void SetValueQuad1D(int N_Points,double *Quad_X,double *Quad_W, Quadrature1D *quad1d);
void FreeQuad1D(Quadrature1D *quad);
#endif

