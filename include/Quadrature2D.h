/*
 * =====================================================================================
 *
 *       Filename:  Quadrature2D.h
 *      
 *      Description:  
 *
 *        Version:  1.0
 *        Created:  2014年10月26日 16时04分03秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Xu Fei (Fly), xufei@lsec.cc.ac.cn
 *        Company:  
 *
 * =====================================================================================
 */
#ifndef __QUADRATURE2D__
#define __QUADRATURE2D__
#include "Mesh3D.h"
#include "Enumerations.h"

typedef struct Quadrature2D_ {
  QuadType2D QuadType;
  int Order;
  int N_Points;
  double *Quad_X;
  double *Quad_Y;
  double *Quad_W;
} Quadrature2D;

void BuildQuad2Face2D(FACEE *face,int N_Points, Quadrature2D *quad2d);
Quadrature2D *BuildQuad2Mesh2D(MESH *mesh, int N_Points);
void Produce_Quad2D(int N_Lines, int N_Points, Quadrature2D *quad2d);
void OutPutQuad2D(Quadrature2D *quad);
void SetValueQuad2D(int N_Points,double *Quad_X,double *Quad_Y,double *Quad_W, Quadrature2D *quad2d);
void FreeQuad2D(Quadrature2D *quad);
#endif
