/*
 * =====================================================================================
 *
 *       Filename:  Base3D.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2014年10月20日 13时31分59秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Xu Fei (Fly), xufei@lsec.cc.ac.cn
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef _BOUNDARY_
#define _BOUNDARY_
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Matrix.h"
#include "Fefunction3D.h"
void BoundaryTreatment(MATRIX *Stiff_Matrix, double *RHS, BoundaryFunction3D *boundaryfunction);
void BoundaryTreatmentEigen(MATRIX *Stiff_Matrix, MATRIX *Mass_Matrix, BoundaryFunction3D *boundaryfunction);
void BoundaryTreatmentEigenMul(MATRIX *Stiff_Matrix, MATRIX *Mass_Matrix, BoundaryFunction3D *boundaryfunction);
Functions3D* BoundFuncMul(int ID_Bound);
void BoundValueMul(double x, double y,double z, int n, double *value);
#endif








