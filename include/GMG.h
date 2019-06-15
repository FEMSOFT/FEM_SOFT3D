/*
 * =====================================================================================
 *
 *       Filename:  GMG.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2015年01月09日 02时46分11秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Xu Fei (Fly), xufei@lsec.cc.ac.cn
 *        Company:  
 *
 * =====================================================================================
 */
#ifndef __GMG__
#define __GMG__
#include "Mesh3D.h"
#include "Enumerations.h"
#include "Matrix.h"
#include "MG.h"
#include "Base3D.h"
#include "Constants.h"
void GMGSolver(MULTILEVEL *multilevel,double *rhs,double *solution,int level,int smooth,int GMGIter);
void GMGIteration(MULTILEVEL *multilevel,double *rhs,double *solution,int level,int smooth,double *return_Residual);
#endif
