/*
 * =====================================================================================
 *
 *       Filename:  Constants.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2014年10月20日 14时09分30秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Xu Fei (Fly), xufei@lsec.cc.ac.cn
 *        Company:  
 *
 * =====================================================================================
 */
#ifndef __CONSTANTS__
#define __CONSTANTS__
#include <stdbool.h>
#include "Mesh3D.h"
#include "Enumerations.h"
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/unistd.h>
#include <string.h>

#ifndef __PI__
#define __PI__
#define PI 3.14159265358979311599796346854
#endif
typedef double Function3D(double, double, double);//???????????????????????????????????
// the prototype of the function for 2D 
// (x,y,n_value,values)
typedef void Functions3D(double, double,double, int, double *);//????????????????????????????
typedef void Functionuser(double, double,double, double *, double *);
typedef void BaseFunct3D(ELEMENT3D *, double,double,double, double, double, double, double *);
// NodalFunct3D(FACEE *face,double *values);                            
typedef void NodalFunct3D(ELEMENT3D *,Functions3D *, double *);
//function for the discrete form for the matrix 
//typedef double DiscreteFormMatrix(double,double,double,int, double *,int, double*,int, double *,int, double*);
typedef double DiscreteFormMatrix(double,double,double,int, double *,int,double*,int, double*,int, double*, int, double*);
typedef double DiscreteFormMatrixFace(int,double,double, double*, double *, double, int, double *,int, double*,int, double *);
typedef double DiscreteFormMatrixLine(int,double,double, double*, double *, double, int, double *,int, double*,int, double *);
//typedef double DiscreteFormRHSVolu(double,double,double,int, double *,int, double *);
typedef double DiscreteFormRHSVolu(double,double,double,int, double *,int, double*,int, double*,int, double*);
typedef double DiscreteFormRHSFace(int ,double ,double ,double ,double *, double *, 
                       double , int , double *,int ,
		       double*);
typedef double DiscreteFormRHSLine(int,double,double,double *,double *, double, int, double *,int, double *);

//define the prototype of the boundary condition 
/**定义区域边界的边界类型的函数： int 表示区域边界的编号，double：表示在区域边界上的参数值 */
typedef void BoundCondFunct3D(int, double, BoundType *);

//define the prototype of the boundary values 
/**给出区域边界的边界函数： */
typedef void BoundValueFunct3D(int ID_Bound, double param, double x, double y, double *value);
typedef Functions3D* BoundaryFunction3D(int ID_Bound); //, Functions3D *boundfun);
double GetTime();
double GetMemory();
//对向量a的 a[left:right]进行排练, 数组长度right-left+1
void QuickSort_Int(int a[], int left, int right);
//对向量a的 a[left:right]进行排练, 数组长度right-left+1
void QuickSort_Double(double a[], int left, int right);

#endif
