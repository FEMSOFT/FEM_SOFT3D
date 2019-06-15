/*
 * =====================================================================================
 *
 *       Filename:  Matrix.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2014年10月20日 13时54分45秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Xu Fei (Fly), xufei@lsec.cc.ac.cn
 *        Company:  
 *
 * =====================================================================================
 */
#ifndef __MATRIX__
#define __MATRIX__
#include "Fefunction3D.h"
#include "Base3D.h"
#include "Enumerations.h"
#include "Constants.h"
#include "DiscreteForm3D.h"

typedef struct MATRIX_ {	
     ///**形成刚度矩阵对应的离散变分形式*/
    DISCRETEFORM3D *DiscreteForm3D;
    /** 下面是存储稀疏矩阵的行压缩形式 */
    /** number rof rows */
    int N_Rows;
    /** number columns */
    int N_Columns;
    /** number of matrix entries */
    int N_Entries;
    /** in which column is the current entry */
    int *KCol;
    /** index in KCol where each row starts */
    int *RowPtr;   
    /** matrix elements in an array */
    double *Entries;
}MATRIX;



/** 下面是一些形成刚度矩阵和对矩阵进行操作的函数*/
/** OutPut the matrix */
void OutPutMatrix(MATRIX *matrix);
/** 根据离散变分形式来初始化矩阵：期间要形成矩阵的结构*/
MATRIX * BuildMatrix(DISCRETEFORM3D *discreteform);
/** 依据试探有限元空间和检验有限元空间产生矩阵的结构 */
// Assemble the matrix 
void AssembleMatrix(MATRIX *Matrix);
//形成在单元上的单刚矩阵
void SubAssembleMatrix3D(int n_bas_test,BASEFUNCTION3D *test_Base,
    int n_bas_anasatz,BASEFUNCTION3D *anasatz_Base,
    int n_test_multiindex,MultiIndex3D *test_MultiIndex,
    int n_anasatz_multiindex,MultiIndex3D *anasatz_MultiIndex,              
    int N_AuxFefunct,Fefunction3D **AuxFefunct,int* N_AuxFefun_MultiIndex,
    MultiIndex3D *AuxFefun_MultiIndex,int N_AuxFefun_Values,double *AuxFefun_Values,int N_UserFunction,Functions3D **UserFunction,
    int N_FeConst,double *FeConst,DiscreteFormMatrix *discreteformvolu, ELEMENT3D *elem3D,double *Test_Values, 
    double *Anasatz_Values,double *Test_MultiIndex_Values,double *Anasatz_MultiIndex_Values,
    double Quad_xi,double Quad_eta, double Quad_zeta, double Quad_W, double *Mat_Entries);
//形成在边界上的单刚矩阵
void SubAssembleMatrix3DLine(int ID_bd,double *normal, double * tangent,double local_direction,int n_bas_test,BASEFUNCTION3D *test_Base,
             int n_bas_anasatz,BASEFUNCTION3D *anasatz_Base,
             int n_test_multiindex,MultiIndex3D *test_MultiIndex,
             int n_anasatz_multiindex,MultiIndex3D *anasatz_MultiIndex,              
             int N_AuxFefunct,Fefunction3D **AuxFefunct,int* N_AuxFefun_MultiIndex,
             MultiIndex3D *AuxFefun_MultiIndex,int N_AuxFefun_Values,double *AuxFefun_Values,
             DiscreteFormMatrixLine *discreteformline, ELEMENT3D *elem3D,double *Test_Values, 
	     double *Anasatz_Values,double *Test_MultiIndex_Values,double *Anasatz_MultiIndex_Values,
	     double Quad_xi,double Quad_eta,double Quad_zeta,double Quad_W, double length,double *Mat_Entries);

void AddSubMatrix(MATRIX *A,double *AA, int n_row,int *row,int n_col,int *col);
//void qqsort(int *s,int l,int r);
//void QuickSort(int *data, int low, int high);
void BubbleSort(int* Data,int low,int high);
int partition(int *data,int low,int high);
void swap(int *a,int *b);

/** 处理边界条件的处理，只针对Matrix中牵涉到的有限元空间进行处理*/
/** 处理椭圆问题的边界条件 */
/** 对于特征值问题的边界处理 */
//void BoundaryTreatment(MATRIX *Stiff_Matrix, double *RHS, BoundValueFunct3D *boundvaluefun);
void BoundaryTreatment(MATRIX *Stiff_Matrix, double *RHS, BoundaryFunction3D *boundaryfunction);
// 对于特征值问题的边界处理
/** 对于特征值问题的边界处理 */
void BoundaryTreatmentEigen(MATRIX *Stiff_Matrix, MATRIX *Mass_Matrix,BoundaryFunction3D *boundaryfunction);
//矩阵加法： A=A+tau*B
void MatrixAdd(MATRIX *A, MATRIX *B, double tau);
//MatVect(MATRIX *A, double *x, double *y);
/** print the matrix */
void Print(MATRIX *matrix);
/** write matrix into file */
int Write(MATRIX *matrix,const char *filename);
void FreeMatrix(MATRIX *Matrix);
//  for mixed element  matrix
void ReAssembleMatrix(MATRIX ***matrix,int nrows, int ncols, MATRIX *matrixb);
void WriteMatrixPS(MATRIX *matrix,const char *filename);
MATRIX* CopyMatrix( MATRIX *A );
MATRIX* EnlargeMatrix( MATRIX *A, int expan, double **vec, double **mat );
MATRIX* ExpandMatrixStruct( MATRIX *A, int expan );
void ExpandMatrix( MATRIX *A, int expan, double **vec, double **mat );

#endif
