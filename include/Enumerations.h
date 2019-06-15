/*
 * =====================================================================================
 *
 *       Filename:  Enumerations.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2014年10月20日 14时14分35秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Xu Fei (Fly), xufei@lsec.cc.ac.cn
 *        Company:  
 *
 * =====================================================================================
 */
#ifndef __ENUMERATIONS__
#define __ENUMERATIONS__

#include <stdbool.h>

//enum boolean {FALSE=0, TRUE=1};
#define TRUE (1)
#define FALSE (0)

#define N_TypeBound 4
typedef enum BoundType { DIRICHLETT=0, NEUMANNN=1, ROBIN=2, NOTBOUNDARY=3 } BoundType;

#define N_FemType  3
typedef enum FemType { C_T_P1_3D, C_T_P2_3D, C_T_P3_3D } FemType;

#define N_MultiIndices3D 10
typedef enum MultiIndex3D  { D000=0,D100,D010,D001,D110,D101,D011,D200,D020,D002 } MultiIndex3D;

#define N_MapType 3
typedef enum MapType {Affine=0, Piola, Actual} MapType; 

#define N_QuadType1D 4
typedef enum QuadType1D { quad_Line1, quad_Line2,quad_Line3,quad_Line4 } QuadType1D;

#define N_QuadType2D 5
typedef enum QuadType2D {quad_Triangle1, quad_Triangle4,quad_Triangle7, quad_Triangle13,quad_Triangle27} QuadType2D;             

#define N_QuadType3D 4
typedef enum QuadType3D { quad_Tetrahedral1, quad_Tetrahedral4,quad_Tetrahedral11,quad_Tetrahedral15 } QuadType3D;
#endif
