/*
 * =====================================================================================
 *
 *       Filename:  MG.h
 *        Created:  2014/12/23 11时29分20秒
 *       Revision:  nonee
 *       Compiler:  gcc
 *         Author:  xufei@lsec.cc.ac.cn
 *        Company:  LSEC
 *
 * =====================================================================================
 */
#ifndef __MGL__
#define __MGL__
#include "Fespace3D.h"
#include "Matrix.h"
#include "Mesh3D.h"

typedef struct MULTILEVEL_ {
  int Num_Levels;
  MESH *mesh;
  Fespace3D *Anasatz_Space;
  Fespace3D *Test_Space;
  MATRIX **Stiff_Matrix;
  MATRIX **Mass_Matrix;
  double **Rhst;
  MATRIX **Prolongs;
  MATRIX **Restricts;
 } MULTILEVEL;
//The idea: 在多重网格的机构中，每层的矩阵都会包含discreteform，但是并不需要我们在每一层网格
//上都存储这些数据。这里关键一点就是每一层在矩阵形成的时候一定要用到当前层网格的有限元空间来存储。
//所以为了节约存储空间，我们可以把之前的网格山除掉。目前还没有考虑网格可以放粗的情况。

//=========================================================================
//=============Some operators for the multilevel methods
//=========================================================================
//MULTILEVEL *BuildMultiLevel(Fespace2D *fespace, MATRIX *stiff_matrix, MATRIX *mass_matrix);
MULTILEVEL *BuildMultiLevel(MESH *mesh,DISCRETEFORM3D *discreteform_stiff, DISCRETEFORM3D *discreteform_mass);
void AddLevelFEM(MULTILEVEL *multilevel);
void AddLevelAFEM(MULTILEVEL *multilevel,int *M,int *num);
void FreeMultiLevel(MULTILEVEL *multilevel);
#endif
