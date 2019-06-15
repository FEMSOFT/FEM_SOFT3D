/*
 * =====================================================================================
 *
 *       Filename:  Fespace3D.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2014年10月20日 13时27分21秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Xu Fei (Fly), xufei@lsec.cc.ac.cn
 *        Company:  
 *
 * =====================================================================================
 */
#ifndef __FESPACE3D__
#define __FESPACE3D__
#include "Mesh3D.h"
#include "Enumerations.h"
#include "Constants.h"
#include "Base3D.h"

// typedef struct BOUNDCOND3D_ {
//     int Num_BDComps;
//     int *ID_BDComp;
//     BoundType *BdTypes;
// }BOUNDCOND3D;

typedef struct Fespace3D_ {
    int DOF_ALL;
    int Num_Volus;
    MESH *Mesh;
    BASEFUNCTION3D *Base;
    //边界函数类型： 应该还需要跟有限元函数结合在一起，
    //有些复杂的情况需要这样处理
    BoundCondFunct3D *BoundCond;//?????????????????????
    //自由度信息
    int *BeginIndex;
    //每个单元上的自由度
    int *GlobalNumbers;    
} Fespace3D;

Fespace3D *BUILDFEMSPACE3D(MESH *mesh, FemType femtype,BoundCondFunct3D *boundcond);
Fespace3D *BUILDFEMSPACE3DBASE(MESH *mesh, BASEFUNCTION3D *base,BoundCondFunct3D *boundcond);
BASEFUNCTION3D *BuildBase3D(FemType femtype);
void OutPutFespace3D(Fespace3D *fesp3D);
BASEFUNCTION3D *BuildBase3D(FemType femtype);
void BuildDOF3D(MESH *mesh,BASEFUNCTION3D *base3D,int * BeginIndex,int *GlobalNumbers);
void FreeFespace3D(Fespace3D *fesp);
void FreeFespace3DL(Fespace3D *fesp);
#endif
