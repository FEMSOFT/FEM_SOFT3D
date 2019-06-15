/*
 * =====================================================================================
 *
 *       Filename:  DiscretForm3D.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2014年10月20日 13时58分55秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Xu Fei (Fly), xufei@lsec.cc.ac.cn
 *        Company:  
 *
 * =====================================================================================
 */
#ifndef __DISCRETEFORM3D__
#define __DISCRETEFORM3D__
#include "Fespace3D.h"
#include "Fefunction3D.h"
#include "Mesh3D.h"
#include "Base3D.h"
#include "Enumerations.h"
#include "Constants.h"
#include "Quadrature3D.h"
#include "Quadrature2D.h"
#include "Quadrature1D.h"


typedef struct DISCRETEFORM3D_ {
	/** the anasatz space */
	/** 试探有限元函数空间 */
    Fespace3D *Anasatz_Space;    
    /** test space */
       /** 检验有限元函数空间 */
    Fespace3D *Test_Space;    
    /** the Multiindex for testspace */
    	/** 形成刚度矩阵的时候需要检验有限元空间的导数个数 */
    int N_Test_MultiIndex;
    	/** 形成刚度矩阵的时候需要检验有限元空间的导数 */
    MultiIndex3D *Test_MultiIndex;    
    /** the Multiindex for anasatzspace */
     /** 形成刚度矩阵的时候需要试探有限元空间的导数个数 */
    int N_Anasatz_MultiIndex;
     /** 形成刚度矩阵的时候需要试探有限元空间的导数 */
    MultiIndex3D *Anasatz_MultiIndex;     
    //解析函数应该不需要，直接在具体的计算的接口给就可以，减少程序的复杂性
    /** Aux Function3D */
    /** 形成刚度矩阵需要的用户自己定义的解析函数*/
    /** 解析函数的个数（或者是返回值的维数） */
    //int Dim_AuxFunct;
    /** 存储解析函数 */
    //Functions3D *AuxFunct;
     /** 存储解析函数在某一点的值 */
    //double *AuxFunct_Values;    
    /** Aux Fefunction */
     /** 形成刚度矩阵所需要的有限元函数对象 */
     /** 需要多少个有限元函数 */
    int N_AuxFeFun;
    /** 存储有限元函数 */
    Fefunction3D **AuxFeFun;
    /** the Multiindex for AuxFefunction3D */
    /** 对这些有限元函数我们需要对他们取什么样的导数值*/
    /** 就是说AuxFefunct中的每个有限元函数需要在数组AuxFefun_MultiIndex 中需要计算的导数个数 */
    //比如 N_AuxFefun_MultiIndex=[3，2 ，1 ]: 表示的是前3个导数值
    //是第0个有限元函数的导数，接下来的2个导数值是第1个有限元函数的导数
    // 最后一个导数是第2个有限元函数的导数 
    int *N_AuxFeFun_MultiIndex;
    /** 对有限元函数需要取值的导数指标 */
    //记录的是一组导数的符号：比如[D00 D10 D01 D11 D00 D10 D01]	
    MultiIndex3D *AuxFeFun_MultiIndex;
    /** 从有限元函数列中取到的在某一点的函数值的个数 */
    /**这个应该可以直接计算出来的，不需要给定：*/
    //具体的值是：N_AuxFefun_Values = sum(N_AuxFefun_MultiIndex*dim_values);
    int N_AuxFeFun_Values;
     /** 用来存储有限元函数值的数组： 维数 N_AuxFefun_Values  */
    double *AuxFeFun_Values;    

    


    //用户函数
    int N_UserFunction;
    Functions3D  **UserFunction;

    int N_FeConst;
    double *FeConst;

    

    /** Discrete Form for the Matrix */
    DiscreteFormMatrix *DiscreteFormVolu;
    /** 定义形成刚度矩阵的变分形式 */
    DiscreteFormMatrixFace *DiscreteFormFace;
    /** 定义在线上的变分形式：积分区域是在每个单元的边上*/
    DiscreteFormMatrixLine *DiscreteFormLine;
    /**定义在边界上的变分形式，积分区域是在单元的边界上 */
    //DiscreteFormMatrixLine *DiscreteFormBD;
    /**定义在点上的变分形式，在配置点方法中可能会用到*/
    DiscreteFormMatrix *DiscreteFormVert;
    /** The information of the quadrature scheme */
    /** 形成刚度矩阵需要的积分点的信息 */
    /**三维的积分信息*/
    Quadrature3D *Quad3D;
    /**二维的积分信息*/
    Quadrature2D *Quad2D;
    /**一维的积分信息*/
    Quadrature1D *Quad1D;

}DISCRETEFORM3D;

//下面是一些操作函数：最重要的是建立离散变分形式
/** 建立离散变分形式 （不需有限元函数）*/
DISCRETEFORM3D *BuildDiscreteForm3D(Fespace3D *anasatz_space,Fespace3D *test_space,int n_test_MultiIndex,
			MultiIndex3D *test_MultiIndex,int n_anasatz_MultiIndex,MultiIndex3D *anasatz_MultiIndex,
			DiscreteFormMatrix *discreteformvolu,int Quad_Points3D,int Quad_Points2D,
			int Quad_Points1D);
/** 建立离散变分形式对象 （需要有限元函数）*/
DISCRETEFORM3D *BuildDiscreteForm3DAuxFeFun(Fespace3D *anasatz_space,Fespace3D *test_space,int n_test_MultiIndex,
					    MultiIndex3D *test_MultiIndex,int n_anasatz_MultiIndex,MultiIndex3D *anasatz_MultiIndex,
					    int n_auxfefun, Fefunction3D **auxfefun, int *n_auxFefun_MultiIndex,
					    MultiIndex3D *auxFefun_MultiIndex,DiscreteFormMatrix *discreteformface,
					    int Quad_Points3D,int Quad_Points2D, int Quad_Points1D);		    
/** 建立矩阵对象 （需要有限元函数和解析函数）*/
DISCRETEFORM3D *BuildDiscreteFormAll3D(Fespace3D *anasatz_space,Fespace3D *test_space,int n_test_MultiIndex,
			               MultiIndex3D *test_MultiIndex,int n_anasatz_MultiIndex,MultiIndex3D 
				       *anasatz_MultiIndex,int n_auxfefun,Fefunction3D **auxfefun,
				       int *N_AuxFefun_MultiIndex,MultiIndex3D *AuxFefun_MultiIndex,
                                       int N_UserFunction, Functions3D **userfunction,
			               DiscreteFormMatrix *discreteformvolu,int Quad_Points3D,
				       int Quad_Points2D,int Quad_Points1D);

DISCRETEFORM3D *BuildDiscreteFormFeConst3D(Fespace3D *anasatz_space,Fespace3D *test_space,int n_test_MultiIndex,
			               MultiIndex3D *test_MultiIndex,int n_anasatz_MultiIndex,MultiIndex3D 
				       *anasatz_MultiIndex,int n_auxfefun,Fefunction3D **auxfefun,
				       int *N_AuxFefun_MultiIndex,MultiIndex3D *AuxFefun_MultiIndex,
                                       int N_UserFunction, Functions3D **userfunction,
			               DiscreteFormMatrix *discreteformvolu,int N_FeConst, double *FeConst,int Quad_Points3D,
				       int Quad_Points2D,int Quad_Points1D);



//释放内存空间
void FreeDiscreteForm3D(DISCRETEFORM3D *discreteformface);
#endif
