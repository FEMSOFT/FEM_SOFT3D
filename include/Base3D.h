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

#ifndef __BASE3D__
#define __BASE3D__
#include "Enumerations.h"
#include "Constants.h"

/** define a type of struct to describe the base functions */
typedef struct BASEFUNCTION3D_ {

    // the number of the base functions
    int Num_Bas;
    // the dimension of the value for each base function
    int Value_Dim;
    // the polynomial degree in the base functions
    int PolynomialDegree;
    // the accuracy of the element
    int Accuracy;
    // define the type of element: Invariant=1 means the element 
    // is invariant under the affine transform, Invariant=0 means the 
    // element is defined on the actual element 
    //bool Invariant;  
    //the maptype: Affine means we have the invariant according to the affine map
    // Piola: another type of mapping used in the H(div) and H(curl)
    // Actual: means we can not use tha mapping in the matrix building process
    MapType Maptype;
    // the distribution of the degree of freedom:
    // DOF[0]: the dof on the vertices
    // DOF[1]: the dof on the lines
    // DOF[2]: the dof on the face
    // DOF[3]: the dof on the volume
    int *DOF;
    // the point coordinations for all the base functions
    double *Nodal_Points;
    // the corresponding based finctions: [D00,D10,D01,D20,D11,D02]
    BaseFunct3D *Base[N_MultiIndices3D];
    // the corresponding nodal function: always needed for interpolation
    NodalFunct3D *NodalF3D;
    
    // some seriers for save temp values????????????????????????????????????????? 
    double *Values_X;
    double *Values_Y;
    double *Values_Z;
    double *Values_XY;
    double *Values_YZ;
    double *Values_XZ;
    
} BASEFUNCTION3D;

BASEFUNCTION3D *BuildBaseFunct3D(int num_bas,int value_dim,int polydeg,
			       int accuracy,MapType maptype, int *dof,
			       double *nodal_points,BaseFunct3D *baseD000,
			       BaseFunct3D *baseD100,BaseFunct3D *baseD010,BaseFunct3D *baseD001,
			       BaseFunct3D *baseD110,BaseFunct3D *baseD101,
			       BaseFunct3D *baseD011,BaseFunct3D *baseD200, BaseFunct3D *baseD020,
			       BaseFunct3D *baseD002,NodalFunct3D *nodal_funct);
void GetBaseValues(BASEFUNCTION3D *Base,MultiIndex3D MultiIndex,
		   ELEMENT3D *elem3D,double Quad_x,double Quad_y,double Quad_z,
		   double Quad_xi,double Quad_eta,double Quad_zeta,double *Values);
void FreeBase3D(BASEFUNCTION3D *Base3D);
#endif

