/*
 * =====================================================================================
 *
 *       Filename:  Base3D.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2014年10月20日 14时33分44秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Xu Fei (Fly), xufei@lsec.cc.ac.cn
 *        Company:  
 *
 * =====================================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include"Base3D.h"
#include "Enumerations.h"
#include "Constants.h"
BASEFUNCTION3D *BuildBaseFunct3D(int num_bas,int value_dim,int polydeg,
			       int accuracy, MapType maptype, int *dof,
			       double *nodal_points,BaseFunct3D *baseD000,
			       BaseFunct3D *baseD100,BaseFunct3D *baseD010,BaseFunct3D *baseD001,
			       BaseFunct3D *baseD110,BaseFunct3D *baseD101,BaseFunct3D *baseD011,
			       BaseFunct3D *baseD200,BaseFunct3D *baseD020,BaseFunct3D *baseD002,
			       NodalFunct3D *nodal_funct3D)
 {
   BASEFUNCTION3D *basefun3D;
   basefun3D = malloc(sizeof(BASEFUNCTION3D));
   basefun3D->Num_Bas = num_bas;
   //printf("Num_Bas=%d\n",num_bas);
   basefun3D->Value_Dim = value_dim;
   basefun3D->PolynomialDegree = polydeg;
   
   basefun3D->Accuracy = accuracy;
   basefun3D->Maptype = maptype;
   basefun3D->DOF = dof;
   basefun3D->Nodal_Points = nodal_points;
   
   basefun3D->Base[D000] = baseD000;
   basefun3D->Base[D100] = baseD100;
   basefun3D->Base[D010] = baseD010;
   basefun3D->Base[D001] = baseD001;
   basefun3D->Base[D110] = baseD110;
   basefun3D->Base[D101] = baseD101;
   basefun3D->Base[D011] = baseD011;
   basefun3D->Base[D200] = baseD200;
   basefun3D->Base[D020] = baseD020;
   basefun3D->Base[D002] = baseD002;
  
   basefun3D->NodalF3D = nodal_funct3D;
   
   // some seriers for save temp values 
   basefun3D->Values_X = malloc(sizeof(double)*num_bas);
   basefun3D->Values_Y = malloc(sizeof(double)*num_bas);
   basefun3D->Values_Z =  malloc(sizeof(double)*num_bas);
   basefun3D->Values_XY =  malloc(sizeof(double)*num_bas);
   basefun3D->Values_XZ =  malloc(sizeof(double)*num_bas);
   basefun3D->Values_YZ =  malloc(sizeof(double)*num_bas);
   
   return basefun3D;
 } 
 

void GetBaseValues(BASEFUNCTION3D *Base,MultiIndex3D MultiIndex,ELEMENT3D *elem3D,double Quad_x,double Quad_y,  double Quad_z,double Quad_xi,double Quad_eta,double Quad_zeta,double *Values)
{
  if(Base->Maptype==Affine)
  {
     int i;
    // get the values on the reference    
     if(MultiIndex==D000)
       Base->Base[MultiIndex](elem3D,Quad_x,Quad_y,Quad_z,Quad_xi,Quad_eta,Quad_zeta,Values);
     
    if((MultiIndex==D100)||(MultiIndex==D010)||(MultiIndex==D001))
    {     
      Base->Base[D100](elem3D,Quad_x,Quad_y,Quad_z,Quad_xi,Quad_eta,Quad_zeta,Base->Values_X);
      Base->Base[D010](elem3D,Quad_x,Quad_y,Quad_z,Quad_xi,Quad_eta,Quad_zeta,Base->Values_Y);
      Base->Base[D001](elem3D,Quad_x,Quad_y,Quad_z,Quad_xi,Quad_eta,Quad_zeta,Base->Values_Z);
      switch(MultiIndex)
      {
	case D100:
	  for(i=0;i<Base->Num_Bas;i++)
	  {
	    Values[i] =   Base->Values_X[i]*elem3D->Inv_Jacobian[0][0] 
	                + Base->Values_Y[i]*elem3D->Inv_Jacobian[1][0] 
	                + Base->Values_Z[i]*elem3D->Inv_Jacobian[2][0];	    
	  }//end for i
	  break;
	case D010:
	  for(i=0;i<Base->Num_Bas;i++)
	  {
	    Values[i] =   Base->Values_X[i]*elem3D->Inv_Jacobian[0][1] 
	                + Base->Values_Y[i]*elem3D->Inv_Jacobian[1][1] 
	                + Base->Values_Z[i]*elem3D->Inv_Jacobian[2][1];	    
	  }//end for i
	  break;
	case D001:
	  for(i=0;i<Base->Num_Bas;i++)
	  {
	    Values[i] =   Base->Values_X[i]*elem3D->Inv_Jacobian[0][2] 
	                + Base->Values_Y[i]*elem3D->Inv_Jacobian[1][2] 
	                + Base->Values_Z[i]*elem3D->Inv_Jacobian[2][2];	    
	  }//end for i
	  break;
         default:
           printf("error occur in Base3D.c\n");
      }//end for switch(MultiIndex)
    }//end for if((MultiIndex==D100)||(MultiIndex==D010)....)
        
    if((MultiIndex==D200)||(MultiIndex==D110)||(MultiIndex==D020)||(MultiIndex==D101)
       ||(MultiIndex==D011)||(MultiIndex==D002))
    {
      Base->Base[D200](elem3D,Quad_x,Quad_y,Quad_z,Quad_xi,Quad_eta,Quad_zeta,Base->Values_X); //H11
      Base->Base[D020](elem3D,Quad_x,Quad_y,Quad_z,Quad_xi,Quad_eta,Quad_zeta,Base->Values_Y); //H11
      Base->Base[D002](elem3D,Quad_x,Quad_y,Quad_z,Quad_xi,Quad_eta,Quad_zeta,Base->Values_Z);  //H22
      Base->Base[D110](elem3D,Quad_x,Quad_y,Quad_z,Quad_xi,Quad_eta,Quad_zeta,Base->Values_XY); //H12, H21
      Base->Base[D101](elem3D,Quad_x,Quad_y,Quad_z,Quad_xi,Quad_eta,Quad_zeta,Base->Values_XZ); //H12, H21
      Base->Base[D011](elem3D,Quad_x,Quad_y,Quad_z,Quad_xi,Quad_eta,Quad_zeta,Base->Values_YZ); //H12, H21
      // double H11,H12,H21,H22;
      switch(MultiIndex)
      {
	case D200:
	  for(i=0;i<Base->Num_Bas;i++)
	  {
	    Values[i] =  Base->Values_X[i] *elem3D->Inv_Jacobian[0][0]*elem3D->Inv_Jacobian[0][0]
	               + Base->Values_XY[i]*elem3D->Inv_Jacobian[0][0]*elem3D->Inv_Jacobian[1][0]
                       + Base->Values_XZ[i]*elem3D->Inv_Jacobian[0][0]*elem3D->Inv_Jacobian[2][0]
                       + Base->Values_XY[i]*elem3D->Inv_Jacobian[0][0]*elem3D->Inv_Jacobian[1][0]
                       + Base->Values_Y[i] *elem3D->Inv_Jacobian[1][0]*elem3D->Inv_Jacobian[1][0]
                       + Base->Values_YZ[i]*elem3D->Inv_Jacobian[1][0]*elem3D->Inv_Jacobian[2][0]
                       + Base->Values_XZ[i]*elem3D->Inv_Jacobian[0][0]*elem3D->Inv_Jacobian[2][0]
                       + Base->Values_YZ[i]*elem3D->Inv_Jacobian[1][0]*elem3D->Inv_Jacobian[2][0]
                       + Base->Values_Z[i] *elem3D->Inv_Jacobian[2][0]*elem3D->Inv_Jacobian[2][0];
	  }//end for i
	  break;
	
	case D020:
	  for(i=0;i<Base->Num_Bas;i++)
	  {
	     Values[i] =  Base->Values_X[i] *elem3D->Inv_Jacobian[0][1]*elem3D->Inv_Jacobian[0][1]
	                + Base->Values_XY[i]*elem3D->Inv_Jacobian[0][1]*elem3D->Inv_Jacobian[1][1]
                        + Base->Values_XZ[i]*elem3D->Inv_Jacobian[0][1]*elem3D->Inv_Jacobian[2][1]
                        + Base->Values_XY[i]*elem3D->Inv_Jacobian[0][1]*elem3D->Inv_Jacobian[1][1]
                        + Base->Values_Y[i] *elem3D->Inv_Jacobian[1][1]*elem3D->Inv_Jacobian[1][1]
                        + Base->Values_YZ[i]*elem3D->Inv_Jacobian[1][1]*elem3D->Inv_Jacobian[2][1]
                        + Base->Values_XZ[i]*elem3D->Inv_Jacobian[0][1]*elem3D->Inv_Jacobian[2][1]
                        + Base->Values_YZ[i]*elem3D->Inv_Jacobian[1][1]*elem3D->Inv_Jacobian[2][1]
                        + Base->Values_Z[i] *elem3D->Inv_Jacobian[2][1]*elem3D->Inv_Jacobian[2][1];
	  }//end for i
	  break;
	
	case D002:
	  for(i=0;i<Base->Num_Bas;i++)
	  {
	      Values[i] =  Base->Values_X[i] *elem3D->Inv_Jacobian[0][2]*elem3D->Inv_Jacobian[0][2]
	                 + Base->Values_XY[i]*elem3D->Inv_Jacobian[0][2]*elem3D->Inv_Jacobian[1][2]
                         + Base->Values_XZ[i]*elem3D->Inv_Jacobian[0][2]*elem3D->Inv_Jacobian[2][2]
                         + Base->Values_XY[i]*elem3D->Inv_Jacobian[0][2]*elem3D->Inv_Jacobian[1][2]
                         + Base->Values_Y[i] *elem3D->Inv_Jacobian[1][2]*elem3D->Inv_Jacobian[1][2]
                         + Base->Values_YZ[i]*elem3D->Inv_Jacobian[1][2]*elem3D->Inv_Jacobian[2][2]
                         + Base->Values_XZ[i]*elem3D->Inv_Jacobian[0][2]*elem3D->Inv_Jacobian[2][2]
                         + Base->Values_YZ[i]*elem3D->Inv_Jacobian[1][2]*elem3D->Inv_Jacobian[2][2]
                         + Base->Values_Z[i] *elem3D->Inv_Jacobian[2][2]*elem3D->Inv_Jacobian[2][2];
	  }//end for i
	  break;	  

        case D110:
	  for(i=0;i<Base->Num_Bas;i++)
	  {
	     Values[i] =  Base->Values_X[i] *elem3D->Inv_Jacobian[0][0]*elem3D->Inv_Jacobian[0][1]
	                + Base->Values_XY[i]*elem3D->Inv_Jacobian[0][1]*elem3D->Inv_Jacobian[1][0]
                        + Base->Values_XZ[i]*elem3D->Inv_Jacobian[0][1]*elem3D->Inv_Jacobian[2][0]
                        + Base->Values_XY[i]*elem3D->Inv_Jacobian[0][0]*elem3D->Inv_Jacobian[1][1]
                        + Base->Values_Y[i] *elem3D->Inv_Jacobian[1][0]*elem3D->Inv_Jacobian[1][1]
                        + Base->Values_YZ[i]*elem3D->Inv_Jacobian[1][1]*elem3D->Inv_Jacobian[2][0]
                        + Base->Values_XZ[i]*elem3D->Inv_Jacobian[0][0]*elem3D->Inv_Jacobian[2][1]
                        + Base->Values_YZ[i]*elem3D->Inv_Jacobian[1][0]*elem3D->Inv_Jacobian[2][1]
                        + Base->Values_Z[i] *elem3D->Inv_Jacobian[2][0]*elem3D->Inv_Jacobian[2][1];


	  }//end for i
	  break;
         case D101:
	  for(i=0;i<Base->Num_Bas;i++)
	  {
             Values[i] =  Base->Values_X[i] *elem3D->Inv_Jacobian[0][0]*elem3D->Inv_Jacobian[0][2]
	                + Base->Values_XY[i]*elem3D->Inv_Jacobian[0][2]*elem3D->Inv_Jacobian[1][0]
                        + Base->Values_XZ[i]*elem3D->Inv_Jacobian[0][2]*elem3D->Inv_Jacobian[2][0]
                        + Base->Values_XY[i]*elem3D->Inv_Jacobian[0][0]*elem3D->Inv_Jacobian[1][2]
                        + Base->Values_Y[i] *elem3D->Inv_Jacobian[1][0]*elem3D->Inv_Jacobian[1][2]
                        + Base->Values_YZ[i]*elem3D->Inv_Jacobian[1][2]*elem3D->Inv_Jacobian[2][0]
                        + Base->Values_XZ[i]*elem3D->Inv_Jacobian[0][0]*elem3D->Inv_Jacobian[2][2]
                        + Base->Values_YZ[i]*elem3D->Inv_Jacobian[1][0]*elem3D->Inv_Jacobian[2][2]
                        + Base->Values_Z[i] *elem3D->Inv_Jacobian[2][0]*elem3D->Inv_Jacobian[2][2];



	  }//end for i
	  break;

         case D011:
	  for(i=0;i<Base->Num_Bas;i++)
	  {
	     Values[i] =  Base->Values_X[i] *elem3D->Inv_Jacobian[0][2]*elem3D->Inv_Jacobian[0][1]
	                + Base->Values_XY[i]*elem3D->Inv_Jacobian[0][1]*elem3D->Inv_Jacobian[1][2]
                        + Base->Values_XZ[i]*elem3D->Inv_Jacobian[0][1]*elem3D->Inv_Jacobian[2][2]
                        + Base->Values_XY[i]*elem3D->Inv_Jacobian[0][2]*elem3D->Inv_Jacobian[1][1]
                        + Base->Values_Y[i] *elem3D->Inv_Jacobian[1][2]*elem3D->Inv_Jacobian[1][1]
                        + Base->Values_YZ[i]*elem3D->Inv_Jacobian[1][1]*elem3D->Inv_Jacobian[2][2]
                        + Base->Values_XZ[i]*elem3D->Inv_Jacobian[0][2]*elem3D->Inv_Jacobian[2][1]
                        + Base->Values_YZ[i]*elem3D->Inv_Jacobian[1][2]*elem3D->Inv_Jacobian[2][1]
                        + Base->Values_Z[i] *elem3D->Inv_Jacobian[2][2]*elem3D->Inv_Jacobian[2][1];



	  }//end for i
	  break;


      }//end for switch(MultiIndex)      
    }//end for if((MultiIndex==D20)||(MultiIndex==D11)||(MultiIndex==D02))
    
  }//end for if(Base->Invariant==TRUE)
  else
  {
    Base->Base[MultiIndex](elem3D,Quad_x,Quad_y,Quad_z,Quad_xi,
	                   Quad_eta, Quad_zeta,Values);
  }//end for else
}
 
void FreeBase3D(BASEFUNCTION3D *Base3D)
 {
   free(Base3D->Values_X);
   free(Base3D->Values_Y);
   free(Base3D->Values_Z);
   free(Base3D->Values_XY);
   free(Base3D->Values_YZ);
   free(Base3D->Values_XZ);
   free(Base3D);
 }
