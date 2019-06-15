/*
 * =====================================================================================
 *
 *       Filename:  Quadrature1D.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2014年10月26日 16时15分07秒
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
#include "Quadrature1D.h"
#include "Mesh3D.h"

Quadrature1D *BuildQuad1D(int N_Points)
{
  Quadrature1D *quad1d;
  quad1d = malloc(1*sizeof(Quadrature1D));
  Produce_Quad1D(N_Points, quad1d); 
  return quad1d;
}

void OutPutQuad1D(Quadrature1D *quad)
{
  printf("Print the information for the quadrature1D:\n");
  int i;
  int N_Points;
  N_Points = quad->N_Points;
  printf("N_Points=%d\n",N_Points);
  for(i=0;i<N_Points;i++)
  {
    printf("X[%d]=%18.15f, W[%d]=%18.15f\n",
	   i+1,quad->Quad_X[i],i+1,quad->Quad_W[i]);
  }
}
void Produce_Quad1D(int N_Points, Quadrature1D *quad1d)
{ 
    switch(N_Points)
    {
      
      case 1:
	quad1d->QuadType = quad_Line1;
	quad1d->Order = 1;
	quad1d->N_Points = 1;
	double Quad_X1[1] = {0.55555555555555555555};
	
	double Quad_W1[1] = {1.0};
	SetValueQuad1D(1,Quad_X1,Quad_W1,quad1d);
	//printf("use the one point quadrature on line [0,1]!\n");
	break;
	
      case 2:
        quad1d->QuadType = quad_Line2;
        quad1d->Order = 3;
        quad1d->N_Points = 2;
        double Quad_X2[2] = {0.21132486540518711774, 0.78867513459481288226};
        double Quad_W2[2] = {0.5, 0.5};
        SetValueQuad1D(2,Quad_X2,Quad_W2,quad1d);
	//printf("use the two point quadrature on line [0,1]!\n");
        break;
        
    case 3:
        quad1d->QuadType = quad_Line3;
        quad1d->Order = 5;
        quad1d->N_Points = 3;
        
        double Quad_X3[3] = {0.11270166537925831149, 0.5, 0.88729833462074168851};
        double Quad_W3[3] = {0.2777777777777778,0.4444444444444444,0.2777777777777778};
        
        SetValueQuad1D(3,Quad_X3,Quad_W3,quad1d);
	//printf("use the three point quadrature on line [0,1]!\n");
        break;
        
        
      case 4:
	quad1d->QuadType = quad_Line4;
	quad1d->Order = 7;
	quad1d->N_Points = 4;
	double Quad_X4[4] = { 0.069431844202973713731097404888715,
                              0.33000947820757187134432797392947,
                              0.66999052179242812865567202607053,
                              0.93056815579702628626890259511129 }; 				
	
	double Quad_W4[4] = { 0.17392742256872692485636378,
                              0.3260725774312730473880606,
                              0.3260725774312730473880606,
                              0.17392742256872692485636378};
	SetValueQuad1D(4,Quad_X4,Quad_W4,quad1d);
	//printf("use the four point quadrature on line [0,1]!\n");
	
	break;
	
	
      default:	  
	printf("This quadrature scheme is not implemented now!\n");
     }//end for switch
}


void SetValueQuad1D(int N_Points,double *Quad_X,double *Quad_W, Quadrature1D *quad1d)
{
  int i;
  quad1d->Quad_X = malloc(N_Points*sizeof(double));
  quad1d->Quad_W = malloc(N_Points*sizeof(double));
  for(i=0;i<N_Points;i++)
  {
    quad1d->Quad_X[i] = Quad_X[i];
    quad1d->Quad_W[i] = Quad_W[i];
  }
}

void FreeQuad1D(Quadrature1D *quad)
{
  //printf("delete the quadrature1D object!\n");
  free(quad->Quad_X);
  free(quad->Quad_W);
  free(quad);
}

