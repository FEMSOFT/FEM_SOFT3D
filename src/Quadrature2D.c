/*
 * =====================================================================================
 *
 *       Filename:  Quadrature2D.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2014年10月26日 16时19分09秒
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
#include "Quadrature2D.h"
#include "Mesh3D.h"
void BuildQuad2Face2D(FACEE *face,int N_Points, Quadrature2D *quad2d)
{
  int Num_Lines = 3;
  Produce_Quad2D(Num_Lines, N_Points, quad2d);
}
Quadrature2D *BuildQuad2Mesh2D(MESH *mesh, int N_Points)
{
  Quadrature2D *quad2d;
  quad2d = malloc(1*sizeof(Quadrature2D));
  FACEE *face;
  face = mesh->Faces[0];
  BuildQuad2Face2D(face,N_Points, quad2d);  
  return quad2d;
}

void OutPutQuad2D(Quadrature2D *quad)
{
  printf("Print the information for the quadrature2D:\n");
  int i;
  int N_Points;
  N_Points = quad->N_Points;
  printf("N_Points=%d\n",N_Points);
  for(i=0;i<N_Points;i++)
  {
    printf("X[%d]=%18.15f, Y[%d]=%18.15f, W[%d]=%18.15f\n",
	   i+1,quad->Quad_X[i],i+1,quad->Quad_Y[i],i+1,quad->Quad_W[i]);
  }
}
void Produce_Quad2D(int N_Lines, int N_Points, Quadrature2D *quad2d)
{ 
  if(N_Lines==3)
  {
    
    switch(N_Points)
    {
      
      case 1:
	quad2d->QuadType = quad_Triangle1;
	quad2d->Order = 1;
	quad2d->N_Points = 1;
	double Quad_X1[1] = {0.33333333333333333333};
	double Quad_Y1[1] = {0.33333333333333333333};
	double Quad_W1[1] = {1.0};
	SetValueQuad2D(1,Quad_X1,Quad_Y1,Quad_W1,quad2d);
	//printf("use the one point quadrature on triangle!\n");
	break;
      case 4:
	quad2d->QuadType = quad_Triangle4;
	quad2d->Order = 4;
	quad2d->N_Points = 4;
	double Quad_X4[4] = {0.07503111022261, 0.17855872826362,
				0.28001991549907, 0.66639024601470 }; 				
	double Quad_Y4[4] = {0.28001991549907, 0.66639024601470, 
			       0.07503111022261, 0.17855872826362};
	quad2d->Quad_Y = Quad_Y4;
	double Quad_W4[4] = {0.18195861825602, 0.31804138174398, 
			      0.18195861825602, 0.31804138174398};
	SetValueQuad2D(4,Quad_X4,Quad_Y4,Quad_W4,quad2d);
	//printf("use the four point quadrature on triangle!\n");
	
	break;
	  
      case 7:
	quad2d->QuadType = quad_Triangle7;
	quad2d->Order = 7;
	quad2d->N_Points = 7;
	double Quad_X7[7] = { 0.333333330000000, 0.059615870000000,
			      0.470142060000000, 0.470142060000000,
			      0.797426990000000, 0.101286510000000,
			      0.101286510000000};  
	double Quad_Y7[7] = { 0.333333330000000, 0.470142060000000,
			      0.059615870000000, 0.470142060000000,
			      0.101286510000000, 0.797426990000000,
			      0.101286510000000};			    
	double Quad_W7[7] = { 0.225000000000000, 0.132394150000000,
			      0.132394150000000, 0.132394150000000,
			      0.125939180000000, 0.125939180000000,
			      0.125939180000000};
	SetValueQuad2D(7,Quad_X7,Quad_Y7,Quad_W7,quad2d);
	//printf("use the seven point quadrature on triangle!\n");

	break;

      case 13:
	quad2d->QuadType = quad_Triangle13;
	quad2d->Order = 13;
	quad2d->N_Points = 13;
	double Quad_X13[13] = {   0.333333333333000,
				  0.260345966079000,
				  0.260345966079000,
				  0.479308067842000,
				  0.065130102902000,
				  0.065130102902000,
				  0.869739794196000,
				  0.312865496005000,
				  0.638444188570000,
				  0.048690315425000,
				  0.638444188570000,
				  0.312865496005000,
				  0.048690315425000};  
	double Quad_Y13[13] = {   0.333333333333000,
				  0.260345966079000,
				  0.479308067842000,
				  0.260345966079000,
				  0.065130102902000,
				  0.869739794196000,
				  0.065130102902000,
				  0.638444188570000,
				  0.048690315425000,
				  0.312865496005000,
				  0.312865496005000,
				  0.048690315425000,
				  0.638444188570000};
	double Quad_W13[13] = {  -0.149570044468000,
				  0.175615257433000,
				  0.175615257433000,
				  0.175615257433000,
				  0.053347235609000,
				  0.053347235609000,
				  0.053347235609000,
				  0.077113760890000,
				  0.077113760890000,
				  0.077113760890000,
				  0.077113760890000,
				  0.077113760890000,
				  0.077113760890000};
	SetValueQuad2D(13,Quad_X13,Quad_Y13,Quad_W13,quad2d);
	//printf("use the 13 point quadrature on triangle!\n");
    
	break;
	  
      case 27:
	quad2d->QuadType = quad_Triangle27;
	quad2d->Order = 27;
	quad2d->N_Points = 27;
	double Quad_X27[27] = {   0.032364948111276,
				  0.032364948111276,
				  0.935270103777448,
				  0.119350912282581,
				  0.119350912282581,
				  0.761298175434837,
				  0.534611048270758,
				  0.534611048270758,
				 -0.069222096541517,
				  0.203309900431282,
				  0.203309900431282,
				  0.593380199137435,
				  0.398969302965855,
				  0.398969302965855,
				  0.202061394068290,
				  0.593201213428213,
				  0.050178138310495,
				  0.356620648261293,
				  0.593201213428213,
				  0.050178138310495,
				  0.356620648261293,
				  0.807489003159792,
				  0.021022016536166,
				  0.171488980304042,
				  0.807489003159792,
				  0.021022016536166,
				  0.171488980304042};  
	double Quad_Y27[27] = {   0.032364948111276,
				  0.935270103777448,
				  0.032364948111276,
				  0.119350912282581,
				  0.761298175434837,
				  0.119350912282581,
				  0.534611048270758,
				 -0.069222096541517,
				  0.534611048270758,
				  0.203309900431282,
				  0.593380199137435,
				  0.203309900431282,
				  0.398969302965855,
				  0.202061394068290,
				  0.398969302965855,
				  0.050178138310495,
				  0.356620648261293,
				  0.593201213428213,
				  0.356620648261293,
				  0.593201213428213,
				  0.050178138310495,
				  0.021022016536166,
				  0.171488980304042,
				  0.807489003159792,
				  0.171488980304042,
				  0.807489003159792,
				  0.021022016536166};
	double Quad_W27[27] = { 0.013659731002678,
				0.013659731002678,
				0.013659731002678,
				0.036184540503418,
				0.036184540503418,
				0.036184540503418,
				0.000927006328961,
				0.000927006328961,
				0.000927006328961,
				0.059322977380774,
				0.059322977380774,
				0.059322977380774,
				0.077149534914813,
				0.077149534914813,
				0.077149534914813,
				0.052337111962204,
				0.052337111962204,
				0.052337111962204,
				0.052337111962204,
				0.052337111962204,
				0.052337111962204,
				0.020707659639141,
				0.020707659639141,
				0.020707659639141,
				0.020707659639141,
				0.020707659639141,
				0.020707659639141};
				
	SetValueQuad2D(27,Quad_X27,Quad_Y27,Quad_W27,quad2d);
	//printf("use the 27 point quadrature on triangle!\n");	  
	break;
	
      default:	  
	printf("This quadrature scheme is not implemented now!\n");
     }//end for switch
     
  }//end for if
}


void SetValueQuad2D(int N_Points,double *Quad_X,double *Quad_Y,double *Quad_W, Quadrature2D *quad2d)
{
  int i;
  quad2d->Quad_X = malloc(N_Points*sizeof(double));
  quad2d->Quad_Y = malloc(N_Points*sizeof(double));
  quad2d->Quad_W = malloc(N_Points*sizeof(double));
  for(i=0;i<N_Points;i++)
  {
    quad2d->Quad_X[i] = Quad_X[i];
    quad2d->Quad_Y[i] = Quad_Y[i];
    quad2d->Quad_W[i] = Quad_W[i];
  }
}

void FreeQuad2D(Quadrature2D *quad)
{
  //printf("delete the quadrature2D object!\n");
  free(quad->Quad_X);
  free(quad->Quad_Y);
  free(quad->Quad_W);
  free(quad);
}
