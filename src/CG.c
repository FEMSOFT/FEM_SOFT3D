#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "CG.h"

void CG(MATRIX *Matrix, double *b, double *x, double accur, int Max_It)
{
  //printf("Beging CG solving!\n");
  //double *x0; 		// initial value vector 
  double *r;		// the residue: r = b - A*x
  double *p;		// the conjugate vector
  double *tmp; 		// temporary vector, tmp = A*p
  double alpha, beta;	// lengths of steps
  double error;		//error estimator
  int niter = 0;
  
  int size;
  size = Matrix->N_Columns;
  r = malloc(size*sizeof(double));
  p = malloc(size*sizeof(double));  
  tmp = malloc(size*sizeof(double));


  MatrixDotVec(Matrix,x,p); //tmp1 = A*x0 
  
  SumVec(b,1.0,p,-1.0,r,size);  ////r = b - tmp1
  memcpy(p,r,size*sizeof(double));
  do{
    MatrixDotVec(Matrix,p,tmp);		//tmp = A*p
    alpha = VecDotVec(r,p,size) / VecDotVec(p,tmp,size);    
    SumVec(x,1.0,p,alpha,x,size);  //x = x0 + tmp1    
    SumVec(r,1.0,tmp,-alpha,r,size);
    error = NormVec(r,size);   //用残量的模来判断误差
    
    if(error<accur)
      break;
    
    beta = -VecDotVec(r,tmp,size) / VecDotVec(p,tmp,size); //beta = -(r,tmp)/(p,tmp)
   
    SumVec(r,1.0,p,beta,p,size);
    niter++;
    
  }while((error >= accur)&&(niter<Max_It));
  
  free(r);free(p);free(tmp);
}
