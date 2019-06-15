/*
 * =====================================================================================
 *
 *       Filename:  GMG.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2015年01月09日 01时59分48秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Xu Fei (Fly), xufei@lsec.cc.ac.cn
 *        Company:  
 *
 * =====================================================================================
 */
#include "GMG.h"
#include "MG.h"
#include "Enumerations.h"
#include "Matrix.h"
#include "Constants.h"
void GMGSolver(MULTILEVEL *multilevel,double *rhs,double *solution,int level,int smooth,int GMGIter)
{
  int i;
  double *residual;
  residual = malloc(sizeof(double));
  printf("Do the %d GMG iterations\n",GMGIter);
  for(i=0;i<GMGIter;++i)
  { 
     printf("begin the %d inner iterations\n",i);
     GMGIteration(multilevel,rhs,solution,level,smooth,residual);
  }
  
  printf("End of the GMG iterations\n");
}

void GMGIteration(MULTILEVEL *multilevel,double *rhs,double *solution,int level,int smooth,double *return_Residual)
{
    int i;
    MATRIX *Matrix,*Prolong,*Restrict;
    double *residual,*tmpval,*tmpval1,*solu;
    Matrix = multilevel->Stiff_Matrix[level];
    if(level==0)//直接求解
      DIRECTSOLVER(Matrix,rhs,solution);
    else//光滑
    {
      Prolong  = multilevel->Prolongs[level-1];
      Restrict = multilevel->Restricts[level-1];
      residual = malloc(sizeof(double)*Restrict->N_Rows);
      solu     = malloc(sizeof(double)*Restrict->N_Rows);
      memset(solu,0.0,sizeof(double)*Restrict->N_Rows);
      tmpval = malloc(sizeof(double)*Prolong->N_Rows);
      memset(tmpval,0.0,sizeof(double)*Prolong->N_Rows);
      tmpval1 = malloc(sizeof(double)*Prolong->N_Rows);//必须加上这个才行，为什么？？？？？？？？？？？？？？？？？


      CG(Matrix,rhs,solution,1e-5,smooth);
      MatrixDotVec(Matrix,solution,tmpval);
      for(i=0;i<Matrix->N_Rows;++i)
       tmpval1[i]=rhs[i]-tmpval[i];
      MatrixDotVec(Restrict,tmpval1,residual);
      double *r_residual;
      r_residual = malloc(sizeof(double));
      GMGIteration(multilevel,residual,solu,level-1,smooth,r_residual);
      MatrixDotVec(Prolong,solu,tmpval);
      for(i=0;i<Matrix->N_Rows;++i)
       solution[i] = solution[i]+tmpval1[i];
      //后光滑
      CG(Matrix,rhs,solution,1e-5,smooth);
      //calculate the residual
      double val;
      val = 0.0;
      MatrixDotVec(Matrix,solution,tmpval);
      for(i=0;i<Matrix->N_Rows;++i)
	  {
		  tmpval[i]=rhs[i]-tmpval[i];
		  val += tmpval[i]*tmpval[i];
	  }
	  return_Residual[0] = sqrt(val);
  
      free(residual);
      free(solu);
      free(tmpval);
    }
}
