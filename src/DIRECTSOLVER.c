/*
 * =====================================================================================
 *
 *       Filename:  DIRECTSOLVER.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2014年10月26日 21时03分04秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Xu Fei (Fly), xufei@lsec.cc.ac.cn
 *        Company:  
 *
 * =====================================================================================
 */
#include <DIRECTSOLVER.h>
#include <umfpack.h>
/**利用umfpack求解线性方程，其具体的实现是用多重阵面法 */
void DIRECTSOLVER(MATRIX *matrix, double *rhs, double *sol)
{
  double t1, t2, t3, t4;
  int ret, i, j, k, l, begin, end;
  double value;
  int N_Eqn;
  int *Row, *kcol;
  double *Values;
  void *Symbolic, *Numeric;
  kcol=malloc(sizeof(int)*(matrix->RowPtr[matrix->N_Rows]));
  for(i=0;i<matrix->RowPtr[matrix->N_Rows];i++)
   {
     kcol[i]=matrix->KCol[i];
   }
  N_Eqn = matrix->N_Rows;
  Row = matrix->RowPtr;
  Values = matrix->Entries;  

  /** 第一步，进行符号登记*/
  ret = umfpack_di_symbolic(N_Eqn, N_Eqn, Row, kcol, Values,&Symbolic, NULL, NULL);
  if (ret!=0)
  {
     printf("error in umfpack_di_symbolic %d\n",ret);
     exit(4711);
  }
  /**第二步：进行数值计算 */
  ret = umfpack_di_numeric(Row, kcol, Values, Symbolic,&Numeric, NULL, NULL);
  umfpack_di_free_symbolic(&Symbolic);
 
  if (ret!=0)
  {
    printf("error in umfpack_di_numeric %d\n",ret);
    exit(4711);
  }
  /** 第三步：进行数值求解 */
  ret = umfpack_di_solve(UMFPACK_At, Row, kcol, Values, sol, rhs, Numeric, NULL, NULL);
  umfpack_di_free_numeric(&Numeric);
  if (ret!=0)
  {
    printf("error in umfpack_di_solve %d\n",ret); 
    exit(4711);
  }
 free(kcol);
}
