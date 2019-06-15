#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <LinearAlg.h>


void MaxVec(double *a,int dim,double *value,int *index)
{
  int i;
  //double value;
  value[0] = a[0];
  index[0] = 0;
  for(i=1;i<dim;i++)
  {
    if(a[i]>value[0])
    {
      value[0] = a[i];
      index[0] = i;
    }
  }
}
void MaxAbsVec(double *a,int dim,double *value,int *index)
{
  int i;
  //double value;
  value[0] = abs(a[0]);
  index[0] = 0;
  for(i=1;i<dim;i++)
  {
    if(abs(a[i])>value[0])
    {
      value[0] = abs(a[i]);
      index[0] = i;
    }
  }
  value[0] = a[index[0]];
}


void NormalizeVec( double *a, int dim )
{
    double value;
    int index;
    MaxAbsVec( a, dim, &value, &index );
    
    int i;
    for( i=0; i<dim; i++ )
    {
        a[i] = a[i]/value;
    }
}

void SumVecSelf(double *a, double alpha, double *b, double beta, int size)
{
  int i;
  for(i=0;i<size;i++)
  {
    b[i] = alpha*a[i] + beta*b[i];
  }
}

void SumVec(double *a, double alpha, double *b, double beta,double *c,int size)
{
  int i;
  for(i=0;i<size;i++)
  {
    c[i] = alpha*a[i] + beta*b[i];
  }
}

// p = alpha*p
void ScalVec(double alpha, double *p, int size)
{
	int i;
	for(i=0;i<size;i++)
	{
	  p[i] = alpha * p[i];
	}	
}

// q = alpha*p
void NumVec(double alpha, double *p, double *q, int size)
{
	int i;
	for(i=0;i<size;i++)
	{
	  q[i] = alpha * p[i];
	}	
}

double
NormVec(double *b, int size)
{
	int i; double norm=0.0;
	for(i=0;i<size;i++)
	{
		norm = norm + b[i]*b[i];
	}
	norm = sqrt(norm);
	
	return norm;
}

void MatrixDotVec(MATRIX *Matrix, double *x, double *r)
{ 
  int N_Rows = Matrix->N_Rows;
  int N_Columns = Matrix->N_Columns;
  int i, j, length, start, end;
  int *RowPtr,*KCol;
  double *Entries, tmp;
  RowPtr = Matrix->RowPtr;
  KCol = Matrix->KCol;
  Entries = Matrix->Entries;
  memset(r,0.0,N_Rows*sizeof(double));
  
  start = RowPtr[0];
  for(i=0;i<N_Rows;i++)
  {   
    end = RowPtr[i+1];
    length = end - start;
    for(j=0;j<length;j++)
    {
      r[i] += Entries[start+j]*x[KCol[start+j]];      
    }
    start = end;
  }
}

 void MatMatVec(MATRIX *MatrixA, MATRIX *MatrixB,double *x, double *r)
{
  double *tmp;
  tmp=malloc(sizeof(double)*MatrixB->N_Rows);
  MatrixDotVec(MatrixB, x, tmp);
  MatrixDotVec(MatrixA, tmp,r);
  free(tmp);
}




void TripleMatrix(MATRIX *matrixa,MATRIX *matrixb,MATRIX *matrixc,int num, double *r)
{

 int i,j;
 double *bb,*cc,*dd;
 bb = malloc(sizeof(double)*matrixc->N_Columns);
 cc = malloc(sizeof(double)*matrixc->N_Rows);
 dd = malloc(sizeof(double)*matrixc->N_Rows);
 memset(bb,0.0,sizeof(double)*matrixc->N_Columns);
 bb[num]=1;
 MatrixDotVec(matrixc,bb,cc); 
 MatrixDotVec(matrixb,cc,dd); 
 MatrixDotVec(matrixa,dd,r); 
free(bb);
free(cc);
free(dd);
}



int power(int a,int b)
{
  int i,sum=1;
  for(i=0;i<b;++i)
      sum*=a;

 return sum;
}


void ArrayArrayArray(double A[4][4],double B[4][4],double C[4][4],double D[4][4])
{

D[0][0] =  C[0][0]*(A[0][0]*B[0][0] + A[1][0]*B[1][0] + A[2][0]*B[2][0] + A[3][0]*B[3][0]) + C[1][0]*(A[0][0]*B[0][1] + A[1][0]*B[1][1] + A[2][0]*B[2][1] + A[3][0]*B[3][1]) 
	     + C[2][0]*(A[0][0]*B[0][2] + A[1][0]*B[1][2] + A[2][0]*B[2][2] + A[3][0]*B[3][2]) + C[3][0]*(A[0][0]*B[0][3] + A[1][0]*B[1][3] + A[2][0]*B[2][3] + A[3][0]*B[3][3]);
D[0][1] =  C[0][1]*(A[0][0]*B[0][0] + A[1][0]*B[1][0] + A[2][0]*B[2][0] + A[3][0]*B[3][0]) + C[1][1]*(A[0][0]*B[0][1] + A[1][0]*B[1][1] + A[2][0]*B[2][1] + A[3][0]*B[3][1]) 
	     + C[2][1]*(A[0][0]*B[0][2] + A[1][0]*B[1][2] + A[2][0]*B[2][2] + A[3][0]*B[3][2]) + C[3][1]*(A[0][0]*B[0][3] + A[1][0]*B[1][3] + A[2][0]*B[2][3] + A[3][0]*B[3][3]); 
D[0][2] =  C[0][2]*(A[0][0]*B[0][0] + A[1][0]*B[1][0] + A[2][0]*B[2][0] + A[3][0]*B[3][0]) + C[1][2]*(A[0][0]*B[0][1] + A[1][0]*B[1][1] + A[2][0]*B[2][1] + A[3][0]*B[3][1]) 
	     + C[2][2]*(A[0][0]*B[0][2] + A[1][0]*B[1][2] + A[2][0]*B[2][2] + A[3][0]*B[3][2]) + C[3][2]*(A[0][0]*B[0][3] + A[1][0]*B[1][3] + A[2][0]*B[2][3] + A[3][0]*B[3][3]); 
D[0][3] =  C[0][3]*(A[0][0]*B[0][0] + A[1][0]*B[1][0] + A[2][0]*B[2][0] + A[3][0]*B[3][0]) + C[1][3]*(A[0][0]*B[0][1] + A[1][0]*B[1][1] + A[2][0]*B[2][1] + A[3][0]*B[3][1]) 
	     + C[2][3]*(A[0][0]*B[0][2] + A[1][0]*B[1][2] + A[2][0]*B[2][2] + A[3][0]*B[3][2]) + C[3][3]*(A[0][0]*B[0][3] + A[1][0]*B[1][3] + A[2][0]*B[2][3] + A[3][0]*B[3][3]);

D[1][0] =  C[0][0]*(A[0][1]*B[0][0] + A[1][1]*B[1][0] + A[2][1]*B[2][0] + A[3][1]*B[3][0]) + C[1][0]*(A[0][1]*B[0][1] + A[1][1]*B[1][1] + A[2][1]*B[2][1] + A[3][1]*B[3][1]) 
         + C[2][0]*(A[0][1]*B[0][2] + A[1][1]*B[1][2] + A[2][1]*B[2][2] + A[3][1]*B[3][2]) + C[3][0]*(A[0][1]*B[0][3] + A[1][1]*B[1][3] + A[2][1]*B[2][3] + A[3][1]*B[3][3]); 
D[1][1] =  C[0][1]*(A[0][1]*B[0][0] + A[1][1]*B[1][0] + A[2][1]*B[2][0] + A[3][1]*B[3][0]) + C[1][1]*(A[0][1]*B[0][1] + A[1][1]*B[1][1] + A[2][1]*B[2][1] + A[3][1]*B[3][1]) 
	     + C[2][1]*(A[0][1]*B[0][2] + A[1][1]*B[1][2] + A[2][1]*B[2][2] + A[3][1]*B[3][2]) + C[3][1]*(A[0][1]*B[0][3] + A[1][1]*B[1][3] + A[2][1]*B[2][3] + A[3][1]*B[3][3]); 
D[1][2] =  C[0][2]*(A[0][1]*B[0][0] + A[1][1]*B[1][0] + A[2][1]*B[2][0] + A[3][1]*B[3][0]) + C[1][2]*(A[0][1]*B[0][1] + A[1][1]*B[1][1] + A[2][1]*B[2][1] + A[3][1]*B[3][1]) 
	     + C[2][2]*(A[0][1]*B[0][2] + A[1][1]*B[1][2] + A[2][1]*B[2][2] + A[3][1]*B[3][2]) + C[3][2]*(A[0][1]*B[0][3] + A[1][1]*B[1][3] + A[2][1]*B[2][3] + A[3][1]*B[3][3]); 
D[1][3] =  C[0][3]*(A[0][1]*B[0][0] + A[1][1]*B[1][0] + A[2][1]*B[2][0] + A[3][1]*B[3][0]) + C[1][3]*(A[0][1]*B[0][1] + A[1][1]*B[1][1] + A[2][1]*B[2][1] + A[3][1]*B[3][1]) 
	     + C[2][3]*(A[0][1]*B[0][2] + A[1][1]*B[1][2] + A[2][1]*B[2][2] + A[3][1]*B[3][2]) + C[3][3]*(A[0][1]*B[0][3] + A[1][1]*B[1][3] + A[2][1]*B[2][3] + A[3][1]*B[3][3]);

D[2][0] =  C[0][0]*(A[0][2]*B[0][0] + A[1][2]*B[1][0] + A[2][2]*B[2][0] + A[3][2]*B[3][0]) + C[1][0]*(A[0][2]*B[0][1] + A[1][2]*B[1][1] + A[2][2]*B[2][1] + A[3][2]*B[3][1]) 
         + C[2][0]*(A[0][2]*B[0][2] + A[1][2]*B[1][2] + A[2][2]*B[2][2] + A[3][2]*B[3][2]) + C[3][0]*(A[0][2]*B[0][3] + A[1][2]*B[1][3] + A[2][2]*B[2][3] + A[3][2]*B[3][3]);
D[2][1] =  C[0][1]*(A[0][2]*B[0][0] + A[1][2]*B[1][0] + A[2][2]*B[2][0] + A[3][2]*B[3][0]) + C[1][1]*(A[0][2]*B[0][1] + A[1][2]*B[1][1] + A[2][2]*B[2][1] + A[3][2]*B[3][1]) 
	     + C[2][1]*(A[0][2]*B[0][2] + A[1][2]*B[1][2] + A[2][2]*B[2][2] + A[3][2]*B[3][2]) + C[3][1]*(A[0][2]*B[0][3] + A[1][2]*B[1][3] + A[2][2]*B[2][3] + A[3][2]*B[3][3]); 
D[2][2] =  C[0][2]*(A[0][2]*B[0][0] + A[1][2]*B[1][0] + A[2][2]*B[2][0] + A[3][2]*B[3][0]) + C[1][2]*(A[0][2]*B[0][1] + A[1][2]*B[1][1] + A[2][2]*B[2][1] + A[3][2]*B[3][1]) 
	     + C[2][2]*(A[0][2]*B[0][2] + A[1][2]*B[1][2] + A[2][2]*B[2][2] + A[3][2]*B[3][2]) + C[3][2]*(A[0][2]*B[0][3] + A[1][2]*B[1][3] + A[2][2]*B[2][3] + A[3][2]*B[3][3]); 
D[2][3] =  C[0][3]*(A[0][2]*B[0][0] + A[1][2]*B[1][0] + A[2][2]*B[2][0] + A[3][2]*B[3][0]) + C[1][3]*(A[0][2]*B[0][1] + A[1][2]*B[1][1] + A[2][2]*B[2][1] + A[3][2]*B[3][1]) 
	     + C[2][3]*(A[0][2]*B[0][2] + A[1][2]*B[1][2] + A[2][2]*B[2][2] + A[3][2]*B[3][2]) + C[3][3]*(A[0][2]*B[0][3] + A[1][2]*B[1][3] + A[2][2]*B[2][3] + A[3][2]*B[3][3]);

D[3][0] =  C[0][0]*(A[0][3]*B[0][0] + A[1][3]*B[1][0] + A[2][3]*B[2][0] + A[3][3]*B[3][0]) + C[1][0]*(A[0][3]*B[0][1] + A[1][3]*B[1][1] + A[2][3]*B[2][1] + A[3][3]*B[3][1]) 
         + C[2][0]*(A[0][3]*B[0][2] + A[1][3]*B[1][2] + A[2][3]*B[2][2] + A[3][3]*B[3][2]) + C[3][0]*(A[0][3]*B[0][3] + A[1][3]*B[1][3] + A[2][3]*B[2][3] + A[3][3]*B[3][3]); 
D[3][1] =  C[0][1]*(A[0][3]*B[0][0] + A[1][3]*B[1][0] + A[2][3]*B[2][0] + A[3][3]*B[3][0]) + C[1][1]*(A[0][3]*B[0][1] + A[1][3]*B[1][1] + A[2][3]*B[2][1] + A[3][3]*B[3][1]) 
	     + C[2][1]*(A[0][3]*B[0][2] + A[1][3]*B[1][2] + A[2][3]*B[2][2] + A[3][3]*B[3][2]) + C[3][1]*(A[0][3]*B[0][3] + A[1][3]*B[1][3] + A[2][3]*B[2][3] + A[3][3]*B[3][3]); 
D[3][2] =  C[0][2]*(A[0][3]*B[0][0] + A[1][3]*B[1][0] + A[2][3]*B[2][0] + A[3][3]*B[3][0]) + C[1][2]*(A[0][3]*B[0][1] + A[1][3]*B[1][1] + A[2][3]*B[2][1] + A[3][3]*B[3][1]) 
	     + C[2][2]*(A[0][3]*B[0][2] + A[1][3]*B[1][2] + A[2][3]*B[2][2] + A[3][3]*B[3][2]) + C[3][2]*(A[0][3]*B[0][3] + A[1][3]*B[1][3] + A[2][3]*B[2][3] + A[3][3]*B[3][3]); 
D[3][3] =  C[0][3]*(A[0][3]*B[0][0] + A[1][3]*B[1][0] + A[2][3]*B[2][0] + A[3][3]*B[3][0]) + C[1][3]*(A[0][3]*B[0][1] + A[1][3]*B[1][1] + A[2][3]*B[2][1] + A[3][3]*B[3][1]) 
	     + C[2][3]*(A[0][3]*B[0][2] + A[1][3]*B[1][2] + A[2][3]*B[2][2] + A[3][3]*B[3][2]) + C[3][3]*(A[0][3]*B[0][3] + A[1][3]*B[1][3] + A[2][3]*B[2][3] + A[3][3]*B[3][3]);

}	

void VecArrayArray(double localval[4], double v1,double v2,double v3,double v4,double A[4][4],double B[4][4])
{
           localval[0] =  v1*(A[0][0]*B[0][0] + A[0][1]*B[1][0] + A[0][2]*B[2][0] + A[0][3]*B[3][0]) 
                         + v2*(A[1][0]*B[0][0] + A[1][1]*B[1][0] + A[1][2]*B[2][0] + A[1][3]*B[3][0]) 
                         + v3*(A[2][0]*B[0][0] + A[2][1]*B[1][0] + A[2][2]*B[2][0] + A[2][3]*B[3][0]) 
                         + v4*(A[3][0]*B[0][0] + A[3][1]*B[1][0] + A[3][2]*B[2][0] + A[3][3]*B[3][0]); 
            localval[1] =  v1*(A[0][0]*B[0][1] + A[0][1]*B[1][1] + A[0][2]*B[2][1] + A[0][3]*B[3][1])
                         + v2*(A[1][0]*B[0][1] + A[1][1]*B[1][1] + A[1][2]*B[2][1] + A[1][3]*B[3][1]) 
                         + v3*(A[2][0]*B[0][1] + A[2][1]*B[1][1] + A[2][2]*B[2][1] + A[2][3]*B[3][1]) 
                         + v4*(A[3][0]*B[0][1] + A[3][1]*B[1][1] + A[3][2]*B[2][1] + A[3][3]*B[3][1]);
            localval[2] =  v1*(A[0][0]*B[0][2] + A[0][1]*B[1][2] + A[0][2]*B[2][2] + A[0][3]*B[3][2]) 
                         + v2*(A[1][0]*B[0][2] + A[1][1]*B[1][2] + A[1][2]*B[2][2] + A[1][3]*B[3][2]) 
                         + v3*(A[2][0]*B[0][2] + A[2][1]*B[1][2] + A[2][2]*B[2][2] + A[2][3]*B[3][2]) 
                         + v4*(A[3][0]*B[0][2] + A[3][1]*B[1][2] + A[3][2]*B[2][2] + A[3][3]*B[3][2]); 
            localval[3] =  v1*(A[0][0]*B[0][3] + A[0][1]*B[1][3] + A[0][2]*B[2][3] + A[0][3]*B[3][3]) 
                         + v2*(A[1][0]*B[0][3] + A[1][1]*B[1][3] + A[1][2]*B[2][3] + A[1][3]*B[3][3])
                         + v3*(A[2][0]*B[0][3] + A[2][1]*B[1][3] + A[2][2]*B[2][3] + A[2][3]*B[3][3]) 
                         + v4*(A[3][0]*B[0][3] + A[3][1]*B[1][3] + A[3][2]*B[2][3] + A[3][3]*B[3][3]);   
 
}


void VecArrayVec(double lastvalstiff[1], double v1,double v2,double v3,double v4,double A[4][4])
{
            lastvalstiff[0] +=v1*(A[0][0]*v1 + A[1][0]*v2 + A[2][0]*v3 + A[3][0]*v4) 
                            + v2*(A[0][1]*v1 + A[1][1]*v2 + A[2][1]*v3 + A[3][1]*v4) 
                            + v3*(A[0][2]*v1 + A[1][2]*v2 + A[2][2]*v3 + A[3][2]*v4) 
                            + v4*(A[0][3]*v1 + A[1][3]*v2 + A[2][3]*v3 + A[3][3]*v4);           
}

void VecArrayVecs(double lastvalstiff[1], double v1,double v2,double v3,double v4,double A[4][4], double w1,double w2,double w3,double w4)
{
            lastvalstiff[0] +=v1*(A[0][0]*w1 + A[1][0]*w2 + A[2][0]*w3 + A[3][0]*w4) 
                            + v2*(A[0][1]*w1 + A[1][1]*w2 + A[2][1]*w3 + A[3][1]*w4) 
                            + v3*(A[0][2]*w1 + A[1][2]*w2 + A[2][2]*w3 + A[3][2]*w4) 
                            + v4*(A[0][3]*w1 + A[1][3]*w2 + A[2][3]*w3 + A[3][3]*w4);           
}


void VecMatVec(double *r,MATRIX *MatrixA,double *c)
{
   double *tmp,t;
   tmp=malloc(sizeof(double)*MatrixA->N_Rows);
   MatrixDotVec(MatrixA, r,tmp);
   t = VecDotVec(r,tmp, MatrixA->N_Rows);
   c[0]=t;
   free(tmp);
   //return 12;
}

void  VecMatVecs(double *r,MATRIX *MatrixA,double *s,double *c)
{
   double *tmp,t;
   tmp=malloc(sizeof(double)*MatrixA->N_Rows);
   MatrixDotVec(MatrixA, s,tmp);
   VecDotVecs(r,tmp, MatrixA->N_Rows,c);
   //c[0]=t;
   //return 12;
}


double
VecDotVec(double *a, double *b, int size)
{
	double value=0;
	int i;
	for(i=0;i<size;i++)
	{
	  value = value + a[i]*b[i];
	}
      	return value;
}



void
VecDotVecs(double *a, double *b, int size,double *c)
{
	double value=0;
	int i;
	for(i=0;i<size;i++)
	{
	  value = value + a[i]*b[i];
	}
        c[0] = value;
       //	return value;
}


void ShowVec(double *p, int size)
{
	printf("\n");
	int i=0;
	for(;i<size;i++)
	{   
		printf("%20.15f\n",*(p++));
	}
	printf("\n");
 
}

void AssignVec(double *p, double *q, int size)
{
  int i;
  for(i=0;i<size;i++)
  {
    p[i] = q[i];   
  }  
}
//do the transpose matrix: B=A^T
void TransposeMatrix(MATRIX *A, MATRIX *B)
{
  int i,j,k, curr; 
  
  B->N_Rows = A->N_Columns;
  B->N_Columns = A->N_Rows;
  B->N_Entries = A->N_Entries;
  int n_rows, n_columns;
  n_rows = B->N_Rows;
  n_columns = B->N_Columns;
  B->RowPtr = malloc((n_rows+1)*sizeof(int));
  int *B_rowptr, *A_rowptr, *A_kcol, *B_kcol;
  B_rowptr = B->RowPtr;
  A_rowptr = A->RowPtr;
  A_kcol = A->KCol;
  memset(B_rowptr,0,(n_rows+1)*sizeof(int));
  
  int start, end;
  for(i=0;i<n_columns;i++)
  {
    start = A_rowptr[i];
    end = A_rowptr[i+1];
    for(curr=start;curr<end;curr++)
    {
      B_rowptr[A_kcol[curr]] += 1;
    }//end for curr
  }//end for i
  int cur_pos, tmp;
  cur_pos = 0;  

  for(i=0;i<n_rows;i++)
  {
    tmp = cur_pos;
    cur_pos += B_rowptr[i];
    B_rowptr[i] = tmp;
  }
   B_rowptr[n_rows] = cur_pos;
   B->KCol = malloc(cur_pos*sizeof(int));
   B->Entries = malloc(cur_pos*sizeof(double));
   B_kcol = B->KCol;
   double *A_entries, *B_entries;
   A_entries = A->Entries;
   B_entries = B->Entries;
   int *pos;
   pos = malloc(n_rows*sizeof(int));
   memset(pos,0,n_rows*sizeof(int));
   for(i=0;i<n_columns;i++ )
   {
     start = A_rowptr[i];
     end = A_rowptr[i+1];
     for(curr=start;curr<end;curr++)
     {
       //get the number of row in B
       j = A_kcol[curr];      
       //get the position in B_kcol
       cur_pos = B_rowptr[j] + pos[j];      
       B_kcol[cur_pos] = i;       
       B_entries[cur_pos] = A_entries[curr];
       pos[j] += 1;
    }//end for curr     
  }//end for i   
  free(pos);
}



void MatMatMat(double *a,double *b,double *d)//a'*b*a
{
 d[0]=a[0]*(a[0]*b[0] + a[12]*b[12] + a[4]*b[4] + a[8]*b[8]) + a[12]*(a[12]*b[15] + a[0]*b[3] + a[8]*b[11] + a[4]*b[7]) + a[8]*(a[0]*b[2] + a[12]*b[14] + a[8]*b[10] + a[4]*b[6]) + a[4]*(a[0]*b[1] + a[12]*b[13] + a[4]*b[5] + a[8]*b[9]);
 d[1]=a[13]*(a[12]*b[15] + a[0]*b[3] + a[8]*b[11] + a[4]*b[7]) + a[1]*(a[0]*b[0] + a[12]*b[12] + a[4]*b[4] + a[8]*b[8]) + a[9]*(a[0]*b[2] + a[12]*b[14] + a[8]*b[10] + a[4]*b[6]) + a[5]*(a[0]*b[1] + a[12]*b[13] + a[4]*b[5] + a[8]*b[9]);
 d[2]=a[10]*(a[0]*b[2] + a[12]*b[14] + a[8]*b[10] + a[4]*b[6]) + a[14]*(a[12]*b[15] + a[0]*b[3] + a[8]*b[11] + a[4]*b[7]) + a[2]*(a[0]*b[0] + a[12]*b[12] + a[4]*b[4] + a[8]*b[8]) + a[6]*(a[0]*b[1] + a[12]*b[13] + a[4]*b[5] + a[8]*b[9]);
 d[3]=a[11]*(a[0]*b[2] + a[12]*b[14] + a[8]*b[10] + a[4]*b[6]) + a[15]*(a[12]*b[15] + a[0]*b[3] + a[8]*b[11] + a[4]*b[7]) + a[3]*(a[0]*b[0] + a[12]*b[12] + a[4]*b[4] + a[8]*b[8]) + a[7]*(a[0]*b[1] + a[12]*b[13] + a[4]*b[5] + a[8]*b[9]);
 d[4]=a[0]*(a[1]*b[0] + a[13]*b[12] + a[5]*b[4] + a[9]*b[8]) + a[12]*(a[13]*b[15] + a[1]*b[3] + a[9]*b[11] + a[5]*b[7]) + a[8]*(a[13]*b[14] + a[1]*b[2] + a[9]*b[10] + a[5]*b[6]) + a[4]*(a[1]*b[1] + a[13]*b[13] + a[5]*b[5] + a[9]*b[9]);
 d[5]=a[13]*(a[13]*b[15] + a[1]*b[3] + a[9]*b[11] + a[5]*b[7]) + a[1]*(a[1]*b[0] + a[13]*b[12] + a[5]*b[4] + a[9]*b[8]) + a[9]*(a[13]*b[14] + a[1]*b[2] + a[9]*b[10] + a[5]*b[6]) + a[5]*(a[1]*b[1] + a[13]*b[13] + a[5]*b[5] + a[9]*b[9]);
 d[6]= a[10]*(a[13]*b[14] + a[1]*b[2] + a[9]*b[10] + a[5]*b[6]) + a[14]*(a[13]*b[15] + a[1]*b[3] + a[9]*b[11] + a[5]*b[7]) + a[2]*(a[1]*b[0] + a[13]*b[12] + a[5]*b[4] + a[9]*b[8]) + a[6]*(a[1]*b[1] + a[13]*b[13] + a[5]*b[5] + a[9]*b[9]);
 d[7]=a[11]*(a[13]*b[14] + a[1]*b[2] + a[9]*b[10] + a[5]*b[6]) + a[15]*(a[13]*b[15] + a[1]*b[3] + a[9]*b[11] + a[5]*b[7]) + a[3]*(a[1]*b[0] + a[13]*b[12] + a[5]*b[4] + a[9]*b[8]) + a[7]*(a[1]*b[1] + a[13]*b[13] + a[5]*b[5] + a[9]*b[9]);
 d[8]=a[0]*(a[2]*b[0] + a[14]*b[12] + a[10]*b[8] + a[6]*b[4]) + a[12]*(a[10]*b[11] + a[14]*b[15] + a[2]*b[3] + a[6]*b[7]) + a[8]*(a[10]*b[10] + a[14]*b[14] + a[2]*b[2] + a[6]*b[6]) + a[4]*(a[14]*b[13] + a[2]*b[1] + a[10]*b[9] + a[6]*b[5]);
 d[9]=a[13]*(a[10]*b[11] + a[14]*b[15] + a[2]*b[3] + a[6]*b[7]) + a[1]*(a[2]*b[0] + a[14]*b[12] + a[10]*b[8] + a[6]*b[4]) + a[9]*(a[10]*b[10] + a[14]*b[14] + a[2]*b[2] + a[6]*b[6]) + a[5]*(a[14]*b[13] + a[2]*b[1] + a[10]*b[9] + a[6]*b[5]);
 d[10]=a[10]*(a[10]*b[10] + a[14]*b[14] + a[2]*b[2] + a[6]*b[6]) + a[14]*(a[10]*b[11] + a[14]*b[15] + a[2]*b[3] + a[6]*b[7]) + a[2]*(a[2]*b[0] + a[14]*b[12] + a[10]*b[8] + a[6]*b[4]) + a[6]*(a[14]*b[13] + a[2]*b[1] + a[10]*b[9] + a[6]*b[5]);
 d[11]=a[11]*(a[10]*b[10] + a[14]*b[14] + a[2]*b[2] + a[6]*b[6]) + a[15]*(a[10]*b[11] + a[14]*b[15] + a[2]*b[3] + a[6]*b[7]) + a[3]*(a[2]*b[0] + a[14]*b[12] + a[10]*b[8] + a[6]*b[4]) + a[7]*(a[14]*b[13] + a[2]*b[1] + a[10]*b[9] + a[6]*b[5]);
 d[12]=a[0]*(a[15]*b[12] + a[3]*b[0] + a[11]*b[8] + a[7]*b[4]) + a[12]*(a[11]*b[11] + a[15]*b[15] + a[3]*b[3] + a[7]*b[7]) + a[8]*(a[11]*b[10] + a[15]*b[14] + a[3]*b[2] + a[7]*b[6]) + a[4]*(a[15]*b[13] + a[3]*b[1] + a[11]*b[9] + a[7]*b[5]);
 d[13]=a[13]*(a[11]*b[11] + a[15]*b[15] + a[3]*b[3] + a[7]*b[7]) + a[1]*(a[15]*b[12] + a[3]*b[0] + a[11]*b[8] + a[7]*b[4]) + a[9]*(a[11]*b[10] + a[15]*b[14] + a[3]*b[2] + a[7]*b[6]) + a[5]*(a[15]*b[13] + a[3]*b[1] + a[11]*b[9] + a[7]*b[5]);
 d[14]=a[10]*(a[11]*b[10] + a[15]*b[14] + a[3]*b[2] + a[7]*b[6]) + a[14]*(a[11]*b[11] + a[15]*b[15] + a[3]*b[3] + a[7]*b[7]) + a[2]*(a[15]*b[12] + a[3]*b[0] + a[11]*b[8] + a[7]*b[4]) + a[6]*(a[15]*b[13] + a[3]*b[1] + a[11]*b[9] + a[7]*b[5]);
 d[15]=a[11]*(a[11]*b[10] + a[15]*b[14] + a[3]*b[2] + a[7]*b[6]) + a[15]*(a[11]*b[11] + a[15]*b[15] + a[3]*b[3] + a[7]*b[7]) + a[3]*(a[15]*b[12] + a[3]*b[0] + a[11]*b[8] + a[7]*b[4]) + a[7]*(a[15]*b[13] + a[3]*b[1] + a[11]*b[9] + a[7]*b[5]);
}


/*
double Minimum(double a, double b)
{
 if(a>b)
    return a;
 else 
    return b;
}
*/



















