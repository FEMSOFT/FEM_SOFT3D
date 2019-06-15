/*
 * =====================================================================================
 *
 *       Filename:  Multilevel.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2014年12月03日 14时35分15秒
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
   #include "Matrix.h"
   #include "Fespace3D.h"
   #include "DiscreteForm3D.h"
   #include "Base3D.h"
   #include "Boundary.h"
   #include "Rhst.h"
   #include "MG.h"
   #include "GMG.h"
   #include "EigenSolver.h"
   #include "Fefunction3D.h"

//Set the matrix entries 
void SetSubMatrix(MATRIX *A,double *AA, int n_row,int *row,int n_col,int *col)
{
    int i,j,k;
    int *RowPtr, *KCol;
    double *Entries;
    RowPtr = A->RowPtr;
    KCol = A->KCol;
    Entries = A->Entries;
    int start, end, length;
    for(i=0;i<n_row;i++)
    {
        start = RowPtr[row[i]];
        end = RowPtr[row[i]+1];
        length = end-start;
        for(j=0;j<n_col;j++)
        {	    
            for(k=0;k<length;k++)
            {
                if(KCol[start+k] == col[j]){
                    Entries[start+k] = AA[i*n_row+j];
		    k = length+1;
		}                
            }//end for k
        }//end for j
    }// end for i
}
//产生相对的矩阵（细单元在大单元中的局部坐标）
void ProduceRelativeElem3D(ELEMENT3D *coarse_elem3D, ELEMENT3D *finer_elem3D, ELEMENT3D *elem3D)
{
  int i;
  double dx, dy,dz;  
  for(i=0;i<4;i++)
  {
    dx = finer_elem3D->Vert_X[i] - coarse_elem3D->Vert_X[0];
    dy = finer_elem3D->Vert_Y[i] - coarse_elem3D->Vert_Y[0];
    dz = finer_elem3D->Vert_Z[i] - coarse_elem3D->Vert_Z[0];
    elem3D->Vert_X[i] = coarse_elem3D->Inv_Jacobian[0][0]*dx + coarse_elem3D->Inv_Jacobian[0][1]*dy
                      + coarse_elem3D->Inv_Jacobian[0][2]*dz;
    elem3D->Vert_Y[i] = coarse_elem3D->Inv_Jacobian[1][0]*dx + coarse_elem3D->Inv_Jacobian[1][1]*dy
                      + coarse_elem3D->Inv_Jacobian[1][2]*dz;
    elem3D->Vert_Z[i] = coarse_elem3D->Inv_Jacobian[2][0]*dx + coarse_elem3D->Inv_Jacobian[2][1]*dy
                      + coarse_elem3D->Inv_Jacobian[2][2]*dz;
  }
}


//产生相对的矩阵（细单元在大单元中的局部坐标）
void ProduceRelativeCoord(ELEMENT3D *coarse_elem3D, double dx, double dy, double dz, double *elem3D)
{
  int i;
    elem3D[0] = coarse_elem3D->Inv_Jacobian[0][0]*dx + coarse_elem3D->Inv_Jacobian[0][1]*dy
                      + coarse_elem3D->Inv_Jacobian[0][2]*dz;
    elem3D[1] = coarse_elem3D->Inv_Jacobian[1][0]*dx + coarse_elem3D->Inv_Jacobian[1][1]*dy
                      + coarse_elem3D->Inv_Jacobian[1][2]*dz;
    elem3D[2] = coarse_elem3D->Inv_Jacobian[2][0]*dx + coarse_elem3D->Inv_Jacobian[2][1]*dy
                      + coarse_elem3D->Inv_Jacobian[2][2]*dz;
}






//Get the local prolong: 目前之提供仿射等价有限元函数的实现
void GetLocalProlong(ELEMENT3D *coarse_elem3D,BASEFUNCTION3D *coarse_base,int n_coarse_base,
		     ELEMENT3D *finer_elem3D, BASEFUNCTION3D *finer_base,int n_finer_base,
		     ELEMENT3D *elem3D,double *Local_Prolong)
{
  //Get Local Matrix
  double *nodal_points = coarse_base->Nodal_Points;
  int i; 
  double Quad_xi, Quad_eta, Quad_zeta, Quad_x, Quad_y, Quad_z;
  Quad_x = 0.0;
  Quad_y = 0.0;
  Quad_z = 0.0;
  //coarse base 里面的每一个基函数在finer 单元上的插值作为Prolong的元素
  if((coarse_base->Maptype==Affine)&&(finer_base->Maptype==Affine))
  {
    ProduceRelativeElem3D(coarse_elem3D,finer_elem3D,elem3D);
    for(i=0;i<n_finer_base;i++)
    {    
     //GetElementCoord(elem3D,nodal_points[3*i],nodal_points[3*i+1],nodal_points[3*i+2],
     //                &Quad_xi,&Quad_eta,&Quad_zeta);  
     //GetBaseValues(coarse_base,D000,elem3D,Quad_x,Quad_y,Quad_z,Quad_xi,Quad_eta,Quad_zeta,
     //	           Local_Prolong+i*n_coarse_base);      
      GetBaseValues(coarse_base,D000,elem3D,Quad_x,Quad_y,Quad_z,elem3D->Vert_X[i],elem3D->Vert_Y[i],
	            elem3D->Vert_Z[i],Local_Prolong+i*n_coarse_base);     
    
    } //end for i
  }//end if
}




void BuildProlongFa(Fespace3D *coarse_space,Fespace3D *finer_space, MATRIX *Prolong)
{
  printf("Produce the interpolation and restriction operators\n");
  int i, j, k, volu_ind, start, end;
  BASEFUNCTION3D *coarse_base, *finer_base;
  coarse_base = coarse_space->Base;
  finer_base = finer_space->Base;
  int n_coarse_base, n_finer_base;
  n_coarse_base = coarse_base->Num_Bas;
  n_finer_base = finer_base->Num_Bas;
  //printf ( "nnnnnnnnnnnnnnnnnnnnnnnnnnn\n" );
  MESH *coarse_mesh, *finer_mesh;
  coarse_mesh = coarse_space->Mesh;
  finer_mesh = finer_space->Mesh;
  VOLU **finer_volus;
  VOLU **coarse_volus;
  coarse_volus = coarse_mesh->Volus;
  finer_volus = finer_mesh->Volus;
  VOLU *coarse_volu, *finer_volu;
  //printf ( "bbbbbbbbbbbbbbbbbbbbbb\n" );
  int *coarse_begin_index, *finer_begin_index, *coarse_global_numbers, *finer_global_numbers;
  coarse_begin_index = coarse_space->BeginIndex;
  coarse_global_numbers = coarse_space->GlobalNumbers;
  
  finer_begin_index = finer_space->BeginIndex;
  finer_global_numbers = finer_space->GlobalNumbers;
  //printf ( "ccccccccccccccccccccccc\n" );
  double *Local_Prolong;
  Local_Prolong = malloc(n_coarse_base*n_finer_base*sizeof(double));
  
  int finer_num_volus;
  finer_num_volus = finer_mesh->Num_Volus_Global;
  //printf ( "ddddddddddddddddddddd\n" );
  //首先建立矩阵的结构（稀疏矩阵）
  //Prolong=malloc(sizeof(MATRIX));
  Prolong->N_Rows = finer_space->DOF_ALL;
  //printf ( "sssssssssssssssss\n" );
  Prolong->N_Columns = coarse_space->DOF_ALL;
  //printf ( "dddddddddddeeeeeeeeeeeee\n" );
  Prolong->RowPtr = malloc((Prolong->N_Rows+1)*sizeof(int));
  //printf ( "eeeeeeeeeeeeeeeeeeeee\n" );
  int n_rows = Prolong->N_Rows;
  //memset(Prolong->RowPtr,0,(n_rows+1)*sizeof(int));
  int *rowptr, *rows, *columns;
  rowptr = Prolong->RowPtr;
  //build the information for RowPtr:
  //printf ( "mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm\n" );
  /*for(volu_ind=0;volu_ind<finer_num_volus;volu_ind++)
  {
    rows = finer_global_numbers + finer_begin_index[volu_ind];
    for(i=0;i<n_finer_base;i++)
    {
      rowptr[rows[i]] = n_coarse_base;      
    }
  }
  printf ( "nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn\n" );
  int cur_pos, tmp;
  cur_pos = 0;
  for(i=0;i<n_rows;i++)
  {
    //tmp = cur_pos;
    cur_pos += rowptr[i];
    rowptr[i] = tmp;    
  }
  */
  for(i=0;i<n_rows+1;i++)
  {
    rowptr[i] = i*n_coarse_base;    
  }

  //rowptr[n_rows] = cur_pos;
  int n_entries;
  n_entries = rowptr[n_rows];
  Prolong->KCol = malloc(n_entries*sizeof(int));
  Prolong->N_Entries = n_entries;
  Prolong->Entries = malloc(n_entries*sizeof(double));
  int *kcol;
  double *entries;
  kcol = Prolong->KCol;
  entries = Prolong->Entries;
  int row_ind, father;

  ELEMENT3D *coarse_elem3D;
  coarse_elem3D = malloc(sizeof(ELEMENT3D));
  InitialElem3D(coarse_mesh,coarse_elem3D);
  ELEMENT3D *finer_elem3D;
  finer_elem3D = malloc(sizeof(ELEMENT3D));
  InitialElem3D(finer_mesh,finer_elem3D);
  //FreeElem3D(coarse_elem3D);
  //printf ( "uuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu\n" );
 for(volu_ind=0;volu_ind<finer_num_volus;volu_ind++)
  {
    finer_volu = finer_volus[volu_ind];
    father = finer_volu->Father;
    rows = finer_global_numbers+finer_begin_index[volu_ind];
    columns = coarse_global_numbers+coarse_begin_index[father];
    for(i=0;i<n_finer_base;i++)
    {
      row_ind = rows[i];
      //printf ( "aaaa\n" );
      start = rowptr[row_ind];
      //printf ( "bbbbb\n" );
      for(j=0;j<n_coarse_base;j++)
      {
        kcol[start+j] = columns[j];
      }
      //对这一行的编号进行排序
      QuickSort_Int(kcol,start,start+n_coarse_base-1);
    }
  }
  //这样就获得了矩阵Prolong的结构  
  //下面来具体形成Prolong矩阵  
  ELEMENT3D *elem3D;
  elem3D = malloc(sizeof(ELEMENT3D));
  InitialElem3D(coarse_mesh,elem3D);
  //printf("Start the interpolation\n");
  for(volu_ind=0;volu_ind<finer_num_volus;volu_ind++)
  {
  finer_volu = finer_volus[volu_ind];
  father = finer_volu->Father;
  coarse_volu = coarse_volus[father];
  GetElement3D(coarse_mesh,father,coarse_elem3D);
  GetElement3D(finer_mesh,volu_ind,finer_elem3D);
  //获得局部的扩展算子
  GetLocalProlong(coarse_elem3D,coarse_base,n_coarse_base,finer_elem3D,finer_base,n_finer_base,
	          elem3D,Local_Prolong);
  //把获得到局部的扩展算子设置到全局的扩张算子
  SetSubMatrix(Prolong,Local_Prolong, n_finer_base,finer_global_numbers+finer_begin_index[volu_ind],
		 n_coarse_base,coarse_global_numbers+coarse_begin_index[father]);	
  }
  free(Local_Prolong);
  FreeElem3D(coarse_elem3D);
  FreeElem3D(finer_elem3D);
  FreeElem3D(elem3D);  
}


void BuildRestrictFa(Fespace3D *coarse_spsace,Fespace3D *finer_space, MATRIX *Restrict)
{
  //printf("Produce the interpolation and restriction operators\n");
  MATRIX *Prolong;
  Prolong=malloc(sizeof(MATRIX));
  int i, j, k, volu_ind, start, end;
  BASEFUNCTION3D *coarse_base, *finer_base;
  coarse_base = coarse_spsace->Base;
  finer_base = finer_space->Base;
  int n_coarse_base, n_finer_base;
  n_coarse_base = coarse_base->Num_Bas;
  n_finer_base = finer_base->Num_Bas;
  
  MESH *coarse_mesh, *finer_mesh;
  coarse_mesh = coarse_spsace->Mesh;
  //ShowMesh(coarse_mesh);
  finer_mesh = finer_space->Mesh;
  VOLU **finer_volus;
  VOLU **coarse_volus;
  coarse_volus = coarse_mesh->Volus;
  finer_volus = finer_mesh->Volus;
  VOLU *coarse_volu, *finer_volu;
  
  int *coarse_begin_index, *finer_begin_index, *coarse_global_numbers, *finer_global_numbers;
  coarse_begin_index = coarse_spsace->BeginIndex;
  coarse_global_numbers = coarse_spsace->GlobalNumbers;
  
  finer_begin_index = finer_space->BeginIndex;
  finer_global_numbers = finer_space->GlobalNumbers;
  
  double *Local_Prolong;
  Local_Prolong = malloc(n_coarse_base*n_finer_base*sizeof(double));
  
  int finer_num_volus;
  finer_num_volus = finer_mesh->Num_Volus_Global;
  //first, we build the struct of the matrix Prolong
  //首先建立矩阵的结构（稀疏矩阵）
  Prolong->N_Rows = finer_space->DOF_ALL;
  Prolong->N_Columns = coarse_spsace->DOF_ALL;
  Prolong->RowPtr = malloc((Prolong->N_Rows+1)*sizeof(int));
  int n_rows = Prolong->N_Rows;
  //memset(Prolong->RowPtr,0,(n_rows+1)*sizeof(int));
  int *rowptr, *rows, *columns;
  rowptr = Prolong->RowPtr;
  //build the information for RowPtr:
  for(i=0;i<n_rows+1;i++)
  {
    rowptr[i] = i*n_coarse_base;    
  }
  
  int n_entries;
  n_entries = rowptr[n_rows];
  Prolong->KCol = malloc(n_entries*sizeof(int));
  Prolong->N_Entries = n_entries;
  Prolong->Entries = malloc(n_entries*sizeof(double));
  int *kcol;
  double *entries;
  kcol = Prolong->KCol;
  entries = Prolong->Entries;
  int row_ind, father;
  ELEMENT3D *coarse_elem3D, *finer_elem3D;
  coarse_elem3D = malloc(sizeof(ELEMENT3D));
  finer_elem3D = malloc(sizeof(ELEMENT3D));
  
  //ShowMesh(finer_mesh);
  InitialElem3D(coarse_mesh,coarse_elem3D);
  InitialElem3D(finer_mesh,finer_elem3D);
  //Now we build the information for KCOL
  for(volu_ind=0;volu_ind<finer_num_volus;volu_ind++)
  {
    finer_volu = finer_volus[volu_ind];
    father = finer_volu->Father;
    rows = finer_global_numbers+finer_begin_index[volu_ind];
    columns = coarse_global_numbers+coarse_begin_index[father];
    for(i=0;i<n_finer_base;i++)
    {
      row_ind = rows[i];
      start = rowptr[row_ind];
      for(j=0;j<n_coarse_base;j++)
      {
        kcol[start+j] = columns[j];
      }
      //对这一行的编号进行排序
      QuickSort_Int(kcol,start,start+n_coarse_base-1);
    }
  }
  //printf("33333333333\n");
  //这样就获得了矩阵Prolong的结构  
  //下面来具体形成Prolong矩阵  
  //iteration on the mesh faces  
  ELEMENT3D *elem3D;
  elem3D = malloc(sizeof(ELEMENT3D));
  elem3D->Num_Verts = 4;
  InitialElem3D(coarse_mesh,elem3D);
  //printf("Start the interpolation\n");
 for(volu_ind=0;volu_ind<finer_num_volus;volu_ind++)
 {
  finer_volu = finer_volus[volu_ind];
  father = finer_volu->Father;
  coarse_volu = coarse_volus[father];
  GetElement3D(coarse_mesh,father,coarse_elem3D);
  GetElement3D(finer_mesh,volu_ind,finer_elem3D);
  //获得局部的扩展算子
  GetLocalProlong(coarse_elem3D,coarse_base,n_coarse_base,finer_elem3D,finer_base,n_finer_base,elem3D,Local_Prolong);
  //把获得到局部的扩展算子设置到全局的扩张算子
  SetSubMatrix(Prolong,Local_Prolong, n_finer_base,finer_global_numbers+finer_begin_index[volu_ind],
		 n_coarse_base,coarse_global_numbers+coarse_begin_index[father]);	
 }  
 //get the restriction matrix according to the prolong matrix
 TransposeMatrix(Prolong, Restrict);
 FreeMatrix(Prolong);
 free(Local_Prolong);
 FreeElem3D(coarse_elem3D);
 FreeElem3D(finer_elem3D);
 FreeElem3D(elem3D);  
}

void BuildProlongAn(Fespace3D *coarse_spsace,Fespace3D *finer_space, MATRIX *Prolong)
{
  //printf("Produce the interpolation and restriction operators\n");
  int i, j, k, volu_ind, start, end;
  BASEFUNCTION3D *coarse_base, *finer_base;
  coarse_base = coarse_spsace->Base;
  finer_base = finer_space->Base;
  int n_coarse_base, n_finer_base;
  n_coarse_base = coarse_base->Num_Bas;
  n_finer_base = finer_base->Num_Bas;
  
  MESH *coarse_mesh, *finer_mesh;
  coarse_mesh = coarse_spsace->Mesh;
  //ShowMesh(coarse_mesh);
  finer_mesh = finer_space->Mesh;
  VOLU **finer_volus;
  VOLU **coarse_volus;
  coarse_volus = coarse_mesh->Volus;
  finer_volus = finer_mesh->Volus;
  VOLU *coarse_volu, *finer_volu;
  
  int *coarse_begin_index, *finer_begin_index, *coarse_global_numbers, *finer_global_numbers;
  coarse_begin_index = coarse_spsace->BeginIndex;
  coarse_global_numbers = coarse_spsace->GlobalNumbers;
  
  finer_begin_index = finer_space->BeginIndex;
  finer_global_numbers = finer_space->GlobalNumbers;
  
  double *Local_Prolong;
  Local_Prolong = malloc(n_coarse_base*n_finer_base*sizeof(double));
  
  int finer_num_volus;
  finer_num_volus = finer_mesh->Num_Volus_Global;
  //first, we build the struct of the matrix Prolong
  //首先建立矩阵的结构（稀疏矩阵）
  Prolong->N_Rows = finer_space->DOF_ALL;
  Prolong->N_Columns = coarse_spsace->DOF_ALL;
  Prolong->RowPtr = malloc((Prolong->N_Rows+1)*sizeof(int));
  int n_rows = Prolong->N_Rows;
  //memset(Prolong->RowPtr,0,(n_rows+1)*sizeof(int));
  int *rowptr, *rows, *columns;
  rowptr = Prolong->RowPtr;
  for(i=0;i<n_rows+1;i++)
  {
    rowptr[i] = i*n_coarse_base;    
  }
  int n_entries;
  n_entries = rowptr[n_rows];
  Prolong->KCol = malloc(n_entries*sizeof(int));
  Prolong->N_Entries = n_entries;
  Prolong->Entries = malloc(n_entries*sizeof(double));
  int *kcol;
  double *entries;
  kcol = Prolong->KCol;
  entries = Prolong->Entries;
  int row_ind, ancestor;
  ELEMENT3D *coarse_elem3D, *finer_elem3D;
  coarse_elem3D = malloc(sizeof(ELEMENT3D));
  finer_elem3D = malloc(sizeof(ELEMENT3D));
  
  //ShowMesh(finer_mesh);
  InitialElem3D(coarse_mesh,coarse_elem3D);
  InitialElem3D(finer_mesh,finer_elem3D);
  //Now we build the information for KCOL
  for(volu_ind=0;volu_ind<finer_num_volus;volu_ind++)
  {
    finer_volu = finer_volus[volu_ind];
    ancestor = finer_volu->Ancestor;
    rows = finer_global_numbers+finer_begin_index[volu_ind];
    columns = coarse_global_numbers+coarse_begin_index[ancestor];
    for(i=0;i<n_finer_base;i++)
    {
      row_ind = rows[i];
      start = rowptr[row_ind];
      for(j=0;j<n_coarse_base;j++)
      {
        kcol[start+j] = columns[j];
      }
      //对这一行的编号进行排序
      QuickSort_Int(kcol,start,start+n_coarse_base-1);
    }
  }
  //printf("33333333333\n");
  //这样就获得了矩阵Prolong的结构  
  //下面来具体形成Prolong矩阵  
  //iteration on the mesh faces  
  ELEMENT3D *elem3D;
  elem3D = malloc(sizeof(ELEMENT3D));
  InitialElem3D(coarse_mesh,elem3D);
 for(volu_ind=0;volu_ind<finer_num_volus;volu_ind++)
 {
  finer_volu = finer_volus[volu_ind];
  ancestor = finer_volu->Ancestor;
  coarse_volu = coarse_volus[ancestor];
  GetElement3D(coarse_mesh,ancestor,coarse_elem3D);
  GetElement3D(finer_mesh,volu_ind,finer_elem3D);
  //获得局部的扩展算子
  GetLocalProlong(coarse_elem3D,coarse_base,n_coarse_base,finer_elem3D,finer_base,n_finer_base,
	          elem3D,Local_Prolong);
  //把获得到局部的扩展算子设置到全局的扩张算子
  SetSubMatrix(Prolong,Local_Prolong, n_finer_base,finer_global_numbers+finer_begin_index[volu_ind],
		 n_coarse_base,coarse_global_numbers+coarse_begin_index[ancestor]);	
 }
 free(Local_Prolong);
 FreeElem3D(coarse_elem3D);
 FreeElem3D(finer_elem3D);
 FreeElem3D(elem3D);  
}


void BuildRestrictAn(Fespace3D *coarse_spsace,Fespace3D *finer_space, MATRIX *Restrict)
{
   //printf ( "qqqqqqqqqqqqqqqqqqqqqqqq\n" );
   //printf("Produce the interpolation and restriction operators\n");
  MATRIX *Prolong;
  Prolong=malloc(sizeof(MATRIX));
  int i, j, k, volu_ind, start, end;
  BASEFUNCTION3D *coarse_base, *finer_base;
  coarse_base = coarse_spsace->Base;
  finer_base = finer_space->Base;
  int n_coarse_base, n_finer_base;
  n_coarse_base = coarse_base->Num_Bas;
  n_finer_base = finer_base->Num_Bas;
  
  MESH *coarse_mesh, *finer_mesh;
  coarse_mesh = coarse_spsace->Mesh;
  finer_mesh = finer_space->Mesh;
  VOLU **finer_volus;
  VOLU **coarse_volus;
  coarse_volus = coarse_mesh->Volus;
  finer_volus = finer_mesh->Volus;
  VOLU *coarse_volu, *finer_volu;
  
  int *coarse_begin_index, *finer_begin_index, *coarse_global_numbers, *finer_global_numbers;
  coarse_begin_index = coarse_spsace->BeginIndex;
  coarse_global_numbers = coarse_spsace->GlobalNumbers;
  
  finer_begin_index = finer_space->BeginIndex;
  finer_global_numbers = finer_space->GlobalNumbers;
  
  double *Local_Prolong;
  Local_Prolong = malloc(n_coarse_base*n_finer_base*sizeof(double));
  
  int finer_num_volus;
  finer_num_volus = finer_mesh->Num_Volus_Global;
  //first, we build the struct of the matrix Prolong
  //首先建立矩阵的结构（稀疏矩阵）
  //printf ( "wwwwwwwwwwwwwwwwwwwwwwwwwwww\n" );
  Prolong->N_Rows = finer_space->DOF_ALL;
  Prolong->N_Columns = coarse_spsace->DOF_ALL;
  Prolong->RowPtr = malloc((Prolong->N_Rows+1)*sizeof(int));
  //printf ( "00000000000000000000\n" );
  int n_rows = Prolong->N_Rows;
  //memset(Prolong->RowPtr,0,(n_rows+1)*sizeof(int));
  int *rowptr, *rows, *columns;
  rowptr = Prolong->RowPtr;
  //printf ( "11111111111111111111\n" );
  //build the information for RowPtr:
  for(i=0;i<n_rows+1;++i)
     rowptr[i]=i*n_coarse_base;
  int n_entries;
  n_entries = rowptr[n_rows];
  Prolong->KCol = malloc(n_entries*sizeof(int));
  Prolong->N_Entries = n_entries;
  Prolong->Entries = malloc(n_entries*sizeof(double));
  int *kcol;
  double *entries;
  kcol = Prolong->KCol;
  entries = Prolong->Entries;
  int row_ind, ancestor;
  ELEMENT3D *coarse_elem3D, *finer_elem3D;
  coarse_elem3D = malloc(sizeof(ELEMENT3D));
  finer_elem3D = malloc(sizeof(ELEMENT3D));
  InitialElem3D(coarse_mesh,coarse_elem3D);
  InitialElem3D(finer_mesh,finer_elem3D);
  //Now we build the information for KCOL
  for(volu_ind=0;volu_ind<finer_num_volus;volu_ind++)
  {
    finer_volu = finer_volus[volu_ind];
    ancestor = finer_volu->Ancestor;
    rows = finer_global_numbers+finer_begin_index[volu_ind];
    columns = coarse_global_numbers+coarse_begin_index[ancestor];
    for(i=0;i<n_finer_base;i++)
    {
      row_ind = rows[i];
      start = rowptr[row_ind];
      for(j=0;j<n_coarse_base;j++)
      {
        kcol[start+j] = columns[j];
      }
      //对这一行的编号进行排序
      QuickSort_Int(kcol,start,start+n_coarse_base-1);
    }
  }
  //printf("33333333333\n");
  //这样就获得了矩阵Prolong的结构  
  //下面来具体形成Prolong矩阵  
  //iteration on the mesh faces  
  ELEMENT3D *elem3D;
  elem3D = malloc(sizeof(ELEMENT3D));
  InitialElem3D(coarse_mesh,elem3D);
 for(volu_ind=0;volu_ind<finer_num_volus;volu_ind++)
 {
  finer_volu = finer_volus[volu_ind];
  ancestor = finer_volu->Ancestor;
  coarse_volu = coarse_volus[ancestor];
  GetElement3D(coarse_mesh,ancestor,coarse_elem3D);
  GetElement3D(finer_mesh,volu_ind,finer_elem3D);
  //获得局部的扩展算子
  GetLocalProlong(coarse_elem3D,coarse_base,n_coarse_base,finer_elem3D,finer_base,n_finer_base,elem3D,Local_Prolong);
  //把获得到局部的扩展算子设置到全局的扩张算子
  SetSubMatrix(Prolong,Local_Prolong, n_finer_base,finer_global_numbers+finer_begin_index[volu_ind],
		 n_coarse_base,coarse_global_numbers+coarse_begin_index[ancestor]);	
 }  
 //get the restriction matrix according to the prolong matrix
 TransposeMatrix(Prolong, Restrict);
 FreeMatrix(Prolong);
 free(Local_Prolong);
 FreeElem3D(coarse_elem3D);
 FreeElem3D(finer_elem3D);
 FreeElem3D(elem3D);  
}




 //在原有矩阵结构的基础上形成矩阵结构
 MATRIX* ExtensionMatrix(MATRIX *matrix,int num)
 {
  MATRIX *matrixa;
  matrixa  = malloc(sizeof(MATRIX));
  matrixa->N_Rows = matrix->N_Rows+num;
  matrixa->N_Columns = matrix->N_Columns+num;
  matrixa->N_Entries = matrix->N_Entries+num*num+2*num*(matrix->N_Rows);
  int *RowPtr,*KCol;
  double *Entries; 
  RowPtr = malloc(sizeof(int)*(matrixa->N_Rows+1));
  KCol = malloc(sizeof(int)*matrixa->N_Entries); 
  Entries = malloc(sizeof(double)*matrixa->N_Entries); 
  int i,j;
  for(i=0;i<matrix->N_Columns+1;++i)
   RowPtr[i] = matrix->RowPtr[i]+i*num;
  for(i=0;i<num;++i)
   RowPtr[matrix->N_Columns+i+1] = RowPtr[matrix->N_Columns+i]+matrix->N_Columns+num;
  //更新列指标
  for(i=0;i<matrix->N_Rows;++i)
  { 
   for(j=0;j<matrix->RowPtr[i+1]-matrix->RowPtr[i];++j)
     {
       KCol[RowPtr[i]+j] = matrix->KCol[matrix->RowPtr[i]+j];
       Entries[RowPtr[i]+j] = matrix->Entries[matrix->RowPtr[i]+j];
     }
   for(j=0;j<num;++j)
      KCol[RowPtr[i+1]-num+j]=matrix->N_Columns+j;     
  }
  for(i=0;i<num;++i)
  {
   for(j=0;j<matrixa->N_Columns;++j)
      KCol[RowPtr[matrix->N_Rows+i]+j]=j;
  }
  matrixa->RowPtr = RowPtr;
  matrixa->KCol = KCol;
  matrixa->Entries = Entries;
  return matrixa;    
 }

//在原有矩阵结构上构建新的矩阵，这个函数和上个比不会复制矩阵元素
MATRIX* ExtensionMatrixSim(MATRIX *matrix,int num)
 {
  MATRIX *matrixa;
  matrixa  = malloc(1*sizeof(MATRIX));
  matrixa->N_Rows = matrix->N_Rows+num;
  matrixa->N_Columns = matrix->N_Columns+num;
  matrixa->N_Entries = matrix->N_Entries+num*num+2*num*(matrix->N_Rows);
  int *RowPtr,*KCol;
  double *Entries; 
  RowPtr = malloc(sizeof(int)*(matrixa->N_Rows+1));
  KCol = malloc(sizeof(int)*matrixa->N_Entries); 
  Entries = malloc(sizeof(double)*matrixa->N_Entries); 
  int i,j;
  for(i=0;i<matrix->N_Columns+1;++i)
   RowPtr[i] = matrix->RowPtr[i]+i*num;
  for(i=0;i<num;++i)
   RowPtr[matrix->N_Columns+i+1] = RowPtr[matrix->N_Columns+i]+matrix->N_Columns+num;
  //更新列指标
  for(i=0;i<matrix->N_Rows;++i)
  { 
   for(j=0;j<matrix->RowPtr[i+1]-matrix->RowPtr[i];++j)
     {
       KCol[RowPtr[i]+j] = matrix->KCol[matrix->RowPtr[i]+j];
     }
   for(j=0;j<num;++j)
      KCol[RowPtr[i+1]-num+j]=matrix->N_Columns+j;     
  }
  for(i=0;i<num;++i)
  {
   for(j=0;j<matrixa->N_Columns;++j)
      KCol[RowPtr[matrix->N_Rows+i]+j]=j;
  }
  matrixa->RowPtr = RowPtr;
  matrixa->KCol = KCol;
  matrixa->Entries = Entries;
  return matrixa;    
 }



//通过两层网格的延拓与限制算子计算任意两层之间的延拓与限制算子
 double* ProlongH2h(MULTILEVEL *multilevel,int num_H,int num_h,double *vector)
 {
  int i,j;
  double *tmp_vector1,*tmp_vector2;
  tmp_vector1 = malloc(sizeof(double)*multilevel->Stiff_Matrix[num_h]->N_Rows);
  tmp_vector2 = malloc(sizeof(double)*multilevel->Stiff_Matrix[num_h]->N_Rows);
  for(i=0;i<multilevel->Stiff_Matrix[num_H]->N_Rows;++i)
    tmp_vector1[i] = vector[i];
  for(i=num_H;i<num_h-1;++i)
   {
    MatrixDotVec(multilevel->Prolongs[i],tmp_vector1,tmp_vector2);
      for(j=0;j<multilevel->Stiff_Matrix[i+1]->N_Rows;++j)
         tmp_vector1[j] = tmp_vector2[j];
   }
  MatrixDotVec(multilevel->Prolongs[num_h-1],tmp_vector1,tmp_vector2);
  free(tmp_vector1);
  return tmp_vector2;
 }


 void ProlongH2h2(MULTILEVEL *multilevel,int num_H,int num_h,double *vector,double *tmp_vector2)
 {
  int i,j;
  double *tmp_vector1;
  tmp_vector1 = malloc(sizeof(double)*multilevel->Stiff_Matrix[num_h]->N_Rows);
  for(i=0;i<multilevel->Stiff_Matrix[num_H]->N_Rows;++i)
    tmp_vector1[i] = vector[i];
  for(i=num_H;i<num_h-1;++i)
   {
    MatrixDotVec(multilevel->Prolongs[i],tmp_vector1,tmp_vector2);
      for(j=0;j<multilevel->Stiff_Matrix[i+1]->N_Rows;++j)
         tmp_vector1[j] = tmp_vector2[j];
   }
  MatrixDotVec(multilevel->Prolongs[num_h-1],tmp_vector1,tmp_vector2);

  free(tmp_vector1);
 }















 double* Restricth2H(MULTILEVEL *multilevel,int num_h,int num_H,double *vector)
 {
  int i,j;
  double *tmp_vector1,*tmp_vector2;
  tmp_vector1 = malloc(sizeof(double)*multilevel->Stiff_Matrix[num_h]->N_Rows);
  tmp_vector2 = malloc(sizeof(double)*multilevel->Stiff_Matrix[num_h]->N_Rows);
  for(i=0;i<multilevel->Stiff_Matrix[num_h]->N_Rows;++i)
    tmp_vector1[i] = vector[i];
  for(i=num_h-1;i>num_H;--i)
   {
    MatrixDotVec(multilevel->Restricts[i],tmp_vector1,tmp_vector2);
      for(j=0;j<multilevel->Stiff_Matrix[i]->N_Rows;++j)
         tmp_vector1[j] = tmp_vector2[j];
   }
  MatrixDotVec(multilevel->Restricts[num_H],tmp_vector1,tmp_vector2);
  free(tmp_vector1);
  
  return tmp_vector2;
 }
 void Restricth2H2(MULTILEVEL *multilevel,int num_h,int num_H,double *vector,double *vector_H)
 {
  int i,j;
  double *tmp_vector1,*tmp_vector2;
  tmp_vector1 = malloc(sizeof(double)*multilevel->Stiff_Matrix[num_h]->N_Rows);
  tmp_vector2 = malloc(sizeof(double)*multilevel->Stiff_Matrix[num_h]->N_Rows);
  for(i=0;i<multilevel->Stiff_Matrix[num_h]->N_Rows;++i)
    tmp_vector1[i] = vector[i];
  for(i=num_h-1;i>num_H;--i)
   {
    MatrixDotVec(multilevel->Restricts[i],tmp_vector1,tmp_vector2);
      for(j=0;j<multilevel->Stiff_Matrix[i]->N_Rows;++j)
         tmp_vector1[j] = tmp_vector2[j];
   }
  MatrixDotVec(multilevel->Restricts[num_H],tmp_vector1,tmp_vector2);

  for(i=0;i<multilevel->Stiff_Matrix[num_H]->N_Rows;++i)
    vector_H[i] = tmp_vector2[i]; 

  free(tmp_vector1);
  free(tmp_vector2);

 }































 MATRIX* BuildMatrixMul(DISCRETEFORM3D *discreteform)
 {
    Fespace3D *fesp1, *fesp2;
    fesp1 = discreteform->Anasatz_Space;
    fesp2 = discreteform->Test_Space; 
    MATRIX *matrix;
    matrix = malloc(1*sizeof(MATRIX));
    matrix->DiscreteForm3D = discreteform;
    
    matrix->N_Rows = fesp2->DOF_ALL+1;
    matrix->N_Columns = fesp1->DOF_ALL+1;
    // then build the matrix structure    
    //here we assume all the FEM spaces are defined on the same mesh
    MESH *mesh;  
    mesh = fesp2->Mesh;
    int num_volus,num_faces,num_verts,num_lines;
    int curr_volu, curr_face, curr_line, i, j, num_vert_face, num_line_face;
    BASEFUNCTION3D *anasatz_base,*test_base;

    num_volus = mesh->Num_Volus_Global;
    num_faces = mesh->Num_Faces_Global;
    num_lines = mesh->Num_Lines_Global;
    num_verts = mesh->Num_Verts_Global;
   
    anasatz_base = fesp1->Base;
    test_base = fesp2->Base;
    int *anasatz_DOF,*test_DOF;
    anasatz_DOF = anasatz_base->DOF;
    test_DOF = test_base->DOF;
    int *N_volu4vert,*N_volu4line,*N_volu4face,*N_face4line,*N_face4vert, *N_line4vert; 
    /** need to be free */
    N_volu4vert = malloc(num_verts*sizeof(int));
    N_volu4line = malloc(num_lines*sizeof(int));
    N_volu4face = malloc(num_faces*sizeof(int));    
    N_face4vert = malloc(num_verts*sizeof(int));
    N_face4line = malloc(num_lines*sizeof(int));
    N_line4vert = malloc(num_verts*sizeof(int));

    memset(N_volu4vert,0,num_verts*sizeof(int));
    memset(N_volu4line,0,num_lines*sizeof(int));
    memset(N_volu4face,0,num_faces*sizeof(int));
    memset(N_face4line,0,num_lines*sizeof(int));
    memset(N_face4vert,0,num_verts*sizeof(int));
    memset(N_line4vert,0,num_verts*sizeof(int));

    VOLU *volu;
    FACEE *face;
    LINE *line;
    int N_local_vert,N_local_line,N_local_face;
    for(curr_volu=0;curr_volu<num_volus;curr_volu++)
    {
        //first give the number of related dof for each vert
        volu = mesh->Volus[curr_volu];
        N_local_vert = 4;
        N_local_line = 6;
	N_local_face = 4;
        for(i=0;i<N_local_vert;i++)
            N_volu4vert[volu->Verts[i]]++;
        for(i=0;i<N_local_line;i++)
            N_volu4line[volu->Lines[i]]++;
        for(i=0;i<N_local_face;i++)
            N_volu4face[volu->Faces[i]]++;
    }
     for(curr_face=0;curr_face<num_faces;curr_face++)
    {
        face = mesh->Faces[curr_face];
        N_local_vert = 3;
        N_local_line = 3;
        for(i=0;i<N_local_vert;i++)
            N_face4vert[face->Verts[i]]++;
        for(i=0;i<N_local_line;i++)
            N_face4line[face->Lines[i]]++;
    }
    for(curr_line=0;curr_line<num_lines;curr_line++)
    {
        line = mesh->Lines[curr_line];
        N_line4vert[line->Verts[0]]++;
        N_line4vert[line->Verts[1]]++;
     }
     int curr_vert;
     int *RowPtr;
     int N_Rows;
     N_Rows = fesp2->DOF_ALL;
     RowPtr = malloc((N_Rows+2)*sizeof(int));
     memset(RowPtr,-1,(N_Rows+2)*sizeof(int));
     RowPtr[0] = 0;
     int curr_row;
     curr_row = 0;
     int vert_dof, line_dof, face_dof,volu_dof;
     int local_volu,local_face,local_line;
     int volu_verts = 4;
     int volu_lines = 6;
     int volu_faces = 4;
// for the vert
     for(i=0;i<test_DOF[0];i++)
     {
       for(curr_vert=0;curr_vert<num_verts;curr_vert++)
       {
	 curr_row ++;
	 vert_dof = 0;
	 // the vert itself
	 vert_dof += anasatz_DOF[0];
	 
	 local_volu = N_volu4vert[curr_vert];
	 vert_dof += local_volu*(anasatz_DOF[3]
		    +(volu_faces-3)*anasatz_DOF[2]);
      	 
	 local_face = N_face4vert[curr_vert];
	 vert_dof += local_face*(anasatz_DOF[2]+anasatz_DOF[1]);
        
	 local_line = N_line4vert[curr_vert];
	 vert_dof += local_line*(anasatz_DOF[0]+anasatz_DOF[1]);
	 RowPtr[curr_row] = RowPtr[curr_row-1] + vert_dof+1;
       }// end for curr_vert
     }//end for i
     //for the lines
     for(i=0;i<test_DOF[1];i++)
     {
       for(curr_line=0;curr_line<num_lines;curr_line++)
       {
	 curr_row++;
	 line_dof = anasatz_DOF[0]*2+anasatz_DOF[1]
	            +N_volu4line[curr_line]*(anasatz_DOF[3]+anasatz_DOF[1]+anasatz_DOF[2]*2)
	            +N_face4line[curr_line]*(anasatz_DOF[2]+anasatz_DOF[1]*2+anasatz_DOF[0]);
	 RowPtr[curr_row] = RowPtr[curr_row-1] + line_dof+1;
       }//end for curr_line       
     }// end for i    

     //for the face
     for(i=0;i<test_DOF[2];i++)
     {
       for(curr_face=0;curr_face<num_faces;curr_face++)
       {
	 curr_row++;
	 face_dof = anasatz_DOF[0]*3+anasatz_DOF[1]*3+anasatz_DOF[2]
	       +N_volu4face[curr_face]*(anasatz_DOF[3]+anasatz_DOF[2]*3+anasatz_DOF[1]*3+anasatz_DOF[0]);
	 RowPtr[curr_row] = RowPtr[curr_row-1] + face_dof+1;
       }
     }
     //for the volu
     for(i=0;i<test_DOF[3];i++)
     {
       for(curr_volu=0;curr_volu<num_volus;curr_volu++)
       {
	 curr_row++;
	 volu_dof = anasatz_DOF[3]+anasatz_DOF[2]*4+anasatz_DOF[1]*6+anasatz_DOF[0]*4;
	 RowPtr[curr_row] = RowPtr[curr_row-1] + volu_dof+1;
       }//end for curr_volu       
     }// end for i   
     curr_row++;
     RowPtr[curr_row] = RowPtr[curr_row-1] + fesp2->DOF_ALL+1;
    int N_Entries;
    N_Entries = RowPtr[curr_row];
    int *KCol;
    KCol = malloc((N_Entries)*sizeof(int));
    memset(KCol,-1,(N_Entries)*sizeof(int));
    double *Entries;
    Entries = malloc(N_Entries*sizeof(double));
    memset(Entries,0.0,N_Entries*sizeof(double));
    int *pos_row;
    /** need to free */
    pos_row = malloc(N_Rows*sizeof(int));
    memset(pos_row,0,N_Rows*sizeof(int));
    int curr_pos_test,curr_pos_anasatz;
    int n_dof_test_volu,n_dof_anasatz_volu;
    // for the face iteration
    int *BeginIndex_test, *GlobalNumber_test;
    int *BeginIndex_anasatz,*GlobalNumber_anasatz;
    BeginIndex_anasatz = fesp1->BeginIndex;
    GlobalNumber_anasatz = fesp1->GlobalNumbers;
    BeginIndex_test = fesp2->BeginIndex;
    GlobalNumber_test = fesp2->GlobalNumbers;
    int curr_index_test,curr_index_anasatz,col_num, curr_entry;
    int k;
    for(curr_volu=0;curr_volu<num_volus;curr_volu++)
    {
        //得到当前单元在GlobalNumber_test中的起始位置
        curr_index_test = BeginIndex_test[curr_volu];
        //当前单元上自由度的个数
        n_dof_test_volu = BeginIndex_test[curr_volu+1]-BeginIndex_test[curr_volu];
        curr_index_anasatz = BeginIndex_anasatz[curr_volu];
        n_dof_anasatz_volu = BeginIndex_anasatz[curr_volu+1]-BeginIndex_anasatz[curr_volu];
        //对单元上自由度进行循环
        for(i=0;i<n_dof_test_volu;i++)
        {
	    //当前单元第i个自由度的全局编号
            curr_pos_test = GlobalNumber_test[curr_index_test+i];
            //KCOL中放置当前单元的第i个自由度的起始位置
	    curr_entry = RowPtr[curr_pos_test];
            //对测试函数空间自由度个数循环
	    for(j=0;j<n_dof_anasatz_volu;j++)
            {   
	        //测试函数空间第i个自由度的全局编号
                col_num = GlobalNumber_anasatz[curr_index_anasatz+j];
                //下面要把测试空间中的这个自由度放进列编号里
		//首先对已经放进去的列编号进行循环，如果发现和当前自由度编号一致的，则循环结束
		for(k=0;k<pos_row[curr_pos_test];k++)
                {
                    if(KCol[curr_entry+k]==col_num)
                        break;
                }
		//如果循环结束不是由break导致，即已有列编号不包括要放进去的编号，则将该编号放入
                if(k==pos_row[curr_pos_test])
                {
                    KCol[curr_entry+pos_row[curr_pos_test]] = col_num;
                    pos_row[curr_pos_test]++;
                }        
             }// j循环结束
        }//单元上自由度循环结束
    }//单元循环结束
    for(i=0;i<fesp2->DOF_ALL+1;i++)
    {
      curr_entry=RowPtr[i+1]-1;
      KCol[curr_entry]= fesp2->DOF_ALL;
    }
    for(i=0;i<fesp2->DOF_ALL;i++)
    {
      curr_entry=RowPtr[fesp2->DOF_ALL];
      KCol[curr_entry+i]= i;
    }
    free(N_volu4line);
    free(N_volu4vert); 
    free(N_line4vert);
    free(N_volu4face);
    free(N_face4line);
    free(N_face4vert); 
    
    free(pos_row);
    matrix->RowPtr = RowPtr;
    matrix->KCol = KCol;
    matrix->Entries = Entries;
    matrix->N_Entries = N_Entries;
    for(i=0;i<N_Rows;i++)
    {
       QuickSort_Int(KCol,RowPtr[i],RowPtr[i+1]-1);
    }
  return matrix;    
}

/**====================================================================*/
// Assemble the matrix 
/**合成刚度矩阵 */
void AssembleMatrixMul(MATRIX *Matrix,Fefunction3D *Fefunct)//matrix是要合成的矩阵，matrixA是细网格刚度矩阵
{
   printf ( "beign assemble multilevel\n" );
  // define some notion
  DISCRETEFORM3D *discreteform;
  Fespace3D *test_space, *anasatz_space,*fesp_finer;
  discreteform = Matrix->DiscreteForm3D;
  test_space = discreteform->Test_Space;
  anasatz_space = discreteform->Anasatz_Space;
  fesp_finer = Fefunct->Fespace;
  //get the mesh: all the fem spaces are defined on the same mesh
  MESH *mesh,*mesh_finer;
  mesh = test_space->Mesh;
  mesh_finer=fesp_finer->Mesh;
  int num_volus, i, j,k;
  num_volus = mesh_finer->Num_Volus_Global;
  Quadrature3D *Quad3D;
  Quad3D = discreteform->Quad3D;
  Quadrature2D *Quad2D;
  Quad2D = discreteform->Quad2D;
  Quadrature1D *Quad1D;
  Quad1D = discreteform->Quad1D;
  ELEMENT3D *elem3D;
  elem3D = malloc(1*sizeof(ELEMENT3D));
  // initializing the element
  InitialElem3D(mesh,elem3D);
  int n_bas_test, n_bas_anasatz,n_bas_finer, N_Mat_Entries, dim_anasatz_value, dim_test_value;
  double *Mat_Entries,*tmp_entries;
  n_bas_test = test_space->Base->Num_Bas;
  dim_test_value = test_space->Base->Value_Dim;
  n_bas_anasatz = anasatz_space->Base->Num_Bas;
  dim_anasatz_value = anasatz_space->Base->Value_Dim;
  n_bas_finer=fesp_finer->Base->Num_Bas;

  N_Mat_Entries = n_bas_test * n_bas_anasatz;
  Mat_Entries = malloc(N_Mat_Entries*sizeof(double));
  tmp_entries = malloc(N_Mat_Entries*sizeof(double));
    
  int *Test_GlobalNumbers, *Anasatz_GlobalNumbers, 
      *Test_BeginIndex, *Anasatz_BeginIndex;
  Test_GlobalNumbers = test_space->GlobalNumbers;
  Test_BeginIndex = test_space->BeginIndex;
  Anasatz_GlobalNumbers = anasatz_space->GlobalNumbers;
  Anasatz_BeginIndex = anasatz_space->BeginIndex;
  
  int n_multiindex_test, n_multiindex_anasatz;
  n_multiindex_test = discreteform->N_Test_MultiIndex;
  n_multiindex_anasatz = discreteform->N_Anasatz_MultiIndex;
  
  double *Test_Values, *Anasatz_Values;
  double *Test_MultiIndex_Values, *Anasatz_MultiIndex_Values;
  
  Test_Values = malloc(n_bas_test*dim_test_value*n_multiindex_test*sizeof(double));
  Anasatz_Values = malloc(n_bas_anasatz*dim_anasatz_value*n_multiindex_anasatz*sizeof(double));
 
  Test_MultiIndex_Values = malloc(n_multiindex_test*dim_test_value*sizeof(double));
  Anasatz_MultiIndex_Values = malloc(n_multiindex_anasatz*dim_anasatz_value*sizeof(double));
  double Quad_xi, Quad_eta,Quad_zeta, Quad_W;
  int curr_volu, kk,ii, Quad_Points;//????????????????????;
  int Quad_Points3D,Quad_Points2D, Quad_Points1D;
  Quad_Points3D = Quad3D->N_Points;
  if(Quad1D)
    Quad_Points1D = Quad1D->N_Points;
  if(Quad2D)
    Quad_Points2D = Quad2D->N_Points;
  int n_test_multiindex,n_anasatz_multiindex;
  MultiIndex3D *test_MultiIndex, *anasatz_MultiIndex;
  //下面两个定义是否重复了？ 没有重复 ？？？？？？？？？？？？？？？？
  n_test_multiindex = discreteform->N_Test_MultiIndex;
  n_anasatz_multiindex = discreteform->N_Anasatz_MultiIndex;
  
  test_MultiIndex = discreteform->Test_MultiIndex;
  anasatz_MultiIndex = discreteform->Anasatz_MultiIndex;
  
  DiscreteFormMatrix *discreteformvolu;
  discreteformvolu = discreteform->DiscreteFormVolu;
  
  DiscreteFormMatrixFace *discreteformface; // *discreteformbd;
  discreteformface = discreteform->DiscreteFormFace;
  
  BASEFUNCTION3D *test_Base, *anasatz_Base,*base_finer;
  test_Base = discreteform->Test_Space->Base;
  anasatz_Base = discreteform->Anasatz_Space->Base;
  base_finer =   fesp_finer->Base;

  int N_AuxFefunct, *AuxFefunPtr, *N_AuxFefun_MultiIndex, N_AuxFefun_Values;
  double *AuxFefun_Values;
  MultiIndex3D *AuxFefun_MultiIndex;
  Fefunction3D **AuxFefunct,*Fefunction,**AuxFefunct_finer;
  
  N_AuxFefunct = discreteform->N_AuxFeFun;
  AuxFefunct = discreteform->AuxFeFun;
  AuxFefunct_finer=malloc(sizeof(Fefunction3D*)*N_AuxFefunct);
  double *value,*value1,*value2,*value3;
  value1=malloc(sizeof(double)*test_space->DOF_ALL);
  value2=malloc(sizeof(double)*fesp_finer->DOF_ALL);
  value3=malloc(sizeof(double)*fesp_finer->DOF_ALL);
  //memset(value3,0.0,fesp_finer->DOF_ALL*sizeof(double));
  MATRIX *matrix_Prolong;
  matrix_Prolong=malloc(sizeof(MATRIX));
  BuildProlongAn(test_space,fesp_finer,matrix_Prolong);
  if(N_AuxFefunct)
  {
    for(i=0;i<N_AuxFefunct;++i)
     {
      value=AuxFefunct[i]->Values;
      MatrixDotVec(matrix_Prolong,value,value3);
      for(j=0;j<fesp_finer->DOF_ALL;++j)
        value3[j]+=value[test_space->DOF_ALL]*Fefunct->Values[j];
      AuxFefunct_finer[i] = BuildFefunction3D(fesp_finer,value3);
     }
    //Normalization(AuxFefunct_finer[0]);
  }
  N_AuxFefun_MultiIndex = discreteform->N_AuxFeFun_MultiIndex;
  AuxFefun_MultiIndex = discreteform->AuxFeFun_MultiIndex;
  N_AuxFefun_Values = discreteform->N_AuxFeFun_Values;
  AuxFefun_Values = discreteform->AuxFeFun_Values;
  //printf ( "!!!!!!!!!55555555555\n" );
  //提出用户函数
  Functions3D **UserFunction;
  int N_UserFunction;
  UserFunction=discreteform->UserFunction;
  N_UserFunction=discreteform->N_UserFunction;

  int N_FeConst;
  N_FeConst=discreteform->N_FeConst;
  double *FeConst;
  FeConst=discreteform->FeConst;

  VOLU **Volus, *volu;
  FACEE **Faces, *face;
  LINE **Lines, *line;
  Volus = mesh->Volus;
  Faces = mesh->Faces;
  Lines = mesh->Lines;  
  int elem_lines, ID_bd, Ancestor;
  
  BoundCondFunct3D *BoundCond;
  BoundCond = anasatz_space->BoundCond;
  BoundType bdtype; 
  
  double *Local_Prolong;
  Local_Prolong = malloc(n_bas_test*n_bas_test*sizeof(double));
  
  ELEMENT3D *coarse_elem3D,*finer_elem3D,*elem;
  coarse_elem3D = malloc(sizeof(ELEMENT3D));
  InitialElem3D(mesh,coarse_elem3D);
  finer_elem3D = malloc(sizeof(ELEMENT3D));
  InitialElem3D(mesh_finer,finer_elem3D);
  elem = malloc(sizeof(ELEMENT3D));
  InitialElem3D(mesh_finer,elem3D);
  /**进行单元循环，得到刚度矩阵 */
  //printf ( "kai shi xi wang ge xun huan\n" );
  for(curr_volu=0;curr_volu<num_volus;curr_volu++)
  { 
    volu = mesh_finer->Volus[curr_volu];
    Ancestor=volu->Ancestor;
    //取出当前的有限元
    GetElement3D(mesh_finer,curr_volu,elem3D);
    memset(Mat_Entries,0.0,N_Mat_Entries*sizeof(double));
    // iteration for the quadrature points
     if(discreteformvolu!=NULL)
     {
       for(i=0;i< Quad_Points3D;i++)
       {
	//get the quadrature information 
	Quad_xi   =  Quad3D->Quad_X[i];
	Quad_eta  =  Quad3D->Quad_Y[i];
	Quad_zeta =  Quad3D->Quad_Z[i];
	Quad_W    =  Quad3D->Quad_W[i];
	SubAssembleMatrix3D(n_bas_test,base_finer,n_bas_anasatz,base_finer,
			    n_test_multiindex,test_MultiIndex,n_anasatz_multiindex,anasatz_MultiIndex, 
		            N_AuxFefunct,AuxFefunct_finer,N_AuxFefun_MultiIndex,AuxFefun_MultiIndex,
		            N_AuxFefun_Values, AuxFefun_Values,N_UserFunction,UserFunction, 
			    N_FeConst,FeConst,discreteformvolu,elem3D,
			    Test_Values,Anasatz_Values,Test_MultiIndex_Values, Anasatz_MultiIndex_Values,
			    Quad_xi,Quad_eta,Quad_zeta,Quad_W,Mat_Entries);
       } //end for i (Quad_Points)       
     } //end for if

    GetElement3D(mesh,Ancestor,coarse_elem3D);
    GetElement3D(mesh_finer,curr_volu,finer_elem3D);
    GetLocalProlong(coarse_elem3D,test_Base,n_bas_test,finer_elem3D,base_finer,n_bas_finer,
	          elem3D,Local_Prolong);
    MatMatMat(Local_Prolong,Mat_Entries,tmp_entries); 
    for(i=0;i<N_Mat_Entries;++i)
    Mat_Entries[i]=tmp_entries[i];
    AddSubMatrix(Matrix,Mat_Entries,n_bas_test,Test_GlobalNumbers+Test_BeginIndex[Ancestor],
		 n_bas_anasatz,Anasatz_GlobalNumbers+Anasatz_BeginIndex[Ancestor]);	
  }//end for curr_volu
/* ------------------------------------------------------------------------------------------------- */
 
  //Fespace3D *fespfine;
  //fespfine=Fefunct->Fespace;
  MATRIX *MatrixR,*MatrixA;
  MatrixR=malloc(sizeof(MATRIX));
  TransposeMatrix(matrix_Prolong,MatrixR);
  //BuildRestrictAn(test_space,fesp_finer,MatrixR);
  DISCRETEFORM3D *discreteform_finer;
  /*
  discreteform_finer = BuildDiscreteFormAll3D(fesp_finer,fesp_finer,n_test_multiindex,test_MultiIndex,
	                         n_anasatz_multiindex, anasatz_MultiIndex,
				 N_AuxFefunct,AuxFefunct_finer,N_AuxFefun_MultiIndex,AuxFefun_MultiIndex,
				 N_UserFunction,UserFunction, discreteformvolu,
				 Quad_Points3D, Quad_Points2D,Quad_Points1D);
  */
  discreteform_finer = BuildDiscreteFormFeConst3D(fesp_finer,fesp_finer,n_test_multiindex,test_MultiIndex,
                                 n_anasatz_multiindex, anasatz_MultiIndex,
                                 N_AuxFefunct,AuxFefunct_finer,N_AuxFefun_MultiIndex,AuxFefun_MultiIndex,
                                 N_UserFunction,UserFunction, discreteformvolu,
                                 N_FeConst,FeConst, Quad_Points3D, Quad_Points2D,Quad_Points1D);
  MatrixA = BuildMatrix(discreteform_finer);
  AssembleMatrix(MatrixA);
  
  double *bb,cc[0];
  bb=malloc(sizeof(double)*test_space->DOF_ALL);
  //printf ( "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n" );
  MatMatVec(MatrixR,MatrixA,Fefunct->Values,bb);
  //printf ( "@@@@@@@@@@@@@@@@@@@@@@@@@@@11111111111111\n" );
  //for(i=0;i<Fefunct->Fespace->DOF_ALL;++i)
  //   Fefunct->Values[i]=2;

  //cc[0]=VecMatVec(Fefunct->Values,MatrixA);//为什么这样不行？？？？？？？？？？？？？？
  //printf ( "cc=%f\n",cc );
  VecMatVec(Fefunct->Values,MatrixA,cc);
  //printf ( "cc=%f\n",cc[0] );

  //getchar();
     //printf ( "%d==%f\n",i,Fefunct->Values[i] );
  //getchar();
  //printf ( "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@22222\n" );
  for(i=0;i<test_space->DOF_ALL;++i)
  {
   Matrix->Entries[Matrix->RowPtr[i+1]-1]=bb[i];
  }
  for(i=0;i<test_space->DOF_ALL;++i)
  {
   Matrix->Entries[Matrix->RowPtr[test_space->DOF_ALL]+i]=bb[i];
  }
   Matrix->Entries[Matrix->RowPtr[test_space->DOF_ALL+1]-1]=cc[0];
   /* ------------------------------------------------------------------------------------------------ */

    /** 释放内存 */
    FreeMatrix(MatrixR);
    FreeMatrix(MatrixA); 
    FreeElem3D(elem3D);
    free(Mat_Entries);
    free(Test_Values);
    free(Anasatz_Values);
    free(Test_MultiIndex_Values);
    free(Anasatz_MultiIndex_Values);
    free(bb);
    FreeDiscreteForm3D(discreteform_finer); 
    //free(AuxFefun_Values);
    //free(normal);
    //free(tangent);
    //free(AuxFunct_Values);
}//end for curr_face  


/**合成刚度矩阵*/
void AssembleMatrixMulSimplify(MATRIX *Matrix,Fefunction3D *Fefunct)//matrix是要合成的矩阵，matrixA是细网格刚度矩阵
{                                                        
  printf ( "beign assemble multilevel\n" );
  int i,j;
  MATRIX *MatrixP;
  MatrixP = malloc(sizeof(MATRIX));
  Fespace3D *test_space, *anasatz_space, *fesp_finer;
  DISCRETEFORM3D *discreteform;
  discreteform = Matrix->DiscreteForm3D;
  test_space = discreteform->Test_Space;
  anasatz_space = discreteform->Anasatz_Space;
  fesp_finer = Fefunct->Fespace;
  BuildProlongAn(test_space, fesp_finer, MatrixP);
  int N_AuxFefunct, *N_AuxFefun_MultiIndex;
  Fefunction3D **AuxFefunct;
  MultiIndex3D *AuxFefun_MultiIndex;
  N_AuxFefun_MultiIndex = discreteform->N_AuxFeFun_MultiIndex;
  AuxFefun_MultiIndex = discreteform->AuxFeFun_MultiIndex;
  N_AuxFefunct = discreteform->N_AuxFeFun;
  AuxFefunct = discreteform->AuxFeFun;
 
  Fefunction3D **AuxFefunct_finer;
  AuxFefunct_finer=malloc(sizeof(Fefunction3D*)*N_AuxFefunct);
 
  double *value,**value3;
  value3=malloc(sizeof(double*)*N_AuxFefunct);
  for(i=0;i<N_AuxFefunct;++i)
    value3[i] = malloc(sizeof(double)*fesp_finer->DOF_ALL);
  if(N_AuxFefunct)
  {
    for(i=0;i<N_AuxFefunct;++i)
     {
      value = AuxFefunct[i]->Values;
      MatrixDotVec(MatrixP,value,value3[i]);
      for(j=0;j<fesp_finer->DOF_ALL;++j)
        value3[i][j]+=value[test_space->DOF_ALL]*Fefunct->Values[j];
      AuxFefunct_finer[i] = BuildFefunction3D(fesp_finer,value3[i]);
     }
  for(j=0;j<fesp_finer->DOF_ALL;++j)
    value3[1][j] = AuxFefunct[1]->Values[j];

 }
   
 
  MATRIX *MatrixR,*MatrixA;
  MatrixR=malloc(sizeof(MATRIX));
  TransposeMatrix(MatrixP,MatrixR);
  DISCRETEFORM3D *discreteform_finer;
  int n_test_multiindex,n_anasatz_multiindex;
  MultiIndex3D *test_MultiIndex, *anasatz_MultiIndex;
  n_test_multiindex = discreteform->N_Test_MultiIndex;
  n_anasatz_multiindex = discreteform->N_Anasatz_MultiIndex;
  test_MultiIndex = discreteform->Test_MultiIndex;
  anasatz_MultiIndex = discreteform->Anasatz_MultiIndex;
  
  DiscreteFormMatrix *discreteformvolu;
  discreteformvolu = discreteform->DiscreteFormVolu;
  
  Functions3D **UserFunction;
  int N_UserFunction, N_FeConst;
  UserFunction=discreteform->UserFunction;
  N_UserFunction=discreteform->N_UserFunction;
  N_FeConst=discreteform->N_FeConst;
  double *FeConst;
  FeConst=discreteform->FeConst;
  discreteform_finer = BuildDiscreteFormFeConst3D(fesp_finer,fesp_finer,n_test_multiindex,test_MultiIndex,
                                 n_anasatz_multiindex, anasatz_MultiIndex,
                                 N_AuxFefunct,AuxFefunct_finer,N_AuxFefun_MultiIndex,AuxFefun_MultiIndex,
                                 N_UserFunction,UserFunction, discreteformvolu,
                                 N_FeConst,FeConst,4,4,4);
  MatrixA = BuildMatrix(discreteform_finer);
  AssembleMatrix(MatrixA);
  
  double *bb,cc[0];
  bb = malloc(sizeof(double)*test_space->DOF_ALL);
  for(i=0;i<test_space->DOF_ALL;++i)
  {
    TripleMatrix(MatrixR,MatrixA,MatrixP,i,bb);
    for(j=0;j<Matrix->RowPtr[i+1]-Matrix->RowPtr[i]-1;++j)
      Matrix->Entries[Matrix->RowPtr[i]+j] = bb[Matrix->KCol[Matrix->RowPtr[i]+j]]; 
  }
  MatMatVec(MatrixR,MatrixA,Fefunct->Values,bb);
  VecMatVec(Fefunct->Values,MatrixA,cc);
  for(i=0;i<test_space->DOF_ALL;++i)
  {
   Matrix->Entries[Matrix->RowPtr[i+1]-1]=bb[i];
  }
  for(i=0;i<test_space->DOF_ALL;++i)
  {
   Matrix->Entries[Matrix->RowPtr[test_space->DOF_ALL]+i]=bb[i];
  }
   Matrix->Entries[Matrix->RowPtr[test_space->DOF_ALL+1]-1]=cc[0];

  /** 释放内存 */
  FreeMatrix(MatrixR);
  FreeMatrix(MatrixP);
  FreeMatrix(MatrixA); 
  free(bb);
  FreeDiscreteForm3D(discreteform_finer); 
}//end of the function 











void AssembleMatrixHinh(MATRIX *Matrix,Fefunction3D *Fefunct)//当利用到细网格单元上的信息时，在细网格上操作
{
   printf ( "beign assemble multilevel\n" );
  // define some notion
  DISCRETEFORM3D *discreteform;
  Fespace3D *test_space, *anasatz_space,*fesp_finer;
  discreteform = Matrix->DiscreteForm3D;
  test_space = discreteform->Test_Space;
  anasatz_space = discreteform->Anasatz_Space;
  fesp_finer = Fefunct->Fespace;
  //get the mesh: all the fem spaces are defined on the same mesh
  MESH *mesh,*mesh_finer;
  mesh = test_space->Mesh;
  mesh_finer=fesp_finer->Mesh;
  int num_volus, i, j,k;
  num_volus = mesh_finer->Num_Volus_Global;
  Quadrature3D *Quad3D;
  Quad3D = discreteform->Quad3D;
  Quadrature2D *Quad2D;
  Quad2D = discreteform->Quad2D;
  Quadrature1D *Quad1D;
  Quad1D = discreteform->Quad1D;
  ELEMENT3D *elem3D;
  elem3D = malloc(1*sizeof(ELEMENT3D));
  // initializing the element
  InitialElem3D(mesh,elem3D);
  int n_bas_test, n_bas_anasatz,n_bas_finer, N_Mat_Entries, dim_anasatz_value, dim_test_value;
  double *Mat_Entries,*tmp_entries;
  n_bas_test = test_space->Base->Num_Bas;
  dim_test_value = test_space->Base->Value_Dim;
  n_bas_anasatz = anasatz_space->Base->Num_Bas;
  dim_anasatz_value = anasatz_space->Base->Value_Dim;
  n_bas_finer=fesp_finer->Base->Num_Bas;

  N_Mat_Entries = n_bas_test * n_bas_anasatz;
  Mat_Entries = malloc(N_Mat_Entries*sizeof(double));
  tmp_entries = malloc(N_Mat_Entries*sizeof(double));
    
  int *Test_GlobalNumbers, *Anasatz_GlobalNumbers, 
      *Test_BeginIndex, *Anasatz_BeginIndex;
  Test_GlobalNumbers = test_space->GlobalNumbers;
  Test_BeginIndex = test_space->BeginIndex;
  Anasatz_GlobalNumbers = anasatz_space->GlobalNumbers;
  Anasatz_BeginIndex = anasatz_space->BeginIndex;
  
  int n_multiindex_test, n_multiindex_anasatz;
  n_multiindex_test = discreteform->N_Test_MultiIndex;
  n_multiindex_anasatz = discreteform->N_Anasatz_MultiIndex;
  
  double *Test_Values, *Anasatz_Values;
  double *Test_MultiIndex_Values, *Anasatz_MultiIndex_Values;
  
  Test_Values = malloc(n_bas_test*dim_test_value*n_multiindex_test*sizeof(double));
  Anasatz_Values = malloc(n_bas_anasatz*dim_anasatz_value*n_multiindex_anasatz*sizeof(double));
 
  Test_MultiIndex_Values = malloc(n_multiindex_test*dim_test_value*sizeof(double));
  Anasatz_MultiIndex_Values = malloc(n_multiindex_anasatz*dim_anasatz_value*sizeof(double));
  double Quad_xi, Quad_eta,Quad_zeta, Quad_W;
  int curr_volu, kk,ii, Quad_Points;//????????????????????;
  int Quad_Points3D,Quad_Points2D, Quad_Points1D;
  Quad_Points3D = Quad3D->N_Points;
  if(Quad1D)
    Quad_Points1D = Quad1D->N_Points;
  if(Quad2D)
    Quad_Points2D = Quad2D->N_Points;
  int n_test_multiindex,n_anasatz_multiindex;
  MultiIndex3D *test_MultiIndex, *anasatz_MultiIndex;
  //下面两个定义是否重复了？ 没有重复 ？？？？？？？？？？？？？？？？
  n_test_multiindex = discreteform->N_Test_MultiIndex;
  n_anasatz_multiindex = discreteform->N_Anasatz_MultiIndex;
  
  test_MultiIndex = discreteform->Test_MultiIndex;
  anasatz_MultiIndex = discreteform->Anasatz_MultiIndex;
  
  DiscreteFormMatrix *discreteformvolu;
  discreteformvolu = discreteform->DiscreteFormVolu;
  
  DiscreteFormMatrixFace *discreteformface; // *discreteformbd;
  discreteformface = discreteform->DiscreteFormFace;
  
  BASEFUNCTION3D *test_Base, *anasatz_Base,*base_finer;
  test_Base = discreteform->Test_Space->Base;
  anasatz_Base = discreteform->Anasatz_Space->Base;
  base_finer =   fesp_finer->Base;

  int N_AuxFefunct, *AuxFefunPtr, *N_AuxFefun_MultiIndex, N_AuxFefun_Values;
  double *AuxFefun_Values;
  MultiIndex3D *AuxFefun_MultiIndex;
  Fefunction3D **AuxFefunct,*Fefunction,**AuxFefunct_finer;
  
  N_AuxFefunct = discreteform->N_AuxFeFun;
  AuxFefunct = discreteform->AuxFeFun;
  AuxFefunct_finer=malloc(sizeof(Fefunction3D*)*N_AuxFefunct);
  double *value,*value1,*value2,*value3;
  value1=malloc(sizeof(double)*test_space->DOF_ALL);
  value2=malloc(sizeof(double)*fesp_finer->DOF_ALL);
  value3=malloc(sizeof(double)*fesp_finer->DOF_ALL);
  memset(value3,0.0,fesp_finer->DOF_ALL*sizeof(double));
  MATRIX *matrix_Prolong;
  matrix_Prolong=malloc(sizeof(MATRIX));
  BuildProlongAn(test_space,fesp_finer,matrix_Prolong);
  if(N_AuxFefunct)
  {
    for(i=0;i<N_AuxFefunct;++i)
     {
      value=AuxFefunct[i]->Values;
      for(k=0;k<test_space->DOF_ALL;++k)
      {
       for(j=0;j<test_space->DOF_ALL;++j)
          value1[j]=0.0;
       value1[k]=value[k];
       MatrixDotVec(matrix_Prolong,value1,value2);
       for(j=0;j<fesp_finer->DOF_ALL;++j)
         value3[j]+=value2[j];
      }       
      for(j=0;j<fesp_finer->DOF_ALL;++j)
        value3[j]+=value[test_space->DOF_ALL]*Fefunct->Values[j];
      AuxFefunct_finer[i] = BuildFefunction3D(fesp_finer,value3);
     }
    Normalization(AuxFefunct_finer[0]);
  }
  N_AuxFefun_MultiIndex = discreteform->N_AuxFeFun_MultiIndex;
  AuxFefun_MultiIndex = discreteform->AuxFeFun_MultiIndex;
  N_AuxFefun_Values = discreteform->N_AuxFeFun_Values;
  AuxFefun_Values = discreteform->AuxFeFun_Values;
  //printf ( "!!!!!!!!!55555555555\n" );
  //提出用户函数
  Functions3D **UserFunction;
  int N_UserFunction;
  UserFunction=discreteform->UserFunction;
  N_UserFunction=discreteform->N_UserFunction;

  int N_FeConst;
  N_FeConst=discreteform->N_FeConst;
  double *FeConst;
  FeConst=discreteform->FeConst;

  VOLU **Volus, *volu;
  FACEE **Faces, *face;
  LINE **Lines, *line;
  Volus = mesh->Volus;
  Faces = mesh->Faces;
  Lines = mesh->Lines;  
  int elem_lines, ID_bd, Ancestor;
  
  BoundCondFunct3D *BoundCond;
  BoundCond = anasatz_space->BoundCond;
  BoundType bdtype; 
  
  double *Local_Prolong;
  Local_Prolong = malloc(n_bas_test*n_bas_test*sizeof(double));
  
  ELEMENT3D *coarse_elem3D,*finer_elem3D,*elem;
  coarse_elem3D = malloc(sizeof(ELEMENT3D));
  InitialElem3D(mesh,coarse_elem3D);
  finer_elem3D = malloc(sizeof(ELEMENT3D));
  InitialElem3D(mesh_finer,finer_elem3D);
  elem = malloc(sizeof(ELEMENT3D));
  InitialElem3D(mesh_finer,elem3D);
  /**进行单元循环，得到刚度矩阵 */
  //printf ( "kai shi xi wang ge xun huan\n" );
  for(curr_volu=0;curr_volu<num_volus;curr_volu++)
  { 
    volu = mesh_finer->Volus[curr_volu];
    Ancestor=volu->Ancestor;
    //取出当前的有限元
    GetElement3D(mesh_finer,curr_volu,elem3D);
    memset(Mat_Entries,0.0,N_Mat_Entries*sizeof(double));
    // iteration for the quadrature points
     if(discreteformvolu!=NULL)
     {
       for(i=0;i< Quad_Points3D;i++)
       {
	//get the quadrature information 
	Quad_xi   =  Quad3D->Quad_X[i];
	Quad_eta  =  Quad3D->Quad_Y[i];
	Quad_zeta =  Quad3D->Quad_Z[i];
	Quad_W    =  Quad3D->Quad_W[i];
	SubAssembleMatrix3D(n_bas_test,base_finer,n_bas_anasatz,base_finer,
			    n_test_multiindex,test_MultiIndex,n_anasatz_multiindex,anasatz_MultiIndex, 
		            N_AuxFefunct,AuxFefunct_finer,N_AuxFefun_MultiIndex,AuxFefun_MultiIndex,
		            N_AuxFefun_Values, AuxFefun_Values,N_UserFunction,UserFunction, 
			    N_FeConst,FeConst,discreteformvolu,elem3D,
			    Test_Values,Anasatz_Values,Test_MultiIndex_Values, Anasatz_MultiIndex_Values,
			    Quad_xi,Quad_eta,Quad_zeta,Quad_W,Mat_Entries);
       } //end for i (Quad_Points)       
     } //end for if

    GetElement3D(mesh,Ancestor,coarse_elem3D);
    GetElement3D(mesh_finer,curr_volu,finer_elem3D);
    GetLocalProlong(coarse_elem3D,test_Base,n_bas_test,finer_elem3D,base_finer,n_bas_finer,
	          elem3D,Local_Prolong);
    MatMatMat(Local_Prolong,Mat_Entries,tmp_entries); 
    for(i=0;i<N_Mat_Entries;++i)
    Mat_Entries[i]=tmp_entries[i];
    AddSubMatrix(Matrix,Mat_Entries,n_bas_test,Test_GlobalNumbers+Test_BeginIndex[Ancestor],
		 n_bas_anasatz,Anasatz_GlobalNumbers+Anasatz_BeginIndex[Ancestor]);	
  }//end for curr_volu
/* ------------------------------------------------------------------------------------------------- */

    /** 释放内存 */
    FreeElem3D(elem3D);
    free(Mat_Entries);
    free(Test_Values);
    free(Anasatz_Values);
    free(Test_MultiIndex_Values);
    free(Anasatz_MultiIndex_Values);
    //free(AuxFefun_Values);
    //free(normal);
    //free(tangent);
    //free(AuxFunct_Values);
}//end for curr_face  




/**合成刚度矩阵 */
void AssembleMatrixMix(MATRIX *Matrix,MATRIX *MatrixB)
{
  // define some notion
  DISCRETEFORM3D *discreteform;
  Fespace3D *test_space, *anasatz_space;
  discreteform = Matrix->DiscreteForm3D;
  test_space = discreteform->Test_Space;
  anasatz_space = discreteform->Anasatz_Space;
  //get the mesh: all the fem spaces are defined on the same mesh
  MESH *mesh;
  mesh = test_space->Mesh;
  int num_volus, i, j;
  num_volus = mesh->Num_Volus_Global;
  Quadrature3D *Quad3D;
  Quad3D = discreteform->Quad3D;
  Quadrature2D *Quad2D;
  Quad2D = discreteform->Quad2D;
  Quadrature1D *Quad1D;
  Quad1D = discreteform->Quad1D;
  ELEMENT3D *elem3D;
  elem3D = malloc(1*sizeof(ELEMENT3D));
  // initializing the element
  InitialElem3D(mesh,elem3D);
  
  int n_bas_test, n_bas_anasatz, N_Mat_Entries, dim_anasatz_value, dim_test_value;
  double *Mat_Entries;
  n_bas_test = test_space->Base->Num_Bas;
  dim_test_value = test_space->Base->Value_Dim;
  n_bas_anasatz = anasatz_space->Base->Num_Bas;
  dim_anasatz_value = anasatz_space->Base->Value_Dim;
  
  N_Mat_Entries = n_bas_test * n_bas_anasatz;
  Mat_Entries = malloc(N_Mat_Entries*sizeof(double));
    
  int *Test_GlobalNumbers, *Anasatz_GlobalNumbers, 
      *Test_BeginIndex, *Anasatz_BeginIndex;
  Test_GlobalNumbers = test_space->GlobalNumbers;
  Test_BeginIndex = test_space->BeginIndex;
  Anasatz_GlobalNumbers = anasatz_space->GlobalNumbers;
  Anasatz_BeginIndex = anasatz_space->BeginIndex;
  
  int n_multiindex_test, n_multiindex_anasatz;
  n_multiindex_test = discreteform->N_Test_MultiIndex;
  n_multiindex_anasatz = discreteform->N_Anasatz_MultiIndex;
  
  double *Test_Values, *Anasatz_Values;
  double *Test_MultiIndex_Values, *Anasatz_MultiIndex_Values;
  
  Test_Values = malloc(n_bas_test*dim_test_value*n_multiindex_test*sizeof(double));
  Anasatz_Values = malloc(n_bas_anasatz*dim_anasatz_value*n_multiindex_anasatz*sizeof(double));
 
  Test_MultiIndex_Values = malloc(n_multiindex_test*dim_test_value*sizeof(double));
  Anasatz_MultiIndex_Values = malloc(n_multiindex_anasatz*dim_anasatz_value*sizeof(double));
  double Quad_xi, Quad_eta,Quad_zeta, Quad_W;
  int curr_volu, kk,ii, Quad_Points;//????????????????????;
  int Quad_Points3D,Quad_Points2D, Quad_Points1D;
  Quad_Points3D = Quad3D->N_Points;
  if(Quad1D)
    Quad_Points1D = Quad1D->N_Points;
  if(Quad2D)
    Quad_Points2D = Quad2D->N_Points;
  int n_test_multiindex,n_anasatz_multiindex;
  MultiIndex3D *test_MultiIndex, *anasatz_MultiIndex;
  //下面两个定义是否重复了？ 没有重复 ？？？？？？？？？？？？？？？？
  n_test_multiindex = discreteform->N_Test_MultiIndex;
  n_anasatz_multiindex = discreteform->N_Anasatz_MultiIndex;
  
  test_MultiIndex = discreteform->Test_MultiIndex;
  anasatz_MultiIndex = discreteform->Anasatz_MultiIndex;
  
  DiscreteFormMatrix *discreteformvolu;
  discreteformvolu = discreteform->DiscreteFormVolu;
  
  DiscreteFormMatrixFace *discreteformface; // *discreteformbd;
  discreteformface = discreteform->DiscreteFormFace;
  
  BASEFUNCTION3D *test_Base, *anasatz_Base;
  test_Base = discreteform->Test_Space->Base;
  anasatz_Base = discreteform->Anasatz_Space->Base;   

  int N_AuxFefunct, *N_AuxFefun_MultiIndex, N_AuxFefun_Values;
  double *AuxFefun_Values;
  MultiIndex3D *AuxFefun_MultiIndex;
  Fefunction3D **AuxFefunct;
  N_AuxFefunct = discreteform->N_AuxFeFun;
  AuxFefunct = discreteform->AuxFeFun;
  N_AuxFefun_MultiIndex = discreteform->N_AuxFeFun_MultiIndex;
  AuxFefun_MultiIndex = discreteform->AuxFeFun_MultiIndex;
  N_AuxFefun_Values = discreteform->N_AuxFeFun_Values;
  AuxFefun_Values = discreteform->AuxFeFun_Values;
  
  //提出用户函数
  Functions3D **UserFunction;
  int N_UserFunction;
  UserFunction=discreteform->UserFunction;
  N_UserFunction=discreteform->N_UserFunction;

  int N_FeConst=discreteform->N_FeConst;
  double *FeConst;
  FeConst=discreteform->FeConst;

  VOLU **Volus, *volu;
  FACEE **Faces, *face;
  LINE **Lines, *line;
  Volus = mesh->Volus;
  Faces = mesh->Faces;
  Lines = mesh->Lines;  
  int elem_lines, ID_bd;
  
  BoundCondFunct3D *BoundCond;
  BoundCond = anasatz_space->BoundCond;
  BoundType bdtype; 
  //iteration of the faces
  /**进行单元循环，得到刚度矩阵 */
  for(curr_volu=0;curr_volu<num_volus;curr_volu++)
  {
    //volu = mesh->Volus[curr_volu];
    //取出当前的有限元
    GetElement3D(mesh,curr_volu,elem3D);
    memset(Mat_Entries,0.0,N_Mat_Entries*sizeof(double));
    // iteration for the quadrature points
    if(discreteformvolu!=NULL)
     {
       //printf ( "volu[%d]=%lf\n",curr_volu,elem3D->Volu );
       for(i=0;i< Quad_Points3D;i++)
       {
	//get the quadrature information 
	Quad_xi   =  Quad3D->Quad_X[i];
	Quad_eta  =  Quad3D->Quad_Y[i];
	Quad_zeta =  Quad3D->Quad_Z[i];
	Quad_W    =  Quad3D->Quad_W[i];
	SubAssembleMatrix3D(n_bas_test,test_Base,n_bas_anasatz,anasatz_Base,
			    n_test_multiindex,test_MultiIndex,n_anasatz_multiindex,anasatz_MultiIndex, 
		            N_AuxFefunct,AuxFefunct,N_AuxFefun_MultiIndex,AuxFefun_MultiIndex,
		            N_AuxFefun_Values, AuxFefun_Values,N_UserFunction,UserFunction,N_FeConst,FeConst, 
			    discreteformvolu,elem3D,
			    Test_Values,Anasatz_Values,Test_MultiIndex_Values, Anasatz_MultiIndex_Values,
			    Quad_xi,Quad_eta,Quad_zeta,Quad_W,Mat_Entries);
       } //end for i (Quad_Points)       
     } //end for if 
       AddSubMatrix(Matrix,Mat_Entries, n_bas_test,Test_GlobalNumbers+Test_BeginIndex[curr_volu],
		 n_bas_anasatz,Anasatz_GlobalNumbers+Anasatz_BeginIndex[curr_volu]);	
  }//end for curr_volu

  double *bb;
  bb=malloc(sizeof(double)*test_space->DOF_ALL);
  MatrixDotVec(MatrixB,AuxFefunct[0]->Values,bb);  
  for(i=0;i<test_space->DOF_ALL;++i)
  {
   Matrix->Entries[Matrix->RowPtr[i+1]-1]=-bb[i];
  }
  for(i=0;i<test_space->DOF_ALL;++i)
  {
   Matrix->Entries[Matrix->RowPtr[test_space->DOF_ALL]+i]=-bb[i];
  }
   Matrix->Entries[Matrix->RowPtr[test_space->DOF_ALL+1]-1]=0;
    
    /** 释放内存 */
    FreeElem3D(elem3D);
    free(Mat_Entries);
    free(Test_Values);
    free(Anasatz_Values);
    free(Test_MultiIndex_Values);
    free(Anasatz_MultiIndex_Values);
    //free(AuxFefun_Values);
    //free(normal);
    //free(tangent);
    //free(AuxFunct_Values);
}//end for curr_face  


//自洽场迭代for BEC
void ScfIteration(DISCRETEFORM3D *discreteformA, DISCRETEFORM3D*discreteformB,int nev,double **vecs,double *evals,double tol, int maxinum)
{
  MATRIX *matrixA,*matrixB;
  matrixA = BuildMatrix(discreteformA);
  matrixB = BuildMatrix(discreteformB);
  AssembleMatrix(matrixA);
  AssembleMatrix(matrixB);
  /**求解特征值问题*/
  //treat the Dirichlet boundary condition
  BoundaryFunction3D *boundaryfunction;
  boundaryfunction = BoundFuncMul;
  //处理边界条件
  BoundaryTreatmentEigen(matrixA, matrixB, boundaryfunction);
  //调用特征值求解的程序
  double tmp_eigen;
  printf("Start Eigenvalue solving!\n");
  DNEigenSolver(matrixA, matrixB, nev, -1, 1.0, evals, vecs,1);
  tmp_eigen=evals[0];
  int i,iter,dof_all;
  dof_all=discreteformA->Test_Space->DOF_ALL;
  Fefunction3D **auxfefun;
  auxfefun=discreteformA->AuxFeFun;
  for(i=0;i<dof_all;++i)
    auxfefun[0]->Values[i]=vecs[0][i];
  Normalization(auxfefun[0]);

  //--------------------------------------开始自洽场迭代-----------------------------------------------  
  for(iter=0;iter<maxinum;++iter)
  {
  printf ( "--------------------------begin the %d-th iteration----------------------\n",iter );
  /** Assemble the matrix */ 
  for(i=0;i<matrixA->N_Entries;++i)
     matrixA->Entries[i]=0.0;
  AssembleMatrix(matrixA);
  /**求解特征值问题*/
  BoundaryFunction3D *boundaryfunction;
  boundaryfunction = BoundFuncMul;
  //处理边界条件
  BoundaryTreatmentEigen(matrixA, matrixB, boundaryfunction);

  //调用特征值求解的程序
  printf("Start Eigenvalue solving!\n");
  DNEigenSolver(matrixA, matrixB, nev, -1, 1.0, evals, vecs,1);
  if(fabs(evals[0]-tmp_eigen)<=tol)
  {
   break;
  }
  else
     tmp_eigen=evals[0];
  for(i=0;i<dof_all;++i)
  {
    auxfefun[0]->Values[i]=vecs[0][i];
    //auxfefun[0]->Values[i]=0.5*auxfefun[0]->Values[i]+0.5*vecs[0][i];
  }
  Normalization(auxfefun[0]);
}//---------------------------------------自洽场迭代结束----------------------------------------
  //Fefunction3D *fef;
  //fef = BuildFefunction3D(discreteformA->Test_Space,vecs[0]);
  //Normalization(fef);
  FreeMatrix(matrixA);
  FreeMatrix(matrixB);
}




//自洽迭代求多元BEC
void ScfIterationCoupleBEC(DISCRETEFORM3D *discreteformA1, DISCRETEFORM3D *discreteformA2, DISCRETEFORM3D*discreteformB,int nev,double *eigenvectors1,double *evals1,double *eigenvectors2,double *evals2,double tol, int maxinum)
{
  /** Assemble the matrix */
  int i,iter,dof_all = discreteformA1->Test_Space->DOF_ALL; 
  Fespace3D *fesp3D;
  fesp3D = discreteformA1->Test_Space; 
  MATRIX *matrixA1,*matrixA2,*matrixB;
  matrixA1 = BuildMatrix(discreteformA1);
  matrixA2 = BuildMatrix(discreteformA2);
  matrixB  = BuildMatrix(discreteformB);
  AssembleMatrix(matrixA1);
  AssembleMatrix(matrixA2);
  AssembleMatrix(matrixB);
  /**求解特征值问题*/
  //treat the Dirichlet boundary condition
  BoundaryFunction3D *boundaryfunction;
  boundaryfunction = BoundFuncMul;
  //处理边界条件
  BoundaryTreatmentEigen(matrixA1, matrixB, boundaryfunction);
  BoundaryTreatmentEigen(matrixA2, matrixB, boundaryfunction);
  //调用特征值求解的程序
  double **vecs1,**vecs2;
  vecs1 = malloc(sizeof(double*));
  vecs1[0] = malloc(sizeof(double)*dof_all);
  vecs2 = malloc(sizeof(double*));
  vecs2[0] = malloc(sizeof(double)*dof_all);
  for(i=0;i<dof_all;++i)
  {
    vecs1[0][i] = 0.0;
    vecs2[0][i] = 0.0;
  }
  printf("Start Eigenvalue solving!\n");
  DNEigenSolver(matrixA1, matrixB, nev, -1, 1.0, evals1, vecs1,1);
  //exit(0);
  DNEigenSolver(matrixA2, matrixB, nev, -1, 1.0, evals2, vecs2,1);

  //exit(0);
  Fefunction3D **auxfefun;
  auxfefun = discreteformA1->AuxFeFun; 
  double tmp_eigen1,tmp_eigen2;

  tmp_eigen1=evals1[0];
  tmp_eigen2=evals2[0];
  for(i=0;i<dof_all;++i)
  {
    auxfefun[0]->Values[i]=vecs1[0][i];
    auxfefun[1]->Values[i]=vecs2[0][i];
  }
  Normalization(auxfefun[0]);
  Normalization(auxfefun[1]);
  //--------------------------------------开始自洽场迭代-----------------------------------------------  
  for(iter=0;iter<maxinum;++iter)
  {
   printf ( "--------------------------begin the %d-th iteration----------------------\n",iter );
  /** Assemble the matrix */ 
   for(i=0;i<matrixA1->N_Entries;++i)
   {
      matrixA1->Entries[i]=0.0;
      matrixA2->Entries[i]=0.0;
   }
   AssembleMatrix(matrixA1);
   AssembleMatrix(matrixA2);
   /**求解特征值问题*/
   //处理边界条件
   BoundaryTreatmentEigen(matrixA1, matrixB, boundaryfunction);
   BoundaryTreatmentEigen(matrixA2, matrixB, boundaryfunction);

   //调用特征值求解的程序
   printf("Start Eigenvalue solving!\n");
   DNEigenSolver(matrixA1, matrixB, nev, -1, 1.0, evals1, vecs1,1);
   DNEigenSolver(matrixA2, matrixB, nev, -1, 1.0, evals2, vecs2,1);
   Normalization_Vector(fesp3D,vecs1[0]);
   Normalization_Vector(fesp3D,vecs2[0]);

   if((fabs(evals1[0]-tmp_eigen1)<=1e-3)&&(fabs(evals2[0]-tmp_eigen2)<=1e-3))
   {
     printf ( "the first eigenvalue error==%f,%f\n", evals1[0]-tmp_eigen1,1e-3);
     printf ( "the second eigenvalue erro =%f,%f\n", evals2[0]-tmp_eigen2,1e-3);
     break;
   }
   else
   {
	  tmp_eigen1=evals1[0];
	  tmp_eigen2=evals2[0];
      printf ( "the first eigenvalue error==%f\n", evals1[0]);
     printf ( "the second eigenvalue erro =%f\n", evals2[0]);
   }
   for(i=0;i<dof_all;++i)
   {
     auxfefun[0]->Values[i]=0.3*auxfefun[0]->Values[i]+0.7*vecs1[0][i];
     auxfefun[1]->Values[i]=0.3*auxfefun[1]->Values[i]+0.7*vecs2[0][i];
   }

   Normalization(auxfefun[0]);
   Normalization(auxfefun[1]);
  }//---------------------------------------自洽场迭代结束----------------------------------------

  for(i=0;i<dof_all;++i)
  {
    eigenvectors1[i] = vecs1[0][i];
    eigenvectors2[i] = vecs2[0][i];
  }

  free(vecs1[0]);
  free(vecs2[0]);
  free(vecs1);
  free(vecs2);


}




//自洽迭代求多元BEC
void ScfIterationCoupleBECGeneral(DISCRETEFORM3D **discreteformA, DISCRETEFORM3D **discreteformB, int nev, double ***vecs, double **evals, int number, double tol, int maxinum, double alpha)
{
  /** Assemble the matrix */
  int i,j,k,iter,dof_all = discreteformA[0]->Test_Space->DOF_ALL; 
  Fespace3D *fesp3D;
  fesp3D = discreteformA[0]->Test_Space; 
  MATRIX **matrixA,*matrixB;
  matrixA = malloc(sizeof(MATRIX*)*nev);
  for(i=0;i<nev;++i)
    matrixA[i] = BuildMatrix(discreteformA[i]);
  matrixB  = BuildMatrix(discreteformB[0]);
  for(i=0;i<nev;++i)
    AssembleMatrix(matrixA[i]);
  AssembleMatrix(matrixB);
  //处理边界条件
  BoundaryFunction3D *boundaryfunction;
  boundaryfunction = BoundFuncMul;
  for(i=0;i<nev;++i)
    BoundaryTreatmentEigen(matrixA[i], matrixB, boundaryfunction);
  //调用特征值求解的程序
  printf("Start Eigenvalue solving!\n");
  for(i=0;i<nev;++i)
    DNEigenSolver(matrixA[i], matrixB, 1, -1, 1.0, evals[i], vecs[i],1);
  Fefunction3D **auxfefun;
  auxfefun = discreteformA[0]->AuxFeFun; 
  double tmp_eigen[nev];

  for(i=0;i<nev;++i)
    tmp_eigen[i]=evals[i][0];
  for(i=0;i<nev;++i)
  {
    for(j=0;j<dof_all;++j)
      auxfefun[i]->Values[j] = alpha*auxfefun[i]->Values[j]+(1-alpha)*vecs[i][0][j];
    Normalization(auxfefun[i]);
  }
  //--------------------------------------开始自洽场迭代-----------------------------------------------  
  for(iter=0;iter<maxinum;++iter)
  {
   printf ( "--------------------------begin the %d-th iteration----------------------\n",iter );
   /** Assemble the matrix */ 
   for(i=0;i<nev;++i)
   {
     for(j=0;j<matrixA[i]->N_Entries;++j)
       matrixA[i]->Entries[j]=0.0;
   }
   for(i=0;i<nev;++i)
      AssembleMatrix(matrixA[i]);
   //处理边界条件
   for(i=0;i<nev;++i)
      BoundaryTreatmentEigen(matrixA[i], matrixB, boundaryfunction);
   //调用特征值求解的程序
   printf("Start Eigenvalue solving!\n");
   for(i=0;i<nev;++i)
     DNEigenSolver(matrixA[i], matrixB, 1, -1, 1.0, evals[i], vecs[i],1);

   if((fabs(evals[0][0]-tmp_eigen[0])<=1e-3)&&(fabs(evals[1][0]-tmp_eigen[1])<=1e-3))
   {
     printf ( "the first  eigenvalue error==%f,%f\n", evals[0][0]-tmp_eigen[0],1e-3);
     printf ( "the second eigenvalue error =%f,%f\n", evals[1][0]-tmp_eigen[1],1e-3);
     break;
   }
   else
   {
     for(i=0;i<nev;++i)
     {
       tmp_eigen[i]=evals[i][0];
       printf ( "the first eigenvalue error==%f\n", evals[i][0]);
     } 
   }
   for(i=0;i<nev;++i)
   {
     for(j=0;j<dof_all;++j)
       auxfefun[i]->Values[j]=alpha*auxfefun[i]->Values[j]+(1-alpha)*vecs[i][0][j];
   }
   for(i=0;i<nev;++i)
     Normalization(auxfefun[i]);
  }//---------------------------------------自洽场迭代结束----------------------------------------

  for(i=0;i<nev;++i)
    FreeMatrix(matrixA[i]);
  FreeMatrix(matrixB);
  free(matrixA);
  

}




/* discreteformA中的auxfefun是迭代过程中用的的迭代fefunction, 参数中的Fefunct是作为基函数不变的 */
void ScfIterationCoupleBECCorrection(DISCRETEFORM3D **discreteformA, DISCRETEFORM3D **discreteformB, RHST **Rhs, int nev, double ***vecs,
                            double **evals, int number, Fefunction3D **Fefunct, double tol, int maxinum, double alpha)
{
  printf("begin coarse mesh correction!\n");
  int i,j,k,iter,dof_all;
  double *quad;
  quad = malloc(sizeof(double));
  dof_all = discreteformA[0]->Test_Space->DOF_ALL+1;
  MATRIX *Prolong;
  Prolong=malloc(sizeof(MATRIX));
  BuildProlongAn(discreteformA[0]->Test_Space, Fefunct[0]->Fespace, Prolong);
  int dof_fine = Fefunct[0]->Fespace->DOF_ALL;
  MATRIX *matrixb, **matrixc, **matrixd,*matrix_tmp;
  matrixc = malloc(sizeof(MATRIX*)*nev);
  matrixd = malloc(sizeof(MATRIX*)*nev);
  matrixb = BuildMatrix(discreteformB[0]);
  AssembleMatrix(matrixb);
  matrix_tmp = BuildMatrix(discreteformA[0]);
  for(i=0;i<nev;++i) 
  {
    matrixd[i] = ExtensionMatrixSim(matrix_tmp,1);
    matrixd[i]->DiscreteForm3D = discreteformA[i];
    AssembleMatrixMulSimplify(matrixd[i],Fefunct[i]);
  }
  for(i=0;i<nev;++i) 
  {
    matrixc[i] = ExtensionMatrix(matrixd[i],1);
    matrixc[i]->DiscreteForm3D = discreteformA[i];
  }
  MATRIX *MatrixR;
  MatrixR=malloc(sizeof(MATRIX));
  TransposeMatrix(Prolong,MatrixR);
  double *bb,*cc;
  cc=malloc(sizeof(double));
  bb=malloc(sizeof(double)*(dof_all-1));

  for(i=0;i<nev;++i) 
  {
    MatMatVec(MatrixR,matrixb,Fefunct[i]->Values,bb);//初始auxfefun中只有最后一个分量非零，即存的值即为Fefunct
    VecMatVec(Fefunct[i]->Values,matrixb,cc);
    for(j=0;j<dof_all-1;++j)
    {
      matrixc[i]->Entries[matrixc[i]->RowPtr[j+1]-1] = -bb[j];
      matrixc[i]->Entries[matrixc[i]->RowPtr[dof_all]+j] = -bb[j];
    }  
    matrixc[i]->Entries[matrixc[i]->RowPtr[dof_all]-1]= -cc[0];
    matrixc[i]->Entries[matrixc[i]->RowPtr[dof_all+1]-2]= -cc[0];
    matrixc[i]->Entries[matrixc[i]->RowPtr[dof_all+1]-1]= 0.0 ;
  }
  //右端项
  for(i=0;i<nev;++i)
    for(j=0;j<Rhs[i]->DOF_ALL;++j)
      Rhs[i]->Entries[j] = 0.0;

  double **rhs;
  rhs = malloc(sizeof(double*)*nev);
  for(i=0;i<nev;++i)
  {
    AssembleRHST(Rhs[i]);
    VecDotVecs(Fefunct[i]->Values,Rhs[i]->Entries,dof_fine,cc);
    MatrixDotVec(MatrixR,Rhs[i]->Entries,bb);
    rhs[i] = malloc(sizeof(double)*(dof_all+1));
    for(j=0;j<dof_all-1;++j)
      rhs[i][j] = bb[j];
    rhs[i][dof_all-1] = cc[0];
    rhs[i][dof_all] = -1.0;
  }

  BoundaryFunction3D *boundaryfunction;
  boundaryfunction = BoundFuncMul;
  for(i=0;i<nev;++i)
    BoundaryTreatment(matrixc[i], rhs[i], boundaryfunction); 
  //调用特征值求解的程序
  double tmp_eigen[nev]; 
  double **sol;
  sol = malloc(nev*sizeof(double*));
  for(i=0;i<nev;++i)
  {
    sol[i] = malloc((dof_all+1)*sizeof(double));
    for(j=0;j<dof_all-1;++j)
      sol[i][j] = 0.0;
    sol[i][dof_all-1] =1.0;
    sol[i][dof_all]=discreteformA[i]->FeConst[i];
    CG(matrixc[i],rhs[i],sol[i],1e-4,1000);
    tmp_eigen[i]=sol[i][dof_all];
  }

  printf("eigenvalue==========:%f,%f\n",sol[0][dof_all],sol[1][dof_all]);
  //--------------------------------------开始自洽场迭代-----------------------------------------------  
  for(iter=0;iter<maxinum;++iter)
  {
    printf ( "--------------------------begin the %d-th coarse_iteration----------------------\n",iter );
    for(i=0;i<nev;++i)
    {
      for(j=0;j<dof_all;++j)
        discreteformA[i]->AuxFeFun[i]->Values[j] = alpha*discreteformA[i]->AuxFeFun[i]->Values[j]+(1-alpha)*sol[i][j];
      discreteformA[i]->FeConst[i] = alpha*discreteformA[i]->FeConst[i]+(1-alpha)*sol[i][dof_all];
    }
    for(i=0;i<Rhs[0]->N_AuxFeFun;++i)
    {
      MatrixDotVec(Prolong,discreteformA[i]->AuxFeFun[i]->Values,Rhs[i]->AuxFeFun[i]->Values);
      for(j=0;j<dof_fine;++j)
        Rhs[i]->AuxFeFun[i]->Values[j] += discreteformA[i]->AuxFeFun[i]->Values[dof_all-1]*Fefunct[i]->Values[j];
    }
    for(i=0;i<nev;++i)
    {
      Quad_Vector(Fefunct[0]->Fespace,Rhs[i]->AuxFeFun[i]->Values, quad);
      for(j=0;j<dof_fine;++j)
        Rhs[i]->AuxFeFun[i]->Values[j] /= quad[0];
      for(j=0;j<dof_all-1;++j)
        discreteformA[i]->AuxFeFun[i]->Values[j] /= quad[0];
    }
    for(i=0;i<nev;++i) 
    {
      for(j=0;j<matrixd[i]->N_Entries;++j)
        matrixd[i]->Entries[j]=0.0;
      AssembleMatrixMulSimplify(matrixd[i],Fefunct[i]);
    }
  
    for(i=0;i<nev;++i)
    {
      for(j=0;j<matrixd[i]->N_Rows;++j)
      {
        for(k=0;k<matrixd[i]->RowPtr[j+1]-matrixd[i]->RowPtr[j];++k)
          matrixc[i]->Entries[matrixc[i]->RowPtr[j]+k] = matrixd[i]->Entries[matrixd[i]->RowPtr[j]+k];
      }
      MatMatVec(MatrixR,matrixb,Rhs[i]->AuxFeFun[i]->Values,bb);//初始auxfefun中只有最后一个分量非零，即存的值即为Fefunct
      VecMatVecs(Fefunct[i]->Values,matrixb,Rhs[i]->AuxFeFun[i]->Values,cc);
      for(j=0;j<dof_all-1;++j)
      {
        matrixc[i]->Entries[matrixc[i]->RowPtr[j+1]-1] = -bb[j];
        matrixc[i]->Entries[matrixc[i]->RowPtr[dof_all]+j] = -bb[j];
      }  
      matrixc[i]->Entries[matrixc[i]->RowPtr[dof_all]-1]= -cc[0];
      matrixc[i]->Entries[matrixc[i]->RowPtr[dof_all+1]-2]= -cc[0];
      matrixc[i]->Entries[matrixc[i]->RowPtr[dof_all+1]-1]= 0.0 ;
    }


    //右端项
    for(i=0;i<nev;++i)
      for(j=0;j<Rhs[0]->DOF_ALL;++j)
        Rhs[i]->Entries[j] = 0.0;

    for(i=0;i<nev;++i)
    {
      AssembleRHST(Rhs[i]);
      VecDotVecs(Fefunct[i]->Values,Rhs[i]->Entries,dof_fine,cc);
      MatrixDotVec(MatrixR,Rhs[i]->Entries,bb);
      for(j=0;j<dof_all-1;++j)
        rhs[i][j] = bb[j];
      rhs[i][dof_all-1] = cc[0];
      rhs[i][dof_all] = -1.0;
    }
    for(i=0;i<nev;++i)
      BoundaryTreatment(matrixc[i], rhs[i], boundaryfunction); 
    //调用特征值求解的程序
    for(i=0;i<nev;++i)
      CG(matrixc[i],rhs[i],sol[i],1e-4,100);
  
    if((fabs(sol[0][dof_all]-tmp_eigen[0])<=tol)&&(fabs(sol[1][dof_all]-tmp_eigen[1])<=tol))
      break;
    else
    {
      for(i=0;i<nev;++i)
        tmp_eigen[i]=sol[i][dof_all];
    }
   printf("here==========================%f,%f\n",tmp_eigen[0],tmp_eigen[1]);
  } //---------------------------------------自洽场迭代结束----------------------------------------
    //下面输出细网格维数的有限元解
    //首先构造延拓算子
  for(i=0;i<nev;++i)
  { 
    MatrixDotVec(Prolong,sol[i],vecs[i][number]);  
    for(j=0;j<dof_fine;++j)
     vecs[i][number][j]+= sol[i][dof_all-1]*Fefunct[i]->Values[j];
    evals[i][number] = sol[i][dof_all];
  }

  FreeMatrix(matrixb);
  FreeMatrix(matrix_tmp);
  FreeMatrix(MatrixR);
  for(i=0;i<nev;++i)
  {
   FreeMatrix(matrixc[i]);
   FreeMatrix(matrixd[i]);
   free(sol[i]);
  }
  free(matrixc);
  free(matrixd);
  free(bb);
  free(cc);
  free(sol);
  free(quad);

}



//自洽场迭代for BEC
void ScfIterationMix( DISCRETEFORM3D *discreteformA, DISCRETEFORM3D*discreteformB, RHST *Rhs, int nev, double **vecs, double *evals,double tol, int maxinum)
{
  MATRIX *matrix,*matrixB;
  matrix  = BuildMatrixMul(discreteformA);
  matrixB = BuildMatrix(discreteformB);
  AssembleMatrix(matrixB);
  AssembleMatrixMix(matrix,matrixB);
  /**建立右端项 */
  AssembleRHSTAllMix(Rhs);
  /**求解线性方程组 */
  double *sol;
  int i,iter,dof_all;
  dof_all=discreteformA->Test_Space->DOF_ALL;
  sol = malloc(sizeof(double)*(dof_all+1));
  for(i=0;i<dof_all;++i)
    sol[i] = vecs[0][i];
  sol[dof_all] = evals[0];
  BoundaryFunction3D *boundaryfunction;
  boundaryfunction = BoundFuncMul;
  //处理边界条件
  BoundaryTreatment(matrix, Rhs->Entries, boundaryfunction);
  //求解
  CG(matrix,Rhs->Entries,sol,tol,200);
  double tmp_eigen;
  tmp_eigen=sol[dof_all];
  discreteformA->FeConst[0]=sol[dof_all];
  Rhs->FeConst[0]=sol[dof_all];
  for(i=0;i<dof_all;++i)
  {
    discreteformA->AuxFeFun[0]->Values[i]=0.5*sol[i]+0.5*discreteformA->AuxFeFun[0]->Values[i];
    Rhs->AuxFeFun[0]->Values[i]=0.5*sol[i]+0.5*Rhs->AuxFeFun[0]->Values[i];
  }
  //--------------------------------------开始自洽场迭代-----------------------------------------------  
  for(iter=0;iter<maxinum;++iter)
  {
    printf ( "--------------------------begin the %d-th iteration----------------------\n",iter );
    for(i=0;i<matrix->N_Entries;++i)
      matrix->Entries[i]=0.0;
    for(i=0;i<Rhs->DOF_ALL;++i)
      Rhs->Entries[i]=0.0;
    AssembleMatrixMix(matrix,matrixB);
    AssembleRHSTAllMix(Rhs);//值是否需要赋0？？？
   //处理边界条件
    BoundaryTreatment(matrix, Rhs->Entries, boundaryfunction);
    //求解
    CG(matrix,Rhs->Entries,sol,tol,200);
    if(fabs(sol[dof_all]-tmp_eigen)<=tol)
    {
      for(i=0;i<dof_all;++i)
         vecs[0][i]=sol[i];
       evals[0] = sol[dof_all];  
      break;
    }
    tmp_eigen=sol[dof_all];
    discreteformA->FeConst[0]=sol[dof_all];
    Rhs->FeConst[0]=sol[dof_all];
    for(i=0;i<dof_all;++i)
    {
      discreteformA->AuxFeFun[0]->Values[i]=0.5*sol[i]+0.5*discreteformA->AuxFeFun[0]->Values[i];
      Rhs->AuxFeFun[0]->Values[i]=0.5*sol[i]+0.5*Rhs->AuxFeFun[0]->Values[i];
    }
  }//---------------------------------------自洽场迭代结束----------------------------------------
  free(sol);
 
}


/* discreteformA中的auxfefun是迭代过程中用的的迭代fefunction, 参数中的Fefunct是作为基函数不变的 */
void ScfIterationMul(DISCRETEFORM3D *discreteformA, DISCRETEFORM3D *discreteformB,int nev,double *vecs,
                     double *evals,Fefunction3D *Fefunct,double tol,int maxinum)
{
  int i,j,k,iter;
  int dof_all;
  dof_all=discreteformA->Test_Space->DOF_ALL+1;
   
  MATRIX *Prolong;
  Prolong=malloc(sizeof(MATRIX));
  BuildProlongAn(discreteformA->Test_Space,Fefunct->Fespace, Prolong);
  Fefunction3D **auxfefun;
  auxfefun=discreteformA->AuxFeFun;
  int dof_fine = Fefunct->Fespace->DOF_ALL;
  double *fef_norm;
  fef_norm=malloc(sizeof(double)*dof_fine);
  Fefunction3D *fef;
  double *norm;
  norm = malloc(sizeof(double)); 
  fef = BuildFefunction3D(Fefunct->Fespace,Fefunct->Values);
  QuadFefunction(fef,norm);
  for(i=0;i<dof_all;++i)
  {
    auxfefun[0]->Values[i]=auxfefun[0]->Values[i]/norm[0];
  }

  MATRIX *matrixa, *matrixb;
  matrixa = BuildMatrixMul(discreteformA); 
  matrixb = BuildMatrixMul(discreteformB);
  //matrixa = ExtensionMatrix(matrix0,1);
  //matrixb = ExtensionMatrix(matrix0,1);
  //matrixa->DiscreteForm3D =  discreteformA;
  //matrixb->DiscreteForm3D =  discreteformB;
  AssembleMatrixMulSimplify(matrixa,Fefunct);
  //OutPutMatrix(matrixa);
  AssembleMatrixMulSimplify(matrixb,Fefunct);
  /**求解特征值问题*/
  //处理边界条件
  BoundaryFunction3D *boundaryfunction;
  boundaryfunction = BoundFuncMul;
  BoundaryTreatmentEigen(matrixa, matrixb, boundaryfunction); 
  //调用特征值求解的程序
  double tmp_eigen; 
  double *sol1;
  sol1 = malloc(nev*dof_all*sizeof(double));
  //memset(sol1,0.0,nev*dof_all*sizeof(double));
  for(i=0;i<dof_all;++i)
    sol1[i] = auxfefun[0]->Values[i];
  double **vecs1;
  vecs1 = malloc(nev*sizeof(double*));
  for(i=0;i<nev;i++)
  {
    vecs1[i] = sol1+i*dof_all;    
  }  
  printf("Start Eigenvalue solving!here!\n");
  DNEigenSolver(matrixa, matrixb, nev, -1, 1.0, evals, vecs1,1);
  tmp_eigen=evals[0];
  //--------------------------------------开始自洽场迭代-----------------------------------------------  
  for(iter=0;iter<maxinum;++iter)
  {
  printf ( "--------------------------begin the %d-th coarse_iteration----------------------\n",iter );
  /** Assemble the matrix */ 
   for(i=0;i<dof_all;++i)
  {
    auxfefun[0]->Values[i]=vecs1[0][i];
    auxfefun[0]->Values[i]=0.01*auxfefun[0]->Values[i]+0.99*vecs1[0][i];
    //auxfefun[0]->Values[i]=0.9*vecs1[0][i]+0.1*auxfefun[0]->Values[i];
  }

  for(i=0;i<matrixa->N_Entries;++i)
     matrixa->Entries[i]=0.0;
  AssembleMatrixMulSimplify(matrixa,Fefunct);
  /**求解特征值问题*/
  BoundaryFunction3D *boundaryfunction;
  boundaryfunction = BoundFuncMul;
  //处理边界条件
  BoundaryTreatmentEigen(matrixa, matrixb, boundaryfunction);
  //调用特征值求解的程序
  printf("Start Eigenvalue solving!\n");
  DNEigenSolver(matrixa, matrixb, nev, -1, 1.0, evals, vecs1,1);
  if(fabs(evals[0]-tmp_eigen)<=tol)
  {
   break;
  }
  else
     tmp_eigen=evals[0];

  printf("eigenvalue=%f\n",evals[0]);  
  //for(i=0;i<dof_all;++i)
  //{
  //  auxfefun[0]->Values[i]=vecs1[0][i];
  //}
  //MatrixDotVec(Prolong, vecs1[0],fef_norm);  
  //for(k=0;k<dof_fine;++k)
 //   fef_norm[k]+= vecs1[0][dof_all-1]*Fefunct->Values[k];
 // fef = BuildFefunction3D(Fefunct->Fespace,fef_norm);
 // QuadFefunction(fef,norm);
}//---------------------------------------自洽场迭代结束----------------------------------------
  //下面输出细网格维数的有限元解
  //首先构造延拓算子
  MatrixDotVec(Prolong,vecs1[0],vecs);  
  for(k=0;k<dof_fine;++k)
   vecs[k]+= vecs1[0][dof_all-1]*Fefunct->Values[k];
  printf ( "eigenvalue=%f\n",evals[0] );
  //fef = BuildFefunction3D(Fefunct->Fespace,vecs);
  //Normalization(fef);
  FreeMatrix(matrixb);
  FreeMatrix(matrixa);
}



/* discreteformA中的auxfefun是迭代过程中用的的迭代fefunction, 参数中的Fefunct是作为基函数不变的 */
void ScfIterationCorrection(DISCRETEFORM3D *discreteformA, DISCRETEFORM3D *discreteformB, RHST *Rhs, int nev, double *vecs,
                            double *evals, Fefunction3D *Fefunct, double tol, int maxinum)
{
  printf("begin coarse mesh correction!\n");
  int i,j,k,iter,dof_all;
  dof_all = discreteformA->Test_Space->DOF_ALL+1;
  MATRIX *Prolong;
  Prolong=malloc(sizeof(MATRIX));
  BuildProlongAn(discreteformA->Test_Space, Fefunct->Fespace, Prolong);
  int dof_fine = Fefunct->Fespace->DOF_ALL;
  MATRIX *matrixa, *matrixb, *matrixc, *matrixd;
  /*
  matrixa = BuildMatrixMul(discreteformA); 
  matrixb = BuildMatrix(discreteformB);
  AssembleMatrixMul(matrixa,Fefunct);
  AssembleMatrix(matrixb);
   matrixc = ExtensionMatrix(matrixa,1);
  matrixc->DiscreteForm3D = discreteformA;
  MATRIX *MatrixR;
  MatrixR=malloc(sizeof(MATRIX));
  TransposeMatrix(Prolong,MatrixR);
  double *bb,*cc;
  cc=malloc(sizeof(double));
  bb=malloc(sizeof(double)*(dof_all-1));
  MatMatVec(MatrixR,matrixb,Fefunct->Values,bb);//初始auxfefun中只有最后一个分量非零，即存的值即为Fefunct
  VecMatVec(Fefunct->Values,matrixb,cc);
  for(i=0;i<dof_all-1;++i)
  {
   matrixc->Entries[matrixc->RowPtr[i+1]-1] = -bb[i];
  }
  for(i=0;i<dof_all-1;++i)
  {
   matrixc->Entries[matrixc->RowPtr[dof_all]+i] = -bb[i];
  }
   matrixc->Entries[matrixc->RowPtr[dof_all]-1]= -cc[0];
   matrixc->Entries[matrixc->RowPtr[dof_all+1]-2]= -cc[0];
   matrixc->Entries[matrixc->RowPtr[dof_all+1]-1]= 0 ;
  */
  Fefunction3D **AuxFefunct_finer,**AuxFefunct_coarse;
  int N_AuxFefunct = Rhs->N_AuxFeFun;
  AuxFefunct_finer=malloc(sizeof(Fefunction3D*)*N_AuxFefunct);
  double *value1;
  value1=malloc(sizeof(double)*dof_fine);
  for(j=0;j<dof_fine;++j)
     value1[j] = Fefunct->Values[j];
  AuxFefunct_finer[0] = BuildFefunction3D(Fefunct->Fespace,value1);
  AuxFefunct_coarse = discreteformA->AuxFeFun;
  matrixb = BuildMatrix(discreteformB);
  AssembleMatrix(matrixb);
 
  matrixd = BuildMatrix(discreteformA);
  AssembleMatrixHinh(matrixd,Fefunct); 
  matrixc = ExtensionMatrix(matrixd,2);
  matrixc->DiscreteForm3D = discreteformA;
  Fespace3D *fespold;
  fespold = discreteformA->Test_Space;  
  discreteformA->AuxFeFun = AuxFefunct_finer;  
  discreteformA->Test_Space = Fefunct->Fespace;  
  discreteformA->Anasatz_Space = Fefunct->Fespace;  
  matrixa = BuildMatrix(discreteformA);
  AssembleMatrix(matrixa);
  MATRIX *MatrixR;
  MatrixR=malloc(sizeof(MATRIX));
  TransposeMatrix(Prolong,MatrixR);
  double *bb,*cc;
  cc=malloc(sizeof(double));
  bb=malloc(sizeof(double)*(dof_all-1));
  MatMatVec(MatrixR,matrixa,Fefunct->Values,bb);//初始auxfefun中只有最后一个分量非零，即存的值即为Fefunct
  VecMatVec(Fefunct->Values,matrixa,cc);
 
  for(i=0;i<dof_all-1;++i)
  {
   matrixc->Entries[matrixc->RowPtr[i+1]-2] = bb[i];
  }
  for(i=0;i<dof_all-1;++i)
  {
   matrixc->Entries[matrixc->RowPtr[dof_all-1]+i] = bb[i];
  }
   matrixc->Entries[matrixc->RowPtr[dof_all]-2]= cc[0];


  MatMatVec(MatrixR,matrixb,Fefunct->Values,bb);//初始auxfefun中只有最后一个分量非零，即存的值即为Fefunct
  VecMatVec(Fefunct->Values,matrixb,cc);
  for(i=0;i<dof_all-1;++i)
  {
   matrixc->Entries[matrixc->RowPtr[i+1]-1] = -bb[i];
  }
  for(i=0;i<dof_all-1;++i)
  {
   matrixc->Entries[matrixc->RowPtr[dof_all]+i] = -bb[i];
  }
   matrixc->Entries[matrixc->RowPtr[dof_all]-1]= -cc[0];
   matrixc->Entries[matrixc->RowPtr[dof_all+1]-2]= -cc[0];
   matrixc->Entries[matrixc->RowPtr[dof_all+1]-1]= 0 ;

 //右端项
  double *rhs;
  /*
  Fefunction3D **AuxFefunct_finer;
  int N_AuxFefunct = Rhs->N_AuxFeFun;
  AuxFefunct_finer=malloc(sizeof(Fefunction3D*)*N_AuxFefunct);
  double *value1;
  value1=malloc(sizeof(double)*Fefunct->DOF_ALL);
  for(j=0;j<Fefunct->Fespace->DOF_ALL;++j)
     value1[j] = Fefunct->Values[j];
  AuxFefunct_finer[0] = BuildFefunction3D(Fefunct->Fespace,value1);
  */
  //Normalization(AuxFefunct_finer[0]);
  Rhs->AuxFeFun = AuxFefunct_finer;
  AssembleRHST(Rhs);
  //for(i=0;i<Rhs->DOF_ALL;++i)
  //  printf("%d,%f\n",i,Fefunct->Values[i]);


  VecDotVecs(Fefunct->Values,Rhs->Entries,dof_fine,cc);
  MatrixDotVec(MatrixR,Rhs->Entries,bb);
  rhs = malloc(sizeof(double)*(dof_all+1));
  for(i=0;i<dof_all-1;++i)
    rhs[i] = bb[i];
  rhs[dof_all-1] = cc[0];
  rhs[dof_all] = -1;
  //matrixc->DiscreteForm3D->Test_Space = fespold;
  //matrixc->DiscreteForm3D->Anasatz_Space = fespold;
  discreteformA->Test_Space = fespold;
  discreteformA->Anasatz_Space = fespold;
  discreteformA->AuxFeFun = AuxFefunct_coarse;

  BoundaryFunction3D *boundaryfunction;
  boundaryfunction = BoundFuncMul;
  BoundaryTreatment(matrixc, rhs, boundaryfunction); 
  //调用特征值求解的程序
  double tmp_eigen; 
  double *sol;
  sol = malloc((dof_all+1)*sizeof(double));
  for(i=0;i<dof_all-1;++i)
   sol[i] = 0.0;
  sol[dof_all-1] =1;
  sol[dof_all]=discreteformA->FeConst[0] ;
  //求解
  //
  
  //for(i=0;i<dof_all+1;++i)
  //  printf("%d,%f\n",i,rhs[i]);

  CG(matrixc,rhs,sol,1e-4,100);
  tmp_eigen=sol[dof_all];

  printf("eigenvalue==i========              :%f\n",sol[dof_all]);
  //--------------------------------------开始自洽场迭代-----------------------------------------------  
  //for(iter=0;iter<maxinum;++iter)
  for(iter=0;iter<maxinum;++iter)
  {
    printf ( "--------------------------begin the %d-th coarse_iteration----------------------\n",iter );
    //discreteformA->AuxFeFun =AuxFefunct_coarse;  
    for(i=0;i<dof_all;++i)
    {  
      discreteformA->AuxFeFun[0]->Values[i] = 0.5*discreteformA->AuxFeFun[0]->Values[i]+0.5*sol[i];
      //Rhs->AuxFeFun[0]->Values[i] =1.0*Rhs->AuxFeFun[0]->Values[i]+0.0*sol[i];
    }
    discreteformA->FeConst[0] =0.5*discreteformA->FeConst[0]+0.5*sol[dof_all];
    Rhs->FeConst[0] = 0.5*Rhs->FeConst[0]+ 0.5*sol[dof_all];

    //value = AuxFefunct[0]->Values;
    MatrixDotVec(Prolong,discreteformA->AuxFeFun[0]->Values,value1);
    for(j=0;j<Fefunct->Fespace->DOF_ALL;++j)
       value1[j]+=discreteformA->AuxFeFun[0]->Values[dof_all-1]*Fefunct->Values[j];
    AuxFefunct_finer[0] = BuildFefunction3D(Fefunct->Fespace, value1);
    //Normalization(AuxFefunct_finer[0]);
    //discreteformA->AuxFeFun =AuxFefunct_finer;
    for(i=0;i<matrixd->N_Entries;++i)
      matrixd->Entries[i]=0.0;
    AssembleMatrixHinh(matrixd,Fefunct); 
    //matrixc = ExtensionMatrix(matrixd,2);//内存有泄漏
    //matrixc->DiscreteForm3D = discreteformA;
    for(i=0;i<matrixd->N_Rows;++i)
    {
      for(j=0;j<matrixd->RowPtr[i+1]-matrixd->RowPtr[i];++j)
       matrixc->Entries[matrixc->RowPtr[i]+j] = matrixd->Entries[matrixd->RowPtr[i]+j];
    }
    discreteformA->AuxFeFun = AuxFefunct_finer;  
    discreteformA->Test_Space = Fefunct->Fespace;  
    discreteformA->Anasatz_Space = Fefunct->Fespace;  
     for(i=0;i<matrixa->N_Entries;++i)
      matrixa->Entries[i]=0.0;
    AssembleMatrix(matrixa);
    MatMatVec(MatrixR,matrixa,Fefunct->Values,bb);//初始auxfefun中只有最后一个分量非零，即存的值即为Fefunct
    VecMatVec(Fefunct->Values,matrixa,cc);
 
    for(i=0;i<dof_all-1;++i)
    {
      matrixc->Entries[matrixc->RowPtr[i+1]-2] = bb[i];
    }
    for(i=0;i<dof_all-1;++i)
    {
     matrixc->Entries[matrixc->RowPtr[dof_all-1]+i] = bb[i];
    }
    matrixc->Entries[matrixc->RowPtr[dof_all]-2]= cc[0];


    MatMatVec(MatrixR,matrixb,value1,bb);//初始auxfefun中只有最后一个分量非零，即存的值即为Fefunct
    VecMatVecs(Fefunct->Values,matrixb,value1,cc);
    for(i=0;i<dof_all-1;++i)
    {
     matrixc->Entries[matrixc->RowPtr[i+1]-1] = -bb[i];
    }
    for(i=0;i<dof_all-1;++i)
    {
     matrixc->Entries[matrixc->RowPtr[dof_all]+i] = -bb[i];
    }
    matrixc->Entries[matrixc->RowPtr[dof_all]-1]= -cc[0];
    matrixc->Entries[matrixc->RowPtr[dof_all+1]-2]= -cc[0];
    matrixc->Entries[matrixc->RowPtr[dof_all+1]-1]= 0 ;

    
    /*
    AssembleMatrixMul(matrixa,Fefunct);
    for(i=0;i<matrixa->N_Rows;++i)
    { 
      for(j=0;j<matrixa->RowPtr[i+1]-matrixa->RowPtr[i];++j)
      {
        matrixc->Entries[matrixc->RowPtr[i]+j] = matrixa->Entries[matrixa->RowPtr[i]+j];
      }
    }
    //for(i=0;i<matrixb->N_Entries;++i)
    //  matrixb->Entries[i]=0.0;
    //AssembleMatrix(matrixb);
    MatMatVec(MatrixR,matrixb,value1,bb);
    VecMatVecs(Fefunct->Values,matrixb,value1,cc);
    for(i=0;i<dof_all-1;++i)
    {
      matrixc->Entries[matrixc->RowPtr[i+1]-1]= -bb[i];
    }
    for(i=0;i<dof_all-1;++i)
    {
      matrixc->Entries[matrixc->RowPtr[dof_all]+i]= -bb[i];
    }
    matrixc->Entries[matrixc->RowPtr[dof_all]-1]= -cc[0];
    matrixc->Entries[matrixc->RowPtr[dof_all+1]-2]=-cc[0];
    matrixc->Entries[matrixc->RowPtr[dof_all+1]-1]= 0 ;
    //OutPutMatrix(matrixc);
    */

    //右端项
    //Rhs->AuxFeFun = AuxFefunct_finer;
    for(i=0;i<Rhs->DOF_ALL;++i)
      Rhs->Entries[i] = 0.0;
    AssembleRHST(Rhs);
    printf("------------------------------------------------\n");
    VecDotVecs(Fefunct->Values,Rhs->Entries,dof_fine,cc);
    MatrixDotVec(MatrixR,Rhs->Entries,bb);
    for(i=0;i<dof_all-1;++i)
     rhs[i] = bb[i];
    rhs[dof_all-1] = cc[0];
    rhs[dof_all] = -1;
    // matrixc->DiscreteForm3D->Test_Space = fespold;
    // matrixc->DiscreteForm3D->Anasatz_Space = fespold;
    discreteformA->Test_Space = fespold;
    discreteformA->Anasatz_Space = fespold;
    discreteformA->AuxFeFun = AuxFefunct_coarse;
    BoundaryFunction3D *boundaryfunction;
    boundaryfunction = BoundFuncMul;
    BoundaryTreatment(matrixc, rhs, boundaryfunction); 
    //调用特征值求解的程序
    double tmp_eigen; 
    //求解
    //
    CG(matrixc,rhs,sol,1e-4,100);
    //调用求解的程序
    CG(matrixc, rhs,sol,tol,100);
    if(fabs(sol[dof_all]-tmp_eigen)<=tol)
      break;
    else
      tmp_eigen = sol[dof_all];

   printf("here=%f\n",tmp_eigen);
  //if(iter==10) 
  //exit(0);
  } //---------------------------------------自洽场迭代结束----------------------------------------
    //下面输出细网格维数的有限元解
    //首先构造延拓算子
    MatrixDotVec(Prolong,sol,vecs);  
    for(k=0;k<dof_fine;++k)
    vecs[k]+= sol[dof_all-1]*Fefunct->Values[k];
    evals[0] = sol[dof_all];
    FreeMatrix(matrixa);
    FreeMatrix(matrixb);
    FreeMatrix(matrixc);
    FreeMatrix(matrixd);
}


void AssembleRHSTAllMix(RHST *Rhs)
{
  printf ( "begin assemble rhsmix\n" );
  Fespace3D *test_space;
  test_space = Rhs->Test_Space;
  MESH *mesh;
  mesh =test_space->Mesh;
  BASEFUNCTION3D *test_base;
  test_base = test_space->Base;
  int n_bas_test, N_Rhs_Entries;
  double *Rhs_Entries;
  n_bas_test = test_base->Num_Bas;  
  N_Rhs_Entries = n_bas_test;
  Rhs_Entries = malloc(N_Rhs_Entries*sizeof(double));
  int dim_test_value;
  dim_test_value = test_base->Value_Dim;
  
  int *Test_GlobalNumbers, *Test_BeginIndex;
  Test_GlobalNumbers = test_space->GlobalNumbers;
  Test_BeginIndex = test_space->BeginIndex;
  
  int num_volus,i,j;
  num_volus = mesh->Num_Volus_Global;
  Quadrature3D *Quad3D;
  Quadrature2D *Quad2D;
  Quadrature1D *Quad1D;

  Quad3D = Rhs->Quad3D;  
  Quad2D = NULL;  
  Quad1D = NULL;
  if(Rhs->Quad1D)
    Quad1D = Rhs->Quad1D;
  if(Rhs->Quad2D)
    Quad2D = Rhs->Quad2D;
  printf ( "rhs33333\n" );
  
  ELEMENT3D *elem3D;
  elem3D = malloc(1*sizeof(ELEMENT3D));
  InitialElem3D(mesh,elem3D); 
  
  //取出右端项所需要计算的导数值
  int n_multiindex_test;
  n_multiindex_test = Rhs->N_Test_MultiIndex;  
  MultiIndex3D *Test_MultiIndex;
  Test_MultiIndex = Rhs->Test_MultiIndex;
  //printf ( "rhs4444444444444\n" );
  double *Test_Values, *Test_MultiIndex_Values;  
  Test_Values = malloc(n_bas_test*dim_test_value*n_multiindex_test*sizeof(double));
  Test_MultiIndex_Values = malloc(n_multiindex_test*dim_test_value*sizeof(double));
  memset(Test_MultiIndex_Values,0.0,n_multiindex_test*dim_test_value*sizeof(double));
  //在体上的离散变分形式函数
  DiscreteFormRHSVolu *discreteformvolu;  
  //在面上的离散变分形式函数
  DiscreteFormRHSFace *discreteformface;  
  //printf ( "rhs5555555555\n" );
  discreteformvolu = Rhs->DiscreteFormVolu;
  discreteformface = Rhs->DiscreteFormFace;
 
  int N_UserFunction;
  Functions3D **UserFunction;
  int N_FeConst;
  double *FeConst;
  N_UserFunction=Rhs->N_UserFunction;
  UserFunction=Rhs->UserFunction;
  N_FeConst=Rhs->N_FeConst;
  FeConst=Rhs->FeConst;

  double Quad_xi, Quad_eta, Quad_zeta, Quad_W;
  int curr_volu, Quad_Points3D, Quad_Points2D, Quad_Points1D;
  Quad_Points3D = Quad3D->N_Points;
  if(Rhs->Quad1D)
    Quad_Points1D = Quad1D->N_Points;
  if(Rhs->Quad2D)
    Quad_Points2D = Quad2D->N_Points;
  //printf ( "rhs666666666666\n" );
  int N_AuxFeFun, *N_AuxFeFun_MultiIndex, N_AuxFeFun_Values;
  MultiIndex3D *AuxFeFun_MultiIndex;
  N_AuxFeFun = Rhs->N_AuxFeFun;
  AuxFeFun_MultiIndex = Rhs->AuxFeFun_MultiIndex;
  N_AuxFeFun_MultiIndex = Rhs->N_AuxFeFun_MultiIndex;
  //printf ( "rhs7777\n" );
  Fefunction3D **AuxFeFuns;
  N_AuxFeFun_Values = Rhs->N_AuxFeFun_Values;
  AuxFeFuns = Rhs->AuxFeFun;
  double *AuxFeFun_Values;
  AuxFeFun_Values = Rhs->AuxFeFun_Values;
  int *Dim_Values_AuxFeFun;
  //printf ( "rhs777778888888888\n" );
  if(N_AuxFeFun)
  {
    Dim_Values_AuxFeFun = malloc(N_AuxFeFun*sizeof(int));
    for(i=0;i<N_AuxFeFun;i++)
      Dim_Values_AuxFeFun[i] = AuxFeFuns[i]->Fespace->Base->Value_Dim;
  }
  //printf ( "rrhs888888888888888\n" );
  double XI,ETA;
  double local_direction;
  int elem_faces, ID_bd;
 
  BoundType bdtype; 
  BoundCondFunct3D *BoundCond;
  BoundCond = test_space->BoundCond;
  
  //N_AuxFeFun_Values = 0;
  //AuxFeFun_Values = NULL;
  printf ( "beign to loop\n" );

  //for(curr_volu=0;curr_volu<1;curr_volu++)
for(curr_volu=0;curr_volu<num_volus;curr_volu++) 
   {
    //volu = Volus[curr_volu];
    GetElement3D(mesh,curr_volu,elem3D);
    memset(Rhs_Entries,0.0,N_Rhs_Entries*sizeof(double));
    //printf ( "2222222222222222222222222\n" );
    if(discreteformvolu!=NULL)
    {
      for(i=0;i<Quad_Points3D;i++)
      {

	Quad_xi = Quad3D->Quad_X[i];
	Quad_eta = Quad3D->Quad_Y[i];
	Quad_zeta = Quad3D->Quad_Z[i];
	Quad_W = Quad3D->Quad_W[i];
	//先获得所有矩阵和右端项在本积分点的值
	//assemble 
	//每个矩阵包含自己的取值空间和试探空间，后面的微分算子也是针对
	//针对这两个空间的。另外可能还有其他的辅助函数空间和相应的微分 
	//算子，最后也会包含辅助函数来定义矩阵   
	//printf ( "33333333333333333333333333333\n" ); 
	SubAssembleRhs3DAll(elem3D,test_base,Test_Values,n_multiindex_test,Test_MultiIndex,
			  Test_MultiIndex_Values, Dim_Values_AuxFeFun,N_AuxFeFun,AuxFeFuns,
			  N_AuxFeFun_MultiIndex,
			  AuxFeFun_MultiIndex, N_AuxFeFun_Values,AuxFeFun_Values,
                          N_UserFunction,UserFunction,N_FeConst,FeConst,
			  Quad_xi,Quad_eta,Quad_zeta,Quad_W,discreteformvolu,Rhs_Entries); 	
      }
    }//end if
    AddSubRhs(Rhs,Rhs_Entries,n_bas_test,Test_GlobalNumbers+Test_BeginIndex[curr_volu]);
  }//单元循环结束 
    //下面补上最后一个元素
    Rhs->Entries[Rhs->DOF_ALL-1]=-1;
    FreeElem3D(elem3D);
    free(Rhs_Entries);
    free(Test_Values);
    free(Test_MultiIndex_Values);
    //free(normal);
    //free(tangent);
    if(N_AuxFeFun)
    {
      free(Dim_Values_AuxFeFun);
    }
}





void NonScfIteration(DISCRETEFORM3D *discreteform,RHST *rhs,double *sol,double tol,int maxinum)
{
  int i,iter;
  MATRIX *matrix;
  matrix = BuildMatrix(discreteform);
  AssembleMatrix(matrix);
  //treat the Dirichlet boundary condition
  BoundaryFunction3D *boundaryfunction;
  boundaryfunction = BoundFuncMul;
  AssembleRHST(rhs);
  //处理边界条件
  BoundaryTreatment(matrix,rhs->Entries, boundaryfunction);
  //调用求解的程序
  DIRECTSOLVER(matrix,rhs->Entries,sol);
  int dof_all;
  dof_all=discreteform->Test_Space->DOF_ALL;
  Fefunction3D **auxfefun;
  auxfefun=rhs->AuxFeFun;
  for(i=0;i<dof_all;++i)
    auxfefun[0]->Values[i]=sol[i];
  Fefunction3D *tmp_fefun;
  double *tmp_value;
  tmp_value=malloc(sizeof(double)*dof_all);
  Fespace3D *fesp3D;
  fesp3D = discreteform->Test_Space;
  for(i=0;i<dof_all;++i)
     tmp_value[i]=sol[i];
  tmp_fefun = BuildFefunction3D(fesp3D,tmp_value);
  MultiIndex3D MuiltiIndex[1] = {D000};
  //--------------------------------------开始自洽场迭代-----------------------------------------------  
  for(iter=0; iter<maxinum; ++iter)
  {
  printf ( "--------------------------begin the %d-th iteration----------------------\n",iter );
  for(i=0;i<matrix->N_Entries;++i)
     matrix->Entries[i]=0.0;
  for(i=0;i<rhs->DOF_ALL;++i)
     rhs->Entries[i]=0.0;
  AssembleMatrix(matrix);
  AssembleRHST(rhs);
  //treat the Dirichlet boundary condition
  BoundaryFunction3D *boundaryfunction;
  boundaryfunction = BoundFuncMul;
  //处理边界条件
  BoundaryTreatment(matrix,rhs->Entries, boundaryfunction);
  DIRECTSOLVER(matrix,rhs->Entries,sol);
   for(i=0;i<dof_all;++i)
     auxfefun[0]->Values[i]=sol[i];
  if(ErrorFefunction(tmp_fefun,auxfefun[0],1,MuiltiIndex)<=1e-4)
    break;
  else
  { 
     for(i=0;i<dof_all;++i)
     tmp_fefun->Values[i]=sol[i];
  }
}//---------------------------------------自洽场迭代结束----------------------------------------
  FreeMatrix(matrix);
}

/*

void NonScfIterationWithMatrix(MATRIX *matrix,RHST *rhs,double *sol,double tol,int maxinum)
{
  int i,iter;
  //treat the Dirichlet boundary condition
  BoundaryFunction3D *boundaryfunction;
  boundaryfunction = BoundFuncMul;
  AssembleRHST(rhs);
  //处理边界条件
  BoundaryTreatment(matrix,rhs->Entries, boundaryfunction);
  //调用求解的程序
  CG(matrix,rhs->Entries,sol,1e-3,100);
  int dof_all;
  dof_all=matrix->N_Rows;
  Fefunction3D **auxfefun;
  auxfefun=rhs->AuxFeFun;
  for(i=0;i<dof_all;++i)
    auxfefun[0]->Values[i]=sol[i];
  Fefunction3D *tmp_fefun;
  double *tmp_value;
  tmp_value=malloc(sizeof(double)*dof_all);
  Fespace3D *fesp3D;
  fesp3D = matrix->DiscreteForm3D->Test_Space;
  for(i=0;i<dof_all;++i)
     tmp_value[i]=sol[i];
  tmp_fefun = BuildFefunction3D(fesp3D,tmp_value);
  MultiIndex3D MuiltiIndex[1] = {D000};
  //--------------------------------------开始自洽场迭代-----------------------------------------------  
  for(iter=0; iter<100; ++iter)
  {
  printf ( "--------------------------begin the %d-th iteration nonscfiter----------------------\n",iter );
  for(i=0;i<rhs->DOF_ALL;++i)
     rhs->Entries[i]=0.0;
  AssembleRHST(rhs);
  //treat the Dirichlet boundary condition
  BoundaryFunction3D *boundaryfunction;
  boundaryfunction = BoundFuncMul;
  //处理边界条件
  BoundaryTreatment(matrix,rhs->Entries, boundaryfunction);
  CG(matrix,rhs->Entries,sol,1e-3,100);
   for(i=0;i<dof_all;++i)
     auxfefun[0]->Values[i]=sol[i];
  if(ErrorFefunction(tmp_fefun,auxfefun[0],1,MuiltiIndex)<=1e-3)
    break;
  else
  { 
     for(i=0;i<dof_all;++i)
     tmp_fefun->Values[i]=sol[i];
  }
}//---------------------------------------自洽场迭代结束----------------------------------------
}









void NonScfCorrection(MATRIX *stiff_matrix, MATRIX *matrix,DISCRETEFORM3D *discreteformH, RHST *rhs, double *result, Fefunction3D *fefuncbase, double tol, int maxinum)
{
  printf ( "begin the coarse mesh correction\n"); 
  int i,iter;
  MATRIX *prolong, *Restrict;
  prolong = malloc(sizeof(MATRIX));
  Restrict = malloc(sizeof(MATRIX));
  BuildProlongAn(discreteformH->Anasatz_Space, fefuncbase->Fespace, prolong);
  TransposeMatrix(prolong, Restrict);
  MATRIX *StiffMatrix;
  StiffMatrix = ExtensionMatrix(matrix,1);
  StiffMatrix->DiscreteForm3D = discreteformH;
  Fefunction3D *fefunc_finer;
  int dof_all = fefuncbase->DOF_ALL;
  int dof_all_H = matrix->N_Rows+1;
  double *val, *Rhs;
  Rhs = malloc(sizeof(double)*dof_all_H);
  //val = malloc(sizeof(double)*dof_all);
  //for(i=0;i<dof_all;++i)
  // val[i] = fefuncbase->Values[i];
  //fefunc_finer = BuildFefunction3D(fefuncbase->Fespace, val);
  fefunc_finer = rhs->AuxFeFun[0];
  for(i=0;i<dof_all;++i)
    fefunc_finer->Values[i] = fefuncbase->Values[i];
  double *bb;
  bb = malloc(sizeof(double)*(dof_all_H-1));
  MatMatVec(Restrict, stiff_matrix, fefuncbase->Values, bb);
  double *cc;
  cc = malloc(sizeof(double));
  VecMatVec(fefuncbase->Values, stiff_matrix, cc);
  for(i=0;i<dof_all_H-1;++i)
  {
   StiffMatrix->Entries[StiffMatrix->RowPtr[i+1]-1]=bb[i];
  }
  for(i=0;i<dof_all_H-1;++i)
  {
   StiffMatrix->Entries[StiffMatrix->RowPtr[dof_all_H-1]+i]=bb[i];
  }
   StiffMatrix->Entries[StiffMatrix->RowPtr[dof_all_H]-1]=cc[0];

  //rhs->AuxFeFun[0] = fefunc_finer;
  AssembleRHST(rhs);
  MatrixDotVec(Restrict, rhs->Entries, Rhs);
  VecDotVecs(rhs->Entries,fefuncbase->Values,dof_all,cc);
  Rhs[dof_all_H-1] = cc[0];
  BoundaryFunction3D *boundaryfunction;
  boundaryfunction = BoundFuncMul;
  //处理边界条件
  //OutPutMatrix(StiffMatrix);
  //for(i=0;i<rhs->DOF_ALL;++i)
  //   printf ( "%d, %f\n", i, rhs->Entries[i] );
  
  BoundaryTreatment(StiffMatrix, Rhs, boundaryfunction);
  //调用求解的程序
  double *sol;
  sol = malloc(sizeof(double)*dof_all_H);
  for(i=0;i<dof_all_H-1;++i)
     sol[i] = 0;
  sol[dof_all_H-1] = 1;
  CG(StiffMatrix,Rhs,sol, 1e-3, 100);
  MatrixDotVec(prolong,sol, rhs->AuxFeFun[0]->Values);
  for(i=0;i<dof_all;++i)
    rhs->AuxFeFun[0]->Values[i]+=sol[dof_all_H-1]*fefuncbase->Values[i]; 
  Fefunction3D *tmp_fefun;
  double *tmp_value;
  tmp_value=malloc(sizeof(double)*dof_all);
  tmp_fefun = BuildFefunction3D(fefuncbase->Fespace,tmp_value);
  MultiIndex3D MuiltiIndex[1] = {D000};
  //--------------------------------------开始自洽场迭代-----------------------------------------------  
  for(iter=0; iter<maxinum; ++iter)
  {
    printf ( "--------------------------begin the %d-th iteration coarse correction----------------------\n",iter );
    for(i=0;i<dof_all;++i)
      tmp_value[i] = rhs->AuxFeFun[0]->Values[i];  
    //fefunc_finer = BuildFefunction3D(fefuncbase->Fespace, val);
    for(i=0;i<dof_all;++i)
     rhs->Entries[i] = 0;
    AssembleRHST(rhs);
    MatrixDotVec(Restrict, rhs->Entries, Rhs);
    VecDotVecs(rhs->Entries,fefuncbase->Values,dof_all,cc);
    Rhs[dof_all_H-1] = cc[0];
    BoundaryFunction3D *boundaryfunction;
    boundaryfunction = BoundFuncMul;
    //处理边界条件
    BoundaryTreatment(StiffMatrix,Rhs, boundaryfunction);
    CG(StiffMatrix,Rhs,sol, 1e-3, 100);
    MatrixDotVec(prolong,sol, rhs->AuxFeFun[0]->Values);
    for(i=0;i<dof_all;++i)
      rhs->AuxFeFun[0]->Values[i]+=sol[dof_all_H-1]*fefuncbase->Values[i]; 
    if(ErrorFefunction(tmp_fefun,fefunc_finer,1,MuiltiIndex)<=1e-5)
      break;
  }//---------------------------------------自洽场迭代结束----------------------------------------
  for(i=0;i<dof_all;++i)
    result[i]=rhs->AuxFeFun[0]->Values[i]; 
  printf ( "end the coarse mesh correction\n"); 
  //内存释放
  free(Rhs);
  //free(val);
  free(sol);
  free(tmp_value);
  free(bb);

}



*/


















