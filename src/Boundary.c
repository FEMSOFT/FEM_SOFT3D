/*
 * =====================================================================================
 *
 *       Filename:  Matrix.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2014年10月20日 17时03分08秒
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
#include <string.h>
#include <math.h>
#include "Matrix.h"
#include "Fefunction3D.h"

/** 处理Dirichlet边界条件 */
void BoundaryTreatment(MATRIX *Stiff_Matrix, double *RHS, BoundaryFunction3D *boundaryfunction)
{
  printf("In the boundary condition!\n");
  int i, j, k, ID_bd, start, end, ID_DOF, vert_dof, line_dof,face_dof, DOF_ALL;
  Fespace3D *anasatz_space, *test_space;
  MESH *mesh;
  mesh = Stiff_Matrix->DiscreteForm3D->Anasatz_Space->Mesh;
  anasatz_space = Stiff_Matrix->DiscreteForm3D->Anasatz_Space;
  DOF_ALL = anasatz_space->DOF_ALL;
  BASEFUNCTION3D *Base;
  Base = anasatz_space->Base;

  BoundCondFunct3D *BoundCond;
  BoundCond = anasatz_space->BoundCond;
  BoundType bdtype;
  
  int Num_Verts, Num_Lines, Num_Faces;  
  VERT **Verts, *vert, *vert0, *vert1;
  LINE **Lines, *line;
  FACEE **Faces, *face;
  VOLU **Volus, *volu;
  Num_Verts = mesh->Num_Verts_Global;
  Num_Lines = mesh->Num_Lines_Global;
  Num_Faces = mesh->Num_Faces_Global;
  
  Verts = mesh->Verts;
  Lines = mesh->Lines;
  Faces = mesh->Faces;
  Volus = mesh->Volus; 
  int *RowPtr_Stiff, *KCol_Stiff;
  double *Entries_Stiff;
  RowPtr_Stiff = Stiff_Matrix->RowPtr;
  KCol_Stiff = Stiff_Matrix->KCol;
  Entries_Stiff = Stiff_Matrix->Entries;  
  
  /** 处理节点上的边界条件 */
  vert_dof = Base->DOF[0];
  line_dof = Base->DOF[1];
  face_dof = Base->DOF[2];
  int dof_base, dof_base_ref;
  double eta, coordx, coordy;
  dof_base = 0;
  dof_base_ref = 0; 
  
  int num_volus,volu_ind; 
  num_volus = mesh->Num_Volus_Global;
  
  int *BeginIndex_anasatz,*GlobalNumber_anasatz;
  BeginIndex_anasatz = anasatz_space->BeginIndex;
  GlobalNumber_anasatz = anasatz_space->GlobalNumbers;
  int curr_index_test,curr_index_anasatz,col_num, curr_entry;
  int *local_dof;
  NodalFunct3D *NodalFun;
  NodalFun = Base->NodalF3D;
  Functions3D *bdvaluefun;
 
  ELEMENT3D *elem3D;
  elem3D = malloc(1*sizeof(ELEMENT3D));
  // initializing the element
  InitialElem3D(mesh,elem3D);

  double *bdvalue;
  int dim_basis,elem_faces, elem_verts, elem_lines;
  dim_basis = Base->Num_Bas;
  bdvalue = malloc(dim_basis*sizeof(double));
  
  int v4f,v1,v2,v3,l1,l2,l3;
  printf ( "----------------------------------------------------\n" );
  for(volu_ind=0;volu_ind<num_volus;volu_ind++)
  {
    //printf ( "%d\n",volu_ind );
    dof_base_ref = 0;
    volu = Volus[volu_ind];
  
    elem_lines = 6;
    elem_verts = 4;
    elem_faces = 4; 
    GetElement3D(mesh,volu_ind,elem3D);
    //printf ( "aaaaaaaaaaaaaaaaaaaaaaaaaaa\n" );
    for(v4f=0;v4f<4;v4f++)
    {//printf ( "face[%d]\n",v4f );
      face = Faces[volu->Faces[v4f]];
      if(face->ID_Boundary<=0)
      continue;      
      BoundCond(face->ID_Boundary,0,&bdtype);
        if(bdtype==DIRICHLETT)
        {
	  v1= face->Verts[0];
	  v2= face->Verts[1];
	  v3= face->Verts[2];
          l1= face->Lines[0];
	  l2= face->Lines[1];
	  l3= face->Lines[2];
	  //首先处理面的三个顶点的自由度
	  for(i=0;i<elem_verts;i++)
	  {//printf ( "@@@@@@@@@@@@@2\n" );
	    // printf ( "%d\n", volu->Verts[i]);
	    // printf ( "%d,%d,%d\n",v1,v2,v3 );
	   if((volu->Verts[i]!=v1)&&(volu->Verts[i]!=v2)&&(volu->Verts[i]!=v3))
	      continue;
	   //printf ( "@@@@@@@@@@@@@@\n" );
	   vert = Verts[volu->Verts[i]];
	   //此时 local_dof 是一个地址
	   local_dof = GlobalNumber_anasatz+BeginIndex_anasatz[volu_ind];
	    for(k=0;k<vert_dof;k++)
	    {
	     ID_DOF = local_dof[k*elem_verts+i];
	     //处理刚度矩阵的相应的行
	     start = RowPtr_Stiff[ID_DOF];
	     end = RowPtr_Stiff[ID_DOF+1];
	      for(j=start;j<end;j++)
	      {
	       Entries_Stiff[j] = 0.0;
	       if(KCol_Stiff[j]==ID_DOF)
	        Entries_Stiff[j] = 1.0;	      		
	      }//end for j	  
	      //右端项该的处理
	     bdvaluefun = boundaryfunction(face->ID_Boundary);
	     NodalFun(elem3D, bdvaluefun, bdvalue);
	     RHS[ID_DOF] = bdvalue[k*elem_verts+i];
           }//end for k	    
          }//end for i      
	  //printf ( "bbbbbbbbbbbbbbbbbbbbbbbbb\n" );
          //第二步处理边上的自由度     
         if(line_dof>0)
	 {
	  for(i=0;i<elem_lines;i++)
          {
           if(volu->Lines[i]!=l1&&volu->Lines[i]!=l2&&volu->Lines[i]!=l3)
	     continue;
	   line = Lines[volu->Lines[i]];
	   local_dof = GlobalNumber_anasatz+BeginIndex_anasatz[volu_ind]+elem_verts*vert_dof;
	    for(k=0;k<line_dof;k++)
	    {
	      ID_DOF = local_dof[k*elem_lines+i];
	      //处理刚度矩阵的相应的行
	      start = RowPtr_Stiff[ID_DOF];
	      end = RowPtr_Stiff[ID_DOF+1];
	      for(j=start;j<end;j++)
	      {
		Entries_Stiff[j] = 0.0;
		if(KCol_Stiff[j]==ID_DOF)
		  Entries_Stiff[j] = 1.0;			
	      }//end for j
	      //右端项该的处理
	      bdvaluefun = boundaryfunction(face->ID_Boundary);//,bdvaluefun);
	      NodalFun(elem3D,bdvaluefun, bdvalue);
	      RHS[ID_DOF] = bdvalue[elem_verts*vert_dof+k*elem_lines+i];	      
	    }//end for k	    
           }//end for i	        
	 } 
	     //最后处理面上的自由度
         if(face_dof>0)
         {
	    printf ( "if face dof >0\n" );
	  local_dof = GlobalNumber_anasatz+BeginIndex_anasatz[volu_ind]+elem_verts*vert_dof
	               +elem_lines*line_dof;
	  for(k=0;k<face_dof;k++)
	  {
	   ID_DOF = local_dof[k*elem_faces+v4f];
	   //处理刚度矩阵的相应的行
	   start = RowPtr_Stiff[ID_DOF];
	   end = RowPtr_Stiff[ID_DOF+1];
	   for(j=start;j<end;j++)
	   {
	    Entries_Stiff[j] = 0.0;
	    if(KCol_Stiff[j]==ID_DOF)
	    Entries_Stiff[j] = 1.0;			
	   }//end for j
	   //右端项该处理
	   bdvaluefun = boundaryfunction(face->ID_Boundary);//,bdvaluefun);
	   NodalFun(elem3D,bdvaluefun, bdvalue);
	   RHS[ID_DOF] = bdvalue[elem_verts*vert_dof+elem_lines*line_dof+k*elem_faces+v4f];	      
	  }//end for k	    
         }//end for if face_dof>0    

	 }//end for bdtype==DIRICHLETT
    }//end for v4f 
    
  }//单元循环结束
  //释放空间
  FreeElem3D(elem3D);
  free(bdvalue);
}

/** 处理特征值问题边界条件 */
void BoundaryTreatmentEigen(MATRIX *Stiff_Matrix, MATRIX *Mass_Matrix, BoundaryFunction3D *boundaryfunction)
{
  printf("In the boundary condition!\n");
  int i, j, k, ID_bd, start, end, ID_DOF, vert_dof, line_dof,face_dof, DOF_ALL;
  Fespace3D *anasatz_space, *test_space;
  MESH *mesh;
  mesh = Stiff_Matrix->DiscreteForm3D->Anasatz_Space->Mesh;
  CheckMesh(mesh);
  anasatz_space = Stiff_Matrix->DiscreteForm3D->Anasatz_Space;
  DOF_ALL = anasatz_space->DOF_ALL;
  BASEFUNCTION3D *Base;
  Base = anasatz_space->Base;

  BoundCondFunct3D *BoundCond;
  BoundCond = anasatz_space->BoundCond;
  BoundType bdtype;
  
  int Num_Verts, Num_Lines, Num_Faces;  
  VERT **Verts, *vert, *vert0, *vert1;
  LINE **Lines, *line;
  FACEE **Faces, *face;
  VOLU **Volus, *volu;

  Num_Verts = mesh->Num_Verts_Global;
  Num_Lines = mesh->Num_Lines_Global;
  Num_Faces = mesh->Num_Faces_Global;
  
  Verts = mesh->Verts;
  Lines = mesh->Lines;
  Faces = mesh->Faces;
  Volus = mesh->Volus; 
  int *RowPtr_Stiff, *KCol_Stiff,*RowPtr_Mass, *KCol_Mass;
  double *Entries_Stiff,*Entries_Mass;
  RowPtr_Stiff  = Stiff_Matrix->RowPtr;
  RowPtr_Mass   = Mass_Matrix->RowPtr;
  KCol_Stiff    = Stiff_Matrix->KCol;
  KCol_Mass    = Mass_Matrix->KCol;
  Entries_Stiff = Stiff_Matrix->Entries;  
  Entries_Mass = Mass_Matrix->Entries;  
  
  /** 处理节点上的边界条件 */
  vert_dof = Base->DOF[0];
  line_dof = Base->DOF[1];
  face_dof = Base->DOF[2];
  int dof_base, dof_base_ref;
  double eta, coordx, coordy;
  dof_base = 0;
  dof_base_ref = 0; 
  
  int num_volus,volu_ind; 
  num_volus = mesh->Num_Volus_Global;
  
  int *BeginIndex_anasatz,*GlobalNumber_anasatz;
  BeginIndex_anasatz = anasatz_space->BeginIndex;
  GlobalNumber_anasatz = anasatz_space->GlobalNumbers;
  int curr_index_test,curr_index_anasatz,col_num, curr_entry;
  int *local_dof;
  NodalFunct3D *NodalFun;
NodalFun = Base->NodalF3D;
  Functions3D *bdvaluefun;
  
  ELEMENT3D *elem3D;
  elem3D = malloc(1*sizeof(ELEMENT3D));
  // initializing the element
  InitialElem3D(mesh,elem3D);
  double *bdvalue;
  int dim_basis,elem_faces, elem_verts, elem_lines;
  dim_basis = Base->Num_Bas;
  bdvalue = malloc(dim_basis*sizeof(double));
 
  int v4f,v1,v2,v3,l1,l2,l3;
  //printf ( "############---%d----#####\n",num_volus );
  //for(i=0;i<Num_Faces;++i)
  //   printf ( "%d\n",Faces[i]->ID_Boundary );
  
  for(volu_ind=0;volu_ind<num_volus;volu_ind++)
  {
  //   printf ( "the----%d-th------\n",volu_ind );
    dof_base_ref = 0;
    volu = Volus[volu_ind];
    elem_lines = 6;
    elem_verts = 4;
    elem_faces = 4;  
    GetElement3D(mesh,volu_ind,elem3D);
 //   printf ( "ffffffffff\n" );
    for(v4f=0;v4f<4;v4f++)
    {//printf ( "face[%d]\n",v4f );
 //      printf ( "aaaaaaaaaaaaacccccccccc\n" );
       face = Faces[volu->Faces[v4f]];
 //      printf ( "bianhao===%d\n",face->ID_Boundary );
      if(face->ID_Boundary<=0)
      continue;      
      BoundCond(face->ID_Boundary,0,&bdtype);
 //     printf ( "2222222222\n" ); 
       if(bdtype==DIRICHLETT)
        {
	  v1= face->Verts[0];
	  v2= face->Verts[1];
	  v3= face->Verts[2];
          l1= face->Lines[0];
	  l2= face->Lines[1];
	  l3= face->Lines[2];
	  //首先处理面的三个顶点的自由度
	  for(i=0;i<elem_verts;i++)
	  {//printf ( "@@@@@@@@@@@@@2\n" );
	    // printf ( "%d\n", volu->Verts[i]);
	    // printf ( "%d,%d,%d\n",v1,v2,v3 );
	   if((volu->Verts[i]!=v1)&&(volu->Verts[i]!=v2)&&(volu->Verts[i]!=v3))
	      continue;
	   //printf ( "@@@@@@@@@@@@@@\n" );
	   vert = Verts[volu->Verts[i]];
	   //此时 local_dof 是一个地址
	   local_dof = GlobalNumber_anasatz+BeginIndex_anasatz[volu_ind];
	    for(k=0;k<vert_dof;k++)
	    {
	     ID_DOF = local_dof[k*elem_verts+i];
	     //处理刚度矩阵的相应的行
	     bdvaluefun = boundaryfunction(face->ID_Boundary);
	     NodalFun(elem3D, bdvaluefun, bdvalue);

	     start = RowPtr_Stiff[ID_DOF];
	     end = RowPtr_Stiff[ID_DOF+1];
	      for(j=start;j<end;j++)
	      {
	       Entries_Stiff[j] = 0.0;
	       Entries_Mass[j]=0.0;
	       if(KCol_Stiff[j]==ID_DOF)
	        Entries_Stiff[j] = 1.0;	
	       if(KCol_Mass[j]==ID_DOF)
	        Entries_Mass[j] = bdvalue[k*elem_verts+i];	
	      }//end for j	  
           }//end for k	    
          }//end for i      
     
          //第二步处理边上的自由度     
         if(line_dof>0)
	 {
	  for(i=0;i<elem_lines;i++)
          {
           if(volu->Lines[i]!=l1&&volu->Lines[i]!=l2&&volu->Lines[i]!=l3)
	     continue;
	   line = Lines[volu->Lines[i]];
	   local_dof = GlobalNumber_anasatz+BeginIndex_anasatz[volu_ind]+elem_verts*vert_dof;
	    for(k=0;k<line_dof;k++)
	    {
	      ID_DOF = local_dof[k*elem_lines+i];
	      //处理刚度矩阵的相应的行
	      bdvaluefun = boundaryfunction(face->ID_Boundary);//,bdvaluefun);
	      NodalFun(elem3D,bdvaluefun, bdvalue);

	      start = RowPtr_Stiff[ID_DOF];
	      end = RowPtr_Stiff[ID_DOF+1];
	      for(j=start;j<end;j++)
	      {
		Entries_Stiff[j] = 0.0;
	        Entries_Mass[j]=0.0;
		if(KCol_Stiff[j]==ID_DOF)
		  Entries_Stiff[j] = 1.0;			
		if(KCol_Mass[j]==ID_DOF)
		  Entries_Mass[j] = bdvalue[elem_verts*vert_dof+k*elem_lines+i];	    		
	      }//end for j
	    }//end for k	    
           }//end for i	        
	 } 
	     //最后处理面上的自由度
         if(face_dof>0)
         {
	    printf ( "if face dof >0\n" );
	  local_dof = GlobalNumber_anasatz+BeginIndex_anasatz[volu_ind]+elem_verts*vert_dof
	               +elem_lines*line_dof;
	  for(k=0;k<face_dof;k++)
	  {
	   ID_DOF = local_dof[k*elem_faces+v4f];
	   //处理刚度矩阵的相应的行
	   bdvaluefun = boundaryfunction(face->ID_Boundary);//,bdvaluefun);
	   NodalFun(elem3D,bdvaluefun, bdvalue);

	   start = RowPtr_Stiff[ID_DOF];
	   end = RowPtr_Stiff[ID_DOF+1];
	   for(j=start;j<end;j++)
	   {
	    Entries_Stiff[j] = 0.0;
	    Entries_Mass[j] = 0.0;
	    if(KCol_Stiff[j]==ID_DOF)
	    Entries_Stiff[j] = 1.0;			
	    if(KCol_Mass[j]==ID_DOF)
	    Entries_Mass[j] = bdvalue[elem_verts*vert_dof+elem_lines*line_dof+k*elem_faces+v4f];	      
	   }//end for j
	  }//end for k	    
         }//end for if face_dof>0    

	 }//end for bdtype==DIRICHLETT
    }//end for v4f 
    
  }//单元循环结束
  //释放空间
  FreeElem3D(elem3D);
  free(bdvalue);
}
/** 处理特征值问题边界条件 */
void BoundaryTreatmentEigenMul(MATRIX *Stiff_Matrix, MATRIX *Mass_Matrix, BoundaryFunction3D *boundaryfunction)
{
  printf("In the boundary condition!\n");
  int i, j, k, ID_bd, start, end, ID_DOF, vert_dof, line_dof,face_dof, DOF_ALL;
  Fespace3D *anasatz_space, *test_space;
  MESH *mesh;
  mesh = Stiff_Matrix->DiscreteForm3D->Anasatz_Space->Mesh;
  CheckMesh(mesh);
  anasatz_space = Stiff_Matrix->DiscreteForm3D->Anasatz_Space;
  DOF_ALL = anasatz_space->DOF_ALL;
  BASEFUNCTION3D *Base;
  Base = anasatz_space->Base;

  BoundCondFunct3D *BoundCond;
  BoundCond = anasatz_space->BoundCond;
  BoundType bdtype;
  
  int Num_Verts, Num_Lines, Num_Faces;  
  VERT **Verts, *vert, *vert0, *vert1;
  LINE **Lines, *line;
  FACEE **Faces, *face;
  VOLU **Volus, *volu;

  Num_Verts = mesh->Num_Verts_Global;
  Num_Lines = mesh->Num_Lines_Global;
  Num_Faces = mesh->Num_Faces_Global;
  
  Verts = mesh->Verts;
  Lines = mesh->Lines;
  Faces = mesh->Faces;
  Volus = mesh->Volus; 
  int *RowPtr_Stiff, *KCol_Stiff,*RowPtr_Mass, *KCol_Mass;
  double *Entries_Stiff,*Entries_Mass;
  RowPtr_Stiff  = Stiff_Matrix->RowPtr;
  RowPtr_Mass   = Mass_Matrix->RowPtr;
  KCol_Stiff    = Stiff_Matrix->KCol;
  KCol_Mass    = Mass_Matrix->KCol;
  Entries_Stiff = Stiff_Matrix->Entries;  
  Entries_Mass = Mass_Matrix->Entries;  
  
  /** 处理节点上的边界条件 */
  vert_dof = Base->DOF[0];
  line_dof = Base->DOF[1];
  face_dof = Base->DOF[2];
  int dof_base, dof_base_ref;
  double eta, coordx, coordy;
  dof_base = 0;
  dof_base_ref = 0; 
  
  int num_volus,volu_ind; 
  num_volus = mesh->Num_Volus_Global;
  
  int *BeginIndex_anasatz,*GlobalNumber_anasatz;
  BeginIndex_anasatz = anasatz_space->BeginIndex;
  GlobalNumber_anasatz = anasatz_space->GlobalNumbers;
  int curr_index_test,curr_index_anasatz,col_num, curr_entry;
  int *local_dof;
  NodalFunct3D *NodalFun;
NodalFun = Base->NodalF3D;
  Functions3D *bdvaluefun;
  
  ELEMENT3D *elem3D;
  elem3D = malloc(1*sizeof(ELEMENT3D));
  // initializing the element
  InitialElem3D(mesh,elem3D);
  double *bdvalue;
  int dim_basis,elem_faces, elem_verts, elem_lines;
  dim_basis = Base->Num_Bas;
  bdvalue = malloc(dim_basis*sizeof(double));
 
  int v4f,v1,v2,v3,l1,l2,l3;
  //printf ( "############---%d----#####\n",num_volus );
  
  for(volu_ind=0;volu_ind<num_volus;volu_ind++)
  {
  //   printf ( "the----%d-th------\n",volu_ind );
    dof_base_ref = 0;
    volu = Volus[volu_ind];
    elem_lines = 6;
    elem_verts = 4;
    elem_faces = 4;  
    GetElement3D(mesh,volu_ind,elem3D);
  //  printf ( "ffffffffff\n" );
    for(v4f=0;v4f<4;v4f++)
    {//printf ( "face[%d]\n",v4f );
  //     printf ( "aaaaaaaaaaaaacccccccccc\n" );
       face = Faces[volu->Faces[v4f]];
  //     printf ( "bianhao===%d\n",face->ID_Boundary );
      if(face->ID_Boundary<=0)
      continue;      
      BoundCond(face->ID_Boundary,0,&bdtype);
  //    printf ( "2222222222\n" ); 
       if(bdtype==DIRICHLETT)
        {
	  v1= face->Verts[0];
	  v2= face->Verts[1];
	  v3= face->Verts[2];
          l1= face->Lines[0];
	  l2= face->Lines[1];
	  l3= face->Lines[2];
	  //首先处理面的三个顶点的自由度
	  for(i=0;i<elem_verts;i++)
	  {//printf ( "@@@@@@@@@@@@@2\n" );
	    // printf ( "%d\n", volu->Verts[i]);
	    // printf ( "%d,%d,%d\n",v1,v2,v3 );
	   if((volu->Verts[i]!=v1)&&(volu->Verts[i]!=v2)&&(volu->Verts[i]!=v3))
	      continue;
	   //printf ( "@@@@@@@@@@@@@@\n" );
	   vert = Verts[volu->Verts[i]];
	   //此时 local_dof 是一个地址
	   local_dof = GlobalNumber_anasatz+BeginIndex_anasatz[volu_ind];
	    for(k=0;k<vert_dof;k++)
	    {
	     ID_DOF = local_dof[k*elem_verts+i];
	     //处理刚度矩阵的相应的行
	     bdvaluefun = boundaryfunction(face->ID_Boundary);
	     NodalFun(elem3D, bdvaluefun, bdvalue);

	     start = RowPtr_Stiff[ID_DOF];
	     end = RowPtr_Stiff[ID_DOF+1];
	      for(j=start;j<end-1;j++)
	      {
	       Entries_Stiff[j] = 0.0;
	       Entries_Mass[j]=0.0;
	       if(KCol_Stiff[j]==ID_DOF)
	        Entries_Stiff[j] = 1.0;	
	       if(KCol_Mass[j]==ID_DOF)
	        Entries_Mass[j] = bdvalue[k*elem_verts+i];	
	      }//end for j	  
           }//end for k	    
          }//end for i      
     
          //第二步处理边上的自由度     
         if(line_dof>0)
	 {
	  for(i=0;i<elem_lines;i++)
          {
           if(volu->Lines[i]!=l1&&volu->Lines[i]!=l2&&volu->Lines[i]!=l3)
	     continue;
	   line = Lines[volu->Lines[i]];
	   local_dof = GlobalNumber_anasatz+BeginIndex_anasatz[volu_ind]+elem_verts*vert_dof;
	    for(k=0;k<line_dof;k++)
	    {
	      ID_DOF = local_dof[k*elem_lines+i];
	      //处理刚度矩阵的相应的行
	      bdvaluefun = boundaryfunction(face->ID_Boundary);//,bdvaluefun);
	      NodalFun(elem3D,bdvaluefun, bdvalue);

	      start = RowPtr_Stiff[ID_DOF];
	      end = RowPtr_Stiff[ID_DOF+1];
	      for(j=start;j<end-1;j++)
	      {
		Entries_Stiff[j] = 0.0;
	        Entries_Mass[j]=0.0;
		if(KCol_Stiff[j]==ID_DOF)
		  Entries_Stiff[j] = 1.0;			
		if(KCol_Mass[j]==ID_DOF)
		  Entries_Mass[j] = bdvalue[elem_verts*vert_dof+k*elem_lines+i];	    		
	      }//end for j
	    }//end for k	    
           }//end for i	        
	 } 
	     //最后处理面上的自由度
         if(face_dof>0)
         {
	    printf ( "if face dof >0\n" );
	  local_dof = GlobalNumber_anasatz+BeginIndex_anasatz[volu_ind]+elem_verts*vert_dof
	               +elem_lines*line_dof;
	  for(k=0;k<face_dof;k++)
	  {
	   ID_DOF = local_dof[k*elem_faces+v4f];
	   //处理刚度矩阵的相应的行
	   bdvaluefun = boundaryfunction(face->ID_Boundary);//,bdvaluefun);
	   NodalFun(elem3D,bdvaluefun, bdvalue);

	   start = RowPtr_Stiff[ID_DOF];
	   end = RowPtr_Stiff[ID_DOF+1];
	   for(j=start;j<end-1;j++)
	   {
	    Entries_Stiff[j] = 0.0;
	    Entries_Mass[j] = 0.0;
	    if(KCol_Stiff[j]==ID_DOF)
	    Entries_Stiff[j] = 1.0;			
	    if(KCol_Mass[j]==ID_DOF)
	    Entries_Mass[j] = bdvalue[elem_verts*vert_dof+elem_lines*line_dof+k*elem_faces+v4f];	      
	   }//end for j
	  }//end for k	    
         }//end for if face_dof>0    

	 }//end for bdtype==DIRICHLETT
    }//end for v4f 
    
  }//单元循环结束
  //释放空间
  FreeElem3D(elem3D);
  free(bdvalue);
}


//边界条件函数
void BoundValueMul(double x, double y,double z, int n, double *value)
{
  /** Example 1 */
  /** Example 2*/
  value[0] = 0;
}


Functions3D * BoundFuncMul(int ID_Bound) //,Functions3D *boundfun)
{
  return  BoundValueMul;
}




