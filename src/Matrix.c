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

/** OutPut the matrix */
void OutPutMatrix(MATRIX *matrix)
{
  int i, j, end;
  printf("========================\n");
  printf("Output the matrix:\n\n"); 
  int N_Rows = matrix->N_Rows;
  printf ( "N_Rows=%d\n",matrix->N_Rows);
  for(i=0;i<N_Rows;i++)
  {
    end = matrix->RowPtr[i+1];
    for(j=matrix->RowPtr[i];j<end;j++)
    {
      printf("a(%d,%d)=",i+1,matrix->KCol[j]+1);
      if(fabs(matrix->Entries[j]) < 1e-12)
	matrix->Entries[j] = 0;
      printf("%18.15f\n",matrix->Entries[j]);
    }
    printf("\n");
   }
}

//anasatz space fesp1,test space fesp2
/** 依据试探有限元空间和检验有限元空间产生矩阵的结构 */
MATRIX *BuildMatrix(DISCRETEFORM3D *discreteform){
	//fesp1: anasatz space, fesp2: test space
    Fespace3D *fesp1, *fesp2;
    fesp1 = discreteform->Anasatz_Space;
    fesp2 = discreteform->Test_Space; 
    MATRIX *matrix;
    matrix = malloc(1*sizeof(MATRIX));
    matrix->DiscreteForm3D = discreteform;
    
    matrix->N_Rows = fesp2->DOF_ALL;
    matrix->N_Columns = fesp1->DOF_ALL;
    // then build the matrix structure    
    //here we assume all the FEM spaces are defined on the same mesh
    MESH *mesh;  
    int num_volus,num_faces,num_verts,num_lines;
    int curr_volu, curr_face, curr_line, i, j, num_vert_face, num_line_face;
    BASEFUNCTION3D *anasatz_base,*test_base;
    
    mesh = fesp2->Mesh;
    
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
	N_local_vert = 4;
	N_local_line = 6;
	N_local_face = 4;
    for(curr_volu=0;curr_volu<num_volus;curr_volu++)
    {
        //first give the number of related dof for each vert
        volu = mesh->Volus[curr_volu];
        for(i=0;i<N_local_vert;i++)
            N_volu4vert[volu->Verts[i]]++;
        for(i=0;i<N_local_line;i++)
            N_volu4line[volu->Lines[i]]++;
        for(i=0;i<N_local_face;i++)
            N_volu4face[volu->Faces[i]]++;
    }
     
     N_local_vert = 3;
     N_local_line = 3;
     for(curr_face=0;curr_face<num_faces;curr_face++)
    {
        face = mesh->Faces[curr_face];
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
     RowPtr = malloc((N_Rows+1)*sizeof(int));
     memset(RowPtr,-1,(N_Rows+1)*sizeof(int));
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
			 vert_dof += local_volu*(anasatz_DOF[3]+(volu_faces-3)*anasatz_DOF[2]);
			 local_face = N_face4vert[curr_vert];
			 vert_dof += local_face*(anasatz_DOF[2]+anasatz_DOF[1]);
			 
			 local_line = N_line4vert[curr_vert];
			 vert_dof += local_line*(anasatz_DOF[0]+anasatz_DOF[1]);
			 RowPtr[curr_row] = RowPtr[curr_row-1] + vert_dof;			 
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
			RowPtr[curr_row] = RowPtr[curr_row-1] + line_dof;  		   
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
			 RowPtr[curr_row] = RowPtr[curr_row-1] + face_dof;			 
		 }//end for curr_face		 
	}//end for i
     //for the volu
     for(i=0;i<test_DOF[3];i++)
     {
		 for(curr_volu=0;curr_volu<num_volus;curr_volu++)
		 {
			 curr_row++;
			 volu_dof = anasatz_DOF[3]+anasatz_DOF[2]*4+anasatz_DOF[1]*6+anasatz_DOF[0]*4;
			 RowPtr[curr_row] = RowPtr[curr_row-1] + volu_dof;			 
		}//end for curr_volu       		 
	}// end for i    

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
    free(N_volu4vert); 
    free(N_line4vert); 
    free(N_volu4line); 
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
    //QuickSort_Int(KCol,RowPtr[N_Rows-1],RowPtr[N_Rows]-1);

  return matrix;    
}

/**====================================================================*/
/**合成刚度矩阵 */
void AssembleMatrix(MATRIX *Matrix)
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
		   //先获得所有矩阵和右端项在本积分点的值
		   //每个矩阵包含自己的取值空间和试探空间，后面的微分算子也是针对
		   //针对这两个空间的。另外可能还有其他的辅助函数空间和相应的微分
			//算子，最后也会包含辅助函数来定义矩阵     
			//printf("mat[6]=%f\n",Mat_Entries[6]);      
			SubAssembleMatrix3D(n_bas_test,test_Base,n_bas_anasatz,anasatz_Base,
						n_test_multiindex,test_MultiIndex,n_anasatz_multiindex,anasatz_MultiIndex, 
							N_AuxFefunct,AuxFefunct,N_AuxFefun_MultiIndex,AuxFefun_MultiIndex,
							N_AuxFefun_Values, AuxFefun_Values,N_UserFunction,UserFunction,N_FeConst,FeConst, 
						discreteformvolu,elem3D,
						Test_Values,Anasatz_Values,Test_MultiIndex_Values, Anasatz_MultiIndex_Values,
						Quad_xi,Quad_eta,Quad_zeta,Quad_W,Mat_Entries);
       } //end for i (Quad_Points)       
     } //end for if 
   //   for(i=0;i<16;++i)
   //   printf ( "MAT_ENTRIES=%f\n",Mat_Entries[i] );
       AddSubMatrix(Matrix,Mat_Entries, n_bas_test,Test_GlobalNumbers+Test_BeginIndex[curr_volu],
		 n_bas_anasatz,Anasatz_GlobalNumbers+Anasatz_BeginIndex[curr_volu]);	
  }//end for curr_volu
    
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
/**====================================================================*/
/**====================================================================*/
/** Sub program */
void SubAssembleMatrix3D(int n_bas_test,BASEFUNCTION3D *test_Base,
             int n_bas_anasatz,BASEFUNCTION3D *anasatz_Base,
             int n_test_multiindex,MultiIndex3D *test_MultiIndex,
             int n_anasatz_multiindex,MultiIndex3D *anasatz_MultiIndex,              
             int N_AuxFefunct,Fefunction3D **AuxFefunct,int* N_AuxFefun_MultiIndex,
             MultiIndex3D *AuxFefun_MultiIndex,int N_AuxFefun_Values,double *AuxFefun_Values,
             int N_UserFunction,Functions3D **UserFunction, int N_FeConst,double *FeConst,
	     DiscreteFormMatrix *discreteformvolu, ELEMENT3D *elem3D,double *Test_Values, 
	     double *Anasatz_Values,double *Test_MultiIndex_Values,double *Anasatz_MultiIndex_Values,
	     double Quad_xi,double Quad_eta,double Quad_zeta,double Quad_W, double *Mat_Entries)
{
    double Quad_x,Quad_y,Quad_z;
    //仔细看看
    GetElementCoord(elem3D,Quad_xi,Quad_eta,Quad_zeta, &Quad_x,&Quad_y,&Quad_z);  
    int i, j, k, m,n;
    int dim_anasatz_value, dim_test_value;
    dim_anasatz_value = anasatz_Base->Value_Dim;
    dim_test_value = test_Base->Value_Dim;
    //compute the values of the test bases
    //对于高维的基函数该如何办
    for(i=0;i<n_test_multiindex;i++)
    {
      //测试函数的值存在Test_Values中
      GetBaseValues(test_Base,test_MultiIndex[i],elem3D,Quad_x,Quad_y,Quad_z,Quad_xi,Quad_eta,Quad_zeta,
	            Test_Values+i*dim_test_value*n_bas_test);
    }
    //compute the values of anasatz bases
    //Base = Matrix->Anasatz_Space->Base;
    //multiindex = Matrix->Anasatz_MultiIndex;
    for(i=0;i<n_anasatz_multiindex;i++)
    {
      //取值函数的值存在Anasatz_Values中
      GetBaseValues(anasatz_Base,anasatz_MultiIndex[i],elem3D,Quad_x,Quad_y,Quad_z,Quad_xi,Quad_eta,
	            Quad_zeta,Anasatz_Values+i*dim_anasatz_value*n_bas_anasatz);
    }
    
    //下面处理辅助有限元函数值的计算
    int curr_pos;
    if(N_AuxFefunct)
    {
      curr_pos = 0;
      for(i=0;i<N_AuxFefunct;i++)
      {
	//求一片导数的函数值
	GetFefunctValues(AuxFefunct[i],elem3D,Quad_x,Quad_y,Quad_z,
			 Quad_xi,Quad_eta,Quad_zeta,N_AuxFefun_MultiIndex[i],
		         AuxFefun_MultiIndex+curr_pos,AuxFefun_Values+curr_pos);
	curr_pos += N_AuxFefun_MultiIndex[i]*AuxFefunct[i]->Fespace->Base->Value_Dim;           
      }//end for i
    }
    else
      AuxFefun_Values = NULL;

    int N_UserFunc;
    N_UserFunc=N_UserFunction;
    double *UserFuncValue;
    UserFuncValue=malloc(sizeof(double)*N_UserFunc);
    if(N_UserFunc)
    {
       for(i=0;i<N_UserFunc;++i)
       {
        UserFunction[i](Quad_x,Quad_y,Quad_z,1,UserFuncValue+i);
       }
    }
    else
    UserFuncValue = NULL;
    //for the assemble
    
    
    for(i=0;i<n_bas_test;i++)
    {
		//先把基函数在积分点所需要的导数值都都装备上
		for(k=0;k<n_test_multiindex;k++)
		{
			for(m=0;m<dim_test_value;m++)
			{
				Test_MultiIndex_Values[k*dim_test_value+m] = Test_Values[k*n_bas_test*dim_test_value+i*dim_test_value+m];				
			}			
		}
		// 处理取值有限元空间的值
               for(j=0;j<n_bas_anasatz;j++)
               {
                 //先把基函数在积分点所需要的导数值都都装备上
                   for(k=0;k<n_anasatz_multiindex;k++)
			{
				for(n=0;n<dim_anasatz_value;n++)
				{
					Anasatz_MultiIndex_Values[k*dim_anasatz_value+n] = Anasatz_Values[k*n_bas_anasatz*dim_anasatz_value
                                                                   +j*dim_anasatz_value+n];					
				}				
			}
			//调用用户自己定义的离散变分函数
			Mat_Entries[i*n_bas_test+j] += discreteformvolu(Quad_x,Quad_y,Quad_z,n_test_multiindex,
									Test_MultiIndex_Values,n_anasatz_multiindex, Anasatz_MultiIndex_Values,
								        N_AuxFefun_Values, AuxFefun_Values,N_UserFunc,UserFuncValue,
								        N_FeConst,FeConst)*elem3D->Volu*Quad_W;
												   //getchar();			
              }//end for j		
    }// end for i 
}//end program
//assemble the sub matrix on the face
/*void SubAssembleMatrix3DFace(int ID_bd,double *normal, double * tangent,double local_direction,int n_bas_test,BASEFUNCTION2D *test_Base,
             int n_bas_anasatz,BASEFUNCTION2D *anasatz_Base,
             int n_test_multiindex,MultiIndex2D *test_MultiIndex,
             int n_anasatz_multiindex,MultiIndex2D *anasatz_MultiIndex,              
             int N_AuxFefunct,Fefunction2D **AuxFefunct,int* N_AuxFefun_MultiIndex,
             MultiIndex2D *AuxFefun_MultiIndex,int N_AuxFefun_Values,double *AuxFefun_Values,
             DiscreteFormMatrixLine *discreteformline, ELEMENT2D *elem2D,double *Test_Values, 
	     double *Anasatz_Values,double *Test_MultiIndex_Values,double *Anasatz_MultiIndex_Values,
	     double Quad_xi,double Quad_eta,double Quad_zeta,double Quad_W, double area,
	     double *Mat_Entries)
{
    double Quad_x,Quad_y,Quad_z;
    GetElementCoord(elem3D,Quad_xi,Quad_eta,Quad_zeta, &Quad_x,&Quad_y,&Quad_z);  
    int i, j, k, m,n;
    int dim_anasatz_value, dim_test_value;
    dim_anasatz_value = anasatz_Base->Value_Dim;
    dim_test_value = test_Base->Value_Dim;
    //compute the values of the test bases
    //对于高维的基函数该如何办
    for(i=0;i<n_test_multiindex;i++)
    {
      //测试函数的值存在Test_Values中
      GetBaseValues(test_Base,test_MultiIndex[i],elem3D,Quad_x,Quad_y,Quad_z,Quad_xi,Quad_eta,Quad_zeta,
		    Test_Values+i*dim_test_value*n_bas_test);
    }
    //compute the values of anasatz bases
    //Base = Matrix->Anasatz_Space->Base;
    //multiindex = Matrix->Anasatz_MultiIndex;
    for(i=0;i<n_anasatz_multiindex;i++)
    {
      //取值函数的值存在Anasatz_Values中
      GetBaseValues(anasatz_Base,anasatz_MultiIndex[i],elem3D,Quad_x,Quad_y,Quad_z,Quad_xi,Quad_eta,
	            Quad_zeta,Anasatz_Values+i*dim_anasatz_value*n_bas_anasatz);

    }
    //下面处理辅助有限元函数值的计算
    int curr_pos;
    if(N_AuxFefunct)
    {
        //N_AuxFefunct = Matrix->N_AuxFefunct;
        //AuxFefun_MultiIndex = Matrix->AuxFefun_MultiIndex;
        //N_AuxFefun_MultiIndex = Matrix->N_AuxFefun_MultiIndex;
        curr_pos = 0;
        for(i=0;i<N_AuxFefunct;i++)
        {
            //AuxFefunct = Matrix->AuxFefunct[i];
            //求一片导数的函数值
            GetFefunctValues(AuxFefunct[i],elem3D,Quad_x,Quad_y,Quad_z,
                             Quad_xi,Quad_eta,Quad_zeta,N_AuxFefun_MultiIndex[i],
                             AuxFefun_MultiIndex+curr_pos,AuxFefun_Values+curr_pos);
            curr_pos += N_AuxFefun_MultiIndex[i]*AuxFefunct[i]->Fespace->Base->Value_Dim;           
        }//end for i
    }
    else
    {
      AuxFefun_Values = NULL;
    }
   
    //for the assemble
    for(i=0;i<n_bas_test;i++)
    {
      // get the values at each element of the test and anasatz with all the multiindices
      //先把基函数在积分点所需要的导数值都都装备上
      for(k=0;k<n_test_multiindex;k++)
      {
	for(m=0;m<dim_test_value;m++)
	{
	  Test_MultiIndex_Values[k*dim_test_value+m] = Test_Values[k*n_bas_test*dim_test_value+i*dim_test_value+m];
	}
      }
      // 处理取值有限元空间的值
        for(j=0;j<n_bas_anasatz;j++)
        {
            // get the values at each element of the test and anasatz with all the multiindices
            //先把基函数在积分点所需要的导数值都都装备上
            //for(k=0;k<n_test_multiindex;k++)
            //    Test_MultiIndex_Values[k] = Test_Values[k*n_bas_test+i];
            for(k=0;k<n_anasatz_multiindex;k++)
	    {
	      for(n=0;n<dim_anasatz_value;n++)
	      {
		Anasatz_MultiIndex_Values[k*dim_anasatz_value+n] = Anasatz_Values[k*n_bas_anasatz*dim_anasatz_value+j*dim_anasatz_value+n];
	      }
	    }
	    //调用用户自己定义的离散变分函数
            Mat_Entries[i*n_bas_test+j] += discreteformline(ID_bd,Quad_x,Quad_y,normal,tangent,local_direction,n_test_multiindex,Test_MultiIndex_Values,n_anasatz_multiindex,
                         Anasatz_MultiIndex_Values,N_AuxFefun_Values,AuxFefun_Values)*area*Quad_W;             
       // printf("Mat_Entries[%d]=%f\n",i*n_bas_test+j,Mat_Entries[i*n_bas_test+j]);
        }//end for j
    }// end for i   
}//end program


//assemble the sub matrix on the line: for DG method
void SubAssembleMatrix2DLine(int ID_bd,double *normal, double * tangent,double local_direction,int n_bas_test,BASEFUNCTION2D *test_Base,
             int n_bas_anasatz,BASEFUNCTION2D *anasatz_Base,
             int n_test_multiindex,MultiIndex2D *test_MultiIndex,
             int n_anasatz_multiindex,MultiIndex2D *anasatz_MultiIndex,              
             int N_AuxFefunct,Fefunction2D **AuxFefunct,int* N_AuxFefun_MultiIndex,
             MultiIndex2D *AuxFefun_MultiIndex,int N_AuxFefun_Values,double *AuxFefun_Values,
             DiscreteFormMatrixLine *discreteformline, ELEMENT2D *elem2D,double *Test_Values, 
	     double *Anasatz_Values,double *Test_MultiIndex_Values,double *Anasatz_MultiIndex_Values,
	     double Quad_xi,double Quad_eta,double Quad_zeta,double Quad_W, double length,
	     double *Mat_Entries)
{
    double Quad_x,Quad_y,Quad_z;
    GetElementCoord(elem3D,Quad_xi,Quad_eta,Quad_zeta, &Quad_x,&Quad_y,&Quad_z);  
    int i, j, k, m,n;
    int dim_anasatz_value, dim_test_value;
    dim_anasatz_value = anasatz_Base->Value_Dim;
    dim_test_value = test_Base->Value_Dim;
    //compute the values of the test bases
    //对于高维的基函数该如何办
    for(i=0;i<n_test_multiindex;i++)
    {
      //测试函数的值存在Test_Values中
      GetBaseValues(test_Base,test_MultiIndex[i],elem3D,Quad_x,Quad_y,Quad_z,Quad_xi,Quad_eta,Quad_zeta,
		    Test_Values+i*dim_test_value*n_bas_test);
    }
    //compute the values of anasatz bases
    //Base = Matrix->Anasatz_Space->Base;
    //multiindex = Matrix->Anasatz_MultiIndex;
    for(i=0;i<n_anasatz_multiindex;i++)
    {
      //取值函数的值存在Anasatz_Values中
      GetBaseValues(anasatz_Base,anasatz_MultiIndex[i],elem3D,Quad_x,Quad_y,Quad_z,Quad_xi,Quad_eta,
	            Quad_zeta,Anasatz_Values+i*dim_anasatz_value*n_bas_anasatz);
    }
    //下面处理辅助有限元函数值的计算
    int curr_pos;
    if(N_AuxFefunct)
    {
        //N_AuxFefunct = Matrix->N_AuxFefunct;
        //AuxFefun_MultiIndex = Matrix->AuxFefun_MultiIndex;
        //N_AuxFefun_MultiIndex = Matrix->N_AuxFefun_MultiIndex;
        curr_pos = 0;
        for(i=0;i<N_AuxFefunct;i++)
        {
            //AuxFefunct = Matrix->AuxFefunct[i];
            //求一片导数的函数值
            GetFefunctValues(AuxFefunct[i],elem3D,Quad_x,Quad_y,Quad_z,
                             Quad_xi,Quad_eta,Quad_zeta,N_AuxFefun_MultiIndex[i],
                             AuxFefun_MultiIndex+curr_pos,AuxFefun_Values+curr_pos);
            curr_pos += N_AuxFefun_MultiIndex[i]*AuxFefunct[i]->Fespace->Base->Value_Dim;           
        }//end for i
    }
    else
    {
      AuxFefun_Values = NULL;
    }
   
    //for the assemble
    for(i=0;i<n_bas_test;i++)
    {
      // get the values at each element of the test and anasatz with all the multiindices
      //先把基函数在积分点所需要的导数值都都装备上
      for(k=0;k<n_test_multiindex;k++)
      {
	for(m=0;m<dim_test_value;m++)
	{
	  Test_MultiIndex_Values[k*dim_test_value+m] = Test_Values[k*n_bas_test*dim_test_value+i*dim_test_value+m];
	}
      }
      // 处理取值有限元空间的值
        for(j=0;j<n_bas_anasatz;j++)
        {
            // get the values at each element of the test and anasatz with all the multiindices
            //先把基函数在积分点所需要的导数值都都装备上
            //for(k=0;k<n_test_multiindex;k++)
            //    Test_MultiIndex_Values[k] = Test_Values[k*n_bas_test+i];
            for(k=0;k<n_anasatz_multiindex;k++)
	    {
	      for(n=0;n<dim_anasatz_value;n++)
	      {
		Anasatz_MultiIndex_Values[k*dim_anasatz_value+n] = Anasatz_Values[k*n_bas_anasatz*dim_anasatz_value+j*dim_anasatz_value+n];
	      }
	    }
	    //调用用户自己定义的离散变分函数
            Mat_Entries[i*n_bas_test+j] += discreteformline(ID_bd,Quad_x,Quad_y,normal,tangent,local_direction,n_test_multiindex,Test_MultiIndex_Values,n_anasatz_multiindex,
                         Anasatz_MultiIndex_Values,N_AuxFefun_Values,AuxFefun_Values)*length*Quad_W;             
       // printf("Mat_Entries[%d]=%f\n",i*n_bas_test+j,Mat_Entries[i*n_bas_test+j]);
        }//end for j
    }// end for i   
}//end program
*/


/**把单刚矩阵加到总刚矩阵上去*/
void AddSubMatrix(MATRIX *A,double *AA, int n_row,int *row,int n_col,int *col)
{
    int i,j, k;
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
                if(KCol[start+k] == col[j])
		{
                    Entries[start+k] += AA[i*n_row+j];
		    k = length+1;
		}                
            }//end for k
        }//end for j
    }// end for i
}

/****************************************************/
//矩阵加法： A=A+tau*B
void MatrixAdd(MATRIX *A, MATRIX *B, double tau)
{
   int i,j,k, N_Entries;

   N_Entries = A->N_Entries;
   double *A_Entries, *B_Entries;
   A_Entries = A->Entries;
   B_Entries = B->Entries;
   
   for(i=0;i<N_Entries;i++)
   {
     A_Entries[i] += tau*B_Entries[i];
   }  
}






void FreeMatrix(MATRIX *Matrix)
{
  //printf("delete the Matrix!\n");
  free(Matrix->KCol);
  free(Matrix->RowPtr);
  free(Matrix->Entries);
  free(Matrix);
}


void WriteMatrixPS(MATRIX *matrix,const char *filename)
{
    FILE *file;
    file = fopen(filename,"w");
    if(!file)
    {
        printf("\ncannot open the matrix.dat!\n");
	exit(0);
    }
    //printf("open the output file successfully!\n");

    int size = 500;
    printf("nrows = %d\n", matrix->N_Rows);
    printf("ncols = %d\n", matrix->N_Columns);
    double dot1 = 20.0/(double)matrix->N_Rows;
    double dot2 = (double)size/1000.0;
    double dot = dot2>dot1 ? dot2 : dot1;
    int lineWidth = 0.5;
    double step = (double)size/(double)(matrix->N_Rows-1);
    int i, j;
    
    fprintf(file,"%%!PS-Adobe-2.0 EPSF-2.0\n");
    fprintf(file,"%%%%BoundingBox: 0 0 %d %d\n", size+20, size+20);
    fprintf(file,"%%%%EndComments\n");
    fprintf(file,"10 10 translate\n");
    fprintf(file,"newpath\n");
    fprintf(file,"%f %f moveto\n", 0-2*dot, 0-2*dot);
    fprintf(file,"0 %f rlineto\n", size+4*dot);
    fprintf(file,"%f 0 rlineto\n", size+4*dot );
    fprintf(file,"0 %f rlineto\n", -size-4*dot);
    fprintf(file,"closepath\n");
    fprintf(file,"%d setlinewidth\n",lineWidth);
    fprintf(file,"stroke\n\n");
    fprintf(file,"/box\n");
    fprintf(file,"{gsave\n");
    fprintf(file," translate\n");
    fprintf(file," newpath\n");
    fprintf(file," 0 0 %f 0 360 arc\n", dot);
    fprintf(file," closepath\n");
    fprintf(file," 0 setgray fill\n");
    fprintf(file," grestore} def\n\n");

    for ( i=0;i<matrix->N_Rows;++i )
    {
	for ( j=matrix->RowPtr[i];j<matrix->RowPtr[i+1];++j )
	{
	    fprintf(file,"%f %f box\n", matrix->KCol[j]*step, size-i*step);
	}
    }
    fprintf(file,"showpage\n");
    fclose(file);
}



void WriteMatrix(MATRIX *matrix,const char *filename)
{
    FILE *file;
    file = fopen(filename,"w");
    if(!file)
    {
        printf("\ncannot open the matrix!\n");
	exit(0);
    }
    int i;
    fprintf(file,  "%d\n%d\n%d\n", matrix->N_Rows,matrix->N_Columns,matrix->N_Entries);
    for(i=0;i<matrix->N_Rows+1;i++)
    {
      fprintf(file,  "%d\n", matrix->RowPtr[i]);
    }
    for(i=0;i<matrix->N_Entries;i++)
    {
      fprintf(file,  "%d\n", matrix->KCol[i]);
    }
    for(i=0;i<matrix->N_Entries;i++)
    {
      fprintf(file,  "%18.15f\n", matrix->Entries[i]);
    }

    fclose(file);
}

//  for mixed element  matrix
void ReAssembleMatrix(MATRIX ***matrix,int nrows, int ncols, MATRIX *matrixb)
{
    int N_Rows=0, N_Cols=0, N_Entries=0;
    int i,j,k,l,tmp_row,tmp_col;
    int r[nrows],c[ncols];//nrows and ncols of every single small matrix
    //get the nrow and ncols of new matrixb
    for(i=0;i<nrows;++i)
    {
       for(j=0;j<ncols;++j)
       {
          if(matrix[i][j]!=NULL)
          {
             tmp_row = matrix[i][j]->N_Rows;
             N_Entries += matrix[i][j]->N_Entries;
          }   
       } 
       N_Rows += tmp_row;// only need to add one time every row
       r[i] = tmp_row;
    }
    for(j=0;j<ncols;++j)
    {
       for(i=0;i<nrows;++i)
       {
          if(matrix[i][j]!=NULL)
          {
             tmp_col = matrix[i][j]->N_Columns;
             break;
          }   
       } 
       N_Cols += tmp_col;
       c[j] = tmp_col;
    }
    
    matrixb->N_Rows = N_Rows;
    matrixb->N_Columns = N_Cols;
    matrixb->N_Entries = N_Entries;
    matrixb->Entries = malloc(sizeof(double)*N_Entries);
    matrixb->KCol = malloc(sizeof(int)*N_Entries);
    matrixb->RowPtr = malloc(sizeof(int)*(N_Rows+1));
    memset(matrixb->RowPtr,0,sizeof(int)*(N_Rows+1));

    // get the RowPtr of the new matrix
    tmp_row = 0;//number of rows
    int tmp_entries = 0, tmp_Entries = 0;
    for(i=0;i<nrows;++i)
    {
       for(j=0;j<ncols;++j)
       {
          if(matrix[i][j]!=NULL)
          {
             for(k=0;k<r[i];++k)
                matrixb->RowPtr[tmp_row+k+1] += matrix[i][j]->RowPtr[k+1];//+tmp_Entries; 
             tmp_entries += matrix[i][j]->N_Entries;
          }   
       }
       for(k=0;k<r[i];++k)
          matrixb->RowPtr[tmp_row+k+1] += tmp_Entries; //only need to plus one time
   
       tmp_row +=  r[i];
       tmp_Entries = tmp_entries; 
    }

    /*
    for(i=0;i<matrixb->N_Rows+1;++i)
       printf("%d\n",matrixb->RowPtr[i]); 
    exit(0);

    */

    tmp_col = 0;// number of rows
    for(i=0;i<nrows;++i)
    {
       for(k=0;k<r[i];++k)
       {
             tmp_row = 0;
             tmp_entries = 0;
             for(j=0;j<ncols;++j)
             {
                if(matrix[i][j]!=NULL)
                {
                   for(l=0;l<matrix[i][j]->RowPtr[k+1]-matrix[i][j]->RowPtr[k];++l)
                   {
                         matrixb->KCol[ matrixb->RowPtr[tmp_col+k]+tmp_entries+l] = matrix[i][j]->KCol[ matrix[i][j]->RowPtr[k]+l ] + tmp_row; 
                      matrixb->Entries[ matrixb->RowPtr[tmp_col+k]+tmp_entries+l] = matrix[i][j]->Entries[ matrix[i][j]->RowPtr[k]+l ]; 
                   }
                   tmp_entries +=  matrix[i][j]->RowPtr[k+1]-matrix[i][j]->RowPtr[k];
                }
                tmp_row += c[j]; //must bu put out of if setence
             } //j 
       } //k
       tmp_col += r[i];
    } //i
    printf("End of the matrix reassemble!\n");
}


//=======================================================youchunguang===============================================================
MATRIX * CopyMatrix( MATRIX *A )
{
    MATRIX *B = malloc( sizeof(MATRIX) );
    int i; 
    B->N_Rows    = A->N_Rows;
    B->N_Columns = A->N_Columns;
    B->N_Entries = A->N_Entries;
    
    B->RowPtr  = malloc( (A->N_Rows+1)*sizeof(int) );
    B->KCol    = malloc( A->N_Entries*sizeof(int) );
    B->Entries = malloc( A->N_Entries*sizeof(double) );
    
    memcpy( B->RowPtr,  A->RowPtr,    (A->N_Rows+1)*sizeof(int) );
    memcpy( B->KCol,    A->KCol,       A->N_Entries*sizeof(int) );
    memcpy( B->Entries, A->Entries,    A->N_Entries*sizeof(double) );
    
    return B;
}


//Only can be used to the symmestic matrix
MATRIX * EnlargeMatrix( MATRIX *A, int expan, double **vec, double **mat )
{
    MATRIX * B = malloc(sizeof(MATRIX));
    int nrow = A->N_Rows + expan;
    int ncol = A->N_Columns + expan;
    int nnz  = A->N_Entries + expan*A->N_Rows
	       + expan*A->N_Columns + expan*expan;
    double *aa = malloc( nnz*sizeof(double) );
    int    *ja = malloc( nnz*sizeof(int) );
    int    *ia = malloc( (nrow+1)*sizeof(int) );

    int i, j, k;
    int curr = 0;
    int start, end, length;
    ia[0] = 0;
    for( i=0; i<A->N_Rows; i++ )
    {
	start = A->RowPtr[i];
	end   = A->RowPtr[i+1];
	length = end - start;
	ia[i+1] = ia[i] + length + expan;
	memcpy( &ja[ia[i]], &A->KCol[start], length*sizeof(int) );
	memcpy( &aa[ia[i]], &A->Entries[start], length*sizeof(double) );
	for( j=0; j<expan; j++ )
	{
	    k = ia[i] + length;
	    ja[k+j] = A->N_Columns+j;
	    aa[k+j] = vec[j][i];
	}
    }

    for( i=0; i<expan; i++ )
    {
	k = A->N_Rows+i;
	ia[k+1] = ia[k] + A->N_Columns + expan;
	for( j=0; j<A->N_Columns; j++ )
	{
	    ja[ia[k]+j] = j;
	    aa[ia[k]+j] = vec[i][j];
	} 
	for( j=0; j<expan; j++ )
	{
	    ja[ia[k]+A->N_Columns+j] = A->N_Columns+j;
	    aa[ia[k]+A->N_Columns+j] = mat[i][j];
	}
    }    
    B->N_Rows    = nrow;
    B->N_Columns = ncol;
    B->N_Entries = nnz;
    
    B->RowPtr  = ia;
    B->KCol    = ja;
    B->Entries = aa;

    return B;
}

MATRIX *ExpandMatrixStruct( MATRIX *A, int expan )
{
    MATRIX * B = malloc(sizeof(MATRIX));
    int nrow = A->N_Rows + expan;
    int ncol = A->N_Columns + expan;
    int nnz  = A->N_Entries + expan*A->N_Rows
	       + expan*A->N_Columns + expan*expan;
    double *aa = malloc( nnz*sizeof(double) );
    int    *ja = malloc( nnz*sizeof(int) );
    int    *ia = malloc( (nrow+1)*sizeof(int) );

    int i, j, k;
    int curr = 0;
    int start, end, length;
    ia[0] = 0;
    for( i=0; i<A->N_Rows; i++ )
    {
	start = A->RowPtr[i];
	end   = A->RowPtr[i+1];
	length = end - start;
	ia[i+1] = ia[i] + length + expan;
	memcpy( &ja[ia[i]], &A->KCol[start], length*sizeof(int) );
	memcpy( &aa[ia[i]], &A->Entries[start], length*sizeof(double) );
	for( j=0; j<expan; j++ )
	{
	    k = ia[i] + length;
	    ja[k+j] = A->N_Columns+j;
	    aa[k+j] = 0;
	}
    }

    for( i=0; i<expan; i++ )
    {
	k = A->N_Rows+i;
	ia[k+1] = ia[k] + A->N_Columns + expan;
	for( j=0; j<A->N_Columns; j++ )
	{
	    ja[ia[k]+j] = j;
	    aa[ia[k]+j] = 0;
	} 
	for( j=0; j<expan; j++ )
	{
	    ja[ia[k]+A->N_Columns+j] = A->N_Columns+j;
	    aa[ia[k]+A->N_Columns+j] = 0;
	}
    }    
    B->N_Rows    = nrow;
    B->N_Columns = ncol;
    B->N_Entries = nnz;
    
    B->RowPtr  = ia;
    B->KCol    = ja;
    B->Entries = aa;

    return B;
}

void ExpandMatrix( MATRIX *A, int expan, double **vec, double **mat )
{
    double *aa = A->Entries;
    int    *ja = A->KCol;
    int    *ia = A->RowPtr;

    int i, j, k;
    for( i=0; i<A->N_Rows-expan; i++ )
    {
	for( j=0; j<expan; j++ )
	{
	    k = ia[i+1] -expan;
	    aa[k+j] = vec[j][i];
	}
    }

    for( i=0; i<expan; i++ )
    {
        k = A->N_Rows - expan + i;
	for( j=0; j<A->N_Columns-expan; j++ )
	{
	    aa[ia[k]+j] = vec[i][j];
	} 
	for( j=0; j<expan; j++ )
	{
	    aa[ia[k]+A->N_Columns-expan+j] = mat[i][j];
	}
    }
}



















