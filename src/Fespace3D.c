/*
 * =====================================================================================
 *
 *       Filename:  Fespace3D.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2014年10月20日 16时17分01秒
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
 
 #include "Fespace3D.h"
 #include "AllBases3D.h"
 #include "Base3D.h"
 #include "Enumerations.h"
 #include "Constants.h"
 #include "C_T_P1_3D.h"

 
 /** build the fespace3D from the mesh, femtype and 
    boundary condition */ 
Fespace3D *BUILDFEMSPACE3D(MESH *mesh, FemType femtype, BoundCondFunct3D *boundcond3D)
{
    Fespace3D *fesp3D;  
    fesp3D = malloc(sizeof(Fespace3D));
    fesp3D->Num_Volus = mesh->Num_Volus_Global;
    fesp3D->Mesh = mesh;
    // build the base function information
    BASEFUNCTION3D *base3D;
    base3D = BuildBase3D(femtype);
    fesp3D->Base = base3D;
	fesp3D->DOF_ALL = (mesh->Num_Verts_Global)*base3D->DOF[0]+ (mesh->Num_Lines_Global)*base3D->DOF[1] +
		      (mesh->Num_Faces_Global)*base3D->DOF[2]+ (mesh->Num_Volus_Global)*base3D->DOF[3] ;
    if(boundcond3D)
        fesp3D->BoundCond = boundcond3D;
    // build the dof structure
    fesp3D->BeginIndex = malloc((mesh->Num_Volus_Global+1)*sizeof(int));
    fesp3D->GlobalNumbers = malloc((mesh->Num_Volus_Global)*(fesp3D->Base->Num_Bas)*sizeof(int));
    BuildDOF3D(mesh,base3D,fesp3D->BeginIndex,fesp3D->GlobalNumbers);
    return fesp3D;  
         
}

Fespace3D *BUILDFEMSPACE3DBASE(MESH *mesh, BASEFUNCTION3D *base3D,BoundCondFunct3D *boundcond3D)
{
    Fespace3D *fesp3D;  	
    fesp3D = malloc(sizeof(Fespace3D));
    fesp3D->Num_Volus = mesh->Num_Volus_Global;
    fesp3D->Mesh = mesh;
    // build the base function information
    //BASEFUNCTION2D *base2D;
    fesp3D->Base = base3D;

    fesp3D->DOF_ALL = (mesh->Num_Verts_Global)*base3D->DOF[0]+ (mesh->Num_Lines_Global)*base3D->DOF[1] +
		      (mesh->Num_Faces_Global)*base3D->DOF[2]+ (mesh->Num_Volus_Global)*base3D->DOF[3] ;
    if(boundcond3D)
        fesp3D->BoundCond = boundcond3D;
    // build the dof structure
    fesp3D->BeginIndex = malloc((mesh->Num_Volus_Global+1)*sizeof(int));
    fesp3D->GlobalNumbers = malloc((mesh->Num_Volus_Global)*(fesp3D->Base->Num_Bas)*sizeof(int));
    BuildDOF3D(mesh,base3D,fesp3D->BeginIndex,fesp3D->GlobalNumbers);
    return fesp3D;  
}


void OutPutFespace3D(Fespace3D *fesp3D)
{
    int num_volus,i,j,dof_all, num_bas;
    MESH *mesh;
    mesh = fesp3D->Mesh;
    num_volus = mesh->Num_Volus_Global;
    dof_all = (mesh->Num_Volus_Global)*(fesp3D->Base->Num_Bas);
    num_bas = fesp3D->Base->Num_Bas;
    printf("The information of the BeginIndex:\n");
    for(i=0;i<=num_volus;i++)
      printf("BeginIndex[%d]=%d\n",i,fesp3D->BeginIndex[i]);
    
    printf("The information of the GlobalNumbers:\n");
    for(i=0;i<num_volus;i++)
    {
      for(j=0;j<num_bas;j++)
	printf("    GNO[%d]=%d",i*num_bas+j,fesp3D->GlobalNumbers[i*num_bas+j]);
      printf("\n");
    }
}

BASEFUNCTION3D *BuildBase3D(FemType femtype)
{
  //BASEFUNCTION2D *basefun2D;

  switch(femtype)
  {
    case C_T_P1_3D:
      return BuildBaseFunct3D(C_T_P1_3D_Num_Bas,C_T_P1_3D_Value_Dim,
				   C_T_P1_3D_Polydeg,C_T_P1_3D_Accuracy,
				   C_T_P1_3D_Maptype, C_T_P1_3D_dof,
				   C_T_P1_3D_nodal_points,C_T_P1_3D_D000,
				   C_T_P1_3D_D100,C_T_P1_3D_D010,C_T_P1_3D_D001,
				   C_T_P1_3D_D110,C_T_P1_3D_D101,C_T_P1_3D_D011,
                                   C_T_P1_3D_D200,C_T_P1_3D_D020,C_T_P1_3D_D002,
				   C_T_P1_3D_Nodal);
      printf("Use the continuous P1 conforming finite element method!\n");
      break;
    
    case C_T_P2_3D:
      return BuildBaseFunct3D(C_T_P2_3D_Num_Bas,C_T_P2_3D_Value_Dim,
				   C_T_P2_3D_Polydeg,C_T_P2_3D_Accuracy,
				   C_T_P2_3D_Maptype, C_T_P2_3D_dof,
				   C_T_P2_3D_nodal_points,C_T_P2_3D_D000,
				   C_T_P2_3D_D100,C_T_P2_3D_D010,C_T_P2_3D_D001,
				   C_T_P2_3D_D110,C_T_P2_3D_D101,C_T_P2_3D_D011,
                                   C_T_P2_3D_D200,C_T_P2_3D_D020,C_T_P2_3D_D002,
				   C_T_P2_3D_Nodal);
      printf("Use the continuous P2 conforming finite element method!\n");
      break;
   /*   
    case C_T_P3_3D:
      return BuildBaseFunct3D(C_T_P3_2D_Num_Bas,C_T_P3_2D_Value_Dim,
				   C_T_P3_3D_Polydeg,C_T_P3_3D_Accuracy,
				   C_T_P3_3D_Maptype, C_T_P3_3D_dof,
				   C_T_P3_3D_nodal_points,C_T_P3_3D_D000,
                                   C_T_P3_3D_D100,C_T_P3_3D_D010,C_T_P3_3D_D001,
				   C_T_P3_3D_D110,C_T_P3_3D_D101,C_T_P3_3D_D011,
                                   C_T_P3_3D_D200,C_T_P3_3D_D020,C_T_P3_3D_D002,
				   C_T_P3_3D_Nodal);
      printf("Use the continuous P3 conforming finite element method!\n");
      break;
   */   
    default:
      printf("No finite element is defined!\n");
  }  
  
  //printf("end of base function building!\n");
  //return basefun2D;
}
// build the dof distribution
void BuildDOF3D(MESH *mesh,BASEFUNCTION3D *base3D,int *BeginIndex,int *GlobalNumbers)
{
  int i,j,k, curr_volu, num_bas;
  int *DOF = base3D->DOF;
  
  num_bas = base3D->Num_Bas;
  int num_vert_volu, num_line_volu,num_face_volu;
  int N_Verts,N_Lines, N_Faces,N_Volus;
  N_Verts = mesh->Num_Verts_Global;
  N_Lines = mesh->Num_Lines_Global;
  N_Faces = mesh->Num_Faces_Global;
  N_Volus = mesh->Num_Volus_Global;
  int curr_index;
  curr_index = 0;
  int base_index_1,base_index_2,base_index_3;
  base_index_1 = DOF[0]*N_Verts;
  base_index_2 = DOF[0]*N_Verts + DOF[1]*N_Lines;
  base_index_3 = DOF[0]*N_Verts + DOF[1]*N_Lines +DOF[2]*N_Faces;
  int Line2vert[6] = {0,0,0,2,3,1};
  
  for(curr_volu=0;curr_volu<N_Volus;curr_volu++)
  {
    num_vert_volu = 4;
    num_line_volu = 6;
    num_face_volu = 4;

    BeginIndex[curr_volu] = curr_index;
    for(i=0;i<DOF[0];i++)
    {
      for(j=0;j<num_vert_volu;j++)
      {
	GlobalNumbers[curr_index] = i*N_Verts + mesh->Volus[curr_volu]->Verts[j];	
	curr_index ++ ;
      }//end for j
    }//end for i
    
    //give the dof on the lines
    for(i=0;i<DOF[1];i++)
    {
      for(j=0;j<num_line_volu;j++)
      {
	if(DOF[1] <= 1)
	{
	  GlobalNumbers[curr_index] = base_index_1 + i*N_Lines + mesh->Volus[curr_volu]->Lines[j]; 
	}
	else//现在先给一个边上自由度排列，有了总体自由度编号后，根据标准单元上基函数的排列，得到每一个单元上局部自由度的排列。边上自由度排列是从边的零点开始，往一点方向排。	   
	{
	  if((mesh->Lines[mesh->Volus[curr_volu]->Lines[j]]->Verts[0])==(mesh->Volus[curr_volu]->Verts[Line2vert[j]]))
	  {
	    GlobalNumbers[curr_index] = base_index_1 + i*N_Lines + mesh->Volus[curr_volu]->Lines[j]; 
	  }
	  else
	  {
	    GlobalNumbers[curr_index] = base_index_1 + (DOF[1]-i-1)*N_Lines + mesh->Volus[curr_volu]->Lines[j]; 
	  }//end for if vert	  
	}//end for DOF
	curr_index++;
      }//end for j
    }//end for i
    // give the dof in the face
    for(i=0;i<DOF[2];i++)//面上自由度的方向？？？？？？？？？？？？？？？
    {
       for(j=0;j<num_face_volu;j++)
       {
        GlobalNumbers[curr_index] = base_index_2 + i*N_Faces + mesh->Volus[curr_volu]->Faces[j];
        curr_index++;
       }
    }//end for i
    //单元上的自由度
    for(i=0;i<DOF[3];i++)
    {
        GlobalNumbers[curr_index] = base_index_3 + i*N_Volus + curr_volu;
        curr_index++;
    }//end for i

  //BeginIndex[curr_volu] = curr_index;
  }//end for curr_volu
  BeginIndex[curr_volu] = curr_index;
  
}

void FreeFespace3D(Fespace3D *fesp)
{
   printf ( "free the fespace!\n" );
  if(fesp->BeginIndex)
    free(fesp->BeginIndex);
  if(fesp->GlobalNumbers)
    free(fesp->GlobalNumbers);  
  //if(fesp->BoundCond)
   // FreeBoundCond2D(fesp->BoundCond);

  if(fesp->Base)
  {
    FreeBase3D(fesp->Base);
    //free(fesp->Base);
  }
  free(fesp);
}
void FreeFespace3DL(Fespace3D *fesp)
{
   printf ( "free the fespaceL!\n" );
  if(fesp->BeginIndex)
    free(fesp->BeginIndex);
  if(fesp->GlobalNumbers)
    free(fesp->GlobalNumbers);  
  //if(fesp->BoundCond)
   // FreeBoundCond2D(fesp->BoundCond);

  if(fesp->Base)
  {
    FreeBase3D(fesp->Base);
    //free(fesp->Base);
  }
}
// void FreeBoundCond3D(BOUNDCOND2D *BoundCond)
// {
//   free(BoundCond->ID_BDComp);
//   free(BoundCond->BdTypes);
//   free(BoundCond);
// }
