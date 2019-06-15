/*
 * =====================================================================================
 *
 *       Filename:  MG.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2014/12/22 03时09分20秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:   xufei@lsec.cc.ac.cn
 *        Company:  
 *    
 *
 * =====================================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "MG.h"
#include "Rhst.h"
#include "Mesh3D.h"
#include "Fespace3D.h"
//=========================================================================
//=============Some operators for the multilevel methods
//=========================================================================

MULTILEVEL *BuildMultiLevel(MESH *Mesh,DISCRETEFORM3D *discreteform_stiff, DISCRETEFORM3D *discreteform_mass)
{
	printf("Build the multilevel object!\n");
	MULTILEVEL *MultiLevel;
	MultiLevel = malloc(sizeof(MULTILEVEL));
	MultiLevel->Num_Levels = 1;
	MultiLevel->mesh = Mesh;
         //MultiLevel->Anasatz_Space = anasatz_space;
	//MultiLevel->Test_Space = test_space;
	MultiLevel->Stiff_Matrix = malloc(sizeof(MATRIX*));
	//MultiLevel->Stiff_Matrix[0] = malloc(sizeof(MATRIX));
	MultiLevel->Stiff_Matrix[0] = BuildMatrix(discreteform_stiff);
	AssembleMatrix(MultiLevel->Stiff_Matrix[0]);

	if(discreteform_mass!=NULL)
	{
		MultiLevel->Mass_Matrix = malloc(sizeof(MATRIX*));
		//MultiLevel->Mass_Matrix[0] = malloc(sizeof(MATRIX));
		MultiLevel->Mass_Matrix[0] = BuildMatrix(discreteform_mass);
		AssembleMatrix(MultiLevel->Mass_Matrix[0]);		
	}
	else
	{
		MultiLevel->Mass_Matrix = NULL;		
	}
        //MultiLevel->Rhst = malloc(sizeof(double*));
        //MultiLevel->Rhst[0] = malloc(sizeof(RHST));
        //MultiLevel->Rhst[0] = rhs;
	//AssembleRHST(rhs);
        //MultiLevel->Rhst[0] = malloc(sizeof(double)*rhs->DOF_ALL);
	//int i;
	//for(i=0;i<rhs->DOF_ALL;++i)
        //MultiLevel->Rhst[0][i]=rhs->Entries[i];

	MultiLevel->Prolongs = NULL;
	MultiLevel->Restricts = NULL;
	return MultiLevel;
}
  //如何定义在现有的Multilevel上增加一层新的网格和相应的矩阵，其实我们只需要矩阵
  //只需要告诉程序需要加密成什么样子的有限元空间，他们之间的转化矩阵，
  //这里包括形成刚度矩阵和质量矩阵

void AddLevelFEM(MULTILEVEL *multilevel)
  {
   printf("Add a level in the multilevel object!\n");
   double t1,t2,t3,t4,t5;
   t1 = GetTime();
   int num_levels,i;
   num_levels = multilevel->Num_Levels;
   MESH *mesh_tmp;
   mesh_tmp = malloc(sizeof(MESH));
   CopyMesh(mesh_tmp,multilevel->mesh);
   //MeshConsist(multilevel->mesh,1);
   MeshBisection(multilevel->mesh,1);

  t2 = GetTime();
   Fespace3D *new_anasatz_fespace, *new_test_space;
   new_test_space = BUILDFEMSPACE3DBASE(multilevel->mesh,multilevel->Stiff_Matrix[0]->DiscreteForm3D->Test_Space->Base,
	                                multilevel->Stiff_Matrix[0]->DiscreteForm3D->Test_Space->BoundCond);
   new_anasatz_fespace = BUILDFEMSPACE3DBASE(multilevel->mesh,multilevel->Stiff_Matrix[0]->DiscreteForm3D->Anasatz_Space->Base,
	                                multilevel->Stiff_Matrix[0]->DiscreteForm3D->Anasatz_Space->BoundCond);
   Fespace3D *old_anasatz_space, *old_test_space;
   old_anasatz_space = multilevel->Stiff_Matrix[num_levels-1]->DiscreteForm3D->Anasatz_Space;
   old_test_space = multilevel->Stiff_Matrix[num_levels-1]->DiscreteForm3D->Test_Space;
   old_anasatz_space->Mesh = mesh_tmp; 
   multilevel->Stiff_Matrix[num_levels-1]->DiscreteForm3D->Anasatz_Space = new_anasatz_fespace;
   multilevel->Stiff_Matrix[num_levels-1]->DiscreteForm3D->Test_Space = new_test_space;
   if(multilevel->Mass_Matrix!=NULL)
   {
    multilevel->Mass_Matrix[num_levels-1]->DiscreteForm3D->Anasatz_Space = new_anasatz_fespace;
    multilevel->Mass_Matrix[num_levels-1]->DiscreteForm3D->Test_Space = new_test_space;
   }
   num_levels += 1;
   multilevel->Stiff_Matrix = realloc(multilevel->Stiff_Matrix,num_levels*sizeof(MATRIX*));
   multilevel->Stiff_Matrix[num_levels-1] = 
                         BuildMatrix(multilevel->Stiff_Matrix[num_levels-2]->DiscreteForm3D);
   AssembleMatrix(multilevel->Stiff_Matrix[num_levels-1]);
   if(multilevel->Mass_Matrix!=NULL)
   {
    multilevel->Mass_Matrix = realloc(multilevel->Mass_Matrix,num_levels*sizeof(MATRIX*));  
    multilevel->Mass_Matrix[num_levels-1] =  
                         BuildMatrix(multilevel->Mass_Matrix[num_levels-2]->DiscreteForm3D);
    AssembleMatrix(multilevel->Mass_Matrix[num_levels-1]);
   }
   t3 = GetTime();
    if(num_levels>1)
   {
	  if(num_levels == 2)
	  {		 
		  printf("Build the interpolation operator!\n");
		  multilevel->Prolongs = malloc((num_levels-1)*sizeof(MATRIX*));
		  multilevel->Restricts = malloc((num_levels-1)*sizeof(MATRIX*));
	  }
	  else
	  {
		  printf("Add an interpolation operator!\n");
		  multilevel->Prolongs = realloc(multilevel->Prolongs, (num_levels-1)*sizeof(MATRIX*));
		  multilevel->Restricts = realloc(multilevel->Restricts, (num_levels-1)*sizeof(MATRIX*));
	  }
	  multilevel->Prolongs[num_levels-2] = malloc(sizeof(MATRIX));
	  multilevel->Restricts[num_levels-2] = malloc(sizeof(MATRIX));
	 
	  BuildProlongFa(old_anasatz_space,new_anasatz_fespace, multilevel->Prolongs[num_levels-2]);
	  //BuildRestrictFa(old_anasatz_space,new_anasatz_fespace,multilevel->Restricts[num_levels-2]);
	  TransposeMatrix(multilevel->Prolongs[num_levels-2],multilevel->Restricts[num_levels-2]);
  }
   //creat a new level 
   //multilevel->Anasatz_Space = new_anasatz_fespace;
   //multilevel->Test_Space = new_test_space;
   multilevel->Num_Levels = num_levels;
   FreeMesh(mesh_tmp);
   //FreeFespace3D(old_anasatz_space);
   //FreeFespace3D(old_test_space);
   printf("End of adding a level!\n");

}




  void AddLevelAFEM(MULTILEVEL *multilevel,int *M,int *num)
  {
   printf("Add a level in the multilevel object!\n");
   int num_levels;
   num_levels = multilevel->Num_Levels;
   MESH *mesh_tmp;
   mesh_tmp = malloc(sizeof(MESH));
   CopyMesh(mesh_tmp,multilevel->mesh);
   MeshRefineBisection(multilevel->mesh,M,num);
   CheckMesh(multilevel->mesh);
   MeshFullInfo4Line(multilevel->mesh);
   MeshRenew(multilevel->mesh);
   MeshFullInfo4Face(multilevel->mesh);
   MeshFullInfo4Num(multilevel->mesh);
   CheckMesh(mesh_tmp);
   Fespace3D *new_anasatz_fespace, *new_test_space;
   new_test_space = BUILDFEMSPACE3DBASE(multilevel->mesh,multilevel->Stiff_Matrix[0]->DiscreteForm3D->Test_Space->Base,
	                                multilevel->Stiff_Matrix[0]->DiscreteForm3D->Test_Space->BoundCond);
   new_anasatz_fespace = BUILDFEMSPACE3DBASE(multilevel->mesh,multilevel->Stiff_Matrix[0]->DiscreteForm3D->Anasatz_Space->Base,
	                                multilevel->Stiff_Matrix[0]->DiscreteForm3D->Anasatz_Space->BoundCond);
   Fespace3D *old_anasatz_space, *old_test_space;
   old_anasatz_space = multilevel->Stiff_Matrix[num_levels-1]->DiscreteForm3D->Anasatz_Space;
   old_test_space = multilevel->Stiff_Matrix[num_levels-1]->DiscreteForm3D->Test_Space;
   old_anasatz_space->Mesh = mesh_tmp; 
   multilevel->Stiff_Matrix[num_levels-1]->DiscreteForm3D->Anasatz_Space = new_anasatz_fespace;
   multilevel->Stiff_Matrix[num_levels-1]->DiscreteForm3D->Test_Space = new_test_space;
   if(multilevel->Mass_Matrix!=NULL)
   {
    multilevel->Mass_Matrix[num_levels-1]->DiscreteForm3D->Anasatz_Space = new_anasatz_fespace;
    multilevel->Mass_Matrix[num_levels-1]->DiscreteForm3D->Test_Space = new_test_space;
   }
   num_levels += 1;
   multilevel->Stiff_Matrix = realloc(multilevel->Stiff_Matrix,num_levels*sizeof(MATRIX*));
   multilevel->Stiff_Matrix[num_levels-1] = 
                         BuildMatrix(multilevel->Stiff_Matrix[num_levels-2]->DiscreteForm3D);
   AssembleMatrix(multilevel->Stiff_Matrix[num_levels-1]);
   if(multilevel->Mass_Matrix!=NULL)
   {
    multilevel->Mass_Matrix = realloc(multilevel->Mass_Matrix,num_levels*sizeof(MATRIX*));  
    multilevel->Mass_Matrix[num_levels-1] =  
                         BuildMatrix(multilevel->Mass_Matrix[num_levels-2]->DiscreteForm3D);
    AssembleMatrix(multilevel->Mass_Matrix[num_levels-1]);
   }
    if(num_levels>1)
   {
	  if(num_levels == 2)
	  {		 
		  printf("Build the interpolation operator!\n");
		  multilevel->Prolongs = malloc((num_levels-1)*sizeof(MATRIX*));
		  multilevel->Restricts = malloc((num_levels-1)*sizeof(MATRIX*));
	  }
	  else
	  {
		  printf("Add an interpolation operator!\n");
		  multilevel->Prolongs = realloc(multilevel->Prolongs, (num_levels-1)*sizeof(MATRIX*));
		  multilevel->Restricts = realloc(multilevel->Restricts, (num_levels-1)*sizeof(MATRIX*));
	  }
	  multilevel->Prolongs[num_levels-2] = malloc(sizeof(MATRIX));
	  multilevel->Restricts[num_levels-2] = malloc(sizeof(MATRIX));
	 
	  BuildProlongFa(old_anasatz_space,new_anasatz_fespace, multilevel->Prolongs[num_levels-2]);
	  //BuildRestrictFa(old_anasatz_space,new_anasatz_fespace,multilevel->Restricts[num_levels-2]);
	  TransposeMatrix(multilevel->Prolongs[num_levels-2],multilevel->Restricts[num_levels-2]);
  }
   multilevel->Num_Levels = num_levels;
   FreeMesh(mesh_tmp);
   printf("End of adding a level!\n");

}
//释放Multievel所占用的空间
void FreeMultiLevel(MULTILEVEL *multilevel)
{
  printf("begin to free Multilevl\n");
  int i, num_levels;
  num_levels = multilevel->Num_Levels;
  //FreeMesh(multilevel->mesh); 
  for(i=0;i<num_levels;i++)
  {
    if(multilevel->Stiff_Matrix[i]->N_Rows>0)
		FreeMatrix(multilevel->Stiff_Matrix[i]);
    if((multilevel->Mass_Matrix!=NULL)&&(multilevel->Mass_Matrix[i]->N_Rows>0))
		FreeMatrix(multilevel->Mass_Matrix[i]);   
    if(num_levels>1&&i<num_levels-1)
      FreeMatrix(multilevel->Prolongs[i]);
    //if(multilevel->Restricts[i]->N_Rows>0)
    if(num_levels>1&&i<num_levels-1)
      FreeMatrix(multilevel->Restricts[i]);   
  }//end for
  if(num_levels > 0)
  {
    free(multilevel->Stiff_Matrix);
    free(multilevel->Prolongs);
    free(multilevel->Restricts);
    if(multilevel->Mass_Matrix!=NULL)
      free(multilevel->Mass_Matrix);
  }
  multilevel->Num_Levels = 0;
  free(multilevel);
 printf("end free multilevel\n");
}
