/*
 * =====================================================================================
 *
 *       Filename:  meshrefineconsistent.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2014年09月20日 01时16分09秒
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
#include <math.h>
#include "Mesh3D.h"
#include "Quadrature3D.h"
#define NVERTS_PER_FACEE 3
#define NLINES_PER_FACEE 3
#define NVERTS_PER_VOLU 4
#define NLINES_PER_VOLU 6
#define NFACEES_PER_VOLU 4
#define NVERTS_PER_LINE 2




  


void
MeshRefineConsistent(MESH *mesh)
{
   int i,j,k;
   int num,NP,NF,NE,NV;
   //加密一次点线面体的个数变化
   int tmp_nverts = mesh->Num_Verts_Global+mesh->Num_Lines_Global;
   int tmp_nedges = mesh->Num_Lines_Global*2+mesh->Num_Faces_Global*3+mesh->Num_Volus_Global;
   int tmp_nfaces = mesh->Num_Faces_Global*4+mesh->Num_Volus_Global*8;
   int tmp_nvolus = mesh->Num_Volus_Global*8;

   NP= mesh->Num_Verts_Global;
   NE= mesh->Num_Lines_Global;
   NF= mesh->Num_Faces_Global;
   NV= mesh->Num_Volus_Global;
   int num_face[4];
   mesh->Verts=(VERT **)realloc(mesh->Verts,tmp_nverts);
   mesh->Lines=(LINE **)realloc(mesh->Lines,tmp_nedges);
   mesh->Faces=(FACEE **)realloc(mesh->Faces,tmp_nfaces);
   mesh->Volus=(VOLU **)realloc(mesh->Volus,tmp_nvolus);
 
   for(i=0;i<NE;++i)
  {
     num=i+NP;
     mesh->Verts[num]=(VERT *)malloc(1);
  }
  for(i=0;i<tmp_nedges-NE;++i)
  {
     num=i+NE;
     mesh->Lines[num]=(LINE *)malloc(1);
  }
  for(i=0;i<tmp_nfaces-NF;++i)
  {
     num=i+NF;
     mesh->Faces[num]=(FACEE *)malloc(1);
  }
  for(i=0;i<tmp_nvolus-NV;++i)
  {
     num=i+NV;
     mesh->Volus[num]=(VOLU *)malloc(1);
  }

   //更新点的信息,旧点编号不变，新点编号为 NP+边的编号。
  for(i=0;i<NE;++i)
  {
     num=i+NP;
     mesh->Verts[num]=(VERT *)malloc(1);
     for(j=0;j<3;++j)
     {
	mesh->Verts[num]->Coord[j]=(mesh->Verts[mesh->Lines[i]->Verts[0]]->Coord[j]
	                           +mesh->Verts[mesh->Lines[i]->Verts[1]]->Coord[j])*1/2;
     }
     mesh->Verts[num]->Num_Lines_Owned = -1;
     mesh->Verts[num]->Lines_Owned = NULL;
     mesh->Verts[num]->Num_Faces_Owned = -1;
     mesh->Verts[num]->Faces_Owned = NULL;
     mesh->Verts[num]->Num_Volus_Owned = -1;
     mesh->Verts[num]->Volus_Owned = NULL;
  }
//更新单元上点的信息。 单元的编号顺序：0点对应的单元编号不变，其它七个排在ＮＶ后面。

  for(i=0;i<NV;++i)
  {
     //七个新单元
     mesh->Volus[NV+i*7+0]->Verts[0]=NP+mesh->Volus[i]->Lines[0];
     mesh->Volus[NV+i*7+0]->Verts[1]=mesh->Volus[i]->Verts[1];
     mesh->Volus[NV+i*7+0]->Verts[2]=NP+mesh->Volus[i]->Lines[3];
     mesh->Volus[NV+i*7+0]->Verts[3]=NP+mesh->Volus[i]->Lines[4];
     mesh->Volus[NV+i*7+0]->Ancestor=mesh->Volus[i]->Ancestor;
     mesh->Volus[NV+i*7+0]->Father=i;


     mesh->Volus[NV+i*7+1]->Verts[0]=NP+mesh->Volus[i]->Lines[1];
     mesh->Volus[NV+i*7+1]->Verts[1]=NP+mesh->Volus[i]->Lines[3];
     mesh->Volus[NV+i*7+1]->Verts[2]=mesh->Volus[i]->Verts[2];
     mesh->Volus[NV+i*7+1]->Verts[3]=NP+mesh->Volus[i]->Lines[5];
     mesh->Volus[NV+i*7+1]->Ancestor=mesh->Volus[i]->Ancestor;
     mesh->Volus[NV+i*7+1]->Father=i;

     mesh->Volus[NV+i*7+2]->Verts[0]=NP+mesh->Volus[i]->Lines[2];
     mesh->Volus[NV+i*7+2]->Verts[1]=NP+mesh->Volus[i]->Lines[4];
     mesh->Volus[NV+i*7+2]->Verts[2]=NP+mesh->Volus[i]->Lines[5];
     mesh->Volus[NV+i*7+2]->Verts[3]=mesh->Volus[i]->Verts[3];
     mesh->Volus[NV+i*7+2]->Ancestor=mesh->Volus[i]->Ancestor;
     mesh->Volus[NV+i*7+2]->Father=i;

     mesh->Volus[NV+i*7+3]->Verts[0]=NP+mesh->Volus[i]->Lines[0];
     mesh->Volus[NV+i*7+3]->Verts[1]=NP+mesh->Volus[i]->Lines[1];
     mesh->Volus[NV+i*7+3]->Verts[2]=NP+mesh->Volus[i]->Lines[2];
     mesh->Volus[NV+i*7+3]->Verts[3]=NP+mesh->Volus[i]->Lines[3];
     mesh->Volus[NV+i*7+3]->Ancestor=mesh->Volus[i]->Ancestor;
     mesh->Volus[NV+i*7+3]->Father=i;

     mesh->Volus[NV+i*7+4]->Verts[0]=NP+mesh->Volus[i]->Lines[0];
     mesh->Volus[NV+i*7+4]->Verts[1]=NP+mesh->Volus[i]->Lines[2];
     mesh->Volus[NV+i*7+4]->Verts[2]=NP+mesh->Volus[i]->Lines[4];
     mesh->Volus[NV+i*7+4]->Verts[3]=NP+mesh->Volus[i]->Lines[3];
     mesh->Volus[NV+i*7+4]->Ancestor=mesh->Volus[i]->Ancestor;
     mesh->Volus[NV+i*7+4]->Father=i;

     mesh->Volus[NV+i*7+5]->Verts[0]=NP+mesh->Volus[i]->Lines[5];
     mesh->Volus[NV+i*7+5]->Verts[1]=NP+mesh->Volus[i]->Lines[1];
     mesh->Volus[NV+i*7+5]->Verts[2]=NP+mesh->Volus[i]->Lines[3];
     mesh->Volus[NV+i*7+5]->Verts[3]=NP+mesh->Volus[i]->Lines[2];
     mesh->Volus[NV+i*7+5]->Ancestor=mesh->Volus[i]->Ancestor;
     mesh->Volus[NV+i*7+5]->Father=i;

     mesh->Volus[NV+i*7+6]->Verts[0]=NP+mesh->Volus[i]->Lines[5];
     mesh->Volus[NV+i*7+6]->Verts[1]=NP+mesh->Volus[i]->Lines[2];
     mesh->Volus[NV+i*7+6]->Verts[2]=NP+mesh->Volus[i]->Lines[3];
     mesh->Volus[NV+i*7+6]->Verts[3]=NP+mesh->Volus[i]->Lines[4];
     mesh->Volus[NV+i*7+6]->Ancestor=mesh->Volus[i]->Ancestor;
     mesh->Volus[NV+i*7+6]->Father=i;
     //旧单元
     mesh->Volus[i]->Verts[1]=mesh->Volus[mesh->Num_Volus_Global+i*7+0]->Verts[0];
     mesh->Volus[i]->Verts[2]=mesh->Volus[mesh->Num_Volus_Global+i*7+1]->Verts[0];
     mesh->Volus[i]->Verts[3]=mesh->Volus[mesh->Num_Volus_Global+i*7+2]->Verts[0];
     mesh->Volus[i]->Father=i;
  }
 //更新面包含的点。首先遍历单元，更新里面的八个面。
 for(i=0;i<NV;++i)
    {
      mesh->Faces[NF*4+i*8+0]->Verts[0]=NP+mesh->Volus[i]->Lines[0];
      mesh->Faces[NF*4+i*8+0]->Verts[1]=NP+mesh->Volus[i]->Lines[1];
      mesh->Faces[NF*4+i*8+0]->Verts[2]=NP+mesh->Volus[i]->Lines[2];
      mesh->Faces[NF*4+i*8+0]->ID_Boundary=0;
      
      mesh->Faces[NF*4+i*8+1]->Verts[0]=NP+mesh->Volus[i]->Lines[0];
      mesh->Faces[NF*4+i*8+1]->Verts[1]=NP+mesh->Volus[i]->Lines[4];
      mesh->Faces[NF*4+i*8+1]->Verts[2]=NP+mesh->Volus[i]->Lines[3];
      mesh->Faces[NF*4+i*8+1]->ID_Boundary=0;
      
      mesh->Faces[NF*4+i*8+2]->Verts[0]=NP+mesh->Volus[i]->Lines[5];
      mesh->Faces[NF*4+i*8+2]->Verts[1]=NP+mesh->Volus[i]->Lines[1];
      mesh->Faces[NF*4+i*8+2]->Verts[2]=NP+mesh->Volus[i]->Lines[3];
      mesh->Faces[NF*4+i*8+2]->ID_Boundary=0;

      mesh->Faces[NF*4+i*8+3]->Verts[0]=NP+mesh->Volus[i]->Lines[5];
      mesh->Faces[NF*4+i*8+3]->Verts[1]=NP+mesh->Volus[i]->Lines[4];
      mesh->Faces[NF*4+i*8+3]->Verts[2]=NP+mesh->Volus[i]->Lines[2];
      mesh->Faces[NF*4+i*8+3]->ID_Boundary=0;

      mesh->Faces[NF*4+i*8+4]->Verts[0]=NP+mesh->Volus[i]->Lines[1];
      mesh->Faces[NF*4+i*8+4]->Verts[1]=NP+mesh->Volus[i]->Lines[2];
      mesh->Faces[NF*4+i*8+4]->Verts[2]=NP+mesh->Volus[i]->Lines[3];
      mesh->Faces[NF*4+i*8+4]->ID_Boundary=0;

      mesh->Faces[NF*4+i*8+5]->Verts[0]=NP+mesh->Volus[i]->Lines[2];
      mesh->Faces[NF*4+i*8+5]->Verts[1]=NP+mesh->Volus[i]->Lines[4];
      mesh->Faces[NF*4+i*8+5]->Verts[2]=NP+mesh->Volus[i]->Lines[3];
      mesh->Faces[NF*4+i*8+5]->ID_Boundary=0;

      mesh->Faces[NF*4+i*8+6]->Verts[0]=NP+mesh->Volus[i]->Lines[0];
      mesh->Faces[NF*4+i*8+6]->Verts[1]=NP+mesh->Volus[i]->Lines[3];
      mesh->Faces[NF*4+i*8+6]->Verts[2]=NP+mesh->Volus[i]->Lines[2];
      mesh->Faces[NF*4+i*8+6]->ID_Boundary=0;

      mesh->Faces[NF*4+i*8+7]->Verts[0]=NP+mesh->Volus[i]->Lines[5];
      mesh->Faces[NF*4+i*8+7]->Verts[1]=NP+mesh->Volus[i]->Lines[2];
      mesh->Faces[NF*4+i*8+7]->Verts[2]=NP+mesh->Volus[i]->Lines[3];
      mesh->Faces[NF*4+i*8+7]->ID_Boundary=0;
    
      mesh->Faces[NF*4+i*8+0]->Num_Volus_Owned=-1;
      mesh->Faces[NF*4+i*8+0]->Volus_Owned=NULL;
     
      mesh->Faces[NF*4+i*8+1]->Num_Volus_Owned=-1;
      mesh->Faces[NF*4+i*8+1]->Volus_Owned=NULL;

      mesh->Faces[NF*4+i*8+2]->Num_Volus_Owned=-1;
      mesh->Faces[NF*4+i*8+2]->Volus_Owned=NULL;

      mesh->Faces[NF*4+i*8+3]->Num_Volus_Owned=-1;
      mesh->Faces[NF*4+i*8+3]->Volus_Owned=NULL;

      mesh->Faces[NF*4+i*8+4]->Num_Volus_Owned=-1;
      mesh->Faces[NF*4+i*8+4]->Volus_Owned=NULL;

      mesh->Faces[NF*4+i*8+5]->Num_Volus_Owned=-1;
      mesh->Faces[NF*4+i*8+5]->Volus_Owned=NULL;

      mesh->Faces[NF*4+i*8+6]->Num_Volus_Owned=-1;
      mesh->Faces[NF*4+i*8+6]->Volus_Owned=NULL;
     
      mesh->Faces[NF*4+i*8+7]->Num_Volus_Owned=-1;
      mesh->Faces[NF*4+i*8+7]->Volus_Owned=NULL;
    } 
//遍历面，更新外面的四个面包含的点的信息。
  for(i=0;i<NF;++i)
   {
     mesh->Faces[NF+i*3+0]->Verts[0]=mesh->Faces[i]->Verts[0];
     mesh->Faces[NF+i*3+0]->Verts[1]=NP+mesh->Faces[i]->Lines[2];
     mesh->Faces[NF+i*3+0]->Verts[2]=NP+mesh->Faces[i]->Lines[1];
     mesh->Faces[NF+i*3+0]->ID_Boundary=mesh->Faces[i]->ID_Boundary;
     
     mesh->Faces[NF+i*3+1]->Verts[0]=NP+mesh->Faces[i]->Lines[2];
     mesh->Faces[NF+i*3+1]->Verts[1]=mesh->Faces[i]->Verts[1];
     mesh->Faces[NF+i*3+1]->Verts[2]=NP+mesh->Faces[i]->Lines[0];
     mesh->Faces[NF+i*3+1]->ID_Boundary=mesh->Faces[i]->ID_Boundary;
   
     mesh->Faces[NF+i*3+2]->Verts[0]=NP+mesh->Faces[i]->Lines[1];
     mesh->Faces[NF+i*3+2]->Verts[1]=NP+mesh->Faces[i]->Lines[0];
     mesh->Faces[NF+i*3+2]->Verts[2]=mesh->Faces[i]->Verts[2];
     mesh->Faces[NF+i*3+2]->ID_Boundary=mesh->Faces[i]->ID_Boundary;
   
     mesh->Faces[i]->Verts[0]=NP+mesh->Faces[i]->Lines[0];
     mesh->Faces[i]->Verts[1]=NP+mesh->Faces[i]->Lines[1];
     mesh->Faces[i]->Verts[2]=NP+mesh->Faces[i]->Lines[2];
      
     mesh->Faces[NF+i*3+0]->Volus_Owned=NULL;
     mesh->Faces[NF+i*3+0]->Num_Volus_Owned=-1;
 
     mesh->Faces[NF+i*3+1]->Volus_Owned=NULL;
     mesh->Faces[NF+i*3+1]->Num_Volus_Owned=-1;

     mesh->Faces[NF+i*3+2]->Volus_Owned=NULL;
     mesh->Faces[NF+i*3+2]->Num_Volus_Owned=-1;
   }
  
//更新单元里面的唯一一条新边的顶点的信息。
  for(i=0;i<NV;++i)
   {
     mesh->Lines[NE*2+NF*3+i]->Verts[0]=mesh->Volus[i]->Lines[2]+NP;    
     mesh->Lines[NE*2+NF*3+i]->Verts[1]=mesh->Volus[i]->Lines[3]+NP;    
    
     mesh->Lines[NE*2+NF*3+i]->Num_Faces_Owned=-1;
     mesh->Lines[NE*2+NF*3+i]->Faces_Owned=NULL;
     mesh->Lines[NE*2+NF*3+i]->Num_Volus_Owned=-1;
     mesh->Lines[NE*2+NF*3+i]->Volus_Owned=NULL;
   }
//更新面上的三条边的顶点信息。
 
  for(i=0;i<NF;++i)
   {
     mesh->Lines[NE*2+i*3+0]->Verts[0]=mesh->Faces[i]->Lines[1]+NP;
     mesh->Lines[NE*2+i*3+0]->Verts[1]=mesh->Faces[i]->Lines[2]+NP;

     mesh->Lines[NE*2+i*3+1]->Verts[0]=mesh->Faces[i]->Lines[0]+NP;
     mesh->Lines[NE*2+i*3+1]->Verts[1]=mesh->Faces[i]->Lines[2]+NP;

     mesh->Lines[NE*2+i*3+2]->Verts[0]=mesh->Faces[i]->Lines[0]+NP;
     mesh->Lines[NE*2+i*3+2]->Verts[1]=mesh->Faces[i]->Lines[1]+NP;
   
     mesh->Lines[NE*2+i*3+0]->Num_Faces_Owned=-1;
     mesh->Lines[NE*2+i*3+0]->Faces_Owned=NULL;
     mesh->Lines[NE*2+i*3+0]->Num_Volus_Owned=-1;
     mesh->Lines[NE*2+i*3+0]->Volus_Owned=NULL;
    
     mesh->Lines[NE*2+i*3+1]->Num_Faces_Owned=-1;
     mesh->Lines[NE*2+i*3+1]->Faces_Owned=NULL;
     mesh->Lines[NE*2+i*3+1]->Num_Volus_Owned=-1;
     mesh->Lines[NE*2+i*3+1]->Volus_Owned=NULL;
   
     mesh->Lines[NE*2+i*3+2]->Num_Faces_Owned=-1;
     mesh->Lines[NE*2+i*3+2]->Faces_Owned=NULL;
     mesh->Lines[NE*2+i*3+2]->Num_Volus_Owned=-1;
     mesh->Lines[NE*2+i*3+2]->Volus_Owned=NULL;
   }
//更新边上的新边点的信息。
  for(i=0;i<NE;++i)
  {
    mesh->Lines[NE+i]->Verts[0]=i+NP;
    mesh->Lines[NE+i]->Verts[1]=mesh->Lines[i]->Verts[1];
    mesh->Lines[i]->Verts[0]=mesh->Lines[i]->Verts[0];
    mesh->Lines[i]->Verts[1]=i+NP;
   
    mesh->Lines[NE+i]->Num_Faces_Owned=-1;
    mesh->Lines[NE+i]->Faces_Owned=NULL;
    mesh->Lines[NE+i]->Num_Volus_Owned=-1;
    mesh->Lines[NE+i]->Volus_Owned=NULL;
 }



//更新边的信息： 首先更新单元包含的边,以及面包含的边
  for(i=0;i<NV;++i)  
  {

    for(j=0;j<4;++j)
      num_face[j]=mesh->Volus[i]->Faces[j];
  
     /* 第一条边 */
   if(mesh->Lines[mesh->Volus[i]->Lines[0]]->Verts[0]==mesh->Volus[NV+i*7+0]->Verts[1])
   {
    mesh->Volus[NV+i*7+0]->Lines[0]=mesh->Volus[i]->Lines[0];
    mesh->Volus[i]->Lines[0]=mesh->Volus[i]->Lines[0]+NE;
   }
   else
   {
    mesh->Volus[NV+i*7+0]->Lines[0]=mesh->Volus[i]->Lines[0]+NE;
    mesh->Volus[i]->Lines[0]=mesh->Volus[i]->Lines[0];
   }
   /* 第二条边 */
   if(mesh->Lines[mesh->Volus[i]->Lines[1]]->Verts[0]==mesh->Volus[NV+i*7+1]->Verts[2])
   {
    mesh->Volus[NV+i*7+1]->Lines[1]=mesh->Volus[i]->Lines[1];
    mesh->Volus[i]->Lines[1]=mesh->Volus[i]->Lines[1]+NE;
   }
   else
   {
    mesh->Volus[NV+i*7+1]->Lines[1]=mesh->Volus[i]->Lines[1]+NE;
    mesh->Volus[i]->Lines[1]=mesh->Volus[i]->Lines[1];
   }
   /*第三条边*/
   if(mesh->Lines[mesh->Volus[i]->Lines[2]]->Verts[0]==mesh->Volus[NV+i*7+2]->Verts[3])
   {
    mesh->Volus[NV+i*7+2]->Lines[2]=mesh->Volus[i]->Lines[2];
    mesh->Volus[i]->Lines[2]=mesh->Volus[i]->Lines[2]+NE;
   }
   else
   {
    mesh->Volus[NV+i*7+2]->Lines[2]=mesh->Volus[i]->Lines[2]+NE;
    mesh->Volus[i]->Lines[2]=mesh->Volus[i]->Lines[2];
   }
   /* 第四条边 */
   if(mesh->Lines[mesh->Volus[i]->Lines[3]]->Verts[0]==mesh->Volus[NV+i*7+0]->Verts[1])
   {
    mesh->Volus[NV+i*7+0]->Lines[3]=mesh->Volus[i]->Lines[3];
    mesh->Volus[NV+i*7+1]->Lines[3]=mesh->Volus[i]->Lines[3]+NE;
   }
   else
   {
    mesh->Volus[NV+i*7+0]->Lines[3]=mesh->Volus[i]->Lines[3]+NE;
    mesh->Volus[NV+i*7+1]->Lines[3]=mesh->Volus[i]->Lines[3];
   }
   /* 第五条边  */
   if(mesh->Lines[mesh->Volus[i]->Lines[4]]->Verts[0]==mesh->Volus[NV+i*7+0]->Verts[1])
   {
    mesh->Volus[NV+i*7+0]->Lines[4]=mesh->Volus[i]->Lines[4];
    mesh->Volus[NV+i*7+2]->Lines[4]=mesh->Volus[i]->Lines[4]+NE;
   }
   else
   {
    mesh->Volus[NV+i*7+0]->Lines[4]=mesh->Volus[i]->Lines[4]+NE;
    mesh->Volus[NV+i*7+2]->Lines[4]=mesh->Volus[i]->Lines[4];
   }
   /* 第六条边 */
   if(mesh->Lines[mesh->Volus[i]->Lines[5]]->Verts[0]==mesh->Volus[NV+i*7+1]->Verts[2])
   {
    mesh->Volus[NV+i*7+1]->Lines[5]=mesh->Volus[i]->Lines[5];
    mesh->Volus[NV+i*7+2]->Lines[5]=mesh->Volus[i]->Lines[5]+NE;
   }
   else
   {
    mesh->Volus[NV+i*7+1]->Lines[5]=mesh->Volus[i]->Lines[5]+NE;
    mesh->Volus[NV+i*7+2]->Lines[5]=mesh->Volus[i]->Lines[5];
   }
   /* 第一个面 */
   if(  (mesh->Faces[mesh->Volus[i]->Faces[0]]->Verts[0]==mesh->Volus[NV+i*7+0]->Verts[3])
   && (mesh->Faces[mesh->Volus[i]->Faces[0]]->Verts[1]==mesh->Volus[NV+i*7+0]->Verts[2]) )
   {
    //第一个面上的边被那些单元包含  
    mesh->Volus[NV+i*7+0]->Lines[5]=NE*2+mesh->Volus[i]->Faces[0]*3+2;
    mesh->Volus[NV+i*7+1]->Lines[4]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Volus[NV+i*7+2]->Lines[3]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Volus[NV+i*7+4]->Lines[5]=NE*2+mesh->Volus[i]->Faces[0]*3+2;
    mesh->Volus[NV+i*7+6]->Lines[5]=NE*2+mesh->Volus[i]->Faces[0]*3+2;
    mesh->Volus[NV+i*7+6]->Lines[2]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Volus[NV+i*7+5]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Volus[NV+i*7+6]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    //第一个面上的边被哪些面包含
    mesh->Faces[NF*4+i*8+1]->Lines[0]=NE*2+mesh->Volus[i]->Faces[0]*3+2;
    mesh->Faces[NF*4+i*8+5]->Lines[0]=NE*2+mesh->Volus[i]->Faces[0]*3+2;
    mesh->Faces[NF*4+i*8+2]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Faces[NF*4+i*8+7]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Faces[NF*4+i*8+3]->Lines[2]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+0]->Lines[0]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+0]->Lines[1]=mesh->Volus[NV+i*7+1]->Lines[3];
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+0]->Lines[2]=mesh->Volus[NV+i*7+1]->Lines[5];
    
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+1]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+1]->Lines[0]=mesh->Volus[NV+i*7+2]->Lines[4];
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+1]->Lines[2]=mesh->Volus[NV+i*7+2]->Lines[5];

    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+2]->Lines[2]=NE*2+mesh->Volus[i]->Faces[0]*3+2;
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+2]->Lines[0]=mesh->Volus[NV+i*7+0]->Lines[4];
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+2]->Lines[1]=mesh->Volus[NV+i*7+0]->Lines[3];

    mesh->Faces[mesh->Volus[i]->Faces[0]]->Lines[0]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Faces[mesh->Volus[i]->Faces[0]]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Faces[mesh->Volus[i]->Faces[0]]->Lines[2]=NE*2+mesh->Volus[i]->Faces[0]*3+2;
    //第一个面上产生的新面被哪些单元包含
    mesh->Volus[NV+i*7+0]->Faces[0]=NF+mesh->Volus[i]->Faces[0]*3+2;
    mesh->Volus[NV+i*7+1]->Faces[0]=NF+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Volus[NV+i*7+2]->Faces[0]=NF+mesh->Volus[i]->Faces[0]*3+1;

   } 
   else if (mesh->Faces[mesh->Volus[i]->Faces[0]]->Verts[0]==mesh->Volus[NV+i*7+0]->Verts[3]
	    &&mesh->Faces[mesh->Volus[i]->Faces[0]]->Verts[2]==mesh->Volus[NV+i*7+0]->Verts[2] )
   {
    mesh->Volus[NV+i*7+0]->Lines[5]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Volus[NV+i*7+1]->Lines[4]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Volus[NV+i*7+2]->Lines[3]=NE*2+mesh->Volus[i]->Faces[0]*3+2;
    mesh->Volus[NV+i*7+4]->Lines[5]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Volus[NV+i*7+6]->Lines[5]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Volus[NV+i*7+6]->Lines[2]=NE*2+mesh->Volus[i]->Faces[0]*3+2;
    mesh->Volus[NV+i*7+5]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Volus[NV+i*7+6]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
 
    mesh->Faces[NF*4+i*8+1]->Lines[0]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Faces[NF*4+i*8+5]->Lines[0]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Faces[NF*4+i*8+2]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Faces[NF*4+i*8+7]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Faces[NF*4+i*8+3]->Lines[2]=NE*2+mesh->Volus[i]->Faces[0]*3+2;

    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+0]->Lines[0]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+0]->Lines[2]=mesh->Volus[NV+i*7+1]->Lines[3];
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+0]->Lines[1]=mesh->Volus[NV+i*7+1]->Lines[5];
    
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+1]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+1]->Lines[2]=mesh->Volus[NV+i*7+0]->Lines[3];
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+1]->Lines[0]=mesh->Volus[NV+i*7+0]->Lines[4];

    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+2]->Lines[2]=NE*2+mesh->Volus[i]->Faces[0]*3+2;
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+2]->Lines[0]=mesh->Volus[NV+i*7+2]->Lines[4];
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+2]->Lines[1]=mesh->Volus[NV+i*7+2]->Lines[5];

    mesh->Faces[mesh->Volus[i]->Faces[0]]->Lines[0]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Faces[mesh->Volus[i]->Faces[0]]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Faces[mesh->Volus[i]->Faces[0]]->Lines[2]=NE*2+mesh->Volus[i]->Faces[0]*3+2;

    mesh->Volus[NV+i*7+0]->Faces[0]=NF+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Volus[NV+i*7+1]->Faces[0]=NF+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Volus[NV+i*7+2]->Faces[0]=NF+mesh->Volus[i]->Faces[0]*3+2;
   }
 
 else if (mesh->Faces[mesh->Volus[i]->Faces[0]]->Verts[1]==mesh->Volus[NV+i*7+0]->Verts[3]
          &&mesh->Faces[mesh->Volus[i]->Faces[0]]->Verts[2]==mesh->Volus[NV+i*7+0]->Verts[2] )
   {
    mesh->Volus[NV+i*7+0]->Lines[5]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Volus[NV+i*7+1]->Lines[4]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Volus[NV+i*7+2]->Lines[3]=NE*2+mesh->Volus[i]->Faces[0]*3+2;
    mesh->Volus[NV+i*7+4]->Lines[5]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Volus[NV+i*7+6]->Lines[5]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Volus[NV+i*7+6]->Lines[2]=NE*2+mesh->Volus[i]->Faces[0]*3+2;
    mesh->Volus[NV+i*7+5]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Volus[NV+i*7+6]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    
    mesh->Faces[NF*4+i*8+1]->Lines[0]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Faces[NF*4+i*8+5]->Lines[0]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Faces[NF*4+i*8+2]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Faces[NF*4+i*8+7]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Faces[NF*4+i*8+3]->Lines[2]=NE*2+mesh->Volus[i]->Faces[0]*3+2;
 
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+0]->Lines[0]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+0]->Lines[2]=mesh->Volus[NV+i*7+0]->Lines[3];
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+0]->Lines[1]=mesh->Volus[NV+i*7+0]->Lines[4];
    
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+1]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+1]->Lines[2]=mesh->Volus[NV+i*7+1]->Lines[3];
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+1]->Lines[0]=mesh->Volus[NV+i*7+1]->Lines[5];

    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+2]->Lines[2]=NE*2+mesh->Volus[i]->Faces[0]*3+2;
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+2]->Lines[0]=mesh->Volus[NV+i*7+2]->Lines[5];
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+2]->Lines[1]=mesh->Volus[NV+i*7+2]->Lines[4];

    mesh->Faces[mesh->Volus[i]->Faces[0]]->Lines[0]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Faces[mesh->Volus[i]->Faces[0]]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Faces[mesh->Volus[i]->Faces[0]]->Lines[2]=NE*2+mesh->Volus[i]->Faces[0]*3+2;

    mesh->Volus[NV+i*7+0]->Faces[0]=NF+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Volus[NV+i*7+1]->Faces[0]=NF+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Volus[NV+i*7+2]->Faces[0]=NF+mesh->Volus[i]->Faces[0]*3+2;


   } 
 else if (mesh->Faces[mesh->Volus[i]->Faces[0]]->Verts[1]==mesh->Volus[NV+i*7+0]->Verts[3]
          &&mesh->Faces[mesh->Volus[i]->Faces[0]]->Verts[0]==mesh->Volus[NV+i*7+0]->Verts[2] )
   {
    mesh->Volus[NV+i*7+0]->Lines[5]=NE*2+mesh->Volus[i]->Faces[0]*3+2;
    mesh->Volus[NV+i*7+1]->Lines[4]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Volus[NV+i*7+2]->Lines[3]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Volus[NV+i*7+4]->Lines[5]=NE*2+mesh->Volus[i]->Faces[0]*3+2;
    mesh->Volus[NV+i*7+6]->Lines[5]=NE*2+mesh->Volus[i]->Faces[0]*3+2;
    mesh->Volus[NV+i*7+6]->Lines[2]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Volus[NV+i*7+5]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Volus[NV+i*7+6]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
 
    mesh->Faces[NF*4+i*8+1]->Lines[0]=NE*2+mesh->Volus[i]->Faces[0]*3+2; 
    mesh->Faces[NF*4+i*8+5]->Lines[0]=NE*2+mesh->Volus[i]->Faces[0]*3+2;
    mesh->Faces[NF*4+i*8+2]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Faces[NF*4+i*8+7]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Faces[NF*4+i*8+3]->Lines[2]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
  
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+0]->Lines[0]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+0]->Lines[2]=mesh->Volus[NV+i*7+2]->Lines[5];
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+0]->Lines[1]=mesh->Volus[NV+i*7+2]->Lines[4];
    
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+1]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+1]->Lines[2]=mesh->Volus[NV+i*7+1]->Lines[5];
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+1]->Lines[0]=mesh->Volus[NV+i*7+1]->Lines[3];

    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+2]->Lines[2]=NE*2+mesh->Volus[i]->Faces[0]*3+2;
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+2]->Lines[0]=mesh->Volus[NV+i*7+0]->Lines[3];
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+2]->Lines[1]=mesh->Volus[NV+i*7+0]->Lines[4];

    mesh->Faces[mesh->Volus[i]->Faces[0]]->Lines[0]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Faces[mesh->Volus[i]->Faces[0]]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Faces[mesh->Volus[i]->Faces[0]]->Lines[2]=NE*2+mesh->Volus[i]->Faces[0]*3+2;
 
    mesh->Volus[NV+i*7+0]->Faces[0]=NF+mesh->Volus[i]->Faces[0]*3+2;
    mesh->Volus[NV+i*7+1]->Faces[0]=NF+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Volus[NV+i*7+2]->Faces[0]=NF+mesh->Volus[i]->Faces[0]*3+0;


  } 
 else if (mesh->Faces[mesh->Volus[i]->Faces[0]]->Verts[2]==mesh->Volus[NV+i*7+0]->Verts[3]
          &&mesh->Faces[mesh->Volus[i]->Faces[0]]->Verts[0]==mesh->Volus[NV+i*7+0]->Verts[2] )
   {
    mesh->Volus[NV+i*7+0]->Lines[5]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Volus[NV+i*7+1]->Lines[4]=NE*2+mesh->Volus[i]->Faces[0]*3+2;
    mesh->Volus[NV+i*7+2]->Lines[3]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Volus[NV+i*7+4]->Lines[5]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Volus[NV+i*7+6]->Lines[5]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Volus[NV+i*7+6]->Lines[2]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Volus[NV+i*7+5]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+2;
    mesh->Volus[NV+i*7+6]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+2;
    
    mesh->Faces[NF*4+i*8+1]->Lines[0]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Faces[NF*4+i*8+5]->Lines[0]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Faces[NF*4+i*8+2]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+2;
    mesh->Faces[NF*4+i*8+7]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+2;
    mesh->Faces[NF*4+i*8+3]->Lines[2]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
  
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+0]->Lines[0]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+0]->Lines[2]=mesh->Volus[NV+i*7+2]->Lines[4];
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+0]->Lines[1]=mesh->Volus[NV+i*7+2]->Lines[5];
    
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+1]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+1]->Lines[2]=mesh->Volus[NV+i*7+0]->Lines[4];
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+1]->Lines[0]=mesh->Volus[NV+i*7+0]->Lines[3];

    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+2]->Lines[2]=NE*2+mesh->Volus[i]->Faces[0]*3+2;
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+2]->Lines[0]=mesh->Volus[NV+i*7+1]->Lines[3];
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+2]->Lines[1]=mesh->Volus[NV+i*7+1]->Lines[5];

    mesh->Faces[mesh->Volus[i]->Faces[0]]->Lines[0]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Faces[mesh->Volus[i]->Faces[0]]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Faces[mesh->Volus[i]->Faces[0]]->Lines[2]=NE*2+mesh->Volus[i]->Faces[0]*3+2;

    mesh->Volus[NV+i*7+0]->Faces[0]=NF+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Volus[NV+i*7+1]->Faces[0]=NF+mesh->Volus[i]->Faces[0]*3+2;
    mesh->Volus[NV+i*7+2]->Faces[0]=NF+mesh->Volus[i]->Faces[0]*3+0;
} 
 //else(mesh->Faces[mesh->Volus[i].Faces[0]]->Verts[2]==mesh->Volus[NV+i*7+0]->Verts[3]&&mesh->Faces[mesh->Volus[i]->Faces[0]]->Verts[1]==mesh->Volus[NV+i*7+0]->Verts[2] )
else
   {
    mesh->Volus[NV+i*7+0]->Lines[5]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Volus[NV+i*7+1]->Lines[4]=NE*2+mesh->Volus[i]->Faces[0]*3+2;
    mesh->Volus[NV+i*7+2]->Lines[3]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Volus[NV+i*7+4]->Lines[5]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Volus[NV+i*7+6]->Lines[5]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Volus[NV+i*7+6]->Lines[2]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Volus[NV+i*7+5]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+2;
    mesh->Volus[NV+i*7+6]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+2;
    
    mesh->Faces[NF*4+i*8+1]->Lines[0]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Faces[NF*4+i*8+5]->Lines[0]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Faces[NF*4+i*8+2]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+2;
    mesh->Faces[NF*4+i*8+7]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+2;
    mesh->Faces[NF*4+i*8+3]->Lines[2]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
   
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+0]->Lines[0]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+0]->Lines[2]=mesh->Volus[NV+i*7+0]->Lines[4];
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+0]->Lines[1]=mesh->Volus[NV+i*7+0]->Lines[3];
    
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+1]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+1]->Lines[2]=mesh->Volus[NV+i*7+2]->Lines[4];
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+1]->Lines[0]=mesh->Volus[NV+i*7+2]->Lines[5];

    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+2]->Lines[2]=NE*2+mesh->Volus[i]->Faces[0]*3+2;
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+2]->Lines[0]=mesh->Volus[NV+i*7+1]->Lines[5];
    mesh->Faces[NF+mesh->Volus[i]->Faces[0]*3+2]->Lines[1]=mesh->Volus[NV+i*7+1]->Lines[3];

    mesh->Faces[mesh->Volus[i]->Faces[0]]->Lines[0]=NE*2+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Faces[mesh->Volus[i]->Faces[0]]->Lines[1]=NE*2+mesh->Volus[i]->Faces[0]*3+1;
    mesh->Faces[mesh->Volus[i]->Faces[0]]->Lines[2]=NE*2+mesh->Volus[i]->Faces[0]*3+2;
  
    mesh->Volus[NV+i*7+0]->Faces[0]=NF+mesh->Volus[i]->Faces[0]*3+0;
    mesh->Volus[NV+i*7+1]->Faces[0]=NF+mesh->Volus[i]->Faces[0]*3+2;
    mesh->Volus[NV+i*7+2]->Faces[0]=NF+mesh->Volus[i]->Faces[0]*3+1;
   }
 

  /* 第二个面 */
   if(mesh->Faces[mesh->Volus[i]->Faces[1]]->Verts[0]==mesh->Volus[NV+i*7+1]->Verts[0]
      &&mesh->Faces[mesh->Volus[i]->Faces[1]]->Verts[1]==mesh->Volus[NV+i*7+1]->Verts[3])
   {
    mesh->Volus[i]->Lines[5]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Volus[NV+i*7+1]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+2;
    mesh->Volus[NV+i*7+2]->Lines[1]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Volus[NV+i*7+3]->Lines[3]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Volus[NV+i*7+5]->Lines[4]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Volus[NV+i*7+5]->Lines[0]=NE*2+mesh->Volus[i]->Faces[1]*3+2;
    mesh->Volus[NV+i*7+6]->Lines[0]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Volus[NV+i*7+5]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
  
    mesh->Faces[NF*4+i*8+2]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+2;
    mesh->Faces[NF*4+i*8+0]->Lines[0]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Faces[NF*4+i*8+4]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Faces[NF*4+i*8+7]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Faces[NF*4+i*8+3]->Lines[1]=NE*2+mesh->Volus[i]->Faces[1]*3+0;

 
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+0]->Lines[0]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+0]->Lines[2]=mesh->Volus[NV+i*7+2]->Lines[2];
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+0]->Lines[1]=mesh->Volus[NV+i*7+2]->Lines[5];
    
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+1]->Lines[1]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+1]->Lines[2]=mesh->Volus[i]->Lines[2];
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+1]->Lines[0]=mesh->Volus[i]->Lines[1];

    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+2]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+2;
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+2]->Lines[0]=mesh->Volus[NV+i*7+1]->Lines[1];
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+2]->Lines[1]=mesh->Volus[NV+i*7+1]->Lines[5];

    mesh->Faces[mesh->Volus[i]->Faces[1]]->Lines[0]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Faces[mesh->Volus[i]->Faces[1]]->Lines[1]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Faces[mesh->Volus[i]->Faces[1]]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+2;

    mesh->Volus[NV+i*7+1]->Faces[1]=NF+mesh->Volus[i]->Faces[1]*3+2;
    mesh->Volus[NV+i*7+2]->Faces[1]=NF+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Volus[i]->Faces[1]=NF+mesh->Volus[i]->Faces[1]*3+1;

 }  
   else if(mesh->Faces[mesh->Volus[i]->Faces[1]]->Verts[1]==mesh->Volus[NV+i*7+1]->Verts[0]
	   &&mesh->Faces[mesh->Volus[i]->Faces[1]]->Verts[0]==mesh->Volus[NV+i*7+1]->Verts[3])
   {
    mesh->Volus[i]->Lines[5]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Volus[NV+i*7+1]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+2;
    mesh->Volus[NV+i*7+2]->Lines[1]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Volus[NV+i*7+3]->Lines[3]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Volus[NV+i*7+5]->Lines[4]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Volus[NV+i*7+5]->Lines[0]=NE*2+mesh->Volus[i]->Faces[1]*3+2;
    mesh->Volus[NV+i*7+6]->Lines[0]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Volus[NV+i*7+5]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
   
    mesh->Faces[NF*4+i*8+2]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+2;
    mesh->Faces[NF*4+i*8+0]->Lines[0]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Faces[NF*4+i*8+4]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Faces[NF*4+i*8+7]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Faces[NF*4+i*8+3]->Lines[1]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
 
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+0]->Lines[0]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+0]->Lines[2]=mesh->Volus[i]->Lines[2];
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+0]->Lines[1]=mesh->Volus[i]->Lines[1];
    
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+1]->Lines[1]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+1]->Lines[2]=mesh->Volus[NV+i*7+2]->Lines[2];
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+1]->Lines[0]=mesh->Volus[NV+i*7+2]->Lines[5];

    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+2]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+2;
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+2]->Lines[0]=mesh->Volus[NV+i*7+1]->Lines[5];
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+2]->Lines[1]=mesh->Volus[NV+i*7+1]->Lines[1];

    mesh->Faces[mesh->Volus[i]->Faces[1]]->Lines[0]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Faces[mesh->Volus[i]->Faces[1]]->Lines[1]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Faces[mesh->Volus[i]->Faces[1]]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+2;

    mesh->Volus[NV+i*7+1]->Faces[1]=NF+mesh->Volus[i]->Faces[1]*3+2;
    mesh->Volus[NV+i*7+2]->Faces[1]=NF+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Volus[i]->Faces[1]=NF+mesh->Volus[i]->Faces[1]*3+0;



}  
   else if(mesh->Faces[mesh->Volus[i]->Faces[1]]->Verts[2]==mesh->Volus[NV+i*7+1]->Verts[0]
	   &&mesh->Faces[mesh->Volus[i]->Faces[1]]->Verts[0]==mesh->Volus[NV+i*7+1]->Verts[3])
   {
    mesh->Volus[i]->Lines[5]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Volus[NV+i*7+1]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Volus[NV+i*7+2]->Lines[1]=NE*2+mesh->Volus[i]->Faces[1]*3+2;
    mesh->Volus[NV+i*7+3]->Lines[3]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Volus[NV+i*7+5]->Lines[4]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Volus[NV+i*7+5]->Lines[0]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Volus[NV+i*7+6]->Lines[0]=NE*2+mesh->Volus[i]->Faces[1]*3+2;
    mesh->Volus[NV+i*7+5]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+2;
  
    mesh->Faces[NF*4+i*8+2]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Faces[NF*4+i*8+0]->Lines[0]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Faces[NF*4+i*8+4]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Faces[NF*4+i*8+7]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+2;
    mesh->Faces[NF*4+i*8+3]->Lines[1]=NE*2+mesh->Volus[i]->Faces[1]*3+2;
 
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+0]->Lines[0]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+0]->Lines[2]=mesh->Volus[i]->Lines[1];
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+0]->Lines[1]=mesh->Volus[i]->Lines[2];
    
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+1]->Lines[1]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+1]->Lines[2]=mesh->Volus[NV+i*7+1]->Lines[1];
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+1]->Lines[0]=mesh->Volus[NV+i*7+1]->Lines[5];

    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+2]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+2;
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+2]->Lines[0]=mesh->Volus[NV+i*7+2]->Lines[5];
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+2]->Lines[1]=mesh->Volus[NV+i*7+2]->Lines[2];

    mesh->Faces[mesh->Volus[i]->Faces[1]]->Lines[0]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Faces[mesh->Volus[i]->Faces[1]]->Lines[1]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Faces[mesh->Volus[i]->Faces[1]]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+2;
 
    mesh->Volus[NV+i*7+1]->Faces[1]=NF+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Volus[NV+i*7+2]->Faces[1]=NF+mesh->Volus[i]->Faces[1]*3+2;
    mesh->Volus[i]->Faces[1]=NF+mesh->Volus[i]->Faces[1]*3+0;



 }  
   else if(mesh->Faces[mesh->Volus[i]->Faces[1]]->Verts[0]==mesh->Volus[NV+i*7+1]->Verts[0]
	   &&mesh->Faces[mesh->Volus[i]->Faces[1]]->Verts[2]==mesh->Volus[NV+i*7+1]->Verts[3])
   {
    mesh->Volus[i]->Lines[5]=NE*2+mesh->Volus[i]->Faces[1]*3+2;
    mesh->Volus[NV+i*7+1]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Volus[NV+i*7+2]->Lines[1]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Volus[NV+i*7+3]->Lines[3]=NE*2+mesh->Volus[i]->Faces[1]*3+2;
    mesh->Volus[NV+i*7+5]->Lines[4]=NE*2+mesh->Volus[i]->Faces[1]*3+2;
    mesh->Volus[NV+i*7+5]->Lines[0]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Volus[NV+i*7+6]->Lines[0]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Volus[NV+i*7+5]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
  
    mesh->Faces[NF*4+i*8+2]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Faces[NF*4+i*8+0]->Lines[0]=NE*2+mesh->Volus[i]->Faces[1]*3+2;
    mesh->Faces[NF*4+i*8+4]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+2;
    mesh->Faces[NF*4+i*8+7]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Faces[NF*4+i*8+3]->Lines[1]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
  
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+0]->Lines[0]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+0]->Lines[2]=mesh->Volus[NV+i*7+2]->Lines[5];
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+0]->Lines[1]=mesh->Volus[NV+i*7+2]->Lines[2];
    
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+1]->Lines[1]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+1]->Lines[2]=mesh->Volus[NV+i*7+1]->Lines[5];
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+1]->Lines[0]=mesh->Volus[NV+i*7+1]->Lines[1];

    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+2]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+2;
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+2]->Lines[0]=mesh->Volus[i]->Lines[1];
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+2]->Lines[1]=mesh->Volus[i]->Lines[2];

    mesh->Faces[mesh->Volus[i]->Faces[1]]->Lines[0]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Faces[mesh->Volus[i]->Faces[1]]->Lines[1]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Faces[mesh->Volus[i]->Faces[1]]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+2;

    mesh->Volus[NV+i*7+1]->Faces[1]=NF+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Volus[NV+i*7+2]->Faces[1]=NF+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Volus[i]->Faces[1]=NF+mesh->Volus[i]->Faces[1]*3+2;


}  
   else if(mesh->Faces[mesh->Volus[i]->Faces[1]]->Verts[1]==mesh->Volus[NV+i*7+1]->Verts[0]
	   &&mesh->Faces[mesh->Volus[i]->Faces[1]]->Verts[2]==mesh->Volus[NV+i*7+1]->Verts[3])
   {
    mesh->Volus[i]->Lines[5]=NE*2+mesh->Volus[i]->Faces[1]*3+2;
    mesh->Volus[NV+i*7+1]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Volus[NV+i*7+2]->Lines[1]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Volus[NV+i*7+3]->Lines[3]=NE*2+mesh->Volus[i]->Faces[1]*3+2;
    mesh->Volus[NV+i*7+5]->Lines[4]=NE*2+mesh->Volus[i]->Faces[1]*3+2;
    mesh->Volus[NV+i*7+5]->Lines[0]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Volus[NV+i*7+6]->Lines[0]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Volus[NV+i*7+5]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+1;

    mesh->Faces[NF*4+i*8+2]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Faces[NF*4+i*8+0]->Lines[0]=NE*2+mesh->Volus[i]->Faces[1]*3+2;
    mesh->Faces[NF*4+i*8+4]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+2;
    mesh->Faces[NF*4+i*8+7]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Faces[NF*4+i*8+3]->Lines[1]=NE*2+mesh->Volus[i]->Faces[1]*3+1;

    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+0]->Lines[0]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+0]->Lines[2]=mesh->Volus[NV+i*7+1]->Lines[5];
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+0]->Lines[1]=mesh->Volus[NV+i*7+1]->Lines[1];
    
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+1]->Lines[1]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+1]->Lines[2]=mesh->Volus[NV+i*7+2]->Lines[5];
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+1]->Lines[0]=mesh->Volus[NV+i*7+2]->Lines[2];

    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+2]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+2;
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+2]->Lines[0]=mesh->Volus[i]->Lines[2];
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+2]->Lines[1]=mesh->Volus[i]->Lines[1];

    mesh->Faces[mesh->Volus[i]->Faces[1]]->Lines[0]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Faces[mesh->Volus[i]->Faces[1]]->Lines[1]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Faces[mesh->Volus[i]->Faces[1]]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+2;

    mesh->Volus[NV+i*7+1]->Faces[1]=NF+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Volus[NV+i*7+2]->Faces[1]=NF+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Volus[i]->Faces[1]=NF+mesh->Volus[i]->Faces[1]*3+2;


 }  
   //else (mesh->Faces[mesh->Volus[i]->Faces[1]]->Verts[2]==mesh->Volus[NV+i*7+1]->Verts[0]&&mesh->Faces[mesh->Volus[i]->Faces[1]]->Verts[1]==mesh->Volus[NV+i*7+1]->Verts[3])
else   
  {
    mesh->Volus[i]->Lines[5]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Volus[NV+i*7+1]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Volus[NV+i*7+2]->Lines[1]=NE*2+mesh->Volus[i]->Faces[1]*3+2;
    mesh->Volus[NV+i*7+3]->Lines[3]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Volus[NV+i*7+5]->Lines[4]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Volus[NV+i*7+5]->Lines[0]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Volus[NV+i*7+6]->Lines[0]=NE*2+mesh->Volus[i]->Faces[1]*3+2;
    mesh->Volus[NV+i*7+5]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+2;
   
    mesh->Faces[NF*4+i*8+2]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Faces[NF*4+i*8+0]->Lines[0]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Faces[NF*4+i*8+4]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Faces[NF*4+i*8+7]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+2;
    mesh->Faces[NF*4+i*8+3]->Lines[1]=NE*2+mesh->Volus[i]->Faces[1]*3+2;

  
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+0]->Lines[0]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+0]->Lines[2]=mesh->Volus[NV+i*7+1]->Lines[1];
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+0]->Lines[1]=mesh->Volus[NV+i*7+1]->Lines[5];
    
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+1]->Lines[1]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+1]->Lines[2]=mesh->Volus[i]->Lines[1];
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+1]->Lines[0]=mesh->Volus[i]->Lines[2];

    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+2]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+2;
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+2]->Lines[0]=mesh->Volus[NV+i*7+2]->Lines[2];
    mesh->Faces[NF+mesh->Volus[i]->Faces[1]*3+2]->Lines[1]=mesh->Volus[NV+i*7+2]->Lines[5];

    mesh->Faces[mesh->Volus[i]->Faces[1]]->Lines[0]=NE*2+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Faces[mesh->Volus[i]->Faces[1]]->Lines[1]=NE*2+mesh->Volus[i]->Faces[1]*3+1;
    mesh->Faces[mesh->Volus[i]->Faces[1]]->Lines[2]=NE*2+mesh->Volus[i]->Faces[1]*3+2;

    mesh->Volus[NV+i*7+1]->Faces[1]=NF+mesh->Volus[i]->Faces[1]*3+0;
    mesh->Volus[NV+i*7+2]->Faces[1]=NF+mesh->Volus[i]->Faces[1]*3+2;
    mesh->Volus[i]->Faces[1]=NF+mesh->Volus[i]->Faces[1]*3+1;
}


  /* 第三个面 */
   if(mesh->Faces[mesh->Volus[i]->Faces[2]]->Verts[0]==mesh->Volus[NV+i*7+0]->Verts[0]
      &&mesh->Faces[mesh->Volus[i]->Faces[2]]->Verts[1]==mesh->Volus[NV+i*7+0]->Verts[3])
   {
    mesh->Volus[i]->Lines[4]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Volus[NV+i*7+0]->Lines[2]=NE*2+mesh->Volus[i]->Faces[2]*3+2;
    mesh->Volus[NV+i*7+2]->Lines[0]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Volus[NV+i*7+4]->Lines[3]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Volus[NV+i*7+6]->Lines[4]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Volus[NV+i*7+4]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+2;
    mesh->Volus[NV+i*7+4]->Lines[0]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Volus[NV+i*7+3]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
  
    mesh->Faces[NF*4+i*8+1]->Lines[2]=NE*2+mesh->Volus[i]->Faces[2]*3+2;
    mesh->Faces[NF*4+i*8+0]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Faces[NF*4+i*8+6]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Faces[NF*4+i*8+3]->Lines[0]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Faces[NF*4+i*8+5]->Lines[2]=NE*2+mesh->Volus[i]->Faces[2]*3+0;

   
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+0]->Lines[0]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+0]->Lines[2]=mesh->Volus[NV+i*7+2]->Lines[2];
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+0]->Lines[1]=mesh->Volus[NV+i*7+2]->Lines[4];
    
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+1]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+1]->Lines[2]=mesh->Volus[i]->Lines[2];
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+1]->Lines[0]=mesh->Volus[i]->Lines[0];

    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+2]->Lines[2]=NE*2+mesh->Volus[i]->Faces[2]*3+2;
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+2]->Lines[0]=mesh->Volus[NV+i*7+0]->Lines[0];
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+2]->Lines[1]=mesh->Volus[NV+i*7+0]->Lines[4];

    mesh->Faces[mesh->Volus[i]->Faces[2]]->Lines[0]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Faces[mesh->Volus[i]->Faces[2]]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Faces[mesh->Volus[i]->Faces[2]]->Lines[2]=NE*2+mesh->Volus[i]->Faces[2]*3+2;

    mesh->Volus[NV+i*7+0]->Faces[2]=NF+mesh->Volus[i]->Faces[2]*3+2;
    mesh->Volus[NV+i*7+2]->Faces[2]=NF+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Volus[i]->Faces[2]=NF+mesh->Volus[i]->Faces[2]*3+1;


 }  
   else if(mesh->Faces[mesh->Volus[i]->Faces[2]]->Verts[0]==mesh->Volus[NV+i*7+0]->Verts[0]
	 &&mesh->Faces[mesh->Volus[i]->Faces[2]]->Verts[2]==mesh->Volus[NV+i*7+0]->Verts[3])
   {
    mesh->Volus[i]->Lines[4]=NE*2+mesh->Volus[i]->Faces[2]*3+2;
    mesh->Volus[NV+i*7+0]->Lines[2]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Volus[NV+i*7+2]->Lines[0]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Volus[NV+i*7+4]->Lines[3]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Volus[NV+i*7+6]->Lines[4]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Volus[NV+i*7+4]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Volus[NV+i*7+4]->Lines[0]=NE*2+mesh->Volus[i]->Faces[2]*3+2;
    mesh->Volus[NV+i*7+3]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+2;
   
    mesh->Faces[NF*4+i*8+1]->Lines[2]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Faces[NF*4+i*8+0]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+2;
    mesh->Faces[NF*4+i*8+6]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+2;
    mesh->Faces[NF*4+i*8+3]->Lines[0]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Faces[NF*4+i*8+5]->Lines[2]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
 
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+0]->Lines[0]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+0]->Lines[2]=mesh->Volus[NV+i*7+2]->Lines[4];
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+0]->Lines[1]=mesh->Volus[NV+i*7+2]->Lines[2];
    
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+1]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+1]->Lines[2]=mesh->Volus[NV+i*7+0]->Lines[4];
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+1]->Lines[0]=mesh->Volus[NV+i*7+0]->Lines[0];

    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+2]->Lines[2]=NE*2+mesh->Volus[i]->Faces[2]*3+2;
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+2]->Lines[0]=mesh->Volus[i]->Lines[0];
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+2]->Lines[1]=mesh->Volus[i]->Lines[2];

    mesh->Faces[mesh->Volus[i]->Faces[2]]->Lines[0]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Faces[mesh->Volus[i]->Faces[2]]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Faces[mesh->Volus[i]->Faces[2]]->Lines[2]=NE*2+mesh->Volus[i]->Faces[2]*3+2;

    mesh->Volus[NV+i*7+0]->Faces[2]=NF+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Volus[NV+i*7+2]->Faces[2]=NF+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Volus[i]->Faces[2]=NF+mesh->Volus[i]->Faces[2]*3+2;



   }  
   else if(mesh->Faces[mesh->Volus[i]->Faces[2]]->Verts[1]==mesh->Volus[NV+i*7+0]->Verts[0]
	 &&mesh->Faces[mesh->Volus[i]->Faces[2]]->Verts[2]==mesh->Volus[NV+i*7+0]->Verts[3])
   {
    mesh->Volus[i]->Lines[4]=NE*2+mesh->Volus[i]->Faces[2]*3+2;
    mesh->Volus[NV+i*7+0]->Lines[2]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Volus[NV+i*7+2]->Lines[0]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Volus[NV+i*7+4]->Lines[3]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Volus[NV+i*7+6]->Lines[4]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Volus[NV+i*7+4]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Volus[NV+i*7+4]->Lines[0]=NE*2+mesh->Volus[i]->Faces[2]*3+2;
    mesh->Volus[NV+i*7+3]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+2;
   
    mesh->Faces[NF*4+i*8+1]->Lines[2]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Faces[NF*4+i*8+0]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+2;
    mesh->Faces[NF*4+i*8+6]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+2;
    mesh->Faces[NF*4+i*8+3]->Lines[0]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Faces[NF*4+i*8+5]->Lines[2]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
 
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+0]->Lines[0]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+0]->Lines[2]=mesh->Volus[NV+i*7+0]->Lines[4];
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+0]->Lines[1]=mesh->Volus[NV+i*7+0]->Lines[0];
    
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+1]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+1]->Lines[2]=mesh->Volus[NV+i*7+2]->Lines[4];
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+1]->Lines[0]=mesh->Volus[NV+i*7+2]->Lines[2];

    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+2]->Lines[2]=NE*2+mesh->Volus[i]->Faces[2]*3+2;
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+2]->Lines[0]=mesh->Volus[i]->Lines[2];
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+2]->Lines[1]=mesh->Volus[i]->Lines[0];

    mesh->Faces[mesh->Volus[i]->Faces[2]]->Lines[0]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Faces[mesh->Volus[i]->Faces[2]]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Faces[mesh->Volus[i]->Faces[2]]->Lines[2]=NE*2+mesh->Volus[i]->Faces[2]*3+2;

    mesh->Volus[NV+i*7+0]->Faces[2]=NF+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Volus[NV+i*7+2]->Faces[2]=NF+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Volus[i]->Faces[2]=NF+mesh->Volus[i]->Faces[2]*3+2;



   }  
   else if(mesh->Faces[mesh->Volus[i]->Faces[2]]->Verts[1]==mesh->Volus[NV+i*7+0]->Verts[0]
	 &&mesh->Faces[mesh->Volus[i]->Faces[2]]->Verts[0]==mesh->Volus[NV+i*7+0]->Verts[3])
   {
    mesh->Volus[i]->Lines[4]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Volus[NV+i*7+0]->Lines[2]=NE*2+mesh->Volus[i]->Faces[2]*3+2;
    mesh->Volus[NV+i*7+2]->Lines[0]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Volus[NV+i*7+4]->Lines[3]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Volus[NV+i*7+6]->Lines[4]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Volus[NV+i*7+4]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+2;
    mesh->Volus[NV+i*7+4]->Lines[0]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Volus[NV+i*7+3]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    
    mesh->Faces[NF*4+i*8+1]->Lines[2]=NE*2+mesh->Volus[i]->Faces[2]*3+2;
    mesh->Faces[NF*4+i*8+0]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Faces[NF*4+i*8+6]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Faces[NF*4+i*8+3]->Lines[0]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Faces[NF*4+i*8+5]->Lines[2]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
 
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+0]->Lines[0]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+0]->Lines[2]=mesh->Volus[i]->Lines[2];
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+0]->Lines[1]=mesh->Volus[i]->Lines[0];
    
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+1]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+1]->Lines[2]=mesh->Volus[NV+i*7+2]->Lines[2];
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+1]->Lines[0]=mesh->Volus[NV+i*7+2]->Lines[4];

    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+2]->Lines[2]=NE*2+mesh->Volus[i]->Faces[2]*3+2;
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+2]->Lines[0]=mesh->Volus[NV+i*7+0]->Lines[4];
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+2]->Lines[1]=mesh->Volus[NV+i*7+0]->Lines[0];

    mesh->Faces[mesh->Volus[i]->Faces[2]]->Lines[0]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Faces[mesh->Volus[i]->Faces[2]]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Faces[mesh->Volus[i]->Faces[2]]->Lines[2]=NE*2+mesh->Volus[i]->Faces[2]*3+2;

    mesh->Volus[NV+i*7+0]->Faces[2]=NF+mesh->Volus[i]->Faces[2]*3+2;
    mesh->Volus[NV+i*7+2]->Faces[2]=NF+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Volus[i]->Faces[2]=NF+mesh->Volus[i]->Faces[2]*3+0;


   }  
   else if(mesh->Faces[mesh->Volus[i]->Faces[2]]->Verts[2]==mesh->Volus[NV+i*7+0]->Verts[0]
	 &&mesh->Faces[mesh->Volus[i]->Faces[2]]->Verts[0]==mesh->Volus[NV+i*7+0]->Verts[3])
   {
    mesh->Volus[i]->Lines[4]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Volus[NV+i*7+0]->Lines[2]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Volus[NV+i*7+2]->Lines[0]=NE*2+mesh->Volus[i]->Faces[2]*3+2;
    mesh->Volus[NV+i*7+4]->Lines[3]=NE*2+mesh->Volus[i]->Faces[2]*3+2;
    mesh->Volus[NV+i*7+6]->Lines[4]=NE*2+mesh->Volus[i]->Faces[2]*3+2;
    mesh->Volus[NV+i*7+4]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Volus[NV+i*7+4]->Lines[0]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Volus[NV+i*7+3]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    
    mesh->Faces[NF*4+i*8+1]->Lines[2]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Faces[NF*4+i*8+0]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Faces[NF*4+i*8+6]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Faces[NF*4+i*8+3]->Lines[0]=NE*2+mesh->Volus[i]->Faces[2]*3+2;
    mesh->Faces[NF*4+i*8+5]->Lines[2]=NE*2+mesh->Volus[i]->Faces[2]*3+2;
 
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+0]->Lines[0]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+0]->Lines[2]=mesh->Volus[i]->Lines[0];
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+0]->Lines[1]=mesh->Volus[i]->Lines[2];
    
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+1]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+1]->Lines[2]=mesh->Volus[NV+i*7+0]->Lines[0];
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+1]->Lines[0]=mesh->Volus[NV+i*7+0]->Lines[4];

    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+2]->Lines[2]=NE*2+mesh->Volus[i]->Faces[2]*3+2;
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+2]->Lines[0]=mesh->Volus[NV+i*7+2]->Lines[4];
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+2]->Lines[1]=mesh->Volus[NV+i*7+2]->Lines[2];

    mesh->Faces[mesh->Volus[i]->Faces[2]]->Lines[0]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Faces[mesh->Volus[i]->Faces[2]]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Faces[mesh->Volus[i]->Faces[2]]->Lines[2]=NE*2+mesh->Volus[i]->Faces[2]*3+2;

    mesh->Volus[NV+i*7+0]->Faces[2]=NF+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Volus[NV+i*7+2]->Faces[2]=NF+mesh->Volus[i]->Faces[2]*3+2;
    mesh->Volus[i]->Faces[2]=NF+mesh->Volus[i]->Faces[2]*3+0;


   }  
   //else (mesh->Faces[mesh->Volus[i]->Faces[2]]->Verts[2]==mesh->Volus[NV+i*7+0]->Verts[0]&&mesh->Faces[mesh->Volus[i]->Faces[1]]->Verts[1]==mesh->Volus[NV+i*7+0]->Verts[3])
else
   {
    mesh->Volus[i]->Lines[4]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Volus[NV+i*7+0]->Lines[2]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Volus[NV+i*7+2]->Lines[0]=NE*2+mesh->Volus[i]->Faces[2]*3+2;
    mesh->Volus[NV+i*7+4]->Lines[3]=NE*2+mesh->Volus[i]->Faces[2]*3+2;
    mesh->Volus[NV+i*7+6]->Lines[4]=NE*2+mesh->Volus[i]->Faces[2]*3+2;
    mesh->Volus[NV+i*7+4]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Volus[NV+i*7+4]->Lines[0]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Volus[NV+i*7+3]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+1;

    mesh->Faces[NF*4+i*8+1]->Lines[2]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Faces[NF*4+i*8+0]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Faces[NF*4+i*8+6]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Faces[NF*4+i*8+3]->Lines[0]=NE*2+mesh->Volus[i]->Faces[2]*3+2;
    mesh->Faces[NF*4+i*8+5]->Lines[2]=NE*2+mesh->Volus[i]->Faces[2]*3+2;
 
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+0]->Lines[0]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+0]->Lines[2]=mesh->Volus[NV+i*7+0]->Lines[0];
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+0]->Lines[1]=mesh->Volus[NV+i*7+0]->Lines[4];
    
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+1]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+1]->Lines[2]=mesh->Volus[i]->Lines[0];
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+1]->Lines[0]=mesh->Volus[i]->Lines[2];

    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+2]->Lines[2]=NE*2+mesh->Volus[i]->Faces[2]*3+2;
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+2]->Lines[0]=mesh->Volus[NV+i*7+2]->Lines[2];
    mesh->Faces[NF+mesh->Volus[i]->Faces[2]*3+2]->Lines[1]=mesh->Volus[NV+i*7+2]->Lines[4];

    mesh->Faces[mesh->Volus[i]->Faces[2]]->Lines[0]=NE*2+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Faces[mesh->Volus[i]->Faces[2]]->Lines[1]=NE*2+mesh->Volus[i]->Faces[2]*3+1;
    mesh->Faces[mesh->Volus[i]->Faces[2]]->Lines[2]=NE*2+mesh->Volus[i]->Faces[2]*3+2;

    mesh->Volus[NV+i*7+0]->Faces[2]=NF+mesh->Volus[i]->Faces[2]*3+0;
    mesh->Volus[NV+i*7+2]->Faces[2]=NF+mesh->Volus[i]->Faces[2]*3+2;
    mesh->Volus[i]->Faces[2]=NF+mesh->Volus[i]->Faces[2]*3+1;



   }  
 /* 第四个面 */
   if(mesh->Faces[mesh->Volus[i]->Faces[3]]->Verts[0]==mesh->Volus[NV+i*7+0]->Verts[0]
     &&mesh->Faces[mesh->Volus[i]->Faces[3]]->Verts[1]==mesh->Volus[NV+i*7+0]->Verts[2])
   {
    mesh->Volus[i]->Lines[3]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Volus[NV+i*7+0]->Lines[1]=NE*2+mesh->Volus[i]->Faces[3]*3+2;
    mesh->Volus[NV+i*7+1]->Lines[0]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Volus[NV+i*7+3]->Lines[4]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Volus[NV+i*7+5]->Lines[3]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Volus[NV+i*7+3]->Lines[0]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Volus[NV+i*7+3]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+2;
    mesh->Volus[NV+i*7+4]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+2;
   
    mesh->Faces[NF*4+i*8+1]->Lines[1]=NE*2+mesh->Volus[i]->Faces[3]*3+2;
    mesh->Faces[NF*4+i*8+0]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Faces[NF*4+i*8+6]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+2;
    mesh->Faces[NF*4+i*8+2]->Lines[0]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Faces[NF*4+i*8+4]->Lines[1]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
 
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+0]->Lines[0]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+0]->Lines[2]=mesh->Volus[NV+i*7+1]->Lines[1];
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+0]->Lines[1]=mesh->Volus[NV+i*7+1]->Lines[3];
    
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+1]->Lines[1]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+1]->Lines[2]=mesh->Volus[i]->Lines[1];
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+1]->Lines[0]=mesh->Volus[i]->Lines[0];

    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+2]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+2;
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+2]->Lines[0]=mesh->Volus[NV+i*7+0]->Lines[0];
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+2]->Lines[1]=mesh->Volus[NV+i*7+0]->Lines[3];

    mesh->Faces[mesh->Volus[i]->Faces[3]]->Lines[0]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Faces[mesh->Volus[i]->Faces[3]]->Lines[1]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Faces[mesh->Volus[i]->Faces[3]]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+2;

    mesh->Volus[NV+i*7+0]->Faces[3]=NF+mesh->Volus[i]->Faces[3]*3+2;
    mesh->Volus[NV+i*7+1]->Faces[3]=NF+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Volus[i]->Faces[3]=NF+mesh->Volus[i]->Faces[3]*3+1;



}  
   else if(mesh->Faces[mesh->Volus[i]->Faces[3]]->Verts[0]==mesh->Volus[NV+i*7+0]->Verts[0]
	 &&mesh->Faces[mesh->Volus[i]->Faces[3]]->Verts[2]==mesh->Volus[NV+i*7+0]->Verts[2])
   {
    mesh->Volus[i]->Lines[3]=NE*2+mesh->Volus[i]->Faces[3]*3+2;
    mesh->Volus[NV+i*7+0]->Lines[1]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Volus[NV+i*7+1]->Lines[0]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Volus[NV+i*7+3]->Lines[4]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Volus[NV+i*7+5]->Lines[3]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Volus[NV+i*7+3]->Lines[0]=NE*2+mesh->Volus[i]->Faces[3]*3+2;
    mesh->Volus[NV+i*7+3]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Volus[NV+i*7+4]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
   
    mesh->Faces[NF*4+i*8+1]->Lines[1]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Faces[NF*4+i*8+0]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+2;
    mesh->Faces[NF*4+i*8+6]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Faces[NF*4+i*8+2]->Lines[0]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Faces[NF*4+i*8+4]->Lines[1]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
 
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+0]->Lines[0]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+0]->Lines[2]=mesh->Volus[NV+i*7+1]->Lines[3];
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+0]->Lines[1]=mesh->Volus[NV+i*7+1]->Lines[1];
    
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+1]->Lines[1]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+1]->Lines[2]=mesh->Volus[NV+i*7+0]->Lines[3];
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+1]->Lines[0]=mesh->Volus[NV+i*7+0]->Lines[0];

    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+2]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+2;
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+2]->Lines[0]=mesh->Volus[i]->Lines[0];
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+2]->Lines[1]=mesh->Volus[i]->Lines[1];

    mesh->Faces[mesh->Volus[i]->Faces[3]]->Lines[0]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Faces[mesh->Volus[i]->Faces[3]]->Lines[1]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Faces[mesh->Volus[i]->Faces[3]]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+2;

    mesh->Volus[NV+i*7+0]->Faces[3]=NF+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Volus[NV+i*7+1]->Faces[3]=NF+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Volus[i]->Faces[3]=NF+mesh->Volus[i]->Faces[3]*3+2;


}  
  else if(mesh->Faces[mesh->Volus[i]->Faces[3]]->Verts[1]==mesh->Volus[NV+i*7+0]->Verts[0]
	&&mesh->Faces[mesh->Volus[i]->Faces[3]]->Verts[2]==mesh->Volus[NV+i*7+0]->Verts[2])
   {
    mesh->Volus[i]->Lines[3]=NE*2+mesh->Volus[i]->Faces[3]*3+2;
    mesh->Volus[NV+i*7+0]->Lines[1]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Volus[NV+i*7+1]->Lines[0]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Volus[NV+i*7+3]->Lines[4]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Volus[NV+i*7+5]->Lines[3]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Volus[NV+i*7+3]->Lines[0]=NE*2+mesh->Volus[i]->Faces[3]*3+2;
    mesh->Volus[NV+i*7+3]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Volus[NV+i*7+4]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
   
    mesh->Faces[NF*4+i*8+1]->Lines[1]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Faces[NF*4+i*8+0]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+2;
    mesh->Faces[NF*4+i*8+6]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Faces[NF*4+i*8+2]->Lines[0]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Faces[NF*4+i*8+4]->Lines[1]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
 
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+0]->Lines[0]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+0]->Lines[2]=mesh->Volus[NV+i*7+0]->Lines[3];
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+0]->Lines[1]=mesh->Volus[NV+i*7+0]->Lines[0];
    
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+1]->Lines[1]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+1]->Lines[2]=mesh->Volus[NV+i*7+1]->Lines[3];
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+1]->Lines[0]=mesh->Volus[NV+i*7+1]->Lines[1];

    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+2]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+2;
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+2]->Lines[0]=mesh->Volus[i]->Lines[1];
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+2]->Lines[1]=mesh->Volus[i]->Lines[0];

    mesh->Faces[mesh->Volus[i]->Faces[3]]->Lines[0]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Faces[mesh->Volus[i]->Faces[3]]->Lines[1]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Faces[mesh->Volus[i]->Faces[3]]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+2;

    mesh->Volus[NV+i*7+0]->Faces[3]=NF+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Volus[NV+i*7+1]->Faces[3]=NF+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Volus[i]->Faces[3]=NF+mesh->Volus[i]->Faces[3]*3+2;


   }  
  else if(mesh->Faces[mesh->Volus[i]->Faces[3]]->Verts[1]==mesh->Volus[NV+i*7+0]->Verts[0]
	&&mesh->Faces[mesh->Volus[i]->Faces[3]]->Verts[0]==mesh->Volus[NV+i*7+0]->Verts[2])
   {
    mesh->Volus[i]->Lines[3]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Volus[NV+i*7+0]->Lines[1]=NE*2+mesh->Volus[i]->Faces[3]*3+2;
    mesh->Volus[NV+i*7+1]->Lines[0]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Volus[NV+i*7+3]->Lines[4]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Volus[NV+i*7+5]->Lines[3]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Volus[NV+i*7+3]->Lines[0]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Volus[NV+i*7+3]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+2;
    mesh->Volus[NV+i*7+4]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+2;
  
    mesh->Faces[NF*4+i*8+1]->Lines[1]=NE*2+mesh->Volus[i]->Faces[3]*3+2;
    mesh->Faces[NF*4+i*8+0]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Faces[NF*4+i*8+6]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+2;
    mesh->Faces[NF*4+i*8+2]->Lines[0]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Faces[NF*4+i*8+4]->Lines[1]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
 
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+0]->Lines[0]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+0]->Lines[2]=mesh->Volus[i]->Lines[1];
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+0]->Lines[1]=mesh->Volus[i]->Lines[0];
    
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+1]->Lines[1]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+1]->Lines[2]=mesh->Volus[NV+i*7+1]->Lines[1];
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+1]->Lines[0]=mesh->Volus[NV+i*7+1]->Lines[3];

    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+2]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+2;
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+2]->Lines[0]=mesh->Volus[NV+i*7+0]->Lines[3];
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+2]->Lines[1]=mesh->Volus[NV+i*7+0]->Lines[0];

    mesh->Faces[mesh->Volus[i]->Faces[3]]->Lines[0]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Faces[mesh->Volus[i]->Faces[3]]->Lines[1]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Faces[mesh->Volus[i]->Faces[3]]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+2;

    mesh->Volus[NV+i*7+0]->Faces[3]=NF+mesh->Volus[i]->Faces[3]*3+2;
    mesh->Volus[NV+i*7+1]->Faces[3]=NF+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Volus[i]->Faces[3]=NF+mesh->Volus[i]->Faces[3]*3+0;

   }  
  else if(mesh->Faces[mesh->Volus[i]->Faces[3]]->Verts[2]==mesh->Volus[NV+i*7+0]->Verts[0]
	&&mesh->Faces[mesh->Volus[i]->Faces[3]]->Verts[0]==mesh->Volus[NV+i*7+0]->Verts[2])
   {
    mesh->Volus[i]->Lines[3]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Volus[NV+i*7+0]->Lines[1]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Volus[NV+i*7+1]->Lines[0]=NE*2+mesh->Volus[i]->Faces[3]*3+2;
    mesh->Volus[NV+i*7+3]->Lines[4]=NE*2+mesh->Volus[i]->Faces[3]*3+2;
    mesh->Volus[NV+i*7+5]->Lines[3]=NE*2+mesh->Volus[i]->Faces[3]*3+2;
    mesh->Volus[NV+i*7+3]->Lines[0]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Volus[NV+i*7+3]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Volus[NV+i*7+4]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
  
    mesh->Faces[NF*4+i*8+1]->Lines[1]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Faces[NF*4+i*8+0]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Faces[NF*4+i*8+6]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Faces[NF*4+i*8+2]->Lines[0]=NE*2+mesh->Volus[i]->Faces[3]*3+2;
    mesh->Faces[NF*4+i*8+4]->Lines[1]=NE*2+mesh->Volus[i]->Faces[3]*3+2;
 
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+0]->Lines[0]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+0]->Lines[2]=mesh->Volus[i]->Lines[0];
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+0]->Lines[1]=mesh->Volus[i]->Lines[1];
    
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+1]->Lines[1]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+1]->Lines[2]=mesh->Volus[NV+i*7+0]->Lines[0];
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+1]->Lines[0]=mesh->Volus[NV+i*7+0]->Lines[3];

    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+2]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+2;
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+2]->Lines[0]=mesh->Volus[NV+i*7+1]->Lines[3];
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+2]->Lines[1]=mesh->Volus[NV+i*7+1]->Lines[1];

    mesh->Faces[mesh->Volus[i]->Faces[3]]->Lines[0]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Faces[mesh->Volus[i]->Faces[3]]->Lines[1]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Faces[mesh->Volus[i]->Faces[3]]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+2;

    mesh->Volus[NV+i*7+0]->Faces[3]=NF+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Volus[NV+i*7+1]->Faces[3]=NF+mesh->Volus[i]->Faces[3]*3+2;
    mesh->Volus[i]->Faces[3]=NF+mesh->Volus[i]->Faces[3]*3+0;



 }  
 // else (mesh->Faces[mesh->Volus[i]->Faces[3]]->Verts[2]==mesh->Volus[NV+i*7+0]->Verts[0]&&mesh->Faces[mesh->Volus[i]->Faces[1]]->Verts[1]==mesh->Volus[NV+i*7+0]->Verts[2])
else
   {
    mesh->Volus[i]->Lines[3]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Volus[NV+i*7+0]->Lines[1]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Volus[NV+i*7+1]->Lines[0]=NE*2+mesh->Volus[i]->Faces[3]*3+2;
    mesh->Volus[NV+i*7+3]->Lines[4]=NE*2+mesh->Volus[i]->Faces[3]*3+2;
    mesh->Volus[NV+i*7+5]->Lines[3]=NE*2+mesh->Volus[i]->Faces[3]*3+2;
    mesh->Volus[NV+i*7+3]->Lines[0]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Volus[NV+i*7+3]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Volus[NV+i*7+4]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
  
  
    mesh->Faces[NF*4+i*8+1]->Lines[1]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Faces[NF*4+i*8+0]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Faces[NF*4+i*8+6]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Faces[NF*4+i*8+2]->Lines[0]=NE*2+mesh->Volus[i]->Faces[3]*3+2;
    mesh->Faces[NF*4+i*8+4]->Lines[1]=NE*2+mesh->Volus[i]->Faces[3]*3+2;
 
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+0]->Lines[0]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+0]->Lines[2]=mesh->Volus[NV+i*7+0]->Lines[0];
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+0]->Lines[1]=mesh->Volus[NV+i*7+0]->Lines[3];
    
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+1]->Lines[1]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+1]->Lines[2]=mesh->Volus[i]->Lines[0];
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+1]->Lines[0]=mesh->Volus[i]->Lines[1];

    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+2]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+2;
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+2]->Lines[0]=mesh->Volus[NV+i*7+1]->Lines[1];
    mesh->Faces[NF+mesh->Volus[i]->Faces[3]*3+2]->Lines[1]=mesh->Volus[NV+i*7+1]->Lines[3];

    mesh->Faces[mesh->Volus[i]->Faces[3]]->Lines[0]=NE*2+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Faces[mesh->Volus[i]->Faces[3]]->Lines[1]=NE*2+mesh->Volus[i]->Faces[3]*3+1;
    mesh->Faces[mesh->Volus[i]->Faces[3]]->Lines[2]=NE*2+mesh->Volus[i]->Faces[3]*3+2;

    mesh->Volus[NV+i*7+0]->Faces[3]=NF+mesh->Volus[i]->Faces[3]*3+0;
    mesh->Volus[NV+i*7+1]->Faces[3]=NF+mesh->Volus[i]->Faces[3]*3+2;
    mesh->Volus[i]->Faces[3]=NF+mesh->Volus[i]->Faces[3]*3+1;
  }  


   //更新最中间的一条边被哪个单元哪个面包含
   mesh->Volus[NV+i*7+3]->Lines[5]=NE*2+NF*3+i;
   mesh->Volus[NV+i*7+4]->Lines[4]=NE*2+NF*3+i;
   mesh->Volus[NV+i*7+5]->Lines[5]=NE*2+NF*3+i;
   mesh->Volus[NV+i*7+6]->Lines[3]=NE*2+NF*3+i;
  
   mesh->Faces[NF*4+i*8+4]->Lines[0]=NE*2+NF*3+i;
   mesh->Faces[NF*4+i*8+5]->Lines[1]=NE*2+NF*3+i;
   mesh->Faces[NF*4+i*8+6]->Lines[0]=NE*2+NF*3+i;
   mesh->Faces[NF*4+i*8+7]->Lines[0]=NE*2+NF*3+i; 
   
   //更新不用判断顺序的面
   mesh->Volus[i]->Faces[0]=NF*4+i*8+0;
   mesh->Volus[NV+i*7+0]->Faces[1]=NF*4+i*8+1;
   mesh->Volus[NV+i*7+1]->Faces[2]=NF*4+i*8+2;
   mesh->Volus[NV+i*7+2]->Faces[3]=NF*4+i*8+3;

   mesh->Volus[NV+i*7+3]->Faces[0]=NF*4+i*8+4;
   mesh->Volus[NV+i*7+3]->Faces[1]=NF*4+i*8+6;
   mesh->Volus[NV+i*7+3]->Faces[2]=num_face[3];
   mesh->Volus[NV+i*7+3]->Faces[3]=NF*4+i*8+0;

   mesh->Volus[NV+i*7+4]->Faces[0]=NF*4+i*8+5;
   mesh->Volus[NV+i*7+4]->Faces[1]=NF*4+i*8+1;
   mesh->Volus[NV+i*7+4]->Faces[2]=NF*4+i*8+6;
   mesh->Volus[NV+i*7+4]->Faces[3]=num_face[2];

   mesh->Volus[NV+i*7+5]->Faces[0]=NF*4+i*8+4;
   mesh->Volus[NV+i*7+5]->Faces[1]=NF*4+i*8+7;
   mesh->Volus[NV+i*7+5]->Faces[2]=num_face[1];
   mesh->Volus[NV+i*7+5]->Faces[3]=NF*4+i*8+2;

   mesh->Volus[NV+i*7+6]->Faces[0]=NF*4+i*8+5;
   mesh->Volus[NV+i*7+6]->Faces[1]=num_face[0];
   mesh->Volus[NV+i*7+6]->Faces[2]=NF*4+i*8+3;
   mesh->Volus[NV+i*7+6]->Faces[3]=NF*4+i*8+7;
}


mesh->Num_Verts_Global=tmp_nverts;
mesh->Num_Lines_Global=tmp_nedges;
mesh->Num_Faces_Global=tmp_nfaces;
mesh->Num_Volus_Global=tmp_nvolus;

return;
}


//调整网格单元使网格单元保持右手法则
void AdjustMeshForRHF(MESH *mesh)
{
  int curr_volu,num_volu;
  VOLU *volu;
  VERT *vert;
  double x[6],y[6],z[6],d1,d2,d3,k1,k2,aa;
  int v0,v1,v2,v3,l0,l1,l2,l3,l4,l5,f0,f1,f2,f3;
  num_volu=mesh->Num_Volus_Global;
  for(curr_volu=0;curr_volu<num_volu;curr_volu++)
  {
    volu=mesh->Volus[curr_volu];
    aa=GetVolume(mesh,volu);
    //printf ( "curr_volu=%d,aa=%f\n",curr_volu,aa );
    if(aa>0)
    continue;
    //把单元点，边，面的信息取出来
    v0=volu->Verts[0];
    v1=volu->Verts[1];
    v2=volu->Verts[2];
    v3=volu->Verts[3];
    l0=volu->Lines[0];
    l1=volu->Lines[1];
    l2=volu->Lines[2];
    l3=volu->Lines[3];
    l4=volu->Lines[4];
    l5=volu->Lines[5];
    f0=volu->Faces[0];
    f1=volu->Faces[1];
    f2=volu->Faces[2];
    f3=volu->Faces[3];
      volu->Verts[0]=v0;
      volu->Verts[1]=v1;
      volu->Verts[2]=v3;
      volu->Verts[3]=v2;

      volu->Lines[0]=l0;
      volu->Lines[1]=l2;
      volu->Lines[2]=l1;
      volu->Lines[3]=l4;
      volu->Lines[4]=l3;
      volu->Lines[5]=l5;
 
      volu->Faces[0]=f0;
      volu->Faces[1]=f1;
      volu->Faces[2]=f3;
      volu->Faces[3]=f2;
  }//end for curr_volu
}




//调整网格单元，使得一致加密时不会出现奇异的边
void AdjustMesh(MESH *mesh)
{
  int curr_volu,num_volu;
  VOLU *volu;
  VERT *vert;
  double x[6],y[6],z[6],d1,d2,d3,k1,k2,k3;
  int v0,v1,v2,v3,l0,l1,l2,l3,l4,l5,f0,f1,f2,f3;
  num_volu=mesh->Num_Volus_Global;
  for(curr_volu=0;curr_volu<num_volu;curr_volu++)
  {
    volu=mesh->Volus[curr_volu];
    //把单元点，边，面的信息取出来
    v0=volu->Verts[0];
    v1=volu->Verts[1];
    v2=volu->Verts[2];
    v3=volu->Verts[3];
    l0=volu->Lines[0];
    l1=volu->Lines[1];
    l2=volu->Lines[2];
    l3=volu->Lines[3];
    l4=volu->Lines[4];
    l5=volu->Lines[5];
    f0=volu->Faces[0];
    f1=volu->Faces[1];
    f2=volu->Faces[2];
    f3=volu->Faces[3];
    //计算六条边中点的坐标
    x[0]=(mesh->Verts[v0]->Coord[0]+mesh->Verts[v1]->Coord[0])/2;
    y[0]=(mesh->Verts[v0]->Coord[1]+mesh->Verts[v1]->Coord[1])/2;
    z[0]=(mesh->Verts[v0]->Coord[2]+mesh->Verts[v1]->Coord[2])/2;

    x[1]=(mesh->Verts[v0]->Coord[0]+mesh->Verts[v2]->Coord[0])/2;
    y[1]=(mesh->Verts[v0]->Coord[1]+mesh->Verts[v2]->Coord[1])/2;
    z[1]=(mesh->Verts[v0]->Coord[2]+mesh->Verts[v2]->Coord[2])/2;

    x[2]=(mesh->Verts[v0]->Coord[0]+mesh->Verts[v3]->Coord[0])/2;
    y[2]=(mesh->Verts[v0]->Coord[1]+mesh->Verts[v3]->Coord[1])/2;
    z[2]=(mesh->Verts[v0]->Coord[2]+mesh->Verts[v3]->Coord[2])/2;

    x[3]=(mesh->Verts[v1]->Coord[0]+mesh->Verts[v2]->Coord[0])/2;
    y[3]=(mesh->Verts[v1]->Coord[1]+mesh->Verts[v2]->Coord[1])/2;
    z[3]=(mesh->Verts[v1]->Coord[2]+mesh->Verts[v2]->Coord[2])/2;

    x[4]=(mesh->Verts[v1]->Coord[0]+mesh->Verts[v3]->Coord[0])/2;
    y[4]=(mesh->Verts[v1]->Coord[1]+mesh->Verts[v3]->Coord[1])/2;
    z[4]=(mesh->Verts[v1]->Coord[2]+mesh->Verts[v3]->Coord[2])/2;
  
    x[5]=(mesh->Verts[v2]->Coord[0]+mesh->Verts[v3]->Coord[0])/2;
    y[5]=(mesh->Verts[v2]->Coord[1]+mesh->Verts[v3]->Coord[1])/2;
    z[5]=(mesh->Verts[v2]->Coord[2]+mesh->Verts[v3]->Coord[2])/2;
    //计算三条中间边的长度
    d1=sqrt((x[0]-x[5])*(x[0]-x[5])+(y[0]-y[5])*(y[0]-y[5])+(z[0]-z[5])*(z[0]-z[5]));
    d2=sqrt((x[1]-x[4])*(x[1]-x[4])+(y[1]-y[4])*(y[1]-y[4])+(z[1]-z[4])*(z[1]-z[4]));
    d3=sqrt((x[2]-x[3])*(x[2]-x[3])+(y[2]-y[3])*(y[2]-y[3])+(z[2]-z[3])*(z[2]-z[3]));
    //根据三条边的长度调整单元局部编号, 因为加密程序中默认使
    //                  用d3这条边，所以d3是最短边时，不用调整
    if((d1<d2)&&(d1<d3))
    {
      volu->Verts[0]=v0;
      volu->Verts[1]=v2;
      volu->Verts[2]=v3;
      volu->Verts[3]=v1;

      volu->Lines[0]=l1;
      volu->Lines[1]=l2;
      volu->Lines[2]=l0;
      volu->Lines[3]=l5;
      volu->Lines[4]=l3;
      volu->Lines[5]=l4;
 
      volu->Faces[0]=f0;
      volu->Faces[1]=f2;
      volu->Faces[2]=f3;
      volu->Faces[3]=f1;
    }
    if((d2<d3)&&(d2<d1))
    {
      volu->Verts[0]=v0;
      volu->Verts[1]=v3;
      volu->Verts[2]=v1;
      volu->Verts[3]=v2;

      volu->Lines[0]=l2;
      volu->Lines[1]=l0;
      volu->Lines[2]=l1;
      volu->Lines[3]=l4;
      volu->Lines[4]=l5;
      volu->Lines[5]=l3;
 
      volu->Faces[0]=f0;
      volu->Faces[1]=f3;
      volu->Faces[2]=f1;
      volu->Faces[3]=f2;
    }
   if((d1<d3)&&(d2<d3)&&(d1==d2))
   {
      if((y[5]-y[0])*(x[5]-x[0])!=0)
       {
        k1=(y[5]-y[0])/(x[5]-x[0]);
       }
      else
       k1=0;
      if( (z[5]-z[0])*(x[5]-x[0])!=0 )
      k2=(z[5]-z[0])/(x[5]-x[0]);
      else
      k2=0;
      if( (z[5]-z[0])*(y[5]-y[0])!=0 )
      k3=(z[5]-z[0])/(y[5]-y[0]);
      else
      k3=0;
      if(k1<0||k2<0||k3<0)
      {
      volu->Verts[0]=v0;
      volu->Verts[1]=v3;
      volu->Verts[2]=v1;
      volu->Verts[3]=v2;

      volu->Lines[0]=l2;
      volu->Lines[1]=l0;
      volu->Lines[2]=l1;
      volu->Lines[3]=l4;
      volu->Lines[4]=l5;
      volu->Lines[5]=l3;
 
      volu->Faces[0]=f0;
      volu->Faces[1]=f3;
      volu->Faces[2]=f1;
      volu->Faces[3]=f2;
    }
    else
    {
      volu->Verts[0]=v0;
      volu->Verts[1]=v2;
      volu->Verts[2]=v3;
      volu->Verts[3]=v1;

      volu->Lines[0]=l1;
      volu->Lines[1]=l2;
      volu->Lines[2]=l0;
      volu->Lines[3]=l5;
      volu->Lines[4]=l3;
      volu->Lines[5]=l4;
 
      volu->Faces[0]=f0;
      volu->Faces[1]=f2;
      volu->Faces[2]=f3;
      volu->Faces[3]=f1;
    }

   } 
 if((d1<d2)&&(d3<d2)&&(d1==d3))
   {
     if((y[2]-y[3])*(x[2]-x[3])!=0)
      k1=(y[2]-y[3])/(x[2]-x[3]);
     else
	k1=0;
     if	( (z[2]-z[3])*(x[2]-x[3])!=0)
     k2=(z[2]-z[3])/(x[2]-x[3]);
     else
	k2=0;
     if ((z[2]-z[3])*(y[2]-y[3]!=0))
     k3=(z[2]-z[3])/(y[2]-y[3]);
     else 
	k3=0;
     if(k1<0||k2<0||k3<0)
      {
      volu->Verts[0]=v0;
      volu->Verts[1]=v2;
      volu->Verts[2]=v3;
      volu->Verts[3]=v1;

      volu->Lines[0]=l1;
      volu->Lines[1]=l2;
      volu->Lines[2]=l0;
      volu->Lines[3]=l5;
      volu->Lines[4]=l3;
      volu->Lines[5]=l4;
 
      volu->Faces[0]=f0;
      volu->Faces[1]=f2;
      volu->Faces[2]=f3;
      volu->Faces[3]=f1;
      }
   } 

 if((d3<d1)&&(d2<d1)&&(d3==d2))
   {
     if((y[2]-y[3])*(x[2]-x[3])!=0)
      k1=(y[2]-y[3])/(x[2]-x[3]);
     else
	k1=0;
     if	( (z[2]-z[3])*(x[2]-x[3])!=0)
     k2=(z[2]-z[3])/(x[2]-x[3]);
     else
	k2=0;
     if ((z[2]-z[3])*(y[2]-y[3]!=0))
     k3=(z[2]-z[3])/(y[2]-y[3]);
     else 
	k3=0;

      if(k1<0||k2<0||k3<0)
      {	 
      volu->Verts[0]=v0;
      volu->Verts[1]=v3;
      volu->Verts[2]=v1;
      volu->Verts[3]=v2;

      volu->Lines[0]=l2;
      volu->Lines[1]=l0;
      volu->Lines[2]=l1;
      volu->Lines[3]=l4;
      volu->Lines[4]=l5;
      volu->Lines[5]=l3;
 
      volu->Faces[0]=f0;
      volu->Faces[1]=f3;
      volu->Faces[2]=f1;
      volu->Faces[3]=f2;
    }
    else
    {
    }

   } 
 if(d1==d2&&d1==d3&&d2==d3)
   {
      if((y[2]-y[3])*(x[2]-x[3])!=0)
      k1=(y[2]-y[3])/(x[2]-x[3]);
     else
	k1=0;
     if	( (z[2]-z[3])*(x[2]-x[3])!=0)
     k2=(z[2]-z[3])/(x[2]-x[3]);
     else
	k2=0;
     if ((z[2]-z[3])*(y[2]-y[3]!=0))
     k3=(z[2]-z[3])/(y[2]-y[3]);
     else 
	k3=0;

      if(k1>=0&&k2>=0&&k3>=0)
      {
      }
      
      if((y[5]-y[0])*(x[5]-x[0])!=0)
       {
        k1=(y[5]-y[0])/(x[5]-x[0]);
       }
      else
       k1=0;
      if( (z[5]-z[0])*(x[5]-x[0])!=0 )
      k2=(z[5]-z[0])/(x[5]-x[0]);
      else
      k2=0;
      if( (z[5]-z[0])*(y[5]-y[0])!=0 )
      k3=(z[5]-z[0])/(y[5]-y[0]);
      else
      k3=0;

     
      if(k1>=0&&k2>=0&&k3>=0)
      {

      volu->Verts[0]=v0;
      volu->Verts[1]=v2;
      volu->Verts[2]=v3;
      volu->Verts[3]=v1;

      volu->Lines[0]=l1;
      volu->Lines[1]=l2;
      volu->Lines[2]=l0;
      volu->Lines[3]=l5;
      volu->Lines[4]=l3;
      volu->Lines[5]=l4;
 
      volu->Faces[0]=f0;
      volu->Faces[1]=f2;
      volu->Faces[2]=f3;
      volu->Faces[3]=f1;

      } 
      else{
      volu->Verts[0]=v0;
      volu->Verts[1]=v3;
      volu->Verts[2]=v1;
      volu->Verts[3]=v2;

      volu->Lines[0]=l2;
      volu->Lines[1]=l0;
      volu->Lines[2]=l1;
      volu->Lines[3]=l4;
      volu->Lines[4]=l5;
      volu->Lines[5]=l3;
 
      volu->Faces[0]=f0;
      volu->Faces[1]=f3;
      volu->Faces[2]=f1;
      volu->Faces[3]=f2;
      } 
   
   }

  }//end for curr_volu
}



int GetLongestLine(double d1,double d2,double d3,double d4,double d5,double d6)
{
   if ((d1>=d2)&&(d1>d3)&&(d1>d4)&&(d1>=d5)&&(d1>=d6))
      return 1;
   else if((d2>=d1)&&(d2>d3)&&(d2>d4)&&(d2>=d5)&&(d2>=d6))
      return 2;
   else if((d3>=d1)&&(d3>=d2)&&(d3>=d4)&&(d3>=d5)&&(d3>=d6))
      return 3;
   else if((d4>=d1)&&(d4>=d2)&&(d4>=d3)&&(d4>=d5)&&(d5>=d6))
      return 4;
   else if((d5>=d1)&&(d5>=d2)&&(d5>d3)&&(d5>d4)&&(d5>=d6))
   {
      return 5;
   }
      else if((d6>=d1)&&(d6>=d2)&&(d6>d3)&&(d6>d4)&&(d6>=d5))
      return 6;
}

//调整单元的简单版本，没有考虑斜率的情况
void AdjustMeshL(MESH *mesh)
{
  int i,curr_volu,num_volu;
  VOLU *volu;
  VERT *vert;
  double x[6],y[6],z[6],d1,d2,d3,d4,d5,d6;
  int v0,v1,v2,v3,l0,l1,l2,l3,l4,l5,f0,f1,f2,f3;
  num_volu=mesh->Num_Volus_Global;
  for(curr_volu=0;curr_volu<num_volu;curr_volu++)
  {
    volu=mesh->Volus[curr_volu];
    //把单元点，边，面的信息取出来
    v0=volu->Verts[0];
    v1=volu->Verts[1];
    v2=volu->Verts[2];
    v3=volu->Verts[3];
    l0=volu->Lines[0];
    l1=volu->Lines[1];
    l2=volu->Lines[2];
    l3=volu->Lines[3];
    l4=volu->Lines[4];
    l5=volu->Lines[5];
    f0=volu->Faces[0];
    f1=volu->Faces[1];
    f2=volu->Faces[2];
    f3=volu->Faces[3];
    //取出四个顶点的坐标
    x[0]=mesh->Verts[v0]->Coord[0];
    y[0]=mesh->Verts[v0]->Coord[1];
    z[0]=mesh->Verts[v0]->Coord[2];

    x[1]=mesh->Verts[v1]->Coord[0];
    y[1]=mesh->Verts[v1]->Coord[1];
    z[1]=mesh->Verts[v1]->Coord[2];

    x[2]=mesh->Verts[v2]->Coord[0];
    y[2]=mesh->Verts[v2]->Coord[1];
    z[2]=mesh->Verts[v2]->Coord[2];

    x[3]=mesh->Verts[v3]->Coord[0];
    y[3]=mesh->Verts[v3]->Coord[1];
    z[3]=mesh->Verts[v3]->Coord[2];
    //计算六条边的长度
    d1=sqrt((x[0]-x[1])*(x[0]-x[1])+(y[0]-y[1])*(y[0]-y[1])+(z[0]-z[1])*(z[0]-z[1]));
    d2=sqrt((x[0]-x[2])*(x[0]-x[2])+(y[0]-y[2])*(y[0]-y[2])+(z[0]-z[2])*(z[0]-z[2]));
    d3=sqrt((x[0]-x[3])*(x[0]-x[3])+(y[0]-y[3])*(y[0]-y[3])+(z[0]-z[3])*(z[0]-z[3]));
    d4=sqrt((x[1]-x[2])*(x[1]-x[2])+(y[1]-y[2])*(y[1]-y[2])+(z[1]-z[2])*(z[1]-z[2]));
    d5=sqrt((x[1]-x[3])*(x[1]-x[3])+(y[1]-y[3])*(y[1]-y[3])+(z[1]-z[3])*(z[1]-z[3]));
    d6=sqrt((x[2]-x[3])*(x[2]-x[3])+(y[2]-y[3])*(y[2]-y[3])+(z[2]-z[3])*(z[2]-z[3]));
    //根据三条边的长度调整单元局部编号, 因为加密程序中默认使
    //                  用d3这条边，所以d3是最短边时，不用调整
    i=GetLongestLine(d1,d2,d3,d4,d5,d6);

    if((i==1)||(i==6))
    {
      volu->Verts[0]=v0;
      volu->Verts[1]=v2;
      volu->Verts[2]=v3;
      volu->Verts[3]=v1;

      volu->Lines[0]=l1;
      volu->Lines[1]=l2;
      volu->Lines[2]=l0;
      volu->Lines[3]=l5;
      volu->Lines[4]=l3;
      volu->Lines[5]=l4;
 
      volu->Faces[0]=f0;
      volu->Faces[1]=f2;
      volu->Faces[2]=f3;
      volu->Faces[3]=f1;
    }
    if((i==2)||(i==5))
    {
      volu->Verts[0]=v0;
      volu->Verts[1]=v3;
      volu->Verts[2]=v1;
      volu->Verts[3]=v2;

      volu->Lines[0]=l2;
      volu->Lines[1]=l0;
      volu->Lines[2]=l1;
      volu->Lines[3]=l4;
      volu->Lines[4]=l5;
      volu->Lines[5]=l3;
 
      volu->Faces[0]=f0;
      volu->Faces[1]=f3;
      volu->Faces[2]=f1;
      volu->Faces[3]=f2;
    } 
  }//end for curr_volu
}

//-----------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------
//-------------------------------------------自适应加密------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------
int MarkJudge(MESH *mesh,VOLU *volu)//判断单元是否需要再次加密
{
  if(mesh->Lines[volu->Lines[0]]->Mark>=0&&mesh->Lines[volu->Lines[0]]->Child[0]<0) {
  return 1;}
  else if (mesh->Lines[volu->Lines[0]]->Mark>=0&&mesh->Lines[volu->Lines[0]]->Child[0]>=0
	&&volu->Verts[0]!=mesh->Lines[volu->Lines[0]]->Child[2]
	&&volu->Verts[1]!=mesh->Lines[volu->Lines[0]]->Child[2]) return 1;
  else if (mesh->Lines[volu->Lines[1]]->Mark>=0&&mesh->Lines[volu->Lines[1]]->Child[0]<0) return 1;
  else if (mesh->Lines[volu->Lines[1]]->Mark>=0&&mesh->Lines[volu->Lines[1]]->Child[0]>=0
	&&volu->Verts[0]!=mesh->Lines[volu->Lines[1]]->Child[2]
	&&volu->Verts[2]!=mesh->Lines[volu->Lines[1]]->Child[2]){
  return 1;
}
  else if (mesh->Lines[volu->Lines[2]]->Mark>=0&&mesh->Lines[volu->Lines[2]]->Child[0]<0) return 1;
  else if (mesh->Lines[volu->Lines[2]]->Mark>=0&&mesh->Lines[volu->Lines[2]]->Child[0]>=0
	&&volu->Verts[0]!=mesh->Lines[volu->Lines[2]]->Child[2]
	&&volu->Verts[3]!=mesh->Lines[volu->Lines[2]]->Child[2]) return 1;
  else if (mesh->Lines[volu->Lines[3]]->Mark>=0&&mesh->Lines[volu->Lines[3]]->Child[0]<0) return 1;
  else if (mesh->Lines[volu->Lines[3]]->Mark>=0&&mesh->Lines[volu->Lines[3]]->Child[0]>=0
	&&volu->Verts[1]!=mesh->Lines[volu->Lines[3]]->Child[2]
	&&volu->Verts[2]!=mesh->Lines[volu->Lines[3]]->Child[2]) return 1;
  else if (mesh->Lines[volu->Lines[4]]->Mark>=0&&mesh->Lines[volu->Lines[4]]->Child[0]<0) return 1; 
  else if (mesh->Lines[volu->Lines[4]]->Mark>=0&&mesh->Lines[volu->Lines[4]]->Child[0]>=0
	&&volu->Verts[1]!=mesh->Lines[volu->Lines[4]]->Child[2]
	&&volu->Verts[3]!=mesh->Lines[volu->Lines[4]]->Child[2]) return 1;
  else if (mesh->Lines[volu->Lines[5]]->Mark>=0&&mesh->Lines[volu->Lines[5]]->Child[0]<0) return 1;
  else if (mesh->Lines[volu->Lines[5]]->Mark>=0&&mesh->Lines[volu->Lines[5]]->Child[0]>=0
	&&volu->Verts[2]!=mesh->Lines[volu->Lines[5]]->Child[2]
	&&volu->Verts[3]!=mesh->Lines[volu->Lines[5]]->Child[2]) return 1;
  else return -1;
}

static int cmp(const double *a,const double *b)
{
    return *(double*)a-*(double*)b>0?0:1;  
}


void
QuickSort(double  *x, int low, int high, size_t size, int *ix, int (*cmp)(const double *a,const double *b))
{
    int i, j, it, flag; size_t t;
    if (low < high) /*要排序的元素起止下标，保证小的放在左边，大的放在右边。这里以下标为low的元素为基准点*/
    {
	i = low; j = high;
	it = ix[0]; t = size*it;
	while (i<j) /*循环扫描*/
	{
	    flag = cmp((x+size*ix[j-low]), (x+t))>0?1:0;
	    while (i<j && flag) /*在右边的只要比基准点大仍放在右边*/
	    {
		j--; /*前移一个位置*/
		flag = cmp((x+size*ix[j-low]), (x+t))>0?1:0;
	    }
	    if (i<j)
	    {
		ix[i-low] = ix[j-low];
		i++; /*后移一个位置，并以此为基准点*/
	    }
	    flag = cmp((x+size*ix[i-low]), (x+t))>0?0:1;
	    while (i<j && flag) /*在左边的只要小于等于基准点仍放在左边*/
	    {
		i++; /*后移一个位置*/
		flag = cmp((x+size*ix[i-low]), (x+t))>0?0:1;
	    }
	    if (i<j)
	    {
		ix[j-low] = ix[i-low];
		j--; /*前移一个位置*/
	    }
	}
	ix[i-low] = it;
	QuickSort(x,low,i-1, size, ix, cmp);   /*对基准点左边的数再执行快速排序*/
	QuickSort(x,i+1,high, size, (ix+i+1-low), cmp);   /*对基准点右边的数再执行快速排序*/
    }
}


void MeshDestroy(MESH *mesh)
{
   int i;
   for(i=0;i<mesh->Num_Lines_Global;i++)
   {
      free(mesh->Lines[i]->Volus_Owned);
   }
   for(i=0;i<mesh->Num_Verts_Global;++i)
   {
      free(mesh->Verts[i]);
   }
   for(i=0;i<mesh->Num_Lines_Global;++i)
   {
      free(mesh->Lines[i]);
   }
   for(i=0;i<mesh->Num_Faces_Global;++i)
   {
      free(mesh->Faces[i]);
   }
   for(i=0;i<mesh->Num_Volus_Global;++i)
   {
      free(mesh->Volus[i]);
   }
    free(mesh->Verts);
    free(mesh->Lines);
    free(mesh->Faces);
    free(mesh->Volus);
    free(mesh);
    printf ( "Mesh Destroy\n" );
}


//读入需要一致加密的初始网格信息
void
ReadMesh(MESH *mesh,char *file)
{
  int i,j;
  FILE *fp = fopen(file, "r");
  fscanf(fp,"%d\n",&mesh->Num_Verts_Global);
  mesh->Verts=(VERT **)malloc(mesh->Num_Verts_Global);
  for(i=0;i<mesh->Num_Verts_Global;++i)
  {
   mesh->Verts[i]=(VERT *)malloc(1);
  }
  for(i=0;i<mesh->Num_Verts_Global;++i)
  {
   fscanf(fp,"%lf %lf %lf\n", &(mesh->Verts[i]->Coord[0]), &(mesh->Verts[i]->Coord[1]),&(mesh->Verts[i]->Coord[2]));
  }
  fscanf(fp,"%d\n",&mesh->Num_Lines_Global);
  mesh->Lines=(LINE **)malloc(mesh->Num_Lines_Global);
  for(i=0;i<mesh->Num_Lines_Global;++i)
  {
   mesh->Lines[i]=(LINE *)malloc(1);
  }
  for(i=0;i<mesh->Num_Lines_Global;++i)
  {
   fscanf(fp,"%d %d\n",&mesh->Lines[i]->Verts[0],&mesh->Lines[i]->Verts[1]);
  }
  fscanf(fp,"%d\n",&mesh->Num_Faces_Global);
  mesh->Faces=(FACEE **)malloc(mesh->Num_Faces_Global);
  for(i=0;i<mesh->Num_Faces_Global;++i)
  {
   mesh->Faces[i]=(FACEE *)malloc(1);
  }
  for(i=0;i<mesh->Num_Faces_Global;++i)
  {
   fscanf(fp,"%d %d %d\n",&mesh->Faces[i]->Verts[0],&mesh->Faces[i]->Verts[1],&mesh->Faces[i]->Verts[2]);
   fscanf(fp,"%d %d %d\n",&mesh->Faces[i]->Lines[0],&mesh->Faces[i]->Lines[1],&mesh->Faces[i]->Lines[2]);
   fscanf(fp,"%d\n",&mesh->Faces[i]->ID_Boundary);
  }
  fscanf(fp,"%d\n",&mesh->Num_Volus_Global);
  mesh->Volus=(VOLU **)malloc(mesh->Num_Volus_Global);
  for(i=0;i<mesh->Num_Volus_Global;++i)
  {
   mesh->Volus[i]=(VOLU *)malloc(1);
  }
  for(i=0;i<mesh->Num_Volus_Global;++i)
  {
   fscanf(fp,"%d %d %d %d\n",&mesh->Volus[i]->Verts[0],&mesh->Volus[i]->Verts[1],&mesh->Volus[i]->Verts[2],&mesh->Volus[i]->Verts[3]);
   fscanf(fp,"%d %d %d %d %d %d\n",&mesh->Volus[i]->Lines[0],&mesh->Volus[i]->Lines[1],&mesh->Volus[i]->Lines[2],
 	&mesh->Volus[i]->Lines[3],&mesh->Volus[i]->Lines[4],&mesh->Volus[i]->Lines[5]);
   fscanf(fp,"%d %d %d %d\n",&mesh->Volus[i]->Faces[0],&mesh->Volus[i]->Faces[1],&mesh->Volus[i]->Faces[2],&mesh->Volus[i]->Faces[3]);
   fscanf(fp,"%d %d\n",&mesh->Volus[i]->Ancestor, &mesh->Volus[i]->Father);
  }
  fclose(fp);
}

//读入需要二分加密的初始网格信息
void
ReadMeshBis(MESH *mesh,char *file)
{
  int i,j;
  FILE *fp = fopen(file, "r");
  fscanf(fp,"%d\n",&mesh->Num_Verts_Global);
  mesh->Verts=malloc(sizeof(VERT *)*mesh->Num_Verts_Global);
  for(i=0;i<mesh->Num_Verts_Global;++i)
  {
   mesh->Verts[i]=malloc(sizeof(VERT));
  }
  for(i=0;i<mesh->Num_Verts_Global;++i)
  {
   fscanf(fp,"%lf %lf %lf\n", &(mesh->Verts[i]->Coord[0]), &(mesh->Verts[i]->Coord[1]),&(mesh->Verts[i]->Coord[2]));
  }
  fscanf(fp,"%d\n",&mesh->Num_Lines_Global);
  mesh->Lines=malloc(sizeof(LINE *)*mesh->Num_Lines_Global);
  for(i=0;i<mesh->Num_Lines_Global;++i)
  {
   mesh->Lines[i]=malloc(sizeof(LINE));
  }
  for(i=0;i<mesh->Num_Lines_Global;++i)
  {
   fscanf(fp,"%d %d\n",&mesh->Lines[i]->Verts[0],&mesh->Lines[i]->Verts[1]);
   fscanf(fp,"%d\n",&mesh->Lines[i]->Mark);
   fscanf(fp,"%d\n",&mesh->Lines[i]->ID_Boundary);
   fscanf(fp,"%d %d %d\n",&mesh->Lines[i]->Child[0],&mesh->Lines[i]->Child[1],&mesh->Lines[i]->Child[2]);
   fscanf(fp,"%d\n",&mesh->Lines[i]->Num_Volus_Owned);
   mesh->Lines[i]->Volus_Owned=(int *)malloc(mesh->Lines[i]->Num_Volus_Owned);
   for(j=0;j<mesh->Lines[i]->Num_Volus_Owned;++j) 
     fscanf(fp,"%d",&mesh->Lines[i]->Volus_Owned[j]); 
   fscanf(fp,"\n");
  }
  fscanf(fp,"%d\n",&mesh->Num_Faces_Global);
  mesh->Faces=malloc(sizeof(FACEE *)*mesh->Num_Faces_Global);
  for(i=0;i<mesh->Num_Faces_Global;++i)
  {
   mesh->Faces[i]=malloc(sizeof(FACEE));
  }
  for(i=0;i<mesh->Num_Faces_Global;++i)
  {
   fscanf(fp,"%d %d %d\n",&mesh->Faces[i]->Verts[0],&mesh->Faces[i]->Verts[1],&mesh->Faces[i]->Verts[2]);
   fscanf(fp,"%d %d %d\n",&mesh->Faces[i]->Lines[0],&mesh->Faces[i]->Lines[1],&mesh->Faces[i]->Lines[2]);
   fscanf(fp,"%d %d %d %d %d\n",&mesh->Faces[i]->Child[0],&mesh->Faces[i]->Child[1],&mesh->Faces[i]->Child[2],&mesh->Faces[i]->Old_Child[0],&mesh->Faces[i]->Old_Child[1]);
   fscanf(fp,"%d\n",&mesh->Faces[i]->ID_Boundary);
  }
  fscanf(fp,"%d\n",&mesh->Num_Volus_Global);
  mesh->Volus=malloc(sizeof(VOLU *)*mesh->Num_Volus_Global);
  for(i=0;i<mesh->Num_Volus_Global;++i)
  {
   mesh->Volus[i]=malloc(sizeof(VOLU));
  }
  for(i=0;i<mesh->Num_Volus_Global;++i)
  {
   fscanf(fp,"%d %d %d %d\n",&mesh->Volus[i]->Verts[0],&mesh->Volus[i]->Verts[1],&mesh->Volus[i]->Verts[2],&mesh->Volus[i]->Verts[3]);
   fscanf(fp,"%d %d %d %d %d %d\n",&mesh->Volus[i]->Lines[0],&mesh->Volus[i]->Lines[1],&mesh->Volus[i]->Lines[2],
	&mesh->Volus[i]->Lines[3],&mesh->Volus[i]->Lines[4],&mesh->Volus[i]->Lines[5]);
   fscanf(fp,"%d %d %d %d\n",&mesh->Volus[i]->Faces[0],&mesh->Volus[i]->Faces[1],&mesh->Volus[i]->Faces[2],&mesh->Volus[i]->Faces[3]);
   fscanf(fp,"%d %d\n",&mesh->Volus[i]->Child[0],&mesh->Volus[i]->Child[1]);
   fscanf(fp,"%d\n",&mesh->Volus[i]->Type);
   fscanf(fp,"%d\n",&mesh->Volus[i]->Mark);
   fscanf(fp,"%d %d\n",&mesh->Volus[i]->Ancestor, &mesh->Volus[i]->Father);
  }
  fclose(fp);
}

//读入需要二分加密的初始网格信息简化版
void
ReadMeshLshape(MESH *mesh,char *file)
{
  int i,neighbours;
  FILE *fp = fopen(file, "r");
  mesh->Verts=(VERT **)malloc(mesh->Num_Verts_Global);
  for(i=0;i<mesh->Num_Verts_Global;++i)
  {
   mesh->Verts[i]=(VERT *)malloc(1);
  }
  for(i=0;i<mesh->Num_Verts_Global;++i)
  {
   fscanf(fp,"%lf %lf %lf\n", &(mesh->Verts[i]->Coord[0]), &(mesh->Verts[i]->Coord[1]),&(mesh->Verts[i]->Coord[2]));
  }
  fscanf(fp,"%d\n",&mesh->Num_Lines_Global);
  mesh->Lines=(LINE **)malloc(mesh->Num_Lines_Global);
  for(i=0;i<mesh->Num_Lines_Global;++i)
  {
   mesh->Lines[i]=(LINE *)malloc(1);
  }
  for(i=0;i<mesh->Num_Lines_Global;++i)
  {
   fscanf(fp,"%d %d\n",&mesh->Lines[i]->Verts[0],&mesh->Lines[i]->Verts[1]);
   fscanf(fp,"\n");
  }
  fscanf(fp,"%d\n",&mesh->Num_Faces_Global);
  mesh->Faces=(FACEE **)malloc(mesh->Num_Faces_Global);
  for(i=0;i<mesh->Num_Faces_Global;++i)
  {
   mesh->Faces[i]=(FACEE *)malloc(1);
  }
  for(i=0;i<mesh->Num_Faces_Global;++i)
  {
   fscanf(fp,"%d %d %d\n",&mesh->Faces[i]->Verts[0],&mesh->Faces[i]->Verts[1],&mesh->Faces[i]->Verts[2]);
   fscanf(fp,"%d %d %d\n",&mesh->Faces[i]->Lines[0],&mesh->Faces[i]->Lines[1],&mesh->Faces[i]->Lines[2]);
   fscanf(fp,"%d\n",&mesh->Faces[i]->ID_Boundary);
  }
  fscanf(fp,"%d\n",&mesh->Num_Volus_Global);
  mesh->Volus=(VOLU **)malloc(mesh->Num_Volus_Global);
  for(i=0;i<mesh->Num_Volus_Global;++i)
  {
   mesh->Volus[i]=(VOLU *)malloc(1);
  }
  for(i=0;i<mesh->Num_Volus_Global;++i)
  {
   fscanf(fp,"%d %d %d %d\n",&mesh->Volus[i]->Verts[0],&mesh->Volus[i]->Verts[1],&mesh->Volus[i]->Verts[2],&mesh->Volus[i]->Verts[3]);
   fscanf(fp,"%d %d %d %d %d %d\n",&mesh->Volus[i]->Lines[0],&mesh->Volus[i]->Lines[1],&mesh->Volus[i]->Lines[2],
	&mesh->Volus[i]->Lines[3],&mesh->Volus[i]->Lines[4],&mesh->Volus[i]->Lines[5]);
   fscanf(fp,"%d %d %d %d\n",&mesh->Volus[i]->Faces[0],&mesh->Volus[i]->Faces[1],&mesh->Volus[i]->Faces[2],&mesh->Volus[i]->Faces[3]);
  }
  fclose(fp);
}

//读入需要二分加密的初始网格信息简化版,凹区域，不包含面的信息，另外子函数补上
void
ReadMeshConcave(MESH *mesh,char *file)
{
  int i,j;
  FILE *fp = fopen(file, "r");
  fscanf(fp,"%d\n",&mesh->Num_Verts_Global);
  mesh->Verts=(VERT **)malloc(mesh->Num_Verts_Global);
  for(i=0;i<mesh->Num_Verts_Global;++i)
  {
   mesh->Verts[i]=(VERT *)malloc(1);
  }
  for(i=0;i<mesh->Num_Verts_Global;++i)
  {
   fscanf(fp,"%lf %lf %lf\n", &(mesh->Verts[i]->Coord[0]), &(mesh->Verts[i]->Coord[1]),&(mesh->Verts[i]->Coord[2]));
  }
  fscanf(fp,"%d\n",&mesh->Num_Lines_Global);
  mesh->Lines=(LINE **)malloc(mesh->Num_Lines_Global);
  for(i=0;i<mesh->Num_Lines_Global;++i)
  {
   mesh->Lines[i]=(LINE *)malloc(1);
  }
  for(i=0;i<mesh->Num_Lines_Global;++i)
  {
   fscanf(fp,"%d %d\n",&mesh->Lines[i]->Verts[0],&mesh->Lines[i]->Verts[1]);
   fscanf(fp,"\n");
  }
  fscanf(fp,"%d\n",&mesh->Num_Volus_Global);
  mesh->Volus=(VOLU **)malloc(mesh->Num_Volus_Global);
  for(i=0;i<mesh->Num_Volus_Global;++i)
  {
   mesh->Volus[i]=(VOLU *)malloc(1);
  }
  for(i=0;i<mesh->Num_Volus_Global;++i)
  {
   fscanf(fp,"%d %d %d %d\n",&mesh->Volus[i]->Verts[0],&mesh->Volus[i]->Verts[1],&mesh->Volus[i]->Verts[2],&mesh->Volus[i]->Verts[3]);
   fscanf(fp,"%d %d %d %d %d %d\n",&mesh->Volus[i]->Lines[0],&mesh->Volus[i]->Lines[1],&mesh->Volus[i]->Lines[2],
	&mesh->Volus[i]->Lines[3],&mesh->Volus[i]->Lines[4],&mesh->Volus[i]->Lines[5]);
   fscanf(fp,"%d %d %d %d\n",&mesh->Volus[i]->Faces[0],&mesh->Volus[i]->Faces[1],&mesh->Volus[i]->Faces[2],&mesh->Volus[i]->Faces[3]);
  }
  fscanf(fp,"%d\n",&mesh->Num_Faces_Global);
  mesh->Faces=(FACEE **)malloc(mesh->Num_Faces_Global);
  for(i=0;i<mesh->Num_Faces_Global;++i)
    mesh->Faces[i]=(FACEE *)malloc(1);
  for(i=0;i<mesh->Num_Faces_Global;++i)
    fscanf(fp,"%d\n",&mesh->Faces[i]->ID_Boundary);
  fclose(fp);
}




void InitialMesh(MESH *mesh)
{
 MeshFullInfo4Line(mesh);
 int i;
 for(i=0;i<mesh->Num_Lines_Global;++i)
  {
    mesh->Lines[i]->Mark = -1;
    mesh->Lines[i]->Child[0] = -1;
    mesh->Lines[i]->Child[1] = -1;
    mesh->Lines[i]->Child[2] = -1;
  }
 for(i=0;i<mesh->Num_Faces_Global;++i)
  {
    mesh->Faces[i]->Child[0] = -1;
    mesh->Faces[i]->Child[1] = -1;
    mesh->Faces[i]->Child[2] = -1;
    mesh->Faces[i]->Old_Child[0] = -1;
    mesh->Faces[i]->Old_Child[1] = -1;
  }
 for(i=0;i<mesh->Num_Volus_Global;++i)
  {
    mesh->Volus[i]->Type = 0;
    mesh->Volus[i]->Mark = -1;
    mesh->Volus[i]->Child[0] = -1;
    mesh->Volus[i]->Child[1] = -1;
    mesh->Volus[i]->Ancestor = i;
    mesh->Volus[i]->Father = -1;
  }

}


void GetFace4Mesh(MESH *mesh)
{
  int i,j,num_volu;
  num_volu = mesh->Num_Volus_Global;
  FACEE *face;
  VOLU *volu;
  for(i=0;i<num_volu;++i)
  {
    volu = mesh->Volus[i];
    for(j=0;j<4;++j)    
    {
       face = mesh->Faces[volu->Faces[j]];
       face->Verts[0] = volu->Verts[(j+1)%4];
       face->Verts[1] = volu->Verts[(j+2)%4];
       face->Verts[2] = volu->Verts[(j+3)%4];
       if(j==0)
       {
         face->Lines[0] = volu->Lines[5];
         face->Lines[1] = volu->Lines[4];
         face->Lines[2] = volu->Lines[3];
       }
       else if (j==1)
       {
         face->Lines[0] = volu->Lines[2];
         face->Lines[1] = volu->Lines[1];
         face->Lines[2] = volu->Lines[5];
       }
       else if (j==2)
       {
         face->Lines[0] = volu->Lines[0];
         face->Lines[1] = volu->Lines[4];
         face->Lines[2] = volu->Lines[2];
       }
       else
       {
         face->Lines[0] = volu->Lines[3];
         face->Lines[1] = volu->Lines[1];
         face->Lines[2] = volu->Lines[0];
       }
    
   }
 }
}



//按vtk格式输出网格信息
void WriteMesh(MESH *mesh, char *file)
{
    FILE *fp = fopen( file,"w");
    int i,j;
    /* 输出体的信息 */
    VOLU *e;
    fprintf (fp, "# vtk DataFile Version 2->0\n" );
    fprintf (fp, "Tetrahedron example\n" );
    fprintf (fp, "ASCII\n" );
    fprintf (fp, "DATASET POLYDATA\n" );
    fprintf (fp, "POINTS %d float\n",mesh->Num_Verts_Global );
    VERT *v;
    for(i=0;i<mesh->Num_Verts_Global;i++)
    {
     v=mesh->Verts[i];
     fprintf (fp,  "%f %f %f\n", v->Coord[0],v->Coord[1],v->Coord[2]);
    }
    /* 输出面的信息 */
    
    fprintf (fp,  "POLYGONS %d %d\n", mesh->Num_Faces_Global, mesh->Num_Faces_Global*4 );
    FACEE *f;
    for(i=0;i<mesh->Num_Faces_Global;i++)
    {
     f=mesh->Faces[i];
     fprintf (fp,  "%d %d %d %d\n",3,mesh->Faces[i]->Verts[0],mesh->Faces[i]->Verts[1],mesh->Faces[i]->Verts[2] );
    }
    fclose(fp);
}



void GetFacesVolusOwned(MESH *mesh)
{
  int i,j;
  FACEE *f;
  for(i=0;i<mesh->Num_Faces_Global;++i)
  {
    f=mesh->Faces[i]; 
    f->Num_Volus_Owned=2;
    f->Volus_Owned =(int *)malloc(2);
    f->Volus_Owned[0] = -1;
    f->Volus_Owned[1] = -1;
  }
  VOLU *v; 
  for(i=0;i<mesh->Num_Volus_Global;++i)
  {
    v=mesh->Volus[i]; 
    for(j=0;j<4;++j)
    {
      f=mesh->Faces[v->Faces[j]];
      if(f->Volus_Owned[0]==-1)
        f->Volus_Owned[0]=i; 
      else
        f->Volus_Owned[1]=i;
    }
  }
}




//输出网格信息供PHG调用
void WriteMesh4PHG(MESH *mesh, char *file)
{
    FILE *fp = fopen( file,"w");
    int i,j;
    /* 输出体的信息 */
    fprintf (fp, "DIM: 3\n" );
    fprintf (fp, "DIM_OF_WORLD: 3\n" );
    fprintf (fp, "number of vertices: %d\n", mesh->Num_Verts_Global);
    fprintf (fp, "number of elements: %d\n", mesh->Num_Volus_Global);
    fprintf (fp, "\n" );
    fprintf (fp, "vertex coordinates:\n" );
    VERT *v;
    for(i=0;i<mesh->Num_Verts_Global;i++)
    {
     v=mesh->Verts[i];
     fprintf (fp,  "%f %f %f\n", v->Coord[0],v->Coord[1],v->Coord[2]);
    }
    /* 输出面的信息 */
    fprintf (fp, "\n" );
    fprintf (fp,  "element vertices:\n");
    VOLU *e;
    for(i=0;i<mesh->Num_Volus_Global;i++)
    {
     e=mesh->Volus[i];
     fprintf (fp,  "%d %d %d %d\n",mesh->Volus[i]->Verts[0],mesh->Volus[i]->Verts[1],mesh->Volus[i]->Verts[2],mesh->Volus[i]->Verts[3] );
    }

    fprintf (fp, "\n" );
    fprintf (fp,  "element boundaries:\n");
    for(i=0;i<mesh->Num_Volus_Global;i++)
    {
     e=mesh->Volus[i];
     fprintf (fp,  "%d %d %d %d\n",mesh->Faces[mesh->Volus[i]->Faces[0]]->ID_Boundary,mesh->Faces[mesh->Volus[i]->Faces[1]]->ID_Boundary,
              mesh->Faces[mesh->Volus[i]->Faces[2]]->ID_Boundary,mesh->Faces[mesh->Volus[i]->Faces[3]]->ID_Boundary);
    }
    fprintf (fp, "\n" );
    fprintf (fp,  "element type:\n");
    for(i=0;i<mesh->Num_Volus_Global;i++)
    {
     fprintf (fp,  "%d\n",0);
    }
    int neigh[4];
    FACEE *f;
    fprintf (fp, "\n" );
    GetFacesVolusOwned(mesh);
    fprintf (fp,  "element neighbours:\n");
    for(i=0;i<mesh->Num_Volus_Global;i++)
    {
     e=mesh->Volus[i];
     for(j=0;j<4;++j)
     {    
       f=mesh->Faces[e->Faces[j]];
       if(f->Volus_Owned[0]==i)
         neigh[j]=f->Volus_Owned[1];
       else
         neigh[j]=f->Volus_Owned[0];
     }
     fprintf (fp, "%d %d %d %d\n",neigh[0],neigh[1],neigh[2],neigh[3]);
    }

    fclose(fp);
}






//读入需要二分加密的初始网格信息
void
WriteMeshBis(MESH *mesh,char *file)
{
  int i,j;
  FILE *fp = fopen(file, "w");
  fprintf(fp,"%d\n",mesh->Num_Verts_Global);
  for(i=0;i<mesh->Num_Verts_Global;++i)
  {
   fprintf(fp,"%lf %lf %lf\n", mesh->Verts[i]->Coord[0], mesh->Verts[i]->Coord[1],mesh->Verts[i]->Coord[2]);
  }
  fprintf(fp,"%d\n",mesh->Num_Lines_Global);
  for(i=0;i<mesh->Num_Lines_Global;++i)
  {
   fprintf(fp,"%d %d\n",mesh->Lines[i]->Verts[0],mesh->Lines[i]->Verts[1]);
   fprintf(fp,"%d\n",-1);
   fprintf(fp,"%d\n",mesh->Lines[i]->ID_Boundary);
   fprintf(fp,"%d %d %d\n",-1,-1,-1);
   fprintf(fp,"%d\n",mesh->Lines[i]->Num_Volus_Owned);
   for(j=0;j<mesh->Lines[i]->Num_Volus_Owned;++j) 
     fprintf(fp,"%d ",mesh->Lines[i]->Volus_Owned[j]); 
   fprintf(fp,"\n");
  }
  fprintf(fp,"%d\n",mesh->Num_Faces_Global);
  for(i=0;i<mesh->Num_Faces_Global;++i)
  {
   fprintf(fp,"%d %d %d\n",mesh->Faces[i]->Verts[0],mesh->Faces[i]->Verts[1],mesh->Faces[i]->Verts[2]);
   fprintf(fp,"%d %d %d\n",mesh->Faces[i]->Lines[0],mesh->Faces[i]->Lines[1],mesh->Faces[i]->Lines[2]);
   fprintf(fp,"%d %d %d %d %d\n",-1,-1,-1,-1,-1);
   fprintf(fp,"%d\n",mesh->Faces[i]->ID_Boundary);
  }
  fprintf(fp,"%d\n",mesh->Num_Volus_Global);
  for(i=0;i<mesh->Num_Volus_Global;++i)
  {
   fprintf(fp,"%d %d %d %d\n",mesh->Volus[i]->Verts[0],mesh->Volus[i]->Verts[1],mesh->Volus[i]->Verts[2],mesh->Volus[i]->Verts[3]);
   fprintf(fp,"%d %d %d %d %d %d\n",mesh->Volus[i]->Lines[0],mesh->Volus[i]->Lines[1],mesh->Volus[i]->Lines[2],
	mesh->Volus[i]->Lines[3],mesh->Volus[i]->Lines[4],mesh->Volus[i]->Lines[5]);
   fprintf(fp,"%d %d %d %d\n",mesh->Volus[i]->Faces[0],mesh->Volus[i]->Faces[1],mesh->Volus[i]->Faces[2],mesh->Volus[i]->Faces[3]);
   fprintf(fp,"%d %d\n",-1,-1);
   fprintf(fp,"%d\n",mesh->Volus[i]->Type);
   fprintf(fp,"%d\n",-1);
   fprintf(fp,"%d %d\n",i, -1);
  }
  fclose(fp);
}









void ShowMesh(MESH *mesh)
{
  int i,j;
  /* 输出体的信息 */
  VOLU *e;
  printf ( "Circle on the volus\n" );
  for(i=0;i<mesh->Num_Volus_Global;i++)
  {
   e=mesh->Volus[i];
   printf ( "verts:%d,%d,%d,%d\n",e->Verts[0],e->Verts[1],e->Verts[2],e->Verts[3]);
   printf ( "lines:%d,%d,%d,%d,%d,%d\n",e->Lines[0],e->Lines[1],e->Lines[2],e->Lines[3],e->Lines[4],e->Lines[5]);
   printf ( "faces:%d,%d,%d,%d\n",e->Faces[0],e->Faces[1],e->Faces[2],e->Faces[3] );
   printf ( "type:%d\n",e->Type);
   printf ( "mark:%d\n",e->Mark);
   printf ( "Child:%d,%d\n",e->Child[0],e->Child[1] );
   printf ( "\n" );	
  }
  /* 输出面的信息 */
  FACEE *f;
  printf ( "Circle on the face\n" );
  for(i=0;i<mesh->Num_Faces_Global;i++)
  {
   f=mesh->Faces[i];
   printf ( "lines: %d,%d,%d\n",mesh->Faces[i]->Lines[0],mesh->Faces[i]->Lines[1],mesh->Faces[i]->Lines[2] );
   printf ( "verts:%d,%d,%d\n",mesh->Faces[i]->Verts[0],mesh->Faces[i]->Verts[1],mesh->Faces[i]->Verts[2] );
   printf ( "ID_Boundary:%d\n",mesh->Faces[i]->ID_Boundary );
   printf ( "Child:%d,%d,%d\n",mesh->Faces[i]->Child[0],mesh->Faces[i]->Child[1],mesh->Faces[i]->Child[2] );
  }
  LINE *ed;
  printf ( "Circle on the line\n" );
  for(i=0;i<mesh->Num_Lines_Global;i++)
  {
   ed=mesh->Lines[i];
   printf ( "verts:%d,%d\n",ed->Verts[0],ed->Verts[1] );
   printf ( "Ownd by %ds simplex:",ed->Num_Volus_Owned );
   for(j=0;j<ed->Num_Volus_Owned;j++) 
     printf(" %d,",ed->Volus_Owned[j]);
   printf ( "\n" );
   printf ( "ID_Boundary:%d\n",ed->ID_Boundary);
   printf ( "mark:%d\n",ed->Mark);
   printf ( "Child:%d,%d\n",ed->Child[0],ed->Child[1],ed->Child[2] );
  }
  VERT *v;
  printf ( "Circle on the verts\n" );
  for(i=0;i<mesh->Num_Verts_Global;i++)
  {
   v=mesh->Verts[i];
   printf ( "Coord:%f,%f,%f\n", v->Coord[0],v->Coord[1],v->Coord[2]);
   printf ( "\n" );
  }
}
//输出一个单元的详细信息
void ShowVolu(MESH *mesh,VOLU *volu)
{
   printf ( "volu type:%d\n",volu->Type );
   printf ( "vert:%d,%d,%d,%d\n",volu->Verts[0],volu->Verts[1],volu->Verts[2],volu->Verts[3] );
   printf ( "line:%d,%d,%d,%d,%d,%d\n",volu->Lines[0],volu->Lines[1],volu->Lines[2],volu->Lines[3] ,volu->Lines[4],volu->Lines[5]);
   printf ( "face:%d,%d,%d,%d\n",volu->Faces[0],volu->Faces[1],volu->Faces[2],volu->Faces[3] );
   printf ( "volu mark and edge mark:%d.%d,%d,%d,%d,%d,%d\n",volu->Mark,mesh->Lines[volu->Lines[0]]->Mark,mesh->Lines[volu->Lines[1]]->Mark,
            mesh->Lines[volu->Lines[2]]->Mark,mesh->Lines[volu->Lines[3]]->Mark,mesh->Lines[volu->Lines[4]]->Mark,mesh->Lines[volu->Lines[5]]->Mark );
   printf ( "three child of l0:%d,%d,%d\n" ,mesh->Lines[volu->Lines[0]]->Child[0],mesh->Lines[volu->Lines[0]]->Child[1],
            mesh->Lines[volu->Lines[0]]->Child[2]);
   printf ( "five child of f2:%d,%d,%d,%d,%d\n",mesh->Faces[volu->Faces[2]]->Child[0] ,mesh->Faces[volu->Faces[2]]->Child[1],
            mesh->Faces[volu->Faces[2]]->Child[2],mesh->Faces[volu->Faces[2]]->Old_Child[0],mesh->Faces[volu->Faces[2]]->Old_Child[1]);
   printf ( "five child of f3:%d,%d,%d,%d,%d\n",mesh->Faces[volu->Faces[3]]->Child[0] ,mesh->Faces[volu->Faces[3]]->Child[1],
            mesh->Faces[volu->Faces[3]]->Child[2],mesh->Faces[volu->Faces[3]]->Old_Child[0],mesh->Faces[volu->Faces[3]]->Old_Child[1]);
}

//网格信息补全，补全边属于的单元
void 
MeshFullInfo4Line(MESH *mesh)
{
   int i, curr_line, curr_volu;
   LINE *line;
   VOLU *volu;
   int Num_Volus_Global=mesh->Num_Volus_Global;
   int Num_Lines_Global=mesh->Num_Lines_Global;
   for(i=0;i<Num_Lines_Global;++i)
   {
    mesh->Lines[i]->Num_Volus_Owned=0;
   }
  //遍历单元，更新边关于ｖｏｌｕｓ的个数
  for(curr_volu=0;curr_volu<Num_Volus_Global;++curr_volu)
  {
   volu=mesh->Volus[curr_volu];
   for(i=0;i<NLINES_PER_VOLU;++i)
   mesh->Lines[volu->Lines[i]]->Num_Volus_Owned+=1;
  }
  //分配内存
  for(i=0;i<Num_Lines_Global;++i)
  {
   if(mesh->Lines[i]->Volus_Owned==NULL)
      mesh->Lines[i]->Volus_Owned=(int *)malloc(mesh->Lines[i]->Num_Volus_Owned);
   else
      mesh->Lines[i]->Volus_Owned=
          (int *)realloc(mesh->Lines[i]->Volus_Owned,mesh->Lines[i]->Num_Volus_Owned);  
   mesh->Lines[i]->Count=0;
  }
  //遍历单元，更边面关于ｖｏｌｕｓ的ｏｗｎ值
 for(curr_volu=0;curr_volu<Num_Volus_Global;++curr_volu)
 {
  volu=mesh->Volus[curr_volu];
  for(i=0;i<NLINES_PER_VOLU;++i)
  {
   line=mesh->Lines[volu->Lines[i]];
   line->Volus_Owned[line->Count++]=curr_volu; 
  }
 }
}
//补全面属于的单元
void 
MeshFullInfo4Face(MESH *mesh)
{
   int i, curr_volu;
   FACEE *face;
   VOLU *volu;
   int Num_Volus_Global = mesh->Num_Volus_Global;
   int Num_Faces_Global = mesh->Num_Faces_Global;
   for(i=0;i<Num_Faces_Global;++i)
   {
    mesh->Faces[i]->Num_Volus_Owned=0;
    mesh->Faces[i]->Volus_Owned=NULL;
   }
   //遍历单元，更新面关于ｖｏｌｕｓ的个数
   for(curr_volu=0;curr_volu<Num_Volus_Global;++curr_volu)
   {
    volu=mesh->Volus[curr_volu];
    for(i=0;i<4;++i)
     mesh->Faces[volu->Faces[i]]->Num_Volus_Owned+=1;
   }
   //分配内存
   for(i=0;i<Num_Faces_Global;++i)
   {
    if(mesh->Faces[i]->Volus_Owned==NULL)
    {
     mesh->Faces[i]->Volus_Owned=(int *)malloc((mesh->Faces[i]->Num_Volus_Owned));
    }
    else
    {
     mesh->Faces[i]->Volus_Owned=
             (int *)realloc(mesh->Faces[i]->Volus_Owned,(mesh->Faces[i]->Num_Volus_Owned));  
    }
    mesh->Faces[i]->Count=0;
   }
   //遍历单元，更边面关于ｖｏｌｕｓ的ｏｗｎ值
   for(curr_volu=0;curr_volu<Num_Volus_Global;++curr_volu)
   {
    volu=mesh->Volus[curr_volu];
    for(i=0;i<4;++i)
    {
     face=mesh->Faces[volu->Faces[i]];
     face->Volus_Owned[face->Count++]=curr_volu; 
    }
   }
}


void 
MeshFullInfo4Num(MESH *mesh)
{
  int i;
  for(i=0;i<mesh->Num_Volus_Global;++i)
     mesh->Volus[i]->Num=i;
  for(i=0;i<mesh->Num_Faces_Global;++i)
     mesh->Faces[i]->Num=i;
  for(i=0;i<mesh->Num_Lines_Global;++i)
     mesh->Lines[i]->Num=i;
  for(i=0;i<mesh->Num_Verts_Global;++i)
     mesh->Verts[i]->Num=i;
}









//网格为下一次加密更新
void MeshRenew(MESH *mesh)
{
   int i;
   int nlines=mesh->Num_Lines_Global;
   int nfaces=mesh->Num_Faces_Global;
   int nvolus=mesh->Num_Volus_Global;
   for(i=0;i<nlines;++i)
   {
     mesh->Lines[i]->Mark=-1;
     mesh->Lines[i]->Child[0]=-1;
     mesh->Lines[i]->Child[1]=-1;
     mesh->Lines[i]->Child[2]=-1;
   }
  for(i=0;i<nfaces;++i)
   {
     mesh->Faces[i]->Old_Child[0]=-1;
     mesh->Faces[i]->Old_Child[1]=-1;
     mesh->Faces[i]->Child[0]=-1;
     mesh->Faces[i]->Child[1]=-1;
     mesh->Faces[i]->Child[2]=-1;
   }

 for(i=0;i<nvolus;++i)
   {
     mesh->Volus[i]->Mark=-1;
     mesh->Volus[i]->Child[0]=-1;
     mesh->Volus[i]->Child[1]=-1;
     mesh->Volus[i]->Trans = -1;
   }

}

//得到一条边在一个单元里的局部编号
int
GetLineNum(LINE *line,int line_num,VOLU *volu)
{
  if(line_num== volu->Lines[0])
   return 0;  
  else if(line_num== volu->Lines[1])
   return 1;  
  else if(line_num== volu->Lines[2])
   return 2;  
  else if(line_num==volu->Lines[3])
   return 3;  
  else if(line_num== volu->Lines[4])
   return 4;  
  else return 5;  
}

//自适应加密时，判断需要加密的单元
void
MarkMesh(MESH *mesh, VOLU *volu, int volu_num, LINE *line, int line_num, int *nline_mark)
{ 
  int k,j,i;
  VOLU *volu_owned;
  if(volu->Mark==-1)//如果单元没被标记
  {//1
    volu->Mark=1;
    if(volu->Type==0){//1->1
       k=GetLineNum(line,line_num,volu);
       switch(k){
            case 0: break;//如果line是单元的第零条边，则不做什么了，如果不是，则它的上一层级的边也需要被加密。
		          // e.g.，如果是单元的第一条边，则第零条边也要加密.
		          //如果是第二条边，则第零条，第一条边需要加密。
	    case 1: mesh->Lines[volu->Lines[0]]->Mark=*nline_mark;
		    (*nline_mark)++;
                    for(i=0;i<mesh->Lines[volu->Lines[0]]->Num_Volus_Owned;++i)
                     {
		      volu_owned=mesh->Volus[mesh->Lines[volu->Lines[0]]->Volus_Owned[i]];	    
		      MarkMesh(mesh,volu_owned,mesh->Lines[volu->Lines[0]]->Volus_Owned[i],
			       mesh->Lines[volu->Lines[0]],volu->Lines[0],
			       nline_mark);
		     }
		    break;
	    case 2: mesh->Lines[volu->Lines[0]]->Mark=(*nline_mark);
		    (*nline_mark)++;
                    for(i=0;i<mesh->Lines[volu->Lines[0]]->Num_Volus_Owned;++i)
                     {
		      volu_owned=mesh->Volus[mesh->Lines[volu->Lines[0]]->Volus_Owned[i]];	    
		      MarkMesh(mesh,volu_owned,mesh->Lines[volu->Lines[0]]->Volus_Owned[i],
			       mesh->Lines[volu->Lines[0]],volu->Lines[0],
			       nline_mark);
		     }
		    if(mesh->Lines[volu->Lines[1]]->Mark<0)//新加的if条件,怕上一步加密完零边后一边被顺带加密了
		    {
		     mesh->Lines[volu->Lines[1]]->Mark=(*nline_mark);
		     (*nline_mark)++;
		     for(i=0;i<mesh->Lines[volu->Lines[1]]->Num_Volus_Owned;++i)
                      {
		       volu_owned=mesh->Volus[mesh->Lines[volu->Lines[1]]->Volus_Owned[i]];	    
		       MarkMesh(mesh,volu_owned,mesh->Lines[volu->Lines[1]]->Volus_Owned[i],
			       mesh->Lines[volu->Lines[1]],volu->Lines[1],
	 		       nline_mark);
	     	      }
		    }
		    break;
	    case 3: mesh->Lines[volu->Lines[0]]->Mark=(*nline_mark);
		    (*nline_mark)++;
                    for(i=0;i<mesh->Lines[volu->Lines[0]]->Num_Volus_Owned;++i)
                     {
		      volu_owned=mesh->Volus[mesh->Lines[volu->Lines[0]]->Volus_Owned[i]];	    
		      MarkMesh(mesh,volu_owned,mesh->Lines[volu->Lines[0]]->Volus_Owned[i],
			       mesh->Lines[volu->Lines[0]],volu->Lines[0],
			       nline_mark);
		     }
		    if(mesh->Lines[volu->Lines[4]]->Mark<0)//新加的if
		    {		  
		     mesh->Lines[volu->Lines[4]]->Mark=(*nline_mark);
		     (*nline_mark)++;
		     for(i=0;i<mesh->Lines[volu->Lines[4]]->Num_Volus_Owned;++i)
                      {
		       volu_owned=mesh->Volus[mesh->Lines[volu->Lines[4]]->Volus_Owned[i]];	    
		       MarkMesh(mesh,volu_owned,mesh->Lines[volu->Lines[4]]->Volus_Owned[i],
			       mesh->Lines[volu->Lines[4]],volu->Lines[4],
			       nline_mark);
		      }
		    }
		    break;
	    case 4: mesh->Lines[volu->Lines[0]]->Mark=(*nline_mark);  
		    (*nline_mark)++;
		    for(i=0;i<mesh->Lines[volu->Lines[0]]->Num_Volus_Owned;++i)
                     {
		      volu_owned=mesh->Volus[mesh->Lines[volu->Lines[0]]->Volus_Owned[i]];	    
		      MarkMesh(mesh,volu_owned,mesh->Lines[volu->Lines[0]]->Volus_Owned[i],
			       mesh->Lines[volu->Lines[0]],volu->Lines[0],
			       nline_mark);
		     }
		    break;
	    case 5: mesh->Lines[volu->Lines[0]]->Mark=(*nline_mark);  
		    (*nline_mark)++;
                    for(i=0;i<mesh->Lines[volu->Lines[0]]->Num_Volus_Owned;++i)
                     {
		      volu_owned=mesh->Volus[mesh->Lines[volu->Lines[0]]->Volus_Owned[i]];	    
		      MarkMesh(mesh,volu_owned,mesh->Lines[volu->Lines[0]]->Volus_Owned[i],
			       mesh->Lines[volu->Lines[0]],volu->Lines[0],
			       nline_mark);
		     }
		    if(mesh->Lines[volu->Lines[1]]->Mark<0)
		    {
		     mesh->Lines[volu->Lines[1]]->Mark=(*nline_mark);  
		     (*nline_mark)++;
                     for(i=0;i<mesh->Lines[volu->Lines[1]]->Num_Volus_Owned;++i)
                      {
		       volu_owned=mesh->Volus[mesh->Lines[volu->Lines[1]]->Volus_Owned[i]];	    
		       MarkMesh(mesh,volu_owned,mesh->Lines[volu->Lines[1]]->Volus_Owned[i],
			       mesh->Lines[volu->Lines[1]],volu->Lines[1],
			       nline_mark);
		      }
		    }
		    if(mesh->Lines[volu->Lines[4]]->Mark<0)
		    {		    
		     mesh->Lines[volu->Lines[4]]->Mark=(*nline_mark);  
		     (*nline_mark)++;
		     for(i=0;i<mesh->Lines[volu->Lines[4]]->Num_Volus_Owned;++i)
                      {
		       volu_owned=mesh->Volus[mesh->Lines[volu->Lines[4]]->Volus_Owned[i]];	    
		       MarkMesh(mesh,volu_owned,mesh->Lines[volu->Lines[4]]->Volus_Owned[i],
			       mesh->Lines[volu->Lines[4]],volu->Lines[4],
			       nline_mark);
		      }
		    }
		    break;
	    default:printf ( "line type is not conforming,type=0,k=%d\n",k ); break;
       }

    }//1->1
    else if(volu->Type==1){//1->2
       k=GetLineNum(line,line_num,volu);
       switch(k){
            case 0: break;
	    case 1: mesh->Lines[volu->Lines[0]]->Mark=(*nline_mark);
		    (*nline_mark)++;
		    for(i=0;i<mesh->Lines[volu->Lines[0]]->Num_Volus_Owned;++i)
                     {
		      volu_owned=mesh->Volus[mesh->Lines[volu->Lines[0]]->Volus_Owned[i]];	    
		      MarkMesh(mesh,volu_owned,mesh->Lines[volu->Lines[0]]->Volus_Owned[i],
			       mesh->Lines[volu->Lines[0]],volu->Lines[0],
			       nline_mark);
		     }
		    break;
	    case 2: mesh->Lines[volu->Lines[0]]->Mark=(*nline_mark);
		    (*nline_mark)++;
                    for(i=0;i<mesh->Lines[volu->Lines[0]]->Num_Volus_Owned;++i)
                     {
		      volu_owned=mesh->Volus[mesh->Lines[volu->Lines[0]]->Volus_Owned[i]];	    
		      MarkMesh(mesh,volu_owned,mesh->Lines[volu->Lines[0]]->Volus_Owned[i],
			       mesh->Lines[volu->Lines[0]],volu->Lines[0],
			       nline_mark);
		     }
		    if(mesh->Lines[volu->Lines[1]]->Mark<0)
		    {
		     mesh->Lines[volu->Lines[1]]->Mark=(*nline_mark);
		     (*nline_mark)++;
		     for(i=0;i<mesh->Lines[volu->Lines[1]]->Num_Volus_Owned;++i)
                      {
		       volu_owned=mesh->Volus[mesh->Lines[volu->Lines[1]]->Volus_Owned[i]];	    
		       MarkMesh(mesh,volu_owned,mesh->Lines[volu->Lines[1]]->Volus_Owned[i],
			       mesh->Lines[volu->Lines[1]],volu->Lines[1],
			       nline_mark);
		      }
		    }
		    break;
	    case 3: mesh->Lines[volu->Lines[0]]->Mark=(*nline_mark);
		    (*nline_mark)++;
		    for(i=0;i<mesh->Lines[volu->Lines[0]]->Num_Volus_Owned;++i)
                     {
		      volu_owned=mesh->Volus[mesh->Lines[volu->Lines[0]]->Volus_Owned[i]];	    
		      MarkMesh(mesh,volu_owned,mesh->Lines[volu->Lines[0]]->Volus_Owned[i],
			       mesh->Lines[volu->Lines[0]],volu->Lines[0],
			       nline_mark);
		     }
		    break;
	    case 4: mesh->Lines[volu->Lines[0]]->Mark=(*nline_mark);
		    (*nline_mark)++;
                    for(i=0;i<mesh->Lines[volu->Lines[0]]->Num_Volus_Owned;++i)
                     {
		      volu_owned=mesh->Volus[mesh->Lines[volu->Lines[0]]->Volus_Owned[i]];	    
		      MarkMesh(mesh,volu_owned,mesh->Lines[volu->Lines[0]]->Volus_Owned[i],
			       mesh->Lines[volu->Lines[0]],volu->Lines[0],
			       nline_mark);
		     }
                    if(mesh->Lines[volu->Lines[3]]->Mark<0)
		    {
		     mesh->Lines[volu->Lines[3]]->Mark=(*nline_mark);
		     (*nline_mark)++;
		     for(i=0;i<mesh->Lines[volu->Lines[3]]->Num_Volus_Owned;++i)
                      {
		       volu_owned=mesh->Volus[mesh->Lines[volu->Lines[3]]->Volus_Owned[i]];	    
		       MarkMesh(mesh,volu_owned,mesh->Lines[volu->Lines[3]]->Volus_Owned[i],
			       mesh->Lines[volu->Lines[3]],volu->Lines[3],
			       nline_mark);
		      }
		    }
		    break;
	    case 5: mesh->Lines[volu->Lines[0]]->Mark=(*nline_mark);
		    (*nline_mark)++;
                    for(i=0;i<mesh->Lines[volu->Lines[0]]->Num_Volus_Owned;++i)
                     {
		      volu_owned=mesh->Volus[mesh->Lines[volu->Lines[0]]->Volus_Owned[i]];	    
		      MarkMesh(mesh,volu_owned,mesh->Lines[volu->Lines[0]]->Volus_Owned[i],
			       mesh->Lines[volu->Lines[0]],volu->Lines[0],
			       nline_mark);
		     }
		    if(mesh->Lines[volu->Lines[1]]->Mark<0)
		    {
                     mesh->Lines[volu->Lines[1]]->Mark=(*nline_mark);
		     (*nline_mark)++;
                     for(i=0;i<mesh->Lines[volu->Lines[1]]->Num_Volus_Owned;++i)
                      {
		       volu_owned=mesh->Volus[mesh->Lines[volu->Lines[1]]->Volus_Owned[i]];	    
		       MarkMesh(mesh,volu_owned,mesh->Lines[volu->Lines[1]]->Volus_Owned[i],
			       mesh->Lines[volu->Lines[1]],volu->Lines[1],
			       nline_mark);
		      }
		    }
                    if(mesh->Lines[volu->Lines[3]]->Mark)
		    {
		     mesh->Lines[volu->Lines[3]]->Mark=(*nline_mark);
		     (*nline_mark)++;
		     for(i=0;i<mesh->Lines[volu->Lines[3]]->Num_Volus_Owned;++i)
                      {
		       volu_owned=mesh->Volus[mesh->Lines[volu->Lines[3]]->Volus_Owned[i]];	    
		       MarkMesh(mesh,volu_owned,mesh->Lines[volu->Lines[3]]->Volus_Owned[i],
			       mesh->Lines[volu->Lines[3]],volu->Lines[3],
			       nline_mark);
		      }
		    }
		    break;
	    default:printf ( "line type is not conforming,type=1,k=%d\n",k ); break;
       }
    }//1->2
    else{//1->3 
       k=GetLineNum(line,line_num,volu);
       switch(k){
            case 0: break;
	    case 1: mesh->Lines[volu->Lines[0]]->Mark=(*nline_mark);
		    (*nline_mark)++;
		    for(i=0;i<mesh->Lines[volu->Lines[0]]->Num_Volus_Owned;++i)
                     {
		      volu_owned=mesh->Volus[mesh->Lines[volu->Lines[0]]->Volus_Owned[i]];	    
		      MarkMesh(mesh,volu_owned,mesh->Lines[volu->Lines[0]]->Volus_Owned[i],
			       mesh->Lines[volu->Lines[0]],volu->Lines[0],
			       nline_mark);
		     }
		    break;
	    case 2: mesh->Lines[volu->Lines[0]]->Mark=(*nline_mark);
		    (*nline_mark)++;
                    for(i=0;i<mesh->Lines[volu->Lines[0]]->Num_Volus_Owned;++i)
                     {
		      volu_owned=mesh->Volus[mesh->Lines[volu->Lines[0]]->Volus_Owned[i]];	    
		      MarkMesh(mesh,volu_owned,mesh->Lines[volu->Lines[0]]->Volus_Owned[i],
			       mesh->Lines[volu->Lines[0]],volu->Lines[0],
			       nline_mark);
		     }
                    if(mesh->Lines[volu->Lines[1]]->Mark<0)
		    {
		     mesh->Lines[volu->Lines[1]]->Mark=(*nline_mark);
		     (*nline_mark)++;
		     for(i=0;i<mesh->Lines[volu->Lines[1]]->Num_Volus_Owned;++i)
                      {
		       volu_owned=mesh->Volus[mesh->Lines[volu->Lines[1]]->Volus_Owned[i]];	    
		       MarkMesh(mesh,volu_owned,mesh->Lines[volu->Lines[1]]->Volus_Owned[i],
			       mesh->Lines[volu->Lines[1]],volu->Lines[1],
			       nline_mark);
		      }
		    }
		    break;
	    case 3: mesh->Lines[volu->Lines[0]]->Mark=(*nline_mark);
		    (*nline_mark)++;
		    for(i=0;i<mesh->Lines[volu->Lines[0]]->Num_Volus_Owned;++i)
                     {
		      volu_owned=mesh->Volus[mesh->Lines[volu->Lines[0]]->Volus_Owned[i]];	    
		      MarkMesh(mesh,volu_owned,mesh->Lines[volu->Lines[0]]->Volus_Owned[i],
			       mesh->Lines[volu->Lines[0]],volu->Lines[0],
			       nline_mark);
		     }
		    break;
	    case 4: mesh->Lines[volu->Lines[0]]->Mark=(*nline_mark);
		    (*nline_mark)++;
                    for(i=0;i<mesh->Lines[volu->Lines[0]]->Num_Volus_Owned;++i)
                     {
		      volu_owned=mesh->Volus[mesh->Lines[volu->Lines[0]]->Volus_Owned[i]];	    
		      MarkMesh(mesh,volu_owned,mesh->Lines[volu->Lines[0]]->Volus_Owned[i],
			       mesh->Lines[volu->Lines[0]],volu->Lines[0],
			       nline_mark);
		     }
                    if(mesh->Lines[volu->Lines[3]]->Mark)
		    {
		     mesh->Lines[volu->Lines[3]]->Mark=(*nline_mark);
		     (*nline_mark)++;
		     for(i=0;i<mesh->Lines[volu->Lines[3]]->Num_Volus_Owned;++i)
                      {
		       volu_owned=mesh->Volus[mesh->Lines[volu->Lines[3]]->Volus_Owned[i]];	    
		        MarkMesh(mesh,volu_owned,mesh->Lines[volu->Lines[3]]->Volus_Owned[i],
			       mesh->Lines[volu->Lines[3]],volu->Lines[3],
			       nline_mark);
		      }
		    }
		    break;
	    case 5: volu->Trans = 1;mesh->Mark = 1; break;
;
	    //case 5: mesh->Mark = 1; break;

	    default:printf ( "line type is not conforming,type=2,k=%d\n",k ); break;
       }

    }//1->3
  }//1
  else//如果单元被标记了
  {//2
    if(volu->Type==0){//2.1
       k=GetLineNum(line,line_num,volu);
       switch(k){//如果单元已经被标记了，当line是第零条边时，则这条边肯定被标记了，如果line是第一条边，则
	         //因为单元已经被标记了，所以第零条边肯定也被标记了，所以不用做什么，如果line是第二条边，
		 //则需要判断一下它是上层级边是否已经加密，如果没有，需要加密。
            case 0: break;
	    case 1: break;
	    case 2: if(mesh->Lines[volu->Lines[1]]->Mark==-1)
		    {
		    mesh->Lines[volu->Lines[1]]->Mark=(*nline_mark);
		    (*nline_mark)++;
		    for(i=0;i<mesh->Lines[volu->Lines[1]]->Num_Volus_Owned;++i)
                     {
		      volu_owned=mesh->Volus[mesh->Lines[volu->Lines[1]]->Volus_Owned[i]];	    
		      MarkMesh(mesh,volu_owned,mesh->Lines[volu->Lines[1]]->Volus_Owned[i],
			       mesh->Lines[volu->Lines[1]],volu->Lines[1],
			       nline_mark);
		     }
		    }
		    break;
	    case 3: if(mesh->Lines[volu->Lines[4]]->Mark==-1)
		    {
		    mesh->Lines[volu->Lines[4]]->Mark=(*nline_mark);  
		    (*nline_mark)++;
		    for(i=0;i<mesh->Lines[volu->Lines[4]]->Num_Volus_Owned;++i)
                     {
		      volu_owned=mesh->Volus[mesh->Lines[volu->Lines[4]]->Volus_Owned[i]];	    
		      MarkMesh(mesh,volu_owned,mesh->Lines[volu->Lines[4]]->Volus_Owned[i],
			       mesh->Lines[volu->Lines[4]],volu->Lines[4],
			       nline_mark);
		     }
		    }
		    break;
	    case 4: break;
	    case 5: if(mesh->Lines[volu->Lines[1]]->Mark==-1)
		    {
		    mesh->Lines[volu->Lines[1]]->Mark=(*nline_mark);  
		    (*nline_mark)++;
                    for(i=0;i<mesh->Lines[volu->Lines[1]]->Num_Volus_Owned;++i)
                     {
		      volu_owned=mesh->Volus[mesh->Lines[volu->Lines[1]]->Volus_Owned[i]];	    
		      MarkMesh(mesh,volu_owned,mesh->Lines[volu->Lines[1]]->Volus_Owned[i],
			       mesh->Lines[volu->Lines[1]],volu->Lines[1],
			       nline_mark);
		     }
		    }
		    if(mesh->Lines[volu->Lines[4]]->Mark==-1)
		    {
		    mesh->Lines[volu->Lines[4]]->Mark=(*nline_mark);  
		    (*nline_mark)++;
		    for(i=0;i<mesh->Lines[volu->Lines[4]]->Num_Volus_Owned;++i)
                     {
		      volu_owned=mesh->Volus[mesh->Lines[volu->Lines[4]]->Volus_Owned[i]];	    
		      MarkMesh(mesh,volu_owned,mesh->Lines[volu->Lines[4]]->Volus_Owned[i],
			       mesh->Lines[volu->Lines[4]],volu->Lines[4],
			       nline_mark);
		     }
		    }
		    break;
	    default:printf ( "line type is not conforming,type=0,biaoji,k=%d\n",k ); break;
       }
    
    }//2.1
    else if(volu->Type==1){//2.2
       k=GetLineNum(line,line_num,volu);
       switch(k){
            case 0: break;
	    case 1: break;
	    case 2: if(mesh->Lines[volu->Lines[1]]->Mark==-1)
		    {
		    mesh->Lines[volu->Lines[1]]->Mark=(*nline_mark);
		    (*nline_mark)++;
		    for(i=0;i<mesh->Lines[volu->Lines[1]]->Num_Volus_Owned;++i)
                     {
		      volu_owned=mesh->Volus[mesh->Lines[volu->Lines[1]]->Volus_Owned[i]];	    
		      MarkMesh(mesh,volu_owned,mesh->Lines[volu->Lines[1]]->Volus_Owned[i],
			       mesh->Lines[volu->Lines[1]],volu->Lines[1],
			       nline_mark);
		     }
		    }
		    break;
	    case 3: break;
	    case 4: if(mesh->Lines[volu->Lines[3]]->Mark==-1)
		    {
                    mesh->Lines[volu->Lines[3]]->Mark=(*nline_mark);
		    (*nline_mark)++;
		    for(i=0;i<mesh->Lines[volu->Lines[3]]->Num_Volus_Owned;++i)
                     {
		      volu_owned=mesh->Volus[mesh->Lines[volu->Lines[3]]->Volus_Owned[i]];	    
		      MarkMesh(mesh,volu_owned,mesh->Lines[volu->Lines[3]]->Volus_Owned[i],
			       mesh->Lines[volu->Lines[3]],volu->Lines[3],
			       nline_mark);
		     }
		    }
		    break;
	    case 5: if(mesh->Lines[volu->Lines[1]]->Mark==-1)
		    {
                    mesh->Lines[volu->Lines[1]]->Mark=(*nline_mark);
		    (*nline_mark)++;
                    for(i=0;i<mesh->Lines[volu->Lines[1]]->Num_Volus_Owned;++i)
                     {
		      volu_owned=mesh->Volus[mesh->Lines[volu->Lines[1]]->Volus_Owned[i]];	    
		      MarkMesh(mesh,volu_owned,mesh->Lines[volu->Lines[1]]->Volus_Owned[i],
			       mesh->Lines[volu->Lines[1]],volu->Lines[1],
			       nline_mark);
		     }
		    }
		    if(mesh->Lines[volu->Lines[3]]->Mark==-1)
		    {
                    mesh->Lines[volu->Lines[3]]->Mark=(*nline_mark);
		    (*nline_mark)++;
		    for(i=0;i<mesh->Lines[volu->Lines[3]]->Num_Volus_Owned;++i)
                     {
		      volu_owned=mesh->Volus[mesh->Lines[volu->Lines[3]]->Volus_Owned[i]];	    
		      MarkMesh(mesh,volu_owned,mesh->Lines[volu->Lines[3]]->Volus_Owned[i],
			       mesh->Lines[volu->Lines[3]],volu->Lines[3],
			       nline_mark);
		     }
		    }
		    break;
	    default:printf ( "line type is not conforming,type=1,biaoji,k=%d\n",k ); break;
       }

    }//2->2
    else{//2.3
       k=GetLineNum(line,line_num,volu);
       switch(k){
            case 0: break;
	    case 1: break;
	    case 2: if(mesh->Lines[volu->Lines[1]]->Mark==-1)//没有改为正数啊
		    {
                       mesh->Lines[volu->Lines[1]]->Mark=(*nline_mark);
		       (*nline_mark)++;
		    for(i=0;i<mesh->Lines[volu->Lines[1]]->Num_Volus_Owned;++i)
                     {
                     		      volu_owned=mesh->Volus[mesh->Lines[volu->Lines[1]]->Volus_Owned[i]];	    
		      MarkMesh(mesh,volu_owned,mesh->Lines[volu->Lines[1]]->Volus_Owned[i],
			       mesh->Lines[volu->Lines[1]],volu->Lines[1],
			       nline_mark);
		     }
		    }
		    break;
	    case 3: break;
	    case 4: if(mesh->Lines[volu->Lines[3]]->Mark==-1)
		    {
		       mesh->Lines[volu->Lines[3]]->Mark=(*nline_mark);
		       (*nline_mark)++;
		    for(i=0;i<mesh->Lines[volu->Lines[3]]->Num_Volus_Owned;++i)
                     {
		      volu_owned=mesh->Volus[mesh->Lines[volu->Lines[3]]->Volus_Owned[i]];	    
		      MarkMesh(mesh,volu_owned,mesh->Lines[volu->Lines[3]]->Volus_Owned[i],
			       mesh->Lines[volu->Lines[3]],volu->Lines[3],
			       nline_mark);
		     }
		    }
		    break;

	    case 5: volu->Trans = 1;mesh->Mark = 1; break;
	    //case 5: mesh->Mark = 1; break;
	    default:printf ( "line type is not conforming,type=2,biaoji,k=%d,volu_num=%d\n",k,volu_num ); 
		    break;
       }


    }//2->3

  }//2
}






//网格二分加密
//
//               0
//                ->                                     网格的六个顶点 0,1,2,3,4,5按右手法则排序。
//               . . .                                  网格六条边按赵定点顺序：第零条边[0]:0-1
//              .   .   .                               [1]:0-2,[2]:0-3,[4]:1-2,[5]:1-3,[6]:2-3。
//             .     .     .                            四个面的编号为对应上方顶点编号：face[0]:1-2-3,
//            .    (2).       .                         face[1]:0-2-3,face[2]:0-3-1,face[3]:0-1-2. 
//         (0).        .         . (1)
//          .           . 3         . 
//         .          .    .           .
//        .    (4).             .(5)      .
//       .     .                     .       .
//      .  .                               .     .
//    1.............................................2           
//                       (3) 
void
MeshRefineBisection(MESH *mesh, int *M, int *num)
{
  printf ( "begin adaptive refine\n" );
  int i,j,s;
  int k=0;
  int tmp_line,tmp_face;
  int nlines_mark=0;
  VOLU *volu,*volu_owned;
  FACEE *face;
  LINE *line;
  //找到所有要加密的单元
  for(i=0;i<mesh->Num_Volus_Global;++i)
    mesh->Volus[i]->Sign=-1;
  for(i=0;i<*num;++i)
    mesh->Volus[M[i]]->Sign=1;
  mesh->Mark=-1;//判断会不会出现type2类型的单元第五条边被标记林的情况，如果变量大于零了，说明出现这种情况了
  for(i=0;i<mesh->Num_Volus_Global;++i)
  {
    volu=mesh->Volus[i];//按误差指示子大小取出单元
    if(volu->Sign>0&&volu->Mark<0)//如果单元是需要被加密的单元，并且在加密过程中没有被标记。
     {
       line=mesh->Lines[volu->Lines[0]];//取出它的第零条边。
       line->Mark=nlines_mark; //第零条边做标记。
       nlines_mark++;          //标记边的总个数加1。
       volu->Mark=1;           //单元做标记
       for(j=0;j<line->Num_Volus_Owned;++j)//对第零条边关联的单元做循环
       {
        volu_owned=mesh->Volus[line->Volus_Owned[j]];//取出第零条边关联的单元
        MarkMesh(mesh,volu_owned,line->Volus_Owned[j],line,volu->Lines[0],&nlines_mark);      	
       }
     }  
  }
  if(mesh->Mark>=0)//出现了这种情况，把网格单元里遇到这种情况的type2类型的单元变为type1类型的
  {
   for(i=0;i<mesh->Num_Volus_Global;++i)
   {
      volu= mesh->Volus[i];
      volu->Mark=-1;
      if(volu->Trans==1)//volu->Trans>0，说明这个单元是type2,且第五条边被标记了
      {
	 volu->Type=1;
         volu->Trans=-1;
      }
   }
  for(i=0;i<mesh->Num_Lines_Global;++i)
     mesh->Lines[i]->Mark = -1;
  mesh->Mark=-1;
  nlines_mark=0;
  for(i=0;i<mesh->Num_Volus_Global;++i)
  {
    volu=mesh->Volus[i];//按误差指示子大小取出单元
    if(volu->Sign>0&&volu->Mark<0)//如果单元是需要被加密的单元，并且在加密过程中没有被标记。
     {
       line=mesh->Lines[volu->Lines[0]];//取出它的第零条边。
       line->Mark=nlines_mark; //第零条边做标记。
       nlines_mark++;          //标记边的总个数加1。
       volu->Mark=1;           //单元做标记
       for(j=0;j<line->Num_Volus_Owned;++j)//对第零条边关联的单元做循环
       {
        volu_owned=mesh->Volus[line->Volus_Owned[j]];//取出第零条边关联的单元
        MarkMesh(mesh,volu_owned,line->Volus_Owned[j],line,volu->Lines[0],&nlines_mark);      	
       }
     }  
  }

}
//没有被加密的单元祖先编号保持不变，父编号变为自身的编号
  for(i=0;i<mesh->Num_Volus_Global;++i)
  {
    volu=mesh->Volus[i];
    if(volu->Mark<0)
       volu->Father=i;
  }


  //首先计算加密后需要的存储空间
  //计算新增加的边数
  int mark_line_on_face=0;
  for(i=0;i<mesh->Num_Faces_Global;++i)
  {   
    face=mesh->Faces[i];
     for(j=0;j<3;++j)
     {
       if(mesh->Lines[face->Lines[j]]->Mark>=0)
       {
	  mark_line_on_face++;
       }
     }
  }     
  int mark_line_on_volu=0;
  int mark_line23_on_volu=0;
  for(i=0;i<mesh->Num_Volus_Global;++i)
  {
     volu=mesh->Volus[i];
     for(j=0;j<5;++j)
     {
	if(mesh->Lines[volu->Lines[j]]->Mark>=0)
	   mark_line_on_volu++;
     }
     if(mesh->Lines[volu->Lines[5]]->Mark>=0)
	   mark_line23_on_volu++;
  }
  int total_vert,total_line, total_face,total_volu;
  total_vert=mesh->Num_Verts_Global+nlines_mark;
  total_line=mesh->Num_Lines_Global+nlines_mark+mark_line_on_face+mark_line23_on_volu;
  total_face=mesh->Num_Faces_Global+mark_line_on_face+mark_line_on_volu+mark_line23_on_volu*3;
  total_volu=mesh->Num_Volus_Global+mark_line_on_volu+mark_line23_on_volu*2;
  
  mesh->Verts=(VERT **)realloc(mesh->Verts,total_vert);
  for(i=mesh->Num_Verts_Global;i<total_vert;++i)
  mesh->Verts[i]=(VERT *)malloc(1);
  
  mesh->Lines=(LINE **)realloc(mesh->Lines,total_line);
  for(i=mesh->Num_Lines_Global;i<total_line;++i)
  mesh->Lines[i]=(LINE *)malloc(1);
  
  mesh->Faces=(FACEE **)realloc(mesh->Faces,total_face);
  for(i=mesh->Num_Faces_Global;i<total_face;++i)
  mesh->Faces[i]=(FACEE *)malloc(1);
  
  mesh->Volus=(VOLU **)realloc(mesh->Volus,total_volu);
  for(i=mesh->Num_Volus_Global;i<total_volu;++i)
  mesh->Volus[i]=(VOLU *)malloc(1);

 //加密程序
 //首先更新边上的新点和新边
 for(i=0;i<mesh->Num_Lines_Global;++i)
 {
   if(mesh->Lines[i]->Mark>=0)
   {
    //新点  
    mesh->Verts[mesh->Num_Verts_Global+mesh->Lines[i]->Mark]->Coord[0]=   
    0.5*(mesh->Verts[mesh->Lines[i]->Verts[0]]->Coord[0]+mesh->Verts[mesh->Lines[i]->Verts[1]]->Coord[0]);
    mesh->Verts[mesh->Num_Verts_Global+mesh->Lines[i]->Mark]->Coord[1]=   
    0.5*(mesh->Verts[mesh->Lines[i]->Verts[0]]->Coord[1]+mesh->Verts[mesh->Lines[i]->Verts[1]]->Coord[1]);
    mesh->Verts[mesh->Num_Verts_Global+mesh->Lines[i]->Mark]->Coord[2]=   
    0.5*(mesh->Verts[mesh->Lines[i]->Verts[0]]->Coord[2]+mesh->Verts[mesh->Lines[i]->Verts[1]]->Coord[2]);
    //两条边
    mesh->Lines[mesh->Num_Lines_Global+mesh->Lines[i]->Mark]->Verts[0]=mesh->Num_Verts_Global
                                                                        +mesh->Lines[i]->Mark; 
    mesh->Lines[mesh->Num_Lines_Global+mesh->Lines[i]->Mark]->Verts[1]=mesh->Lines[i]->Verts[1];
    mesh->Lines[mesh->Num_Lines_Global+mesh->Lines[i]->Mark]->Mark=-1;
    mesh->Lines[mesh->Num_Lines_Global+mesh->Lines[i]->Mark]->Child[0]=-1;
    mesh->Lines[mesh->Num_Lines_Global+mesh->Lines[i]->Mark]->Child[1]=-1;
    mesh->Lines[mesh->Num_Lines_Global+mesh->Lines[i]->Mark]->Child[2]=-1;
    mesh->Lines[mesh->Num_Lines_Global+mesh->Lines[i]->Mark]->ID_Boundary=mesh->Lines[i]->ID_Boundary;
    mesh->Lines[mesh->Num_Lines_Global+mesh->Lines[i]->Mark]->Num_Volus_Owned=-1;
    mesh->Lines[mesh->Num_Lines_Global+mesh->Lines[i]->Mark]->Volus_Owned=NULL;
 
    mesh->Lines[i]->Verts[0]=mesh->Lines[i]->Verts[0];
    mesh->Lines[i]->Verts[1]=mesh->Num_Verts_Global+mesh->Lines[i]->Mark;
    mesh->Lines[i]->Child[0]=i;
    mesh->Lines[i]->Child[1]=mesh->Num_Lines_Global+mesh->Lines[i]->Mark;
    mesh->Lines[i]->Child[2]=mesh->Num_Verts_Global+mesh->Lines[i]->Mark;
   }
 }
 
  int Iter_volu=mesh->Num_Volus_Global;
  int idx_vert=mesh->Num_Verts_Global+nlines_mark;
  int idx_line=mesh->Num_Lines_Global+nlines_mark;
  int idx_face=mesh->Num_Faces_Global;
  int idx_volu=mesh->Num_Volus_Global;
  int v0,v1,v2,v3,l0,l1,l2,l3,l4,l5,f0,f1,f2,f3;
  printf ( "begin three times iteration\n" );
  for(i=0;i<3;++i)
 {//加密三次循环
    for(j=0;j<Iter_volu;++j)//Iter_volu是需要循环的单元个数，每一次循环后都要加上新产生的单元个数
    {//单元遍历
       volu=mesh->Volus[j];
       v0=volu->Verts[0];
       v1=volu->Verts[1];
       v2=volu->Verts[2];
       v3=volu->Verts[3];
       l0=volu->Lines[0];
       l1=volu->Lines[1];
       l2=volu->Lines[2];
       l3=volu->Lines[3];
       l4=volu->Lines[4];
       l5=volu->Lines[5];
       f0=volu->Faces[0];
       f1=volu->Faces[1];
       f2=volu->Faces[2];
       f3=volu->Faces[3];

      if(volu->Mark<0) continue;
      if(volu->Type==0)
      {//type0
          //首先是两个相邻面都需要加密的情况   
          if (  (mesh->Faces[volu->Faces[3]]->Child[0]<0&&mesh->Faces[volu->Faces[2]]->Child[0]<0)
	      ||(mesh->Faces[volu->Faces[3]]->Child[0]<0&&mesh->Faces[volu->Faces[2]]->Child[0]>=0
		  &&(mesh->Faces[volu->Faces[2]]->Child[2]==volu->Lines[0]||mesh->Faces[volu->Faces[2]]
	         ->Child[2]==volu->Lines[2]||mesh->Faces[volu->Faces[2]]->Child[2]==volu->Lines[4]))
	      ||(mesh->Faces[volu->Faces[2]]->Child[0]<0&&mesh->Faces[volu->Faces[3]]->Child[0]>=0
	         &&(mesh->Faces[volu->Faces[3]]->Child[2]==volu->Lines[0]||mesh->Faces[volu->Faces[3]]
		    ->Child[2]==volu->Lines[1]||mesh->Faces[volu->Faces[3]]->Child[2]==volu->Lines[3]) )
	      ||(mesh->Faces[volu->Faces[3]]->Child[0]>=0&&(mesh->Faces[volu->Faces[3]]->Child[2]==
		    volu->Lines[0]||mesh->Faces[volu->Faces[3]]
		    ->Child[2]==volu->Lines[1]||mesh->Faces[volu->Faces[3]]->Child[2]==volu->Lines[3]) 
		 &&mesh->Faces[volu->Faces[2]]->Child[0]>=0&&(mesh->Faces[volu->Faces[2]]->Child[2]==
		    volu->Lines[0]||mesh->Faces[volu->Faces[2]]
	         ->Child[2]==volu->Lines[2]||mesh->Faces[volu->Faces[2]]->Child[2]==volu->Lines[4])
		 )
	     )   
		  //if(mesh->Faces[volu->Faces[3]]->Child[0]<0
	   //&&mesh->Faces[volu->Faces[2]]->Child[0]<0)//第一种情况：两个面都没加密 
	  {//1
	   //更新面上两条新边的信息 
	   mesh->Lines[idx_line]->Verts[1]= mesh->Lines[l0]->Child[2];
	   mesh->Lines[idx_line]->Mark=-1;
	   mesh->Lines[idx_line]->ID_Boundary=mesh->Faces[f3]->ID_Boundary;
	   mesh->Lines[idx_line]->Num_Volus_Owned=-1;
	   mesh->Lines[idx_line]->Volus_Owned=NULL;
           mesh->Lines[idx_line]->Child[0]=-1;
	   mesh->Lines[idx_line]->Child[1]=-1;
	   mesh->Lines[idx_line]->Child[2]=-1;
	   mesh->Lines[idx_line++]->Verts[0]=v2;
	   mesh->Lines[idx_line]->Verts[1]= mesh->Lines[l0]->Child[2];
	   mesh->Lines[idx_line]->Mark=-1;
	   mesh->Lines[idx_line]->Num_Volus_Owned=-1;
	   mesh->Lines[idx_line]->Volus_Owned=NULL;
	   mesh->Lines[idx_line]->ID_Boundary=mesh->Faces[f2]->ID_Boundary;
	   mesh->Lines[idx_line]->Child[0]=-1;
	   mesh->Lines[idx_line]->Child[1]=-1;
	   mesh->Lines[idx_line]->Child[2]=-1;
           mesh->Lines[idx_line++]->Verts[0]=v3;
           //更新第一个面上的信息
	   mesh->Faces[idx_face]->Verts[0]=v1;
           mesh->Faces[idx_face]->Verts[1]=v2;
           mesh->Faces[idx_face]->Verts[2]=mesh->Lines[l0]->Child[2];

	   mesh->Faces[idx_face]->Child[0]=-1;
           mesh->Faces[idx_face]->Child[1]=-1;
           mesh->Faces[idx_face]->Child[2]=-1;
           mesh->Faces[idx_face]->Old_Child[0]=-1;
           mesh->Faces[idx_face]->Old_Child[1]=-1;

	   mesh->Faces[idx_face]->ID_Boundary=mesh->Faces[f3]->ID_Boundary;

           mesh->Faces[idx_face]->Lines[0]=idx_line-2;
           mesh->Faces[idx_face]->Lines[1]=
	   v0==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
           mesh->Faces[idx_face++]->Lines[2]=l3;

	   mesh->Faces[f3]->Verts[0]=volu->Verts[0];
           mesh->Faces[f3]->Verts[1]=mesh->Lines[l0]->Child[2];
           mesh->Faces[f3]->Verts[2]=v2;
         
	   mesh->Faces[f3]->Lines[0]=idx_line-2;
           mesh->Faces[f3]->Lines[1]=l1;
           mesh->Faces[f3]->Lines[2]=
           v0==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];

	   mesh->Faces[f3]->Old_Child[0]=mesh->Faces[f3]->Child[1];
	   mesh->Faces[f3]->Old_Child[1]=mesh->Faces[f3]->Child[2];
	   mesh->Faces[f3]->Child[0]=f3;
           mesh->Faces[f3]->Child[1]=idx_face-1;
           mesh->Faces[f3]->Child[2]=idx_line-2;
	   //更新第二个面上的信息
           mesh->Faces[idx_face]->Verts[0]=v1;
           mesh->Faces[idx_face]->Verts[1]=v3;
           mesh->Faces[idx_face]->Verts[2]=mesh->Lines[l0]->Child[2];

	   mesh->Faces[idx_face]->Child[0]=-1;
           mesh->Faces[idx_face]->Child[1]=-1;
           mesh->Faces[idx_face]->Child[2]=-1;
      	   mesh->Faces[idx_face]->Old_Child[0]=-1;
           mesh->Faces[idx_face]->Old_Child[1]=-1;
           mesh->Faces[idx_face]->ID_Boundary=mesh->Faces[f2]->ID_Boundary;//?????????????????

           mesh->Faces[idx_face]->Lines[0]=idx_line-1;
           mesh->Faces[idx_face]->Lines[1]=
           v0==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
           mesh->Faces[idx_face++]->Lines[2]=l4;

	   mesh->Faces[f2]->Verts[0]=v0;
           mesh->Faces[f2]->Verts[1]=mesh->Lines[l0]->Child[2];
           mesh->Faces[f2]->Verts[2]=volu->Verts[3];
         
           mesh->Faces[f2]->Old_Child[0]=mesh->Faces[f2]->Child[1];
           mesh->Faces[f2]->Old_Child[1]=mesh->Faces[f2]->Child[2];
	   mesh->Faces[f2]->Child[0]=f2;
           mesh->Faces[f2]->Child[1]=idx_face-1;
           mesh->Faces[f2]->Child[2]=idx_line-1;

	   mesh->Faces[f2]->Lines[0]=idx_line-1;
           mesh->Faces[f2]->Lines[1]=l2;
           mesh->Faces[f2]->Lines[2]=
           v0==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];

           //更新中间的新面的信息
	   mesh->Faces[idx_face]->Verts[0]=mesh->Lines[l0]->Child[2];
           mesh->Faces[idx_face]->Verts[1]=v2;
           mesh->Faces[idx_face]->Verts[2]=v3;
         
	   mesh->Faces[idx_face]->Child[0]=-1;
           mesh->Faces[idx_face]->Child[1]=-1;
           mesh->Faces[idx_face]->Child[2]=-1;
           mesh->Faces[idx_face]->Old_Child[0]=-1;
           mesh->Faces[idx_face]->Old_Child[1]=-1;

	   mesh->Faces[idx_face]->ID_Boundary=0;

	   mesh->Faces[idx_face]->Lines[0]=l5;
           mesh->Faces[idx_face]->Lines[1]=idx_line-1;
           mesh->Faces[idx_face++]->Lines[2]=idx_line-2;
	      //更新单元信息
	   mesh->Volus[idx_volu]->Verts[0]=v1;
	   mesh->Volus[idx_volu]->Verts[1]=v3;
	   mesh->Volus[idx_volu]->Verts[2]=v2;
           mesh->Volus[idx_volu]->Verts[3]=mesh->Lines[l0]->Child[2];
           
	   mesh->Volus[idx_volu]->Lines[0]=l4;
	   mesh->Volus[idx_volu]->Lines[1]=l3;
	   mesh->Volus[idx_volu]->Lines[2]=
	   volu->Verts[0]==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];//yong volu->verts shi fou you wen ti????
 	   mesh->Volus[idx_volu]->Lines[3]=l5;
           mesh->Volus[idx_volu]->Lines[4]=idx_line-1;
 	   mesh->Volus[idx_volu]->Lines[5]=idx_line-2;

	   if(i==0)
	   {
            mesh->Volus[idx_volu]->Father=j;
            mesh->Volus[j]->Father=j;
	   }
           else
            mesh->Volus[idx_volu]->Father=mesh->Volus[j]->Father;
           
	   mesh->Volus[idx_volu]->Ancestor=mesh->Volus[j]->Ancestor;

           mesh->Volus[idx_volu]->Type=1;

	   mesh->Volus[idx_volu]->Faces[0]=idx_face-1;
	   mesh->Volus[idx_volu]->Faces[1]=idx_face-3;
	   mesh->Volus[idx_volu]->Faces[2]=idx_face-2;
           mesh->Volus[idx_volu++]->Faces[3]=f0;

	   mesh->Volus[j]->Verts[0]=v0;
	   mesh->Volus[j]->Verts[1]=v2;
	   mesh->Volus[j]->Verts[2]=v3;
	   mesh->Volus[j]->Verts[3]=mesh->Lines[l0]->Child[2];


	   tmp_line=volu->Verts[0]==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];
           mesh->Volus[j]->Lines[0]=l1;
	   mesh->Volus[j]->Lines[1]=l2;
	   mesh->Volus[j]->Lines[2]=tmp_line;
	   mesh->Volus[j]->Lines[3]=l5;
           mesh->Volus[j]->Lines[4]=idx_line-2;
 	   mesh->Volus[j]->Lines[5]=idx_line-1;

	   ////tmp_face=volu->Faces[1];
	   //tmp_face=f1;
	   mesh->Volus[j]->Faces[0]=idx_face-1;
	   mesh->Volus[j]->Faces[1]=f2;
	   mesh->Volus[j]->Faces[2]=f3;
	   mesh->Volus[j]->Faces[3]=f1;

           mesh->Volus[j]->Type=1;

	   mesh->Volus[j]->Child[0]=j;
	   mesh->Volus[j]->Child[1]=idx_volu-1;

	  
	   if(MarkJudge(mesh,mesh->Volus[j])<0)//有问题？？？
	   {
            mesh->Volus[j]->Mark=-1;
	   }
           if(MarkJudge(mesh,mesh->Volus[idx_volu-1])>=0) 
	   {
            mesh->Volus[idx_volu-1]->Mark=1;
	   }
	   else
	   mesh->Volus[idx_volu-1]->Mark=-1;

	  }//1
           //-----------------------------------------------------------------------------------
           //第二种情况：第一个面没加密，第二个面加密
	 else if(  (mesh->Faces[volu->Faces[3]]->Child[0]<0
		     &&mesh->Faces[volu->Faces[2]]->Child[0]>=0
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[0]
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[2]
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[4])
		  ||(mesh->Faces[volu->Faces[3]]->Child[0]>=0
		     &&(mesh->Faces[volu->Faces[3]]->Child[2]==
		    volu->Lines[0]||mesh->Faces[volu->Faces[3]]
		    ->Child[2]==volu->Lines[1]||mesh->Faces[volu->Faces[3]]->Child[2]==volu->Lines[3]) 
		     &&mesh->Faces[volu->Faces[2]]->Child[0]>=0
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[0]
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[2]
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[4]
		    )
	       )

	   //else if(mesh->Faces[volu->Faces[3]]->Child[0]<0&&mesh->Faces[volu->Faces[2]]->Child[0]>=0)
	   {//2 
               mesh->Lines[idx_line]->Verts[1]=mesh->Lines[l0]->Child[2];
	       mesh->Lines[idx_line]->Child[0]=-1;
	       mesh->Lines[idx_line]->Child[1]=-1;
               mesh->Lines[idx_line]->Child[2]=-1;
	       mesh->Lines[idx_line]->Num_Volus_Owned=-1;
	       mesh->Lines[idx_line]->Volus_Owned=NULL;
               mesh->Lines[idx_line]->Mark=-1;
	       mesh->Lines[idx_line]->ID_Boundary=mesh->Faces[f3]->ID_Boundary;
	       mesh->Lines[idx_line++]->Verts[0]=v2;
	       //更新第一个面上的新面
               mesh->Faces[f3]->Verts[0]=v0;
               mesh->Faces[f3]->Verts[1]=mesh->Lines[l0]->Child[2];
               mesh->Faces[f3]->Verts[2]=v2;
             
	       mesh->Faces[f3]->Old_Child[0]=mesh->Faces[f3]->Child[1];
               mesh->Faces[f3]->Old_Child[1]=mesh->Faces[f3]->Child[2];
               mesh->Faces[f3]->Child[0]=f3;
               mesh->Faces[f3]->Child[1]=idx_face;
               mesh->Faces[f3]->Child[2]=idx_line-1;
	       mesh->Faces[f3]->Lines[0]=idx_line-1;
               mesh->Faces[f3]->Lines[1]=l1;
               mesh->Faces[f3]->Lines[2]=
	       volu->Verts[0]==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];

               mesh->Faces[idx_face]->Verts[0]=volu->Verts[1];
               mesh->Faces[idx_face]->Verts[1]=volu->Verts[2];
               mesh->Faces[idx_face]->Verts[2]=mesh->Lines[l0]->Child[2];

               mesh->Faces[idx_face]->Child[0]=-1;
               mesh->Faces[idx_face]->Child[1]=-1;
               mesh->Faces[idx_face]->Child[2]=-1;

	       mesh->Faces[idx_face]->Old_Child[0]=-1;
               mesh->Faces[idx_face]->Old_Child[1]=-1;

               mesh->Faces[idx_face]->ID_Boundary=mesh->Faces[f3]->ID_Boundary;

               mesh->Faces[idx_face]->Lines[0]=idx_line-1;
               mesh->Faces[idx_face]->Lines[1]=
	       v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
               mesh->Faces[idx_face++]->Lines[2]=l3;
	       //更新中间新面的信息
               mesh->Faces[idx_face]->Verts[0]=mesh->Lines[l0]->Child[2];
               mesh->Faces[idx_face]->Verts[1]=v2;
               mesh->Faces[idx_face]->Verts[2]=v3;
              
               mesh->Faces[idx_face]->Child[0]=-1;
               mesh->Faces[idx_face]->Child[1]=-1;
               mesh->Faces[idx_face]->Child[2]=-1;

       	       mesh->Faces[idx_face]->Old_Child[0]=-1;
               mesh->Faces[idx_face]->Old_Child[1]=-1;

               mesh->Faces[idx_face]->ID_Boundary=0;

	       mesh->Faces[idx_face]->Lines[0]=l5;
               if(   (mesh->Faces[f2]->Old_Child[1]>=0)&&(mesh->Faces[f2]->Old_Child[1]!=l0)&&
		     (mesh->Faces[f2]->Old_Child[1]!=l2)&&(mesh->Faces[f2]->Old_Child[1]!=l4)   )
	       {
	         mesh->Faces[idx_face]->Lines[1]=mesh->Faces[f2]->Old_Child[1];
	       } 
	       else
	       {
	         mesh->Faces[idx_face]->Lines[1]=mesh->Faces[f2]->Child[2];
	       } 
	       //mesh->Faces[idx_face]->Lines[1]=mesh->Faces[f2]->Child[2];
               mesh->Faces[idx_face++]->Lines[2]=idx_line-1;
          
	       //更新两个单元的信息
               mesh->Volus[idx_volu]->Verts[0]=v1;
               mesh->Volus[idx_volu]->Verts[1]=v3;
               mesh->Volus[idx_volu]->Verts[2]=v2;
               mesh->Volus[idx_volu]->Verts[3]=mesh->Lines[l0]->Child[2];
               
	       mesh->Volus[idx_volu]->Lines[0]=volu->Lines[4];
               mesh->Volus[idx_volu]->Lines[1]=volu->Lines[3];
               mesh->Volus[idx_volu]->Lines[2]=
	       volu->Verts[0]==mesh->Lines[volu->Lines[0]]->Verts[0]?
	       mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
               mesh->Volus[idx_volu]->Lines[3]=l5;
               /*
	       if(mesh->Faces[f2]->Old_Child[1]<0)
               mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f2]->Child[2];
	       else
               mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f2]->Old_Child[1];
	       */
	       mesh->Volus[idx_volu]->Lines[5]=idx_line-1;
               
	       if(   (mesh->Faces[f2]->Old_Child[1]>=0)&&(mesh->Faces[f2]->Old_Child[1]!=l0)&&
		     (mesh->Faces[f2]->Old_Child[1]!=l2)&&(mesh->Faces[f2]->Old_Child[1]!=l4)   )
               mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f2]->Old_Child[1];
	       else
               mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f2]->Child[2];

	       //mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f2]->Child[2];

	       if(i==0)
	       {
                mesh->Volus[idx_volu]->Father=j;
                mesh->Volus[j]->Father=j;
	       }
	       else
	       mesh->Volus[idx_volu]->Father=mesh->Volus[j]->Father;
               mesh->Volus[idx_volu]->Ancestor=mesh->Volus[j]->Ancestor;
               mesh->Volus[idx_volu]->Type=1;

               mesh->Volus[idx_volu]->Faces[0]=idx_face-1;
               mesh->Volus[idx_volu]->Faces[1]=idx_face-2;
               /*
	       if(mesh->Faces[f2]->Old_Child[1]<0)
	       {
	       mesh->Volus[idx_volu]->Faces[2]=
               v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	       mesh->Faces[f2]->Child[1]:mesh->Faces[f2]->Child[0];
	       }
	       else
	      {
	       mesh->Volus[idx_volu]->Faces[2]=
               v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	       mesh->Faces[f2]->Old_Child[0]:mesh->Faces[f2]->Child[0];
	       }
	       */

	      if(   (mesh->Faces[f2]->Old_Child[1]>=0)&&(mesh->Faces[f2]->Old_Child[1]!=l0)&&
		     (mesh->Faces[f2]->Old_Child[1]!=l2)&&(mesh->Faces[f2]->Old_Child[1]!=l4)   )
	       {
	        mesh->Volus[idx_volu]->Faces[2]=
                v1==mesh->Faces[mesh->Faces[f2]->Old_Child[0]]->Verts[0]?
	        mesh->Faces[f2]->Old_Child[0]:mesh->Faces[f2]->Child[0];
	       }
	      else
	       {
	        mesh->Volus[idx_volu]->Faces[2]=
                v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	        mesh->Faces[f2]->Child[1]:mesh->Faces[f2]->Child[0];
	       }

               //mesh->Volus[idx_volu]->Faces[2]=
               //v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	       //mesh->Faces[f2]->Child[1]:mesh->Faces[f2]->Child[0];

	       mesh->Volus[idx_volu++]->Faces[3]=f0;
	       ////mesh->Volus[idx_volu++]->Faces[3]=volu->Faces[0];

	       mesh->Volus[j]->Verts[0]=v0;
               mesh->Volus[j]->Verts[1]=v2;
               mesh->Volus[j]->Verts[2]=v3;
               mesh->Volus[j]->Verts[3]=mesh->Lines[l0]->Child[2];

	       tmp_line= v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];
               mesh->Volus[j]->Lines[0]=l1;
               mesh->Volus[j]->Lines[1]=l2;
               mesh->Volus[j]->Lines[2]=tmp_line;
               mesh->Volus[j]->Lines[3]=l5;
               mesh->Volus[j]->Lines[4]=idx_line-1;
               /*
	      
	       if(mesh->Faces[f2]->Old_Child[1]<0)//gaid
               mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f2]->Child[2];
	       else
               mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f2]->Old_Child[1];
               */
               if(   (mesh->Faces[f2]->Old_Child[1]>=0)&&(mesh->Faces[f2]->Old_Child[1]!=l0)&&
		     (mesh->Faces[f2]->Old_Child[1]!=l2)&&(mesh->Faces[f2]->Old_Child[1]!=l4)   )
               mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Old_Child[1];
	       else
               mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Child[2];
	      

	       ////tmp_face=volu->Faces[1];
	       tmp_face=f1;
	       mesh->Volus[j]->Faces[0]=idx_face-1;
	       /*
	       if(mesh->Faces[f2]->Old_Child[1]<0)
	       {
	       mesh->Volus[j]->Faces[1]=
	       v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	       mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Child[1];
               mesh->Volus[j]->Faces[2]=f3;
               mesh->Volus[j]->Faces[3]=tmp_face;
	       }
	       else
	       {
	       mesh->Volus[j]->Faces[1]=
               v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	       mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Old_Child[0];
               mesh->Volus[j]->Faces[2]=f3;
               mesh->Volus[j]->Faces[3]=tmp_face;
	       }
	       */
	      
	       if(   (mesh->Faces[f2]->Old_Child[1]>=0)&&(mesh->Faces[f2]->Old_Child[1]!=l0)&&
		     (mesh->Faces[f2]->Old_Child[1]!=l2)&&(mesh->Faces[f2]->Old_Child[1]!=l4)   )
	       {
	        mesh->Volus[j]->Faces[1]=
                v1==mesh->Faces[mesh->Faces[f2]->Old_Child[0]]->Verts[0]?
	        mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Old_Child[0];
	       }
	      else
	       {
	        mesh->Volus[j]->Faces[1]=
                v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	        mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Child[1];
	       }
              

	       mesh->Volus[j]->Faces[2]=f3;
               mesh->Volus[j]->Faces[3]=tmp_face;

	       mesh->Volus[j]->Type=1;

               mesh->Volus[j]->Child[0]=j;
	       mesh->Volus[j]->Child[1]=idx_volu-1;


                 //ShowVolu(mesh, volu);
	       //ReConstruct(mesh->Volus[idx_volu-1]);
               //ReConstruct(mesh->Volus[idx_volu-2]);
           if(MarkJudge(mesh,mesh->Volus[j])<0)
	   {
	      mesh->Volus[j]->Mark=-1; 
	   }
	   
           if(MarkJudge(mesh,mesh->Volus[idx_volu-1])>=0) 
	   {
              mesh->Volus[idx_volu-1]->Mark=1;
	   }
           else
	   mesh->Volus[idx_volu-1]->Mark=-1;

        }//2
         //---------------------------------------------------------------------------------------------
	  else if ( (mesh->Faces[volu->Faces[3]]->Child[0]>=0
		     &&mesh->Faces[volu->Faces[3]]->Child[2]!=volu->Lines[0]
		     &&mesh->Faces[volu->Faces[3]]->Child[2]!=volu->Lines[1]
		     &&mesh->Faces[volu->Faces[3]]->Child[2]!=volu->Lines[3]
		     )
		  &&(mesh->Faces[volu->Faces[2]]->Child[0]>=0
		    &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[0]
		    &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[2]
		    &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[4]
		    )
		  ) 
	  //else if(mesh->Faces[volu->Faces[3]]->Child[0]>=0&&mesh->Faces[volu->Faces[2]]->Child[0]>=0)
	  //第三种情况：两个面都加密 
	  {//3
	   mesh->Faces[idx_face]->Verts[0]=mesh->Lines[l0]->Child[2];
           mesh->Faces[idx_face]->Verts[1]=v2;
           mesh->Faces[idx_face]->Verts[2]=v3;
         
	   mesh->Faces[idx_face]->Child[0]=-1;
           mesh->Faces[idx_face]->Child[1]=-1;
           mesh->Faces[idx_face]->Child[2]=-1;
 	   mesh->Faces[idx_face]->Old_Child[0]=-1;
           mesh->Faces[idx_face]->Old_Child[1]=-1;

           mesh->Faces[idx_face]->ID_Boundary=0;

	   mesh->Faces[idx_face]->Lines[0]=volu->Lines[5];
           if(   (mesh->Faces[f2]->Old_Child[1]>=0)&&(mesh->Faces[f2]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f2]->Old_Child[1]!=l2)&&(mesh->Faces[f2]->Old_Child[1]!=l4)   )
             {
	       mesh->Faces[idx_face]->Lines[1]=mesh->Faces[f2]->Old_Child[1];
	     }
	   else
	     {
	       mesh->Faces[idx_face]->Lines[1]=mesh->Faces[f2]->Child[2];
	     }
	   if(   (mesh->Faces[f3]->Old_Child[1]>=0)&&(mesh->Faces[f3]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f3]->Old_Child[1]!=l1)&&(mesh->Faces[f3]->Old_Child[1]!=l3)   )
             {
               mesh->Faces[idx_face++]->Lines[2]=mesh->Faces[f3]->Old_Child[1];
	     }
	   else
	     {
               mesh->Faces[idx_face++]->Lines[2]=mesh->Faces[f3]->Child[2];
	     }
	   //mesh->Faces[idx_face]->Lines[1]=mesh->Faces[f2]->Child[2];
           //mesh->Faces[idx_face++]->Lines[2]=mesh->Faces[f3]->Child[2];

           //更新单元信息
	   mesh->Volus[idx_volu]->Verts[0]=v1;
	   mesh->Volus[idx_volu]->Verts[1]=v3;
	   mesh->Volus[idx_volu]->Verts[2]=v2;
           mesh->Volus[idx_volu]->Verts[3]=mesh->Lines[l0]->Child[2];

 	   mesh->Volus[idx_volu]->Lines[0]=l4;
           mesh->Volus[idx_volu]->Lines[1]=l3;
           mesh->Volus[idx_volu]->Lines[2]=
           v0==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
	   mesh->Volus[idx_volu]->Lines[3]=l5;
	   //mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[volu->Faces[2]]->Child[2];
	   //mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[volu->Faces[3]]->Child[2];

           /*
           if(mesh->Faces[volu->Faces[2]]->Old_Child[0]<0)    
	      mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f2]->Child[2];
	   else
              mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f2]->Old_Child[1];

           if(mesh->Faces[f3]->Old_Child[1]<0)//gaid     
	      mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f3]->Child[2];
	   else
              mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f3]->Old_Child[1];
           */
	   if(   (mesh->Faces[f2]->Old_Child[1]>=0)&&(mesh->Faces[f2]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f2]->Old_Child[1]!=l2)&&(mesh->Faces[f2]->Old_Child[1]!=l4)   )
             {
	         mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f2]->Old_Child[1];
	     }
	   else
	     {
	         mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f2]->Child[2];
	     }
	    if(   (mesh->Faces[f3]->Old_Child[1]>=0)&&(mesh->Faces[f3]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f3]->Old_Child[1]!=l1)&&(mesh->Faces[f3]->Old_Child[1]!=l3)   )
             {
	         mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f3]->Old_Child[1];
	     }
	   else
	     {
	         mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f3]->Child[2];
	     }

	   //mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f2]->Child[2];
	   //mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f3]->Child[2];

	   if(i==0)
	   {
            mesh->Volus[idx_volu]->Father=j;
            mesh->Volus[j]->Father=j;
	   }
           else
           mesh->Volus[idx_volu]->Father=mesh->Volus[j]->Father;
           
	   mesh->Volus[idx_volu]->Ancestor=mesh->Volus[j]->Ancestor;

           mesh->Volus[idx_volu]->Type=1;

	   mesh->Volus[idx_volu]->Faces[0]=idx_face-1;
	   /*
	   if(mesh->Faces[f3]->Old_Child[0]<0)
	   mesh->Volus[idx_volu]->Faces[1]=
           volu->Verts[0]==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	   mesh->Faces[f3]->Child[1]:mesh->Faces[f3]->Child[0];
	   else
	   mesh->Volus[idx_volu]->Faces[1]=
           v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	   mesh->Faces[f3]->Old_Child[0]:mesh->Faces[f3]->Child[0];
           */
   	   if(   (mesh->Faces[f3]->Old_Child[1]>=0)&&(mesh->Faces[f3]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f3]->Old_Child[1]!=l1)&&(mesh->Faces[f3]->Old_Child[1]!=l3)   )
	   {
	      mesh->Volus[idx_volu]->Faces[1]=
              v1==mesh->Faces[mesh->Faces[f3]->Old_Child[0]]->Verts[0]?
	      mesh->Faces[f3]->Old_Child[0]:mesh->Faces[f3]->Child[0];
	   }
	   else
	   {
              mesh->Volus[idx_volu]->Faces[1]=
              v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	      mesh->Faces[f3]->Child[1]:mesh->Faces[f3]->Child[0];
	   }
           /*
	   if(mesh->Faces[f2]->Old_Child[0]<0)
	   mesh->Volus[idx_volu]->Faces[2]=
           v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	   mesh->Faces[f2]->Child[1]:mesh->Faces[f2]->Child[0];
	   else
           mesh->Volus[idx_volu]->Faces[2]=
           v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	   mesh->Faces[f2]->Old_Child[0]:mesh->Faces[f2]->Child[0];
           */
           if(   (mesh->Faces[f2]->Old_Child[1]>=0)&&(mesh->Faces[f2]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f2]->Old_Child[1]!=l2)&&(mesh->Faces[f2]->Old_Child[1]!=l4)   )
	   {
	      mesh->Volus[idx_volu]->Faces[2]=
              v1==mesh->Faces[mesh->Faces[f2]->Old_Child[0]]->Verts[0]?
	      mesh->Faces[f2]->Old_Child[0]:mesh->Faces[f2]->Child[0];
	   }
	   else
	   {
              mesh->Volus[idx_volu]->Faces[2]=
              v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	      mesh->Faces[f2]->Child[1]:mesh->Faces[f2]->Child[0];
	   }
	   mesh->Volus[idx_volu++]->Faces[3]=f0;
           
	   mesh->Volus[j]->Verts[0]=v0;
	   mesh->Volus[j]->Verts[1]=v2;
	   mesh->Volus[j]->Verts[2]=v3;
           mesh->Volus[j]->Verts[3]=mesh->Lines[l0]->Child[2];
           tmp_line=volu->Verts[0]==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];
	   //mesh->Volus[j]->Lines[4]=mesh->Faces[volu->Faces[3]]->Child[2];
	   mesh->Volus[j]->Lines[3]=l5;
           //mesh->Volus[j]->Lines[5]=mesh->Faces[volu->Faces[2]]->Child[2];
	   mesh->Volus[j]->Lines[0]=l1;
	   mesh->Volus[j]->Lines[1]=l2;
           mesh->Volus[j]->Lines[2]=tmp_line;
	   /*   
	   if(mesh->Faces[volu->Faces[2]]->Old_Child[1]<0)//gaid    
	     mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Child[2];
   	   else
             mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Old_Child[1];
	   if(mesh->Faces[volu->Faces[3]]->Old_Child[1]<0)      
	     mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Child[2];
	   else
             mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Old_Child[1];
           */
	    
	   if(   (mesh->Faces[f2]->Old_Child[1]>=0)&&(mesh->Faces[f2]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f2]->Old_Child[1]!=l2)&&(mesh->Faces[f2]->Old_Child[1]!=l4)   )
             {
	         mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Old_Child[1];
	     }
	   else
	     {
	         mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Child[2];
	     }
	    if(   (mesh->Faces[f3]->Old_Child[1]>=0)&&(mesh->Faces[f3]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f3]->Old_Child[1]!=l1)&&(mesh->Faces[f3]->Old_Child[1]!=l3)   )
             {
	         mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Old_Child[1];
	     }
	   else
	     {
	         mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Child[2];
	     }
	   //mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Child[2];
	   //mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Child[2];

	   tmp_face=f1;
	   mesh->Volus[j]->Faces[0]=idx_face-1;
	   /*
	   if(mesh->Faces[f2]->Old_Child[1]<0)
	   mesh->Volus[j]->Faces[1]=
           v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	   mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Child[1];
	   else
           mesh->Volus[j]->Faces[1]=
           volu->Verts[0]==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	   mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Old_Child[0];
           if(mesh->Faces[f3]->Child[1]<0)
	   mesh->Volus[j]->Faces[2]=
           v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	   mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Child[1];
	   else
           mesh->Volus[j]->Faces[2]=
           v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	   mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Old_Child[0];
           */
           mesh->Volus[j]->Faces[3]=tmp_face;
	   if(   (mesh->Faces[f2]->Old_Child[1]>=0)&&(mesh->Faces[f2]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f2]->Old_Child[1]!=l2)&&(mesh->Faces[f2]->Old_Child[1]!=l4)   )
	   {
             mesh->Volus[j]->Faces[1]=
             v0==mesh->Faces[mesh->Faces[f2]->Old_Child[0]]->Verts[0]?
	     mesh->Faces[f2]->Old_Child[0]:mesh->Faces[f2]->Child[0];
	   } 
	   else
	   {
             mesh->Volus[j]->Faces[1]=
             v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	     mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Child[1];
	   }

	   if(   (mesh->Faces[f3]->Old_Child[1]>=0)&&(mesh->Faces[f3]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f3]->Old_Child[1]!=l1)&&(mesh->Faces[f3]->Old_Child[1]!=l3)   )
	   {
	     mesh->Volus[j]->Faces[2]=
             v0==mesh->Faces[mesh->Faces[f3]->Old_Child[0]]->Verts[0]?
	     mesh->Faces[f3]->Old_Child[0]:mesh->Faces[f3]->Child[0];
	   }  
	   else
	   {
	     mesh->Volus[j]->Faces[2]=
             v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	     mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Child[1];
	   }

	   mesh->Volus[j]->Type=1;

	   mesh->Volus[j]->Child[0]=j;
	   mesh->Volus[j]->Child[1]=idx_volu-1;

 	   if(MarkJudge(mesh,mesh->Volus[j])<0)
	   {
	    mesh->Volus[j]->Mark=-1;
	   }
           if(MarkJudge(mesh,mesh->Volus[idx_volu-1])>=0) 
	   {
            mesh->Volus[idx_volu-1]->Mark=1;
	   }
           else
	   mesh->Volus[idx_volu-1]->Mark=-1;

	   }//3
	   //-------------------------------------------------------------------------------------------
           //第四种情况：第一个面(012)加密，第二个面没加密 
	   else
	   {//4
	       mesh->Lines[idx_line]->Verts[1]=mesh->Lines[l0]->Child[2];
	       //printf ( "h440\n" );
	       mesh->Lines[idx_line]->Child[0]=-1;
	       mesh->Lines[idx_line]->Child[1]=-1;
               mesh->Lines[idx_line]->Child[2]=-1;
	       mesh->Lines[idx_line]->Num_Volus_Owned=-1;
	       mesh->Lines[idx_line]->Volus_Owned=NULL;
               mesh->Lines[idx_line]->Mark=-1;
	       mesh->Lines[idx_line]->ID_Boundary=mesh->Faces[f2]->ID_Boundary;
	       mesh->Lines[idx_line++]->Verts[0]=v3;
	       //printf ( "h441\n" );   	       
	       //更新第二个面上的新面
               mesh->Faces[f2]->Verts[0]=v0;
               mesh->Faces[f2]->Verts[1]=mesh->Lines[l0]->Child[2];
               mesh->Faces[f2]->Verts[2]=v3;
             
               mesh->Faces[f2]->Old_Child[0]=mesh->Faces[f2]->Child[1];
               mesh->Faces[f2]->Old_Child[1]=mesh->Faces[f2]->Child[2];
               mesh->Faces[f2]->Child[0]=f2;
               mesh->Faces[f2]->Child[1]=idx_face;
               mesh->Faces[f2]->Child[2]=idx_line-1;

	       mesh->Faces[f2]->Lines[0]=idx_line-1;
               mesh->Faces[f2]->Lines[1]=l2;
               mesh->Faces[f2]->Lines[2]=
	       volu->Verts[0]==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];
	       //printf ( "h442\n" );
               mesh->Faces[idx_face]->Verts[0]=v1;
               mesh->Faces[idx_face]->Verts[1]=v3;
               mesh->Faces[idx_face]->Verts[2]=mesh->Lines[l0]->Child[2];

               mesh->Faces[idx_face]->Child[0]=-1;
               mesh->Faces[idx_face]->Child[1]=-1;
               mesh->Faces[idx_face]->Child[2]=-1;

	       mesh->Faces[idx_face]->Old_Child[0]=-1;
               mesh->Faces[idx_face]->Old_Child[1]=-1;

               mesh->Faces[idx_face]->ID_Boundary=mesh->Faces[f2]->ID_Boundary;

               mesh->Faces[idx_face]->Lines[0]=idx_line-1;
               mesh->Faces[idx_face]->Lines[1]=
	       v0==mesh->Lines[volu->Lines[0]]->Verts[0]?
	       mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
               mesh->Faces[idx_face++]->Lines[2]=l4;

	       //更新中间新面的信息
               mesh->Faces[idx_face]->Verts[0]=mesh->Lines[l0]->Child[2];
               mesh->Faces[idx_face]->Verts[1]=v2;
               mesh->Faces[idx_face]->Verts[2]=v3;
              
               mesh->Faces[idx_face]->Child[0]=-1;
               mesh->Faces[idx_face]->Child[1]=-1;
               mesh->Faces[idx_face]->Child[2]=-1;

	       mesh->Faces[idx_face]->Old_Child[0]=-1;
               mesh->Faces[idx_face]->Old_Child[1]=-1;

               mesh->Faces[idx_face]->ID_Boundary=0;

	       mesh->Faces[idx_face]->Lines[0]=l5;
               mesh->Faces[idx_face]->Lines[1]=idx_line-1;
               if(   (mesh->Faces[f3]->Old_Child[1]>=0)&&(mesh->Faces[f3]->Old_Child[1]!=l0)&&
	             (mesh->Faces[f3]->Old_Child[1]!=l1)&&(mesh->Faces[f3]->Old_Child[1]!=l3)   )
	       {
	         mesh->Faces[idx_face++]->Lines[2]=mesh->Faces[f3]->Old_Child[1];
	       }
	       else
	       {
	         mesh->Faces[idx_face++]->Lines[2]=mesh->Faces[f3]->Child[2];
	       } 
	       //更新两个单元的信息
               mesh->Volus[idx_volu]->Verts[0]=v1;
               mesh->Volus[idx_volu]->Verts[1]=v3;
               mesh->Volus[idx_volu]->Verts[2]=v2;
               mesh->Volus[idx_volu]->Verts[3]=mesh->Lines[l0]->Child[2];

               mesh->Volus[idx_volu]->Lines[0]=l4;
               mesh->Volus[idx_volu]->Lines[1]=l3;
	       mesh->Volus[idx_volu]->Lines[2]=
	       v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
               
	       mesh->Volus[idx_volu]->Lines[3]=l5;
               mesh->Volus[idx_volu]->Lines[4]=idx_line-1;
               /*  
	       if(mesh->Faces[f3]->Old_Child[1]<0)//gaid
	       mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f3]->Child[2];
               else
               mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f3]->Old_Child[1];
               */
	       if(   (mesh->Faces[f3]->Old_Child[1]>=0)&&(mesh->Faces[f3]->Old_Child[1]!=l0)&&
	             (mesh->Faces[f3]->Old_Child[1]!=l1)&&(mesh->Faces[f3]->Old_Child[1]!=l3)   )
	       {
                  mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f3]->Old_Child[1];
	       }
	       else
	       {
	          mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f3]->Child[2];
	       }
	       //mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f3]->Child[2];
    	       
	       
	       if(i==0)
   	       {
                mesh->Volus[idx_volu]->Father=j;
                mesh->Volus[j]->Father=j;
	       }
               else
               mesh->Volus[idx_volu]->Father=mesh->Volus[j]->Father;
               mesh->Volus[idx_volu]->Ancestor=mesh->Volus[j]->Ancestor;

               mesh->Volus[idx_volu]->Type=1;

               mesh->Volus[idx_volu]->Faces[0]=idx_face-1;
               /*
	       if(mesh->Faces[f3]->Old_Child[1]<0)
               mesh->Volus[idx_volu]->Faces[1]=
               v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	       mesh->Faces[f3]->Child[1]:mesh->Faces[f3]->Child[0];
	       else
               mesh->Volus[idx_volu]->Faces[1]=
               v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	       mesh->Faces[f3]->Old_Child[0]:mesh->Faces[f3]->Child[0];
               */
               
	       if(   (mesh->Faces[f3]->Old_Child[1]>=0)&&(mesh->Faces[f3]->Old_Child[1]!=l0)&&
	             (mesh->Faces[f3]->Old_Child[1]!=l1)&&(mesh->Faces[f3]->Old_Child[1]!=l3)   )
	       {
	          mesh->Volus[idx_volu]->Faces[1]=
                  v1==mesh->Faces[mesh->Faces[f3]->Old_Child[0]]->Verts[0]?
	          mesh->Faces[f3]->Old_Child[0]:mesh->Faces[f3]->Child[0];
	       }
	       else
	       {
		  mesh->Volus[idx_volu]->Faces[1]=
                  v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	          mesh->Faces[f3]->Child[1]:mesh->Faces[f3]->Child[0];
	       }
	       mesh->Volus[idx_volu]->Faces[2]=idx_face-2;
               mesh->Volus[idx_volu++]->Faces[3]=f0;

	       mesh->Volus[j]->Verts[0]=v0;
               mesh->Volus[j]->Verts[1]=v2;
               mesh->Volus[j]->Verts[2]=v3;
               mesh->Volus[j]->Verts[3]=mesh->Lines[l0]->Child[2];

	       tmp_line= v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];
               mesh->Volus[j]->Lines[0]=l1;
               mesh->Volus[j]->Lines[1]=l2;
               mesh->Volus[j]->Lines[2]=tmp_line;
               mesh->Volus[j]->Lines[3]=l5;
	       /*
	       if(mesh->Faces[f3]->Old_Child[0]<0)
               mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Child[2];
	       else
               mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Old_Child[1];
               */

	       if(   (mesh->Faces[f3]->Old_Child[1]>=0)&&(mesh->Faces[f3]->Old_Child[1]!=l0)&&
	             (mesh->Faces[f3]->Old_Child[1]!=l1)&&(mesh->Faces[f3]->Old_Child[1]!=l3)   )
	       {
                 mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Old_Child[1];
	       }
	       else
	       {
                 mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Child[2];
	       }
               mesh->Volus[j]->Lines[5]=idx_line-1;

	       ////tmp_face=volu->Faces[1];
	       tmp_face=f1;
               mesh->Volus[j]->Faces[0]=idx_face-1;
               
               mesh->Volus[j]->Faces[1]=f2;
	       ////mesh->Volus[j]->Faces[1]=volu->Faces[2];
	       /*
	       if( mesh->Faces[f3]->Old_Child[0]<0)
               mesh->Volus[j]->Faces[2]=
	       v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	       mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Child[1];
               else
               mesh->Volus[j]->Faces[2]=
	       v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	       mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Old_Child[0];
               */
	       if(   (mesh->Faces[f3]->Old_Child[1]>=0)&&(mesh->Faces[f3]->Old_Child[1]!=l0)&&
	             (mesh->Faces[f3]->Old_Child[1]!=l1)&&(mesh->Faces[f3]->Old_Child[1]!=l3)   )
	       {
	         mesh->Volus[j]->Faces[2]=
	         v0==mesh->Faces[mesh->Faces[f3]->Old_Child[0]]->Verts[0]?
	         mesh->Faces[f3]->Old_Child[0]:mesh->Faces[f3]->Child[0];
	       }
	       else
	       {
	         mesh->Volus[j]->Faces[2]=
	         v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	         mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Child[1];
	       }

	       mesh->Volus[j]->Faces[3]=tmp_face;

               mesh->Volus[j]->Type=1;

	       mesh->Volus[j]->Child[0]=j;
	       mesh->Volus[j]->Child[1]=idx_volu-1;


	       //ReConstruct(mesh->Volus[idx_volu-1]);
               //ReConstruct(mesh->Volus[idx_volu-2]);
           
	       if(MarkJudge(mesh,mesh->Volus[j])<0)
	       {
	        mesh->Volus[j]->Mark=-1;
	       }
               if(MarkJudge(mesh,mesh->Volus[idx_volu-1])>=0) 
	       {
                mesh->Volus[idx_volu-1]->Mark=1;
	       }
               else
	        mesh->Volus[idx_volu-1]->Mark=-1;
             }//4
      }//type0
      
      else if(volu->Type==1)
       {//type1
          if (  (mesh->Faces[volu->Faces[3]]->Child[0]<0&&mesh->Faces[volu->Faces[2]]->Child[0]<0)
	      ||(mesh->Faces[volu->Faces[3]]->Child[0]<0&&mesh->Faces[volu->Faces[2]]->Child[0]>=0
		  &&(mesh->Faces[volu->Faces[2]]->Child[2]==volu->Lines[0]||mesh->Faces[volu->Faces[2]]
	         ->Child[2]==volu->Lines[2]||mesh->Faces[volu->Faces[2]]->Child[2]==volu->Lines[4]))
	      ||(mesh->Faces[volu->Faces[2]]->Child[0]<0&&mesh->Faces[volu->Faces[3]]->Child[0]>=0
	         &&(mesh->Faces[volu->Faces[3]]->Child[2]==volu->Lines[0]||mesh->Faces[volu->Faces[3]]
		    ->Child[2]==volu->Lines[1]||mesh->Faces[volu->Faces[3]]->Child[2]==volu->Lines[3]) )
	      ||(mesh->Faces[volu->Faces[3]]->Child[0]>=0&&(mesh->Faces[volu->Faces[3]]->Child[2]==
		    volu->Lines[0]||mesh->Faces[volu->Faces[3]]
		    ->Child[2]==volu->Lines[1]||mesh->Faces[volu->Faces[3]]->Child[2]==volu->Lines[3]) 
		 &&mesh->Faces[volu->Faces[2]]->Child[0]>=0&&(mesh->Faces[volu->Faces[2]]->Child[2]==
		    volu->Lines[0]||mesh->Faces[volu->Faces[2]]
	         ->Child[2]==volu->Lines[2]||mesh->Faces[volu->Faces[2]]->Child[2]==volu->Lines[4])
		 )
	     )   
          {//1
	   //更新面上两条新边的信息 
	   //printf ( "1-11\n" );
	   mesh->Lines[idx_line]->Verts[1]= mesh->Lines[l0]->Child[2];
	   mesh->Lines[idx_line]->Num_Volus_Owned=-1;
	   mesh->Lines[idx_line]->Volus_Owned=NULL;
           mesh->Lines[idx_line]->Mark=-1;
	   mesh->Lines[idx_line]->ID_Boundary=mesh->Faces[f3]->ID_Boundary;
	   mesh->Lines[idx_line]->Child[0]=-1;
	   mesh->Lines[idx_line]->Child[1]=-1;
	   mesh->Lines[idx_line]->Child[2]=-1;
	   mesh->Lines[idx_line++]->Verts[0]=v2;
	   //printf ( "----------------------1\n" );
	   mesh->Lines[idx_line]->Verts[1]= mesh->Lines[l0]->Child[2];
	   mesh->Lines[idx_line]->Num_Volus_Owned=-1;
	   mesh->Lines[idx_line]->Volus_Owned=NULL;
           mesh->Lines[idx_line]->Mark=-1;
	   mesh->Lines[idx_line]->ID_Boundary=mesh->Faces[f2]->ID_Boundary;
	   mesh->Lines[idx_line]->Child[0]=-1;
	   mesh->Lines[idx_line]->Child[1]=-1;
	   mesh->Lines[idx_line]->Child[2]=-1;
           mesh->Lines[idx_line++]->Verts[0]=v3;
           //更新第一个面上的信息
           mesh->Faces[idx_face]->Verts[0]=v1;
           mesh->Faces[idx_face]->Verts[1]=v2;
           mesh->Faces[idx_face]->Verts[2]=mesh->Lines[l0]->Child[2];

	   mesh->Faces[idx_face]->Old_Child[0]=-1;
           mesh->Faces[idx_face]->Old_Child[1]=-1;

	   mesh->Faces[idx_face]->Child[0]=-1;
           mesh->Faces[idx_face]->Child[1]=-1;
           mesh->Faces[idx_face]->Child[2]=-1;
           mesh->Faces[idx_face]->ID_Boundary=mesh->Faces[f3]->ID_Boundary;

           mesh->Faces[idx_face]->Lines[0]=idx_line-2;
           mesh->Faces[idx_face]->Lines[1]=
	   volu->Verts[0]==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
           mesh->Faces[idx_face++]->Lines[2]=l3;
	   //printf ( "---------------------------------2\n" );
	   mesh->Faces[f3]->Verts[0]=v0;
           mesh->Faces[f3]->Verts[1]=mesh->Lines[l0]->Child[2];
           mesh->Faces[f3]->Verts[2]=v2;
        
           mesh->Faces[f3]->Old_Child[0]=mesh->Faces[f3]->Child[1];
	   mesh->Faces[f3]->Old_Child[1]=mesh->Faces[f3]->Child[2];
	   ////mesh->Faces[f3]->Child[0]=volu->Faces[3];
	   mesh->Faces[f3]->Child[0]=f3;
	   mesh->Faces[f3]->Child[1]=idx_face-1;
           mesh->Faces[f3]->Child[2]=idx_line-2;

	   mesh->Faces[f3]->Lines[0]=idx_line-2;
           mesh->Faces[f3]->Lines[1]=l1;
           mesh->Faces[f3]->Lines[2]=
           volu->Verts[0]==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];
	   //printf ( "-------------------------3\n" );
           //更新第二个面上的信息
           mesh->Faces[idx_face]->Verts[0]=v1;
           mesh->Faces[idx_face]->Verts[1]=v3;
           mesh->Faces[idx_face]->Verts[2]=mesh->Lines[l0]->Child[2];

           mesh->Faces[idx_face]->Old_Child[0]=-1;
           mesh->Faces[idx_face]->Old_Child[1]=-1;

	   mesh->Faces[idx_face]->Child[0]=-1;
           mesh->Faces[idx_face]->Child[1]=-1;
           mesh->Faces[idx_face]->Child[2]=-1;
           mesh->Faces[idx_face]->ID_Boundary=mesh->Faces[f2]->ID_Boundary;

           mesh->Faces[idx_face]->Lines[0]=idx_line-1;
           mesh->Faces[idx_face]->Lines[1]=
           volu->Verts[0]==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
           mesh->Faces[idx_face++]->Lines[2]=volu->Lines[4];
	   //printf ( "--------------------------4\n" );
	   mesh->Faces[f2]->Verts[0]=v0;
           mesh->Faces[f2]->Verts[1]=mesh->Lines[l0]->Child[2];
           mesh->Faces[f2]->Verts[2]=v3;
        
           mesh->Faces[f2]->Old_Child[0]= mesh->Faces[f2]->Child[1];
           mesh->Faces[f2]->Old_Child[1]= mesh->Faces[f2]->Child[2];
	   mesh->Faces[f2]->Child[0]=f2;
           mesh->Faces[f2]->Child[1]=idx_face-1;
           mesh->Faces[f2]->Child[2]=idx_line-1;

	   mesh->Faces[f2]->Lines[0]=idx_line-1;
           mesh->Faces[f2]->Lines[1]=l2;
           mesh->Faces[f2]->Lines[2]=
           v0==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];
	   //printf ( "-----------------------------5\n" );
           //更新中间的新面的信息
	   mesh->Faces[idx_face]->Verts[0]=mesh->Lines[l0]->Child[2];
           mesh->Faces[idx_face]->Verts[1]=v2;
           mesh->Faces[idx_face]->Verts[2]=v3;
	   
	   mesh->Faces[idx_face]->Old_Child[0]=-1;
           mesh->Faces[idx_face]->Old_Child[1]=-1;

	   mesh->Faces[idx_face]->Child[0]=-1;
           mesh->Faces[idx_face]->Child[1]=-1;
           mesh->Faces[idx_face]->Child[2]=-1;
         
	   mesh->Faces[idx_face]->ID_Boundary=0;

	   mesh->Faces[idx_face]->Lines[0]=l5;
           mesh->Faces[idx_face]->Lines[1]=idx_line-1;
           mesh->Faces[idx_face++]->Lines[2]=idx_line-2;

	   //printf ( "-----------------6\n" );
	      //更新单元信息
	   mesh->Volus[idx_volu]->Verts[0]=v1;
	   mesh->Volus[idx_volu]->Verts[1]=v2;
	   mesh->Volus[idx_volu]->Verts[2]=v3;
           mesh->Volus[idx_volu]->Verts[3]=mesh->Lines[l0]->Child[2];
           
	   mesh->Volus[idx_volu]->Lines[0]=l3;
	   mesh->Volus[idx_volu]->Lines[1]=l4;
	   mesh->Volus[idx_volu]->Lines[2]=
	   v0==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
 	   mesh->Volus[idx_volu]->Lines[3]=l5;
           mesh->Volus[idx_volu]->Lines[4]=idx_line-2;
 	   mesh->Volus[idx_volu]->Lines[5]=idx_line-1;

	   if(i==0)
   	   {
            mesh->Volus[idx_volu]->Father=j;
            mesh->Volus[j]->Father=j;
           }
           else
           mesh->Volus[idx_volu]->Father=mesh->Volus[j]->Father;
           mesh->Volus[idx_volu]->Ancestor=mesh->Volus[j]->Ancestor;


           mesh->Volus[idx_volu]->Type=2;

	   mesh->Volus[idx_volu]->Faces[0]=idx_face-1;
	   mesh->Volus[idx_volu]->Faces[1]=idx_face-2;
	   mesh->Volus[idx_volu]->Faces[2]=idx_face-3;
           ////mesh->Volus[idx_volu++]->Faces[3]=volu->Faces[0];
           mesh->Volus[idx_volu++]->Faces[3]=f0;
	   //printf ( "--------------------7\n" );
	   mesh->Volus[j]->Verts[0]=v0;
	   mesh->Volus[j]->Verts[1]=v2;
	   mesh->Volus[j]->Verts[2]=v3;
	   mesh->Volus[j]->Verts[3]=mesh->Lines[l0]->Child[2];

	   tmp_line=volu->Verts[0]==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];
           mesh->Volus[j]->Lines[0]=l1;
	   mesh->Volus[j]->Lines[1]=l2;
	   mesh->Volus[j]->Lines[2]=tmp_line;
 	   mesh->Volus[j]->Lines[3]=l5;
           mesh->Volus[j]->Lines[4]=idx_line-2;
 	   mesh->Volus[j]->Lines[5]=idx_line-1;

	   //tmp_face=volu->Faces[1];
	   mesh->Volus[j]->Faces[0]=idx_face-1;
	   mesh->Volus[j]->Faces[1]=f2;
	   mesh->Volus[j]->Faces[2]=f3;
	   mesh->Volus[j]->Faces[3]=f1;

           mesh->Volus[j]->Type=2;

	   mesh->Volus[j]->Child[0]=j;
	   mesh->Volus[j]->Child[1]=idx_volu-1;
	   //printf ( "--------------------------------8\n" );
	   /*
	   printf ( "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",l0,l1,l2,l3,l4,l5,v0,v1,v2,v3 );
	   printf ( "%d,%d,%d,%d,%d,%d\n",mesh->Lines[l1]->Mark, mesh->Lines[l1]->Verts[0],mesh->Lines[l1]->Verts[1],mesh->Lines[l1]->Child[0], mesh->Lines[l1]->Child[1], mesh->Lines[l1]->Child[2]);
	   printf ( "%d,%d,%d,%d,%d,%d\n",mesh->Lines[l0]->Mark ,mesh->Lines[l1]->Mark,mesh->Lines[l2]->Mark,mesh->Lines[l3]->Mark,mesh->Lines[l4]->Mark,mesh->Lines[l5]->Mark);
	   printf ( "%d,%d,%d,%d\n",mesh->Num_Lines_Global,l0,mesh->Lines[l0]->Child[0],mesh->Lines[l0]->Child[1] );
	   printf ( "lines===%d\n", mesh->Volus[j]->Lines[2]);
           printf ( "mark:%d,%d\n", mesh->Lines[mesh->Volus[j]->Lines[2]]->Mark,mesh->Lines[mesh->Volus[j]->Lines[2]]->Child[0]);
           */
	   //ReConstruct(mesh->Volus[idx_volu-1]);
           //ReConstruct(mesh->Volus[MarkSet[i]]);
           if(MarkJudge(mesh,mesh->Volus[j])<0)
	   {
            mesh->Volus[j]->Mark=-1;
	   }
	   //printf ( "-----------------------89\n" );
           if(MarkJudge(mesh,mesh->Volus[idx_volu-1])>=0) 
	   {
            mesh->Volus[idx_volu-1]->Mark=1;
	    //printf ( "-------aa\n" );
	   }
           else
	   mesh->Volus[idx_volu-1]->Mark=-1;
	   //printf ( "----------------------------9\n" );

         }//1
           //-----------------------------------------------------------------------------------
           //第二种情况：第一个面没加密，第二个面加密 
	    else if(  (mesh->Faces[volu->Faces[3]]->Child[0]<0
		     &&mesh->Faces[volu->Faces[2]]->Child[0]>=0
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[0]
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[2]
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[4])
		  ||(mesh->Faces[volu->Faces[3]]->Child[0]>=0
		     &&(mesh->Faces[volu->Faces[3]]->Child[2]==
		    volu->Lines[0]||mesh->Faces[volu->Faces[3]]
		    ->Child[2]==volu->Lines[1]||mesh->Faces[volu->Faces[3]]->Child[2]==volu->Lines[3]) 
		     &&mesh->Faces[volu->Faces[2]]->Child[0]>=0
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[0]
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[2]
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[4]
		    )
	       )
	   {//2
	       //更新第一个面上的边	   
	       //printf ( "1-22\n" );
               mesh->Lines[idx_line]->Verts[1]=mesh->Lines[l0]->Child[2];
	       mesh->Lines[idx_line]->Child[0]=-1;
	       mesh->Lines[idx_line]->Child[1]=-1;
               mesh->Lines[idx_line]->Child[2]=-1;
	       mesh->Lines[idx_line]->Num_Volus_Owned=-1;
	       mesh->Lines[idx_line]->Volus_Owned=NULL;
               mesh->Lines[idx_line]->Mark=-1;
	       mesh->Lines[idx_line]->ID_Boundary=mesh->Faces[volu->Faces[3]]->ID_Boundary;
	       mesh->Lines[idx_line++]->Verts[0]=v2;
	       	       
	       //更新第一个面上的新面
               mesh->Faces[f3]->Verts[0]=v0;
               mesh->Faces[f3]->Verts[1]=mesh->Lines[l0]->Child[2];
               mesh->Faces[f3]->Verts[2]=v2;
            
	       mesh->Faces[f3]->Old_Child[0]=mesh->Faces[f3]->Child[1];
               mesh->Faces[f3]->Old_Child[1]=mesh->Faces[f3]->Child[2];
  
               mesh->Faces[f3]->Child[0]=f3;
               mesh->Faces[f3]->Child[1]=idx_face;
               mesh->Faces[f3]->Child[2]=idx_line-1;

	       mesh->Faces[f3]->Lines[0]=idx_line-1;
               mesh->Faces[f3]->Lines[1]=l1;
               mesh->Faces[f3]->Lines[2]=
	       v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];

               mesh->Faces[idx_face]->Verts[0]=v1;
               mesh->Faces[idx_face]->Verts[1]=v2;
               mesh->Faces[idx_face]->Verts[2]=mesh->Lines[l0]->Child[2];

               mesh->Faces[idx_face]->Child[0]=-1;
               mesh->Faces[idx_face]->Child[1]=-1;
               mesh->Faces[idx_face]->Child[2]=-1;

	       mesh->Faces[idx_face]->Old_Child[0]=-1;
               mesh->Faces[idx_face]->Old_Child[1]=-1;

               mesh->Faces[idx_face]->ID_Boundary=mesh->Faces[f3]->ID_Boundary;

               mesh->Faces[idx_face]->Lines[0]=idx_line-1;
               mesh->Faces[idx_face]->Lines[1]=
	       v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
               mesh->Faces[idx_face++]->Lines[2]=l3;

	       //更新中间新面的信息
               mesh->Faces[idx_face]->Verts[0]=mesh->Lines[l0]->Child[2];
               mesh->Faces[idx_face]->Verts[1]=v2;
               mesh->Faces[idx_face]->Verts[2]=v3;
              
               mesh->Faces[idx_face]->Child[0]=-1;
               mesh->Faces[idx_face]->Child[1]=-1;
               mesh->Faces[idx_face]->Child[2]=-1;

	       mesh->Faces[idx_face]->Old_Child[0]=-1;
               mesh->Faces[idx_face]->Old_Child[1]=-1;

               mesh->Faces[idx_face]->ID_Boundary=0;

	       mesh->Faces[idx_face]->Lines[0]=l5;
               if(   (mesh->Faces[f2]->Old_Child[1]>=0)&&(mesh->Faces[f2]->Old_Child[1]!=l0)&&
	             (mesh->Faces[f2]->Old_Child[1]!=l2)&&(mesh->Faces[f2]->Old_Child[1]!=l4)   )
	       {
	         mesh->Faces[idx_face]->Lines[1]=mesh->Faces[f2]->Old_Child[1];
	       } 
	       else
	       {
	         mesh->Faces[idx_face]->Lines[1]=mesh->Faces[f2]->Child[2];
	       }
               mesh->Faces[idx_face++]->Lines[2]=idx_line-1;
          
	       //更新两个单元的信息
               mesh->Volus[idx_volu]->Verts[0]=v1;
               mesh->Volus[idx_volu]->Verts[1]=v2;
               mesh->Volus[idx_volu]->Verts[2]=v3;
               mesh->Volus[idx_volu]->Verts[3]=mesh->Lines[l0]->Child[2];

               mesh->Volus[idx_volu]->Lines[0]=l3;
               mesh->Volus[idx_volu]->Lines[1]=l4;
               mesh->Volus[idx_volu]->Lines[2]=
	       v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
               mesh->Volus[idx_volu]->Lines[3]=l5;
               /*     if(mesh->Faces[f2]->Old_Child[1]<0)
	       mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f2]->Child[2];
	       else
               mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f2]->Old_Child[1];
               */
	       if(   (mesh->Faces[f2]->Old_Child[1]>=0)&&(mesh->Faces[f2]->Old_Child[1]!=l0)&&
	             (mesh->Faces[f2]->Old_Child[1]!=l2)&&(mesh->Faces[f2]->Old_Child[1]!=l4)   )
	       {
	         mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f2]->Old_Child[1];
	       }
	       else
	       {
	         mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f2]->Child[2];
	       }
	       mesh->Volus[idx_volu]->Lines[4]=idx_line-1;

               if(i==0)
   	       {
                mesh->Volus[idx_volu]->Father=j;
                mesh->Volus[j]->Father=j;
               }
               else
               mesh->Volus[idx_volu]->Father=mesh->Volus[j]->Father;
               mesh->Volus[idx_volu]->Ancestor=mesh->Volus[j]->Ancestor;
               
               
               mesh->Volus[idx_volu]->Type=2;
               
               mesh->Volus[idx_volu]->Faces[0]=idx_face-1;
               mesh->Volus[idx_volu]->Faces[2]=idx_face-2;
	       /*
	       if( mesh->Faces[volu->Faces[2]]->Old_Child[1]<0)
               mesh->Volus[idx_volu]->Faces[1]=
               v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	       mesh->Faces[f2]->Child[1]:mesh->Faces[f2]->Child[0];
	       else
	       mesh->Volus[idx_volu]->Faces[1]=
               v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	       mesh->Faces[f2]->Old_Child[0]:mesh->Faces[f2]->Child[0];
               */
               if(   (mesh->Faces[f2]->Old_Child[1]>=0)&&(mesh->Faces[f2]->Old_Child[1]!=l0)&&
	             (mesh->Faces[f2]->Old_Child[1]!=l2)&&(mesh->Faces[f2]->Old_Child[1]!=l4)   )
	       {

		  mesh->Volus[idx_volu]->Faces[1]=
                  v1==mesh->Faces[mesh->Faces[f2]->Old_Child[0]]->Verts[0]?
	          mesh->Faces[f2]->Old_Child[0]:mesh->Faces[f2]->Child[0];
	       }
	       else
	       {
                  mesh->Volus[idx_volu]->Faces[1]=
                  v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	          mesh->Faces[f2]->Child[1]:mesh->Faces[f2]->Child[0];
	       }
	       mesh->Volus[idx_volu++]->Faces[3]=f0;

	       mesh->Volus[j]->Verts[0]=v0;
               mesh->Volus[j]->Verts[1]=v2;
               mesh->Volus[j]->Verts[2]=v3;
               mesh->Volus[j]->Verts[3]=mesh->Lines[l0]->Child[2];

	       tmp_line=v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];
               mesh->Volus[j]->Lines[0]=l1;
               mesh->Volus[j]->Lines[1]=l2;
               mesh->Volus[j]->Lines[2]=tmp_line;
               mesh->Volus[j]->Lines[3]=l5;
               mesh->Volus[j]->Lines[4]=idx_line-1;
	       /* if(mesh->Faces[f2]->Old_Child[1]<0)
               mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Child[2];
	       else
               mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Old_Child[1];
               */

	       if(   (mesh->Faces[f2]->Old_Child[1]>=0)&&(mesh->Faces[f2]->Old_Child[1]!=l0)&&
	             (mesh->Faces[f2]->Old_Child[1]!=l2)&&(mesh->Faces[f2]->Old_Child[1]!=l4)   )
	       {
                 mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Old_Child[1];
	       }
	       else
	       {
                 mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Child[2];
	       }

	       tmp_face=f1;
               mesh->Volus[j]->Faces[0]=idx_face-1;
	       /*
	       if(mesh->Faces[f2]->Old_Child[1]<0)
               mesh->Volus[j]->Faces[1]=
	       v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	       mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Child[1];
	       else
	       mesh->Volus[j]->Faces[1]=
	       v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	       mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Old_Child[0];
               */
	       if(   (mesh->Faces[f2]->Old_Child[1]>=0)&&(mesh->Faces[f2]->Old_Child[1]!=l0)&&
	             (mesh->Faces[f2]->Old_Child[1]!=l2)&&(mesh->Faces[f2]->Old_Child[1]!=l4)   )
	       {
	         mesh->Volus[j]->Faces[1]=
	         v0==mesh->Faces[mesh->Faces[f2]->Old_Child[0]]->Verts[0]?
	         mesh->Faces[f2]->Old_Child[0]:mesh->Faces[f2]->Child[0];
	       }
	       else
	       {
	         mesh->Volus[j]->Faces[1]=
	         v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	         mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Child[1];
	       }
	       mesh->Volus[j]->Faces[2]=f3;
               mesh->Volus[j]->Faces[3]=tmp_face;

               mesh->Volus[j]->Type=2;
               mesh->Volus[j]->Child[0]=j;
	       mesh->Volus[j]->Child[1]=idx_volu-1;

	       //ReConstruct(mesh->Volus[idx_volu-1]);
               //ReConstruct(mesh->Volus[idx_volu-2]);
              if(MarkJudge(mesh,mesh->Volus[j])<0)
	      {
	       mesh->Volus[j]->Mark=-1; 
	      }
              if(MarkJudge(mesh,mesh->Volus[idx_volu-1])>=0) 
	      {
               mesh->Volus[idx_volu-1]->Mark=1;
	      }
              else
	       mesh->Volus[idx_volu-1]->Mark=-1;
           }//2
           //---------------------------------------------------------------------------------------------
	   else if ( (mesh->Faces[volu->Faces[3]]->Child[0]>=0
		     &&mesh->Faces[volu->Faces[3]]->Child[2]!=volu->Lines[0]
		     &&mesh->Faces[volu->Faces[3]]->Child[2]!=volu->Lines[1]
		     &&mesh->Faces[volu->Faces[3]]->Child[2]!=volu->Lines[3]
		     )
		  &&(mesh->Faces[volu->Faces[2]]->Child[0]>=0
		    &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[0]
		    &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[2]
		    &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[4]
		    )
		  ) 
	   {//3
           //更新中间的新面的信息
           mesh->Faces[idx_face]->Verts[0]=mesh->Lines[l0]->Child[2];
           mesh->Faces[idx_face]->Verts[1]=v2;
           mesh->Faces[idx_face]->Verts[2]=v3;
         
	   mesh->Faces[idx_face]->Old_Child[0]=-1;
           mesh->Faces[idx_face]->Old_Child[1]=-1;

	   mesh->Faces[idx_face]->Child[0]=-1;
           mesh->Faces[idx_face]->Child[1]=-1;
           mesh->Faces[idx_face]->Child[2]=-1;
 
           mesh->Faces[idx_face]->ID_Boundary=0;

	   mesh->Faces[idx_face]->Lines[0]=l5;
           
      	   if(   (mesh->Faces[f2]->Old_Child[1]>=0)&&(mesh->Faces[f2]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f2]->Old_Child[1]!=l2)&&(mesh->Faces[f2]->Old_Child[1]!=l4)   )
	   {
	     mesh->Faces[idx_face]->Lines[1]=mesh->Faces[f2]->Old_Child[1];
	   }
	   else
	   {
	     mesh->Faces[idx_face]->Lines[1]=mesh->Faces[f2]->Child[2];
	   }
	   if(   (mesh->Faces[f3]->Old_Child[1]>=0)&&(mesh->Faces[f3]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f3]->Old_Child[1]!=l1)&&(mesh->Faces[f3]->Old_Child[1]!=l3)   )
	   {
             mesh->Faces[idx_face++]->Lines[2]=mesh->Faces[f3]->Old_Child[1];
	   }
	   else
	   {
             mesh->Faces[idx_face++]->Lines[2]=mesh->Faces[f3]->Child[2];
	   }
	   //mesh->Faces[idx_face]->Lines[1]=mesh->Faces[f2]->Child[2];
           //mesh->Faces[idx_face++]->Lines[2]=mesh->Faces[f3]->Child[2];

           //更新单元信息
	   mesh->Volus[idx_volu]->Verts[0]=v1;
	   mesh->Volus[idx_volu]->Verts[1]=v2;
	   mesh->Volus[idx_volu]->Verts[2]=v3;
           mesh->Volus[idx_volu]->Verts[3]=mesh->Lines[l0]->Child[2];

 	   mesh->Volus[idx_volu]->Lines[0]=l3;
           mesh->Volus[idx_volu]->Lines[1]=l4;
           mesh->Volus[idx_volu]->Lines[2]=
           v0==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
	   mesh->Volus[idx_volu]->Lines[3]=l5;
	   /*
	   if(mesh->Faces[f2]->Old_Child[1]<0)
	   mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f2]->Child[2];
	   else
	   mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f2]->Old_Child[1];
           if(mesh->Faces[f3]->Old_Child[1]<0)
	   mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f3]->Child[2];
           else
           mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f3]->Old_Child[1];
           */
	   if(   (mesh->Faces[f2]->Old_Child[1]>=0)&&(mesh->Faces[f2]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f2]->Old_Child[1]!=l2)&&(mesh->Faces[f2]->Old_Child[1]!=l4)   )
	   {
	     mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f2]->Old_Child[1];
	   }
	   else
	   {
	     mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f2]->Child[2];
	   }
	   if(   (mesh->Faces[f3]->Old_Child[1]>=0)&&(mesh->Faces[f3]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f3]->Old_Child[1]!=l1)&&(mesh->Faces[f3]->Old_Child[1]!=l3)   )
	   {
	     mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f3]->Old_Child[1];
	   }
	   else
	   {
	     mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f3]->Child[2];
	   }

	   if(i==0)
   	       {
                mesh->Volus[idx_volu]->Father=j;
                mesh->Volus[j]->Father=j;
               }
           else
           mesh->Volus[idx_volu]->Father=mesh->Volus[j]->Father;
           mesh->Volus[idx_volu]->Ancestor=mesh->Volus[j]->Ancestor;


           mesh->Volus[idx_volu]->Type=2;

	   mesh->Volus[idx_volu]->Faces[0]=idx_face-1;
	   /*
	   if(mesh->Faces[f3]->Old_Child[0]<0)
	   mesh->Volus[idx_volu]->Faces[2]=
           v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	   mesh->Faces[f3]->Child[1]:mesh->Faces[f3]->Child[0];
	   else
           mesh->Volus[idx_volu]->Faces[2]=
           v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	   mesh->Faces[f3]->Old_Child[0]:mesh->Faces[f3]->Child[0];
	   if(mesh->Faces[f2]->Old_Child[0]<0)
	   mesh->Volus[idx_volu]->Faces[1]=
           v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	   mesh->Faces[f2]->Child[1]:mesh->Faces[f2]->Child[0];
	   else
	   mesh->Volus[idx_volu]->Faces[1]=
           v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	   mesh->Faces[f2]->Old_Child[0]:mesh->Faces[f2]->Child[0];
           */
           if(   (mesh->Faces[f2]->Old_Child[1]>=0)&&(mesh->Faces[f2]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f2]->Old_Child[1]!=l2)&&(mesh->Faces[f2]->Old_Child[1]!=l4)   )
	   {
	     mesh->Volus[idx_volu]->Faces[1]=
             v1==mesh->Faces[mesh->Faces[f2]->Old_Child[0]]->Verts[0]?
	     mesh->Faces[f2]->Old_Child[0]:mesh->Faces[f2]->Child[0];
	   }
	   else
	   {
	     mesh->Volus[idx_volu]->Faces[1]=
             v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	     mesh->Faces[f2]->Child[1]:mesh->Faces[f2]->Child[0];
	   }
	   if(   (mesh->Faces[f3]->Old_Child[1]>=0)&&(mesh->Faces[f3]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f3]->Old_Child[1]!=l1)&&(mesh->Faces[f3]->Old_Child[1]!=l3)   )
	   {
	     mesh->Volus[idx_volu]->Faces[2]=
             v1==mesh->Faces[mesh->Faces[f3]->Old_Child[0]]->Verts[0]?
	     mesh->Faces[f3]->Old_Child[0]:mesh->Faces[f3]->Child[0];
	   }
	   else
	   {
	     mesh->Volus[idx_volu]->Faces[2]=
             v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	     mesh->Faces[f3]->Child[1]:mesh->Faces[f3]->Child[0];
	   }

	   mesh->Volus[idx_volu++]->Faces[3]=f0;

           mesh->Volus[j]->Verts[0]=v0;
	   mesh->Volus[j]->Verts[1]=v2;
	   mesh->Volus[j]->Verts[2]=v3;
           mesh->Volus[j]->Verts[3]=mesh->Lines[l0]->Child[2];

	   tmp_line=v0==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];
	   mesh->Volus[j]->Lines[0]=l1;
	   mesh->Volus[j]->Lines[1]=l2;
           mesh->Volus[j]->Lines[2]=tmp_line;
	   mesh->Volus[j]->Lines[3]=l5;
           /*if(mesh->Faces[f3]->Old_Child[1]<0)
	   mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Child[2];
	   else
   	   mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Old_Child[1];
	   if(mesh->Faces[f2]->Old_Child[1]<0)
	   else
           mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Old_Child[1];
           */

	   if(   (mesh->Faces[f2]->Old_Child[1]>=0)&&(mesh->Faces[f2]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f2]->Old_Child[1]!=l2)&&(mesh->Faces[f2]->Old_Child[1]!=l4)   )
	   {
	     mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Old_Child[1];
	   }
	   else
	   {
	     mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Child[2];
	   }
	   if(   (mesh->Faces[f3]->Old_Child[1]>=0)&&(mesh->Faces[f3]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f3]->Old_Child[1]!=l1)&&(mesh->Faces[f3]->Old_Child[1]!=l3)   )
	   {
	     mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Old_Child[1];
	   }
	   else
	   {
	     mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Child[2];
	   }

	   //mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Child[2];
           //mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Child[2];

	   tmp_face=f1;
	   mesh->Volus[j]->Faces[0]=idx_face-1;
	   /*
	   if(mesh->Faces[f2]->Old_Child[1]<0)
	   mesh->Volus[j]->Faces[1]=
           v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	   mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Child[1];
	   else
	   mesh->Volus[j]->Faces[1]=
           v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	   mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Old_Child[0];
           if(mesh->Faces[f3]->Old_Child[1]<0)
	   mesh->Volus[j]->Faces[2]=
           v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	   mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Child[1];
	   else
	   mesh->Volus[j]->Faces[2]=
           v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	   mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Old_Child[0];
           */
	   if(   (mesh->Faces[f2]->Old_Child[1]>=0)&&(mesh->Faces[f2]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f2]->Old_Child[1]!=l2)&&(mesh->Faces[f2]->Old_Child[1]!=l4)   )
	   {
	     mesh->Volus[j]->Faces[1]=
             v0==mesh->Faces[mesh->Faces[f2]->Old_Child[0]]->Verts[0]?
	     mesh->Faces[f2]->Old_Child[0]:mesh->Faces[f2]->Child[0];
	   }
	   else
	   {
	     mesh->Volus[j]->Faces[1]=
             v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	     mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Child[1];
	   }
	   if(   (mesh->Faces[f3]->Old_Child[1]>=0)&&(mesh->Faces[f3]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f3]->Old_Child[1]!=l1)&&(mesh->Faces[f3]->Old_Child[1]!=l3)   )
	   {
	     mesh->Volus[j]->Faces[2]=
             v0==mesh->Faces[mesh->Faces[f3]->Old_Child[0]]->Verts[0]?
	     mesh->Faces[f3]->Old_Child[0]:mesh->Faces[f3]->Child[0];
	   }
	   else
	   {
	     mesh->Volus[j]->Faces[2]=
             v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	     mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Child[1];
	   }
	   mesh->Volus[j]->Faces[3]=tmp_face;

           mesh->Volus[j]->Type=2;

	   mesh->Volus[j]->Child[0]=j;
	   mesh->Volus[j]->Child[1]=idx_volu-1;


           //ReConstruct(mesh->Volus[idx_volu-1]);
           //ReConstruct(mesh->Volus[idx_volu-2]);
 	   if(MarkJudge(mesh,mesh->Volus[j])<0)
	   {
	    mesh->Volus[j]->Mark=-1;
	   }
           if(MarkJudge(mesh,mesh->Volus[idx_volu-1])>=0) 
	   {
            mesh->Volus[idx_volu-1]->Mark=1;
	   }
          else
	   mesh->Volus[idx_volu-1]->Mark=-1;


	   }//3
	   //-------------------------------------------------------------------------------------------
           //第四种情况：第一个面(012)加密，第二个面没加密 
	   else
	   {//4
	      //更新第二个面上的边	   
	   //printf ( "1-44\n" );
               mesh->Lines[idx_line]->Verts[1]=mesh->Lines[l0]->Child[2];
	       mesh->Lines[idx_line]->Child[0]=-1;
	       mesh->Lines[idx_line]->Child[1]=-1;
               mesh->Lines[idx_line]->Child[2]=-1;
	       mesh->Lines[idx_line]->Num_Volus_Owned=-1;
	       mesh->Lines[idx_line]->Volus_Owned=NULL;
               mesh->Lines[idx_line]->Mark=-1;
	       mesh->Lines[idx_line]->ID_Boundary=mesh->Faces[volu->Faces[2]]->ID_Boundary;
	       mesh->Lines[idx_line++]->Verts[0]=v3;
	       	       
	       //更新第二个面上的新面
               mesh->Faces[volu->Faces[2]]->Verts[0]=v0;
               mesh->Faces[volu->Faces[2]]->Verts[1]=mesh->Lines[l0]->Child[2];
               mesh->Faces[volu->Faces[2]]->Verts[2]=v3;
              
	       mesh->Faces[volu->Faces[2]]->Lines[0]=idx_line-1;
               mesh->Faces[volu->Faces[2]]->Lines[1]=l2;
               mesh->Faces[volu->Faces[2]]->Lines[2]=
	       volu->Verts[0]==mesh->Lines[volu->Lines[0]]->Verts[0]?
	       mesh->Lines[volu->Lines[0]]->Child[0]:mesh->Lines[l0]->Child[1];

               mesh->Faces[idx_face]->Verts[0]=v1;
               mesh->Faces[idx_face]->Verts[1]=v3;
               mesh->Faces[idx_face]->Verts[2]=mesh->Lines[l0]->Child[2];

	       mesh->Faces[idx_face]->Old_Child[0]=-1;
               mesh->Faces[idx_face]->Old_Child[1]=-1;

               mesh->Faces[idx_face]->Child[0]=-1;
               mesh->Faces[idx_face]->Child[1]=-1;
               mesh->Faces[idx_face]->Child[2]=-1;

               mesh->Faces[idx_face]->ID_Boundary=mesh->Faces[volu->Faces[2]]->ID_Boundary;

               mesh->Faces[idx_face]->Lines[0]=idx_line-1;
               mesh->Faces[idx_face]->Lines[1]=
	       v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
               mesh->Faces[idx_face++]->Lines[2]=l4;

	       mesh->Faces[f2]->Old_Child[0]=mesh->Faces[f2]->Child[1];
	       mesh->Faces[f2]->Old_Child[1]=mesh->Faces[f2]->Child[2];

	       mesh->Faces[f2]->Child[0]=f2;
	       mesh->Faces[f2]->Child[1]=idx_face-1;
	       mesh->Faces[f2]->Child[2]=idx_line-1;
	      
	       //更新中间新面的信息
               mesh->Faces[idx_face]->Verts[0]=mesh->Lines[l0]->Child[2];
               mesh->Faces[idx_face]->Verts[1]=v2;
               mesh->Faces[idx_face]->Verts[2]=v3;
              
	       mesh->Faces[idx_face]->Old_Child[0]=-1;
               mesh->Faces[idx_face]->Old_Child[1]=-1;

               mesh->Faces[idx_face]->Child[0]=-1;
               mesh->Faces[idx_face]->Child[1]=-1;
               mesh->Faces[idx_face]->Child[2]=-1;

               mesh->Faces[idx_face]->ID_Boundary=0;

	       mesh->Faces[idx_face]->Lines[0]=l5;
               mesh->Faces[idx_face]->Lines[1]=idx_line-1;
               if(   (mesh->Faces[f3]->Old_Child[1]>=0)&&(mesh->Faces[f3]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f3]->Old_Child[1]!=l1)&&(mesh->Faces[f3]->Old_Child[1]!=l3)   )
	       {
	         mesh->Faces[idx_face++]->Lines[2]=mesh->Faces[f3]->Old_Child[1];
	       } 
	       else
	       {
	         mesh->Faces[idx_face++]->Lines[2]=mesh->Faces[f3]->Child[2];
	       }

	       //更新两个单元的信息
               mesh->Volus[idx_volu]->Verts[0]=v1;
               mesh->Volus[idx_volu]->Verts[1]=v2;
               mesh->Volus[idx_volu]->Verts[2]=v3;
               mesh->Volus[idx_volu]->Verts[3]=mesh->Lines[l0]->Child[2];

               mesh->Volus[idx_volu]->Lines[0]=l3;
               mesh->Volus[idx_volu]->Lines[1]=l4;
               mesh->Volus[idx_volu]->Lines[2]=
	       v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
               mesh->Volus[idx_volu]->Lines[3]=l5;
               mesh->Volus[idx_volu]->Lines[5]=idx_line-1;
	       /*  if(mesh->Faces[volu->Faces[3]]->Old_Child[1]<0)
               mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f3]->Child[2];
               else
               mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f3]->Old_Child[1];
               */

    	       if(   (mesh->Faces[f3]->Old_Child[1]>=0)&&(mesh->Faces[f3]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f3]->Old_Child[1]!=l1)&&(mesh->Faces[f3]->Old_Child[1]!=l3)   )
	       {
                 mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f3]->Old_Child[1];
	       }
	       else
	       {
                 mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f3]->Child[2];
	       }
	       if(i==0)
   	       {
                mesh->Volus[idx_volu]->Father=j;
                mesh->Volus[j]->Father=j;
               }
               else
               mesh->Volus[idx_volu]->Father=mesh->Volus[j]->Father;
               mesh->Volus[idx_volu]->Ancestor=mesh->Volus[j]->Ancestor;

               mesh->Volus[idx_volu]->Type=2;

               mesh->Volus[idx_volu]->Faces[0]=idx_face-1;
	       /*
	       if(mesh->Faces[f3]->Old_Child[0]<0)
               mesh->Volus[idx_volu]->Faces[2]=
               v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	       mesh->Faces[f3]->Child[1]:mesh->Faces[f3]->Child[0];
	       else
	       mesh->Volus[idx_volu]->Faces[2]=
               v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	       mesh->Faces[f3]->Old_Child[0]:mesh->Faces[f3]->Child[0];
               */
    	       if(   (mesh->Faces[f3]->Old_Child[1]>=0)&&(mesh->Faces[f3]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f3]->Old_Child[1]!=l1)&&(mesh->Faces[f3]->Old_Child[1]!=l3)   )
	       {
		  mesh->Volus[idx_volu]->Faces[2]=
                  v1==mesh->Faces[mesh->Faces[f3]->Old_Child[0]]->Verts[0]?
	          mesh->Faces[f3]->Old_Child[0]:mesh->Faces[f3]->Child[0];
	       }
	       else
	       {
		  mesh->Volus[idx_volu]->Faces[2]=
                  v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	          mesh->Faces[f3]->Child[1]:mesh->Faces[f3]->Child[0];
	       }		  
               mesh->Volus[idx_volu]->Faces[1]=idx_face-2;
               ////mesh->Volus[idx_volu++]->Faces[3]=volu->Faces[0];
               mesh->Volus[idx_volu++]->Faces[3]=f0;
	       
	       mesh->Volus[j]->Verts[0]=v0;
               mesh->Volus[j]->Verts[1]=v2;
               mesh->Volus[j]->Verts[2]=v3;
               mesh->Volus[j]->Verts[3]=mesh->Lines[l0]->Child[2];

	       tmp_line=v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];
               mesh->Volus[j]->Lines[0]=l1;
               mesh->Volus[j]->Lines[1]=l2;
               mesh->Volus[j]->Lines[2]=tmp_line;
               mesh->Volus[j]->Lines[3]=l5;
	       /*if(mesh->Faces[f3]->Old_Child[1]<0)
               mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Child[2];
               else
               mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Old_Child[1];
               */
    	       if(   (mesh->Faces[f3]->Old_Child[1]>=0)&&(mesh->Faces[f3]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f3]->Old_Child[1]!=l1)&&(mesh->Faces[f3]->Old_Child[1]!=l3)   )
	       {
                 mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Old_Child[1];
	       }
	       else
	       {
                 mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Child[2];
	       }		  
               mesh->Volus[j]->Lines[5]=idx_line-1;

	       ////tmp_face=volu->Faces[1];
	       tmp_face=f1;
               mesh->Volus[j]->Faces[0]=idx_face-1;
               mesh->Volus[j]->Faces[1]=f2;
	       /*
	       if(mesh->Faces[f3]->Old_Child[1]<0)
               mesh->Volus[j]->Faces[2]=
	       v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	       mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Child[1];
	       else
               mesh->Volus[j]->Faces[2]=
	       v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	       mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Old_Child[0];
               */
     	       if(   (mesh->Faces[f3]->Old_Child[1]>=0)&&(mesh->Faces[f3]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f3]->Old_Child[1]!=l1)&&(mesh->Faces[f3]->Old_Child[1]!=l3)   )
	       {
		  mesh->Volus[j]->Faces[2]=
	          v0==mesh->Faces[mesh->Faces[f3]->Old_Child[0]]->Verts[0]?
	          mesh->Faces[f3]->Old_Child[0]:mesh->Faces[f3]->Child[0];
	       }
	       else
	       {
		  mesh->Volus[j]->Faces[2]=
	          v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	          mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Child[1];
	       }		  
	       mesh->Volus[j]->Faces[3]=tmp_face;

               mesh->Volus[j]->Type=2;

	       mesh->Volus[j]->Child[0]=j;
	       mesh->Volus[j]->Child[1]=idx_volu-1;


	       //ReConstruct(mesh->Volus[idx_volu-1]);
               //ReConstruct(mesh->Volus[idx_volu-2]);
           
	       if(MarkJudge(mesh,mesh->Volus[j])<0)
	       {
	        mesh->Volus[j]->Mark=-1;
	       }
               if(MarkJudge(mesh,mesh->Volus[idx_volu-1])>=0) 
	       {
                mesh->Volus[idx_volu-1]->Mark=1;
	       }
               else
	        mesh->Volus[idx_volu-1]->Mark=-1;

             }//4

       }//type1
       else 
       {//type2
          if (  (mesh->Faces[volu->Faces[3]]->Child[0]<0&&mesh->Faces[volu->Faces[2]]->Child[0]<0)
	      ||(mesh->Faces[volu->Faces[3]]->Child[0]<0&&mesh->Faces[volu->Faces[2]]->Child[0]>=0
		  &&(mesh->Faces[volu->Faces[2]]->Child[2]==volu->Lines[0]||mesh->Faces[volu->Faces[2]]
	         ->Child[2]==volu->Lines[2]||mesh->Faces[volu->Faces[2]]->Child[2]==volu->Lines[4]))
	      ||(mesh->Faces[volu->Faces[2]]->Child[0]<0&&mesh->Faces[volu->Faces[3]]->Child[0]>=0
	         &&(mesh->Faces[volu->Faces[3]]->Child[2]==volu->Lines[0]||mesh->Faces[volu->Faces[3]]
		    ->Child[2]==volu->Lines[1]||mesh->Faces[volu->Faces[3]]->Child[2]==volu->Lines[3]) )
	      ||(mesh->Faces[volu->Faces[3]]->Child[0]>=0&&(mesh->Faces[volu->Faces[3]]->Child[2]==
		    volu->Lines[0]||mesh->Faces[volu->Faces[3]]
		    ->Child[2]==volu->Lines[1]||mesh->Faces[volu->Faces[3]]->Child[2]==volu->Lines[3]) 
		 &&mesh->Faces[volu->Faces[2]]->Child[0]>=0&&(mesh->Faces[volu->Faces[2]]->Child[2]==
		    volu->Lines[0]||mesh->Faces[volu->Faces[2]]
	         ->Child[2]==volu->Lines[2]||mesh->Faces[volu->Faces[2]]->Child[2]==volu->Lines[4])
		 )
	     )   
          {//1
	   //更新面上两条新边的信息
	   mesh->Lines[idx_line]->Verts[1]= mesh->Lines[l0]->Child[2];
           mesh->Lines[idx_line]->Mark=-1;
	   mesh->Lines[idx_line]->Num_Volus_Owned=-1;
	   mesh->Lines[idx_line]->Volus_Owned=NULL;
           mesh->Lines[idx_line]->Child[0]=-1;
	   mesh->Lines[idx_line]->Child[1]=-1;
	   mesh->Lines[idx_line]->Child[2]=-1;
	   //printf ( "ccc\n" );
	   mesh->Lines[idx_line]->ID_Boundary=mesh->Faces[f3]->ID_Boundary;
	   mesh->Lines[idx_line++]->Verts[0]=v2;
	   //printf ( "xxxxx\n" );
	   mesh->Lines[idx_line]->Verts[1]= mesh->Lines[l0]->Child[2];
           mesh->Lines[idx_line]->Mark=-1;
	   mesh->Lines[idx_line]->Num_Volus_Owned=-1;
	   mesh->Lines[idx_line]->Volus_Owned=NULL;
	   mesh->Lines[idx_line]->Child[0]=-1;
	   mesh->Lines[idx_line]->Child[1]=-1;
	   mesh->Lines[idx_line]->Child[2]=-1;
           mesh->Lines[idx_line]->ID_Boundary=mesh->Faces[f2]->ID_Boundary;
	   mesh->Lines[idx_line++]->Verts[0]=v3;
	   //printf ( "aaa\n" );
	   //更新第一个面上的信息
           mesh->Faces[idx_face]->Verts[0]=v1;
           mesh->Faces[idx_face]->Verts[1]=v2;
           mesh->Faces[idx_face]->Verts[2]=mesh->Lines[l0]->Child[2];

	   mesh->Faces[idx_face]->Old_Child[0]=-1;
           mesh->Faces[idx_face]->Old_Child[1]=-1;

	   mesh->Faces[idx_face]->Child[0]=-1;
           mesh->Faces[idx_face]->Child[1]=-1;
           mesh->Faces[idx_face]->Child[2]=-1;
           mesh->Faces[idx_face]->ID_Boundary=mesh->Faces[f3]->ID_Boundary;

           mesh->Faces[idx_face]->Lines[0]=idx_line-2;
           mesh->Faces[idx_face]->Lines[1]=
	   v0==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
           mesh->Faces[idx_face++]->Lines[2]=l3;

	   mesh->Faces[f3]->Verts[0]=v0;
           mesh->Faces[f3]->Verts[1]=mesh->Lines[l0]->Child[2];
           mesh->Faces[f3]->Verts[2]=v2;
         
	   mesh->Faces[f3]->Lines[0]=idx_line-2;
           mesh->Faces[f3]->Lines[1]=l1;
           mesh->Faces[f3]->Lines[2]=
           v0==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];
           
	   mesh->Faces[f3]->Old_Child[0]=mesh->Faces[f3]->Child[1];
	   mesh->Faces[f3]->Old_Child[1]=mesh->Faces[f3]->Child[2];
	   mesh->Faces[f3]->Child[0]=f3;
           mesh->Faces[f3]->Child[1]=idx_face-1;
           mesh->Faces[f3]->Child[2]=idx_line-2;
	   //printf ( "bbb\n" );
           //更新第二个面上的信息
           mesh->Faces[idx_face]->Verts[0]=v1;
           mesh->Faces[idx_face]->Verts[1]=v3;
           mesh->Faces[idx_face]->Verts[2]=mesh->Lines[l0]->Child[2];

           mesh->Faces[idx_face]->Old_Child[0]=-1;
           mesh->Faces[idx_face]->Old_Child[1]=-1;

	   mesh->Faces[idx_face]->Child[0]=-1;
           mesh->Faces[idx_face]->Child[1]=-1;
           mesh->Faces[idx_face]->Child[2]=-1;
           mesh->Faces[idx_face]->ID_Boundary=mesh->Faces[f2]->ID_Boundary;

           mesh->Faces[idx_face]->Lines[0]=idx_line-1;
           mesh->Faces[idx_face]->Lines[1]=
           v0==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
           mesh->Faces[idx_face++]->Lines[2]=l4;

	   mesh->Faces[f2]->Verts[0]=v0;
           mesh->Faces[f2]->Verts[1]=mesh->Lines[l0]->Child[2];
           mesh->Faces[f2]->Verts[2]=v3;
         

	   mesh->Faces[f2]->Lines[0]=idx_line-1;
           mesh->Faces[f2]->Lines[1]=l2;
           mesh->Faces[f2]->Lines[2]=
           v0==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];

	   mesh->Faces[f2]->Old_Child[0]=mesh->Faces[f2]->Child[1];
           mesh->Faces[f2]->Old_Child[1]=mesh->Faces[f2]->Child[2];

	   mesh->Faces[f2]->Child[0]=f2;
           mesh->Faces[f2]->Child[1]=idx_face-1;
           mesh->Faces[f2]->Child[2]=idx_line-1;

	   //printf ( "ccc\n" );
           //更新中间的新面的信息
	   mesh->Faces[idx_face]->Verts[0]=mesh->Lines[l0]->Child[2];
           mesh->Faces[idx_face]->Verts[1]=v2;
           mesh->Faces[idx_face]->Verts[2]=v3;
         
	   mesh->Faces[idx_face]->Old_Child[0]=-1;
           mesh->Faces[idx_face]->Old_Child[1]=-1;

	   mesh->Faces[idx_face]->Child[0]=-1;
           mesh->Faces[idx_face]->Child[1]=-1;
           mesh->Faces[idx_face]->Child[2]=-1;
         
	   mesh->Faces[idx_face]->ID_Boundary=0;

	   mesh->Faces[idx_face]->Lines[0]=l5;
           mesh->Faces[idx_face]->Lines[1]=idx_line-1;
           mesh->Faces[idx_face++]->Lines[2]=idx_line-2;

	   //printf ( "dddd\n" );
	      //更新单元信息
	   mesh->Volus[idx_volu]->Verts[0]=v1;
	   mesh->Volus[idx_volu]->Verts[1]=v2;
	   mesh->Volus[idx_volu]->Verts[2]=v3;
           mesh->Volus[idx_volu]->Verts[3]=mesh->Lines[l0]->Child[2];
           
	   mesh->Volus[idx_volu]->Lines[0]=l3;
	   mesh->Volus[idx_volu]->Lines[1]=l4;
	   mesh->Volus[idx_volu]->Lines[2]=
	   v0==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
 	   mesh->Volus[idx_volu]->Lines[3]=l5;
           mesh->Volus[idx_volu]->Lines[4]=idx_line-2;
 	   mesh->Volus[idx_volu]->Lines[5]=idx_line-1;

	   if(i==0)
   	       {
                mesh->Volus[idx_volu]->Father=j;
                mesh->Volus[j]->Father=j;
               }
           else
           mesh->Volus[idx_volu]->Father=mesh->Volus[j]->Father;
           mesh->Volus[idx_volu]->Ancestor=mesh->Volus[j]->Ancestor;


           mesh->Volus[idx_volu]->Type=0;

	   mesh->Volus[idx_volu]->Faces[0]=idx_face-1;
	   mesh->Volus[idx_volu]->Faces[1]=idx_face-2;
	   mesh->Volus[idx_volu]->Faces[2]=idx_face-3;
           mesh->Volus[idx_volu++]->Faces[3]=f0;
	   ////mesh->Volus[idx_volu++]->Faces[3]=volu->Faces[0];
	   //printf ( "eee\n" );
	   mesh->Volus[j]->Verts[0]=v0;
	   mesh->Volus[j]->Verts[1]=v2;
	   mesh->Volus[j]->Verts[2]=v3;
	   mesh->Volus[j]->Verts[3]=mesh->Lines[l0]->Child[2];

	   tmp_line=v0==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];
           mesh->Volus[j]->Lines[0]=l1;
	   mesh->Volus[j]->Lines[1]=l2;
	   mesh->Volus[j]->Lines[2]=tmp_line;
 	   mesh->Volus[j]->Lines[3]=l5;
           mesh->Volus[j]->Lines[4]=idx_line-2;
 	   mesh->Volus[j]->Lines[5]=idx_line-1;

	   tmp_face=f1;
	   mesh->Volus[j]->Faces[0]=idx_face-1;
	   mesh->Volus[j]->Faces[1]=f2;
	   mesh->Volus[j]->Faces[2]=f3;
	   mesh->Volus[j]->Faces[3]=tmp_face;

           mesh->Volus[j]->Type=0;

	   mesh->Volus[j]->Child[0]=j;
	   mesh->Volus[j]->Child[1]=idx_volu-1;

           //ReConstruct(mesh->Volus[idx_volu-1]);
           //ReConstruct(mesh->Volus[MarkSet[i]]);
           if(MarkJudge(mesh,mesh->Volus[j])<0)
	   {
            mesh->Volus[j]->Mark=-1;
	   }
           if(MarkJudge(mesh,mesh->Volus[idx_volu-1])>=0) 
	   {
            mesh->Volus[idx_volu-1]->Mark=1;
	   }
           else
	    mesh->Volus[idx_volu-1]->Mark=-1;
      }//1
           //-----------------------------------------------------------------------------------
           //第二种情况：第一个面没加密，第二个面加密 
	 else if(  (mesh->Faces[volu->Faces[3]]->Child[0]<0
		     &&mesh->Faces[volu->Faces[2]]->Child[0]>=0
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[0]
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[2]
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[4])
		  ||(mesh->Faces[volu->Faces[3]]->Child[0]>=0
		     &&(mesh->Faces[volu->Faces[3]]->Child[2]==
		    volu->Lines[0]||mesh->Faces[volu->Faces[3]]
		    ->Child[2]==volu->Lines[1]||mesh->Faces[volu->Faces[3]]->Child[2]==volu->Lines[3]) 
		     &&mesh->Faces[volu->Faces[2]]->Child[0]>=0
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[0]
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[2]
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[4]
		    )
	       )

	   {//2
	         //更新第一个面上的边	   
	       mesh->Lines[idx_line]->Verts[1]=mesh->Lines[l0]->Child[2];
	       mesh->Lines[idx_line]->Child[0]=-1;
	       mesh->Lines[idx_line]->Child[1]=-1;
               mesh->Lines[idx_line]->Child[2]=-1;
	       mesh->Lines[idx_line]->Num_Volus_Owned=-1;
	       mesh->Lines[idx_line]->Volus_Owned=NULL;
               mesh->Lines[idx_line]->Mark=-1;
	       mesh->Lines[idx_line]->ID_Boundary=mesh->Faces[f3]->ID_Boundary;
	       mesh->Lines[idx_line++]->Verts[0]=v2;
	       	       
	       //更新第一个面上的新面
               mesh->Faces[f3]->Verts[0]=v0;
               mesh->Faces[f3]->Verts[1]=mesh->Lines[l0]->Child[2];
               mesh->Faces[f3]->Verts[2]=v2;
              
	       mesh->Faces[f3]->Lines[0]=idx_line-1;
               mesh->Faces[f3]->Lines[1]=l1;
               mesh->Faces[f3]->Lines[2]=
	       volu->Verts[0]==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];

               mesh->Faces[idx_face]->Verts[0]=v1;
               mesh->Faces[idx_face]->Verts[1]=v2;
               mesh->Faces[idx_face]->Verts[2]=mesh->Lines[l0]->Child[2];

	       mesh->Faces[idx_face]->Old_Child[0]=-1;
               mesh->Faces[idx_face]->Old_Child[1]=-1;
 
               mesh->Faces[idx_face]->Child[0]=-1;
               mesh->Faces[idx_face]->Child[1]=-1;
               mesh->Faces[idx_face]->Child[2]=-1;

               mesh->Faces[idx_face]->ID_Boundary=mesh->Faces[f3]->ID_Boundary;

               mesh->Faces[idx_face]->Lines[0]=idx_line-1;
               mesh->Faces[idx_face]->Lines[1]=
	       v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
               mesh->Faces[idx_face++]->Lines[2]=l3;

	       mesh->Faces[f3]->Old_Child[0]=mesh->Faces[f3]->Child[1];
	       mesh->Faces[f3]->Old_Child[1]=mesh->Faces[f3]->Child[2];
	       //mesh->Faces[f3]->Child[0]=volu->Faces[3];
	       mesh->Faces[f3]->Child[0]=f3;
	       mesh->Faces[f3]->Child[1]=idx_face-1;
	       mesh->Faces[f3]->Child[2]=idx_line-1;
	      
	       //更新中间新面的信息
               mesh->Faces[idx_face]->Verts[0]=mesh->Lines[l0]->Child[2];
               mesh->Faces[idx_face]->Verts[1]=v2;
               mesh->Faces[idx_face]->Verts[2]=v3;
             
	       mesh->Faces[idx_face]->Old_Child[0]=-1;
               mesh->Faces[idx_face]->Old_Child[1]=-1;
 
               mesh->Faces[idx_face]->Child[0]=-1;
               mesh->Faces[idx_face]->Child[1]=-1;
               mesh->Faces[idx_face]->Child[2]=-1;

               mesh->Faces[idx_face]->ID_Boundary=0;

	       mesh->Faces[idx_face]->Lines[0]=l5;
               if(   (mesh->Faces[f2]->Old_Child[1]>=0)&&(mesh->Faces[f2]->Old_Child[1]!=l0)&&
	             (mesh->Faces[f2]->Old_Child[1]!=l2)&&(mesh->Faces[f2]->Old_Child[1]!=l4)   )
	       {
	         mesh->Faces[idx_face]->Lines[1]=mesh->Faces[f2]->Old_Child[1];
	       }
	       else
	       {
	         mesh->Faces[idx_face]->Lines[1]=mesh->Faces[f2]->Child[2];
	       }		  
               mesh->Faces[idx_face++]->Lines[2]=idx_line-1;
          
	       //更新两个单元的信息
               mesh->Volus[idx_volu]->Verts[0]=v1;
               mesh->Volus[idx_volu]->Verts[1]=v2;
               mesh->Volus[idx_volu]->Verts[2]=v3;
               mesh->Volus[idx_volu]->Verts[3]=mesh->Lines[l0]->Child[2];

               mesh->Volus[idx_volu]->Lines[0]=l3;
               mesh->Volus[idx_volu]->Lines[1]=l4;
               mesh->Volus[idx_volu]->Lines[2]=
	       v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
               mesh->Volus[idx_volu]->Lines[3]=l5;
	       /*if(mesh->Faces[f2]->Old_Child[1]<0)
               mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f2]->Child[2];
	       else
               mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f2]->Old_Child[1];
	       */
	       mesh->Volus[idx_volu]->Lines[4]=idx_line-1;
               if(   (mesh->Faces[f2]->Old_Child[1]>=0)&&(mesh->Faces[f2]->Old_Child[1]!=l0)&&
	             (mesh->Faces[f2]->Old_Child[1]!=l2)&&(mesh->Faces[f2]->Old_Child[1]!=l4)   )
		{
	          mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f2]->Old_Child[1];
		}
	       else
		{
	          mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f2]->Child[2];
		}
	       if(i==0)
   	       {
                mesh->Volus[idx_volu]->Father=j;
                mesh->Volus[j]->Father=j;
               }
               else
               mesh->Volus[idx_volu]->Father=mesh->Volus[j]->Father;
               mesh->Volus[idx_volu]->Ancestor=mesh->Volus[j]->Ancestor;

               mesh->Volus[idx_volu]->Type=0;

               mesh->Volus[idx_volu]->Faces[0]=idx_face-1;
               mesh->Volus[idx_volu]->Faces[2]=idx_face-2;
	       /*
	       if(mesh->Faces[volu->Faces[2]]->Old_Child[0]<0)
               mesh->Volus[idx_volu]->Faces[1]=
               v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	       mesh->Faces[f2]->Child[1]:mesh->Faces[f2]->Child[0];
	       else
	        mesh->Volus[idx_volu]->Faces[1]=
               v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	       mesh->Faces[f2]->Old_Child[0]:mesh->Faces[f2]->Child[0];
               */
               if(   (mesh->Faces[f2]->Old_Child[1]>=0)&&(mesh->Faces[f2]->Old_Child[1]!=l0)&&
	             (mesh->Faces[f2]->Old_Child[1]!=l2)&&(mesh->Faces[f2]->Old_Child[1]!=l4)   )
	       {
		  mesh->Volus[idx_volu]->Faces[1]=
                  v1==mesh->Faces[mesh->Faces[f2]->Old_Child[0]]->Verts[0]?
	          mesh->Faces[f2]->Old_Child[0]:mesh->Faces[f2]->Child[0];
	       }
	       else
	       {
		  mesh->Volus[idx_volu]->Faces[1]=
                  v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	          mesh->Faces[f2]->Child[1]:mesh->Faces[f2]->Child[0];
	       }		  

	       mesh->Volus[idx_volu++]->Faces[3]=f0;
	       //mesh->Volus[idx_volu++]->Faces[3]=volu->Faces[0];

	       mesh->Volus[j]->Verts[0]=v0;
               mesh->Volus[j]->Verts[1]=v2;
               mesh->Volus[j]->Verts[2]=v3;
               mesh->Volus[j]->Verts[3]=mesh->Lines[l0]->Child[2];

	       tmp_line=v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];
               mesh->Volus[j]->Lines[0]=l1;
               mesh->Volus[j]->Lines[1]=l2;
               mesh->Volus[j]->Lines[2]=tmp_line;
               mesh->Volus[j]->Lines[3]=l5;
               mesh->Volus[j]->Lines[4]=idx_line-1;
	       /*if(mesh->Faces[f2]->Old_Child[1]<0)
               mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Child[2];
	       else
                mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Old_Child[1];
               */
               if(   (mesh->Faces[f2]->Old_Child[1]>=0)&&(mesh->Faces[f2]->Old_Child[1]!=l0)&&
	             (mesh->Faces[f2]->Old_Child[1]!=l2)&&(mesh->Faces[f2]->Old_Child[1]!=l4)   )
	       {
	         mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Old_Child[1];
	       }
	       else
	       {
	         mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Child[2];
	       }		  
	       ////tmp_face=volu->Faces[1];
	       tmp_face=f1;
	       mesh->Volus[j]->Faces[0]=idx_face-1;
	       /*
	       if(mesh->Faces[f2]->Old_Child[1]<0)
               mesh->Volus[j]->Faces[1]=
	       v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	       mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Child[1];
	       else
	       mesh->Volus[j]->Faces[1]=
	       v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	       mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Old_Child[0];
               */
               if(   (mesh->Faces[f2]->Old_Child[1]>=0)&&(mesh->Faces[f2]->Old_Child[1]!=l0)&&
	             (mesh->Faces[f2]->Old_Child[1]!=l2)&&(mesh->Faces[f2]->Old_Child[1]!=l4)   )
	       {
		  mesh->Volus[j]->Faces[1]=
	          v0==mesh->Faces[mesh->Faces[f2]->Old_Child[0]]->Verts[0]?
	          mesh->Faces[f2]->Old_Child[0]:mesh->Faces[f2]->Child[0];
	       }
	       else
	       {
		  mesh->Volus[j]->Faces[1]=
	          v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	          mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Child[1];
	       }		  

	       mesh->Volus[j]->Faces[2]=f3;
               mesh->Volus[j]->Faces[3]=tmp_face;

               mesh->Volus[j]->Type=0;

               mesh->Volus[j]->Child[0]=j;
	       mesh->Volus[j]->Child[1]=idx_volu-1;


	       //ReConstruct(mesh->Volus[idx_volu-1]);
               //ReConstruct(mesh->Volus[idx_volu-2]);
              if(MarkJudge(mesh,mesh->Volus[j])<0)
	      {
	       mesh->Volus[j]->Mark=-1; 
	      }
              if(MarkJudge(mesh,mesh->Volus[idx_volu-1])>=0) 
	      {
               mesh->Volus[idx_volu-1]->Mark=1;
	      }
              else
	       mesh->Volus[idx_volu-1]->Mark=-1;

        }//2
         //---------------------------------------------------------------------------------------------
         else if ( (mesh->Faces[volu->Faces[3]]->Child[0]>=0
		     &&mesh->Faces[volu->Faces[3]]->Child[2]!=volu->Lines[0]
		     &&mesh->Faces[volu->Faces[3]]->Child[2]!=volu->Lines[1]
		     &&mesh->Faces[volu->Faces[3]]->Child[2]!=volu->Lines[3]
		     )
		  &&(mesh->Faces[volu->Faces[2]]->Child[0]>=0
		    &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[0]
		    &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[2]
		    &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[4]
		    )
		  ) 

	   {//3
	    //更新中间的新面的信息
	   mesh->Faces[idx_face]->Verts[0]=mesh->Lines[l0]->Child[2];
           mesh->Faces[idx_face]->Verts[1]=v2;
           mesh->Faces[idx_face]->Verts[2]=v3;
        
	   mesh->Faces[idx_face]->Old_Child[0]=-1;
           mesh->Faces[idx_face]->Old_Child[1]=-1;
 
	   mesh->Faces[idx_face]->Child[0]=-1;
           mesh->Faces[idx_face]->Child[1]=-1;
           mesh->Faces[idx_face]->Child[2]=-1;
 
           mesh->Faces[idx_face]->ID_Boundary=0;

	   mesh->Faces[idx_face]->Lines[0]=l5;

	   if(   (mesh->Faces[f2]->Old_Child[1]>=0)&&(mesh->Faces[f2]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f2]->Old_Child[1]!=l2)&&(mesh->Faces[f2]->Old_Child[1]!=l4)   )
	   {
             mesh->Faces[idx_face]->Lines[1]=mesh->Faces[f2]->Old_Child[1];
	   }
	   else
	   {
             mesh->Faces[idx_face]->Lines[1]=mesh->Faces[f2]->Child[2];
	   }
	   if(   (mesh->Faces[f3]->Old_Child[1]>=0)&&(mesh->Faces[f3]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f3]->Old_Child[1]!=l1)&&(mesh->Faces[f3]->Old_Child[1]!=l3)   )
	   {
             mesh->Faces[idx_face++]->Lines[2]=mesh->Faces[f3]->Old_Child[1];
	   }
	   else
	   {
             mesh->Faces[idx_face++]->Lines[2]=mesh->Faces[f3]->Child[2];
	   }


           //更新单元信息
           mesh->Volus[idx_volu]->Verts[0]=v1;
           mesh->Volus[idx_volu]->Verts[1]=v2;
           mesh->Volus[idx_volu]->Verts[2]=v3;
           mesh->Volus[idx_volu]->Verts[3]=mesh->Lines[l0]->Child[2];

 	   
	   mesh->Volus[idx_volu]->Lines[0]=l3;
           mesh->Volus[idx_volu]->Lines[1]=l4;
           mesh->Volus[idx_volu]->Lines[2]=
           v0==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
	   mesh->Volus[idx_volu]->Lines[3]=l5;
	   /*if(mesh->Faces[volu->Faces[2]]->Old_Child[1]<0)
	   mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f2]->Child[2];
	   else
           mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f2]->Old_Child[1];
           if(mesh->Faces[f3]->Old_Child[1]<0)
            mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f3]->Child[2];
	   else
	   mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f3]->Old_Child[1];
           */
	   if(   (mesh->Faces[f2]->Old_Child[1]>=0)&&(mesh->Faces[f2]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f2]->Old_Child[1]!=l2)&&(mesh->Faces[f2]->Old_Child[1]!=l4)   )
	   {
	     mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f2]->Old_Child[1];
	   }
	   else
	   {
	     mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f2]->Child[2];
	   }
	   if(   (mesh->Faces[f3]->Old_Child[1]>=0)&&(mesh->Faces[f3]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f3]->Old_Child[1]!=l1)&&(mesh->Faces[f3]->Old_Child[1]!=l3)   )
	   {
	     mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f3]->Old_Child[1];
	   }
	   else
	   {
	     mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f3]->Child[2];
	   }
	   //mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f2]->Child[2];
           //mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f3]->Child[2];

	   if(i==0)
   	       {
                mesh->Volus[idx_volu]->Father=j;
                mesh->Volus[j]->Father=j;
               }
           else
           mesh->Volus[idx_volu]->Father=mesh->Volus[j]->Father;
           mesh->Volus[idx_volu]->Ancestor=mesh->Volus[j]->Ancestor;

           mesh->Volus[idx_volu]->Type=0;

	   mesh->Volus[idx_volu]->Faces[0]=idx_face-1;
	   /*
	   if(mesh->Faces[f3]->Old_Child[0]<0)
	   mesh->Volus[idx_volu]->Faces[2]=
           v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	   mesh->Faces[f3]->Child[1]:mesh->Faces[f3]->Child[0];
	   else
	   mesh->Volus[idx_volu]->Faces[2]=
           v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	   mesh->Faces[f3]->Old_Child[0]:mesh->Faces[f3]->Child[0];
	   if(mesh->Faces[f2]->Old_Child[0]<0)
	   mesh->Volus[idx_volu]->Faces[1]=
           v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	   mesh->Faces[f2]->Child[1]:mesh->Faces[f2]->Child[0];
	   else
           mesh->Volus[idx_volu]->Faces[1]=
           v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	   mesh->Faces[f2]->Old_Child[0]:mesh->Faces[f2]->Child[0];
           */
           if(   (mesh->Faces[f2]->Old_Child[1]>=0)&&(mesh->Faces[f2]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f2]->Old_Child[1]!=l2)&&(mesh->Faces[f2]->Old_Child[1]!=l4)   )
	   {
	     mesh->Volus[idx_volu]->Faces[1]=
             v1==mesh->Faces[mesh->Faces[f2]->Old_Child[0]]->Verts[0]?
	     mesh->Faces[f2]->Old_Child[0]:mesh->Faces[f2]->Child[0];
	   }
	   else
	   {
	     mesh->Volus[idx_volu]->Faces[1]=
             v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	     mesh->Faces[f2]->Child[1]:mesh->Faces[f2]->Child[0];
	   }
	   if(   (mesh->Faces[f3]->Old_Child[1]>=0)&&(mesh->Faces[f3]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f3]->Old_Child[1]!=l1)&&(mesh->Faces[f3]->Old_Child[1]!=l3)   )
	   {
	     mesh->Volus[idx_volu]->Faces[2]=
             v1==mesh->Faces[mesh->Faces[f3]->Old_Child[0]]->Verts[0]?
	     mesh->Faces[f3]->Old_Child[0]:mesh->Faces[f3]->Child[0];
	   }
	   else
	   {
	     mesh->Volus[idx_volu]->Faces[2]=
             v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	     mesh->Faces[f3]->Child[1]:mesh->Faces[f3]->Child[0];
	   }

	   //mesh->Volus[idx_volu]->Faces[2]=
           //v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	   //mesh->Faces[f3]->Child[1]:mesh->Faces[f3]->Child[0];
	   //mesh->Volus[idx_volu]->Faces[1]=
           //v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	   //mesh->Faces[f2]->Child[1]:mesh->Faces[f2]->Child[0];

	   mesh->Volus[idx_volu++]->Faces[3]=f0;

           mesh->Volus[j]->Verts[0]=v0;
	   mesh->Volus[j]->Verts[1]=v2;
	   mesh->Volus[j]->Verts[2]=v3;
           mesh->Volus[j]->Verts[3]=mesh->Lines[l0]->Child[2];
	   //printf ( "test identity-=%d,%d,%d,%d\n", volu->Verts[0],volu->Verts[2],volu->Verts[3],mesh->Lines[volu->Lines[0]]->Child[2]);
	  
	   //printf ( "test identity-=%d,%d,%d,%d\n", mesh->Volus[j]->Verts[0],mesh->Volus[j]->Verts[1],mesh->Volus[j]->Verts[2],mesh->Volus[j]->Verts[3]);
	   tmp_line=v0==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];
	   mesh->Volus[j]->Lines[0]=l1;
	   mesh->Volus[j]->Lines[1]=l2;
           mesh->Volus[j]->Lines[2]=tmp_line;
	   mesh->Volus[j]->Lines[3]=l5;
      	   if(   (mesh->Faces[f2]->Old_Child[1]>=0)&&(mesh->Faces[f2]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f2]->Old_Child[1]!=l2)&&(mesh->Faces[f2]->Old_Child[1]!=l4)   )
	   {
	     mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Old_Child[1];
	   }
	   else
	   {
	     mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Child[2];
	   }
	   if(   (mesh->Faces[f3]->Old_Child[1]>=0)&&(mesh->Faces[f3]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f3]->Old_Child[1]!=l1)&&(mesh->Faces[f3]->Old_Child[1]!=l3)   )
	   {
	     mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Old_Child[1];
	   }
	   else
	   {
	     mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Child[2];
	   }

	   //mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Child[2];
           //mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Child[2];

	   ////tmp_face=volu->Faces[1];
	   tmp_face=f1;
	   mesh->Volus[j]->Faces[0]=idx_face-1;
	   /*
	   if(mesh->Faces[f2]->Old_Child[1]<0)
	   mesh->Volus[j]->Faces[1]=
           v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	   mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Child[1];
	   else
           mesh->Volus[j]->Faces[1]=
           v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	   mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Old_Child[0];
           if(mesh->Faces[f3]->Old_Child[1]<0)
           mesh->Volus[j]->Faces[2]=
           v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	   mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Child[1];
	   else
	   mesh->Volus[j]->Faces[2]=
           v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	   mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Old_Child[0];
           */
           if(   (mesh->Faces[f2]->Old_Child[1]>=0)&&(mesh->Faces[f2]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f2]->Old_Child[1]!=l2)&&(mesh->Faces[f2]->Old_Child[1]!=l4)   )
	   {
	     mesh->Volus[j]->Faces[1]=
             v0==mesh->Faces[mesh->Faces[f2]->Old_Child[0]]->Verts[0]?
	     mesh->Faces[f2]->Old_Child[0]:mesh->Faces[f2]->Child[0];
	   }
	   else
	   {
	     mesh->Volus[j]->Faces[1]=
             v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	     mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Child[1];
	   }
	   if(   (mesh->Faces[f3]->Old_Child[1]>=0)&&(mesh->Faces[f3]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f3]->Old_Child[1]!=l1)&&(mesh->Faces[f3]->Old_Child[1]!=l3)   )
	   {
	     mesh->Volus[j]->Faces[2]=
             v0==mesh->Faces[mesh->Faces[f3]->Old_Child[0]]->Verts[0]?
	     mesh->Faces[f3]->Old_Child[0]:mesh->Faces[f3]->Child[0];
	   }
	   else
	   {
	     mesh->Volus[j]->Faces[2]=
             v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	     mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Child[1];
	   }

	   
	   //mesh->Volus[j]->Faces[1]=
           //v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	   //mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Child[1];
           //mesh->Volus[j]->Faces[2]=
           //v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	   //mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Child[1];

	   mesh->Volus[j]->Faces[3]=tmp_face;

           mesh->Volus[j]->Type=0;

	   mesh->Volus[j]->Child[0]=j;
	   mesh->Volus[j]->Child[1]=idx_volu-1;


 	   if(MarkJudge(mesh,mesh->Volus[j])<0)
	   {
	    mesh->Volus[j]->Mark=-1;
	   }
	   if(MarkJudge(mesh,mesh->Volus[idx_volu-1])>=0) 
	   {
            mesh->Volus[idx_volu-1]->Mark=1;
	   }
          else
	   mesh->Volus[idx_volu-1]->Mark=-1;

	   }//3
	   //-------------------------------------------------------------------------------------------
           //第四种情况：第一个面(012)加密，第二个面没加密 
	   else
	   {//4
	       //更新第二个面上的边	   
	       mesh->Lines[idx_line]->Verts[1]=mesh->Lines[l0]->Child[2];
	       mesh->Lines[idx_line]->Child[0]=-1;
	       mesh->Lines[idx_line]->Child[1]=-1;
               mesh->Lines[idx_line]->Child[2]=-1;
	       mesh->Lines[idx_line]->Num_Volus_Owned=-1;
	       mesh->Lines[idx_line]->Volus_Owned=NULL;
               mesh->Lines[idx_line]->Mark=-1;
	       mesh->Lines[idx_line]->ID_Boundary=mesh->Faces[f2]->ID_Boundary;
	       mesh->Lines[idx_line++]->Verts[0]=v3;
	       	       
	       //更新第二个面上的新面
               mesh->Faces[f2]->Verts[0]=v0;
               mesh->Faces[f2]->Verts[1]=mesh->Lines[l0]->Child[2];
               mesh->Faces[f2]->Verts[2]=v3;
              
	       mesh->Faces[f2]->Lines[0]=idx_line-1;
               mesh->Faces[f2]->Lines[1]=volu->Lines[2];
               mesh->Faces[f2]->Lines[2]=
	       v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];

               mesh->Faces[idx_face]->Verts[0]=v1;
               mesh->Faces[idx_face]->Verts[1]=v3;
               mesh->Faces[idx_face]->Verts[2]=mesh->Lines[l0]->Child[2];

	       mesh->Faces[idx_face]->Old_Child[0]=-1;
               mesh->Faces[idx_face]->Old_Child[1]=-1;

               mesh->Faces[idx_face]->Child[0]=-1;
               mesh->Faces[idx_face]->Child[1]=-1;
               mesh->Faces[idx_face]->Child[2]=-1;

               mesh->Faces[idx_face]->ID_Boundary=mesh->Faces[volu->Faces[2]]->ID_Boundary;

               mesh->Faces[idx_face]->Lines[0]=idx_line-1;
               mesh->Faces[idx_face]->Lines[1]=
	       v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
               mesh->Faces[idx_face++]->Lines[2]=l4;

	       mesh->Faces[f2]->Old_Child[0]=mesh->Faces[f2]->Child[1];
	       mesh->Faces[f2]->Old_Child[1]=mesh->Faces[f2]->Child[2];

	       mesh->Faces[f2]->Child[0]=f2;
	       mesh->Faces[f2]->Child[1]=idx_face-1;
	       mesh->Faces[f2]->Child[2]=idx_line-1;
	      
	       //更新中间新面的信息
               mesh->Faces[idx_face]->Verts[0]=mesh->Lines[l0]->Child[2];
               mesh->Faces[idx_face]->Verts[1]=v2;
               mesh->Faces[idx_face]->Verts[2]=v3;
              
	       mesh->Faces[idx_face]->Old_Child[0]=-1;
               mesh->Faces[idx_face]->Old_Child[1]=-1;
               mesh->Faces[idx_face]->Child[0]=-1;
               mesh->Faces[idx_face]->Child[1]=-1;
               mesh->Faces[idx_face]->Child[2]=-1;

               mesh->Faces[idx_face]->ID_Boundary=0;

	       mesh->Faces[idx_face]->Lines[0]=l5;
               mesh->Faces[idx_face]->Lines[1]=idx_line-1;

	       if(   (mesh->Faces[f3]->Old_Child[1]>=0)&&(mesh->Faces[f3]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f3]->Old_Child[1]!=l1)&&(mesh->Faces[f3]->Old_Child[1]!=l3)   )
	       {
                 mesh->Faces[idx_face++]->Lines[2]=mesh->Faces[f3]->Old_Child[1];
	       }
	       else
	       {
                 mesh->Faces[idx_face++]->Lines[2]=mesh->Faces[f3]->Child[2];
	       }

	       //更新两个单元的信息
               mesh->Volus[idx_volu]->Verts[0]=v1;
               mesh->Volus[idx_volu]->Verts[1]=v2;
               mesh->Volus[idx_volu]->Verts[2]=v3;
               mesh->Volus[idx_volu]->Verts[3]=mesh->Lines[l0]->Child[2];

               mesh->Volus[idx_volu]->Lines[0]=l3;
               mesh->Volus[idx_volu]->Lines[1]=l4;
               mesh->Volus[idx_volu]->Lines[2]=
	       v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
               mesh->Volus[idx_volu]->Lines[3]=l5;
               mesh->Volus[idx_volu]->Lines[5]=idx_line-1;
	       /* if(mesh->Faces[f3]->Old_Child[1]<0)
               mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f3]->Child[2];
               else
	       mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f3]->Old_Child[1];
               */

	       if(   (mesh->Faces[f3]->Old_Child[1]>=0)&&(mesh->Faces[f3]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f3]->Old_Child[1]!=l1)&&(mesh->Faces[f3]->Old_Child[1]!=l3)   )
	       {
                 mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f3]->Old_Child[1];
	       }
	       else
	       {
                 mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f3]->Child[2];
	       }

               //mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f3]->Child[2];
	       if(i==0)
   	       {
                mesh->Volus[idx_volu]->Father=j;
                mesh->Volus[j]->Father=j;
               }
               else
               mesh->Volus[idx_volu]->Father=mesh->Volus[j]->Father;
               mesh->Volus[idx_volu]->Ancestor=mesh->Volus[j]->Ancestor;

               mesh->Volus[idx_volu]->Type=0;

               mesh->Volus[idx_volu]->Faces[0]=idx_face-1;
	       /*
	       if(mesh->Faces[f3]->Old_Child[0]<0)
               mesh->Volus[idx_volu]->Faces[2]=
               v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	       mesh->Faces[f3]->Child[1]:mesh->Faces[f3]->Child[0];
               else
	       mesh->Volus[idx_volu]->Faces[2]=
               v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	       mesh->Faces[f3]->Old_Child[0]:mesh->Faces[f3]->Child[0];
               */
               if(   (mesh->Faces[f3]->Old_Child[1]>=0)&&(mesh->Faces[f3]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f3]->Old_Child[1]!=l1)&&(mesh->Faces[f3]->Old_Child[1]!=l3)   )
	       {
		  mesh->Volus[idx_volu]->Faces[2]=
                  v1==mesh->Faces[mesh->Faces[f3]->Old_Child[0]]->Verts[0]?
	          mesh->Faces[f3]->Old_Child[0]:mesh->Faces[f3]->Child[0];
	       }
	       else
	       {
		  mesh->Volus[idx_volu]->Faces[2]=
                  v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	          mesh->Faces[f3]->Child[1]:mesh->Faces[f3]->Child[0];
	       }	
	       //mesh->Volus[idx_volu]->Faces[2]=
               //v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	       //mesh->Faces[f3]->Child[1]:mesh->Faces[f3]->Child[0];
	       mesh->Volus[idx_volu]->Faces[1]=idx_face-2;
               mesh->Volus[idx_volu++]->Faces[3]=f0;


	       mesh->Volus[j]->Verts[0]=v0;
               mesh->Volus[j]->Verts[1]=v2;
               mesh->Volus[j]->Verts[2]=v3;
               mesh->Volus[j]->Verts[3]=mesh->Lines[l0]->Child[2];

               tmp_line= v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];
	       mesh->Volus[j]->Lines[0]=l1;
               mesh->Volus[j]->Lines[1]=l2;
               mesh->Volus[j]->Lines[2]=tmp_line;
               mesh->Volus[j]->Lines[3]=l5;
	       /*if(mesh->Faces[f3]->Old_Child[1]<0)
               mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Child[2];
	       else
               mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Old_Child[1];
               */

	       if(   (mesh->Faces[f3]->Old_Child[1]>=0)&&(mesh->Faces[f3]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f3]->Old_Child[1]!=l1)&&(mesh->Faces[f3]->Old_Child[1]!=l3)   )
	       {
                 mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Old_Child[1];
	       }
	       else
	       {
                 mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Child[2];
	       }	
               //mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Child[2];
               mesh->Volus[j]->Lines[5]=idx_line-1;

	       //tmp_face=volu->Faces[1];
	       tmp_face=f1;
               mesh->Volus[j]->Faces[0]=idx_face-1;
               mesh->Volus[j]->Faces[1]=f2;
               /*  
	       if(mesh->Faces[f3]->Old_Child[1]<0)
               mesh->Volus[j]->Faces[2]=
	       v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	       mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Child[1];
	       else
	       mesh->Volus[j]->Faces[2]=
	       v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	       mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Old_Child[0];
               */
               if(   (mesh->Faces[f3]->Old_Child[1]>=0)&&(mesh->Faces[f3]->Old_Child[1]!=l0)&&
	         (mesh->Faces[f3]->Old_Child[1]!=l1)&&(mesh->Faces[f3]->Old_Child[1]!=l3)   )
	       {
		  mesh->Volus[j]->Faces[2]=
	          v0==mesh->Faces[mesh->Faces[f3]->Old_Child[0]]->Verts[0]?
	          mesh->Faces[f3]->Old_Child[0]:mesh->Faces[f3]->Child[0];
	       }
	       else
	       {
		  mesh->Volus[j]->Faces[2]=
	          v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	          mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Child[1];
	       }		
	       //mesh->Volus[j]->Faces[2]=
	       //v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	       //mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Child[1];

               mesh->Volus[j]->Faces[3]=tmp_face;

               mesh->Volus[j]->Type=0;


	       mesh->Volus[j]->Child[0]=j;
	       mesh->Volus[j]->Child[1]=idx_volu-1;


	       //ReConstruct(mesh->Volus[idx_volu-1]);
               //ReConstruct(mesh->Volus[idx_volu-2]);
           
	       if(MarkJudge(mesh,mesh->Volus[j])<0)
	       {
	        mesh->Volus[j]->Mark=-1;
	       }
               if(MarkJudge(mesh,mesh->Volus[idx_volu-1])>=0) 
	       {
                mesh->Volus[idx_volu-1]->Mark=1;
	       }
               else
	        mesh->Volus[idx_volu-1]->Mark=-1;
             }//4
         }//type2

     }//单元遍历
     printf ( "neibu jiesu\n" );
     Iter_volu=idx_volu;
     printf ( "itervolu==%d\n",idx_volu );
  }//加密三次循环
  printf ( "the end of refine\n" );
  mesh->Num_Verts_Global=total_vert;
  mesh->Num_Lines_Global=total_line;
  mesh->Num_Faces_Global=total_face;
  mesh->Num_Volus_Global=total_volu;
   
}//MeshRefineBisection结束






//这个函数是二分加密函数 
 void
 MeshRefineUniformBisectionNew1(MESH *mesh)
 { 
  int i,j,s;
  int k=0;
  int tmp_line,tmp_face;
  int nlines_mark=0;
  VOLU *volu,*volu_owned;
  FACEE *face;
  LINE *line;
  //找到所有要加密的单元
  for(i=0;i<mesh->Num_Volus_Global;++i)
  {
    mesh->Volus[i]->Sign=1;
    mesh->Volus[i]->Mark=1;
  }
  for(i=0;i<mesh->Num_Lines_Global;++i)
  {
    line=mesh->Lines[i];
    line->Mark=i; 
  }
  nlines_mark=mesh->Num_Lines_Global;
  //首先计算加密后需要的存储空间
  //计算新增加的边数
  int mark_line_on_face=0;
  for(i=0;i<mesh->Num_Faces_Global;++i)
  {   
    face=mesh->Faces[i];
     for(j=0;j<3;++j)
     {
       if(mesh->Lines[face->Lines[j]]->Mark>=0)
       {
	  mark_line_on_face++;
       }
     }
  }     
  int mark_line_on_volu=0;
  int mark_line23_on_volu=0;
  for(i=0;i<mesh->Num_Volus_Global;++i)
  {
     volu=mesh->Volus[i];
     for(j=0;j<5;++j)
     {
	if(mesh->Lines[volu->Lines[j]]->Mark>=0)
	   mark_line_on_volu++;
     }
     if(mesh->Lines[volu->Lines[5]]->Mark>=0)
	   mark_line23_on_volu++;
  }
  int total_vert,total_line, total_face,total_volu;
  total_vert=mesh->Num_Verts_Global+nlines_mark;
  total_line=mesh->Num_Lines_Global+nlines_mark+mark_line_on_face+mark_line23_on_volu;
  total_face=mesh->Num_Faces_Global+mark_line_on_face+mark_line_on_volu+mark_line23_on_volu*3;
  total_volu=mesh->Num_Volus_Global+mark_line_on_volu+mark_line23_on_volu*2;
  mesh->Verts=(VERT **)realloc(mesh->Verts,total_vert);
  for(i=mesh->Num_Verts_Global;i<total_vert;++i)
  mesh->Verts[i]=(VERT *)malloc(1);
  
  mesh->Lines=(LINE **)realloc(mesh->Lines,total_line);
  for(i=mesh->Num_Lines_Global;i<total_line;++i)
  mesh->Lines[i]=(LINE *)malloc(1);
  
  mesh->Faces=(FACEE **)realloc(mesh->Faces,total_face);
  for(i=mesh->Num_Faces_Global;i<total_face;++i)
  mesh->Faces[i]=(FACEE *)malloc(1);
  
  mesh->Volus=(VOLU **)realloc(mesh->Volus,total_volu);
  for(i=mesh->Num_Volus_Global;i<total_volu;++i)
  mesh->Volus[i]=(VOLU *)malloc(1);
  
 //加密程序
 
 //首先更新边上的新点和新边
 for(i=0;i<mesh->Num_Lines_Global;++i)
 {
   if(mesh->Lines[i]->Mark>=0)
   {
    //新点  
    mesh->Verts[mesh->Num_Verts_Global+mesh->Lines[i]->Mark]->Coord[0]=   
    0.5*(mesh->Verts[mesh->Lines[i]->Verts[0]]->Coord[0]+mesh->Verts[mesh->Lines[i]->Verts[1]]->Coord[0]);
    mesh->Verts[mesh->Num_Verts_Global+mesh->Lines[i]->Mark]->Coord[1]=   
    0.5*(mesh->Verts[mesh->Lines[i]->Verts[0]]->Coord[1]+mesh->Verts[mesh->Lines[i]->Verts[1]]->Coord[1]);
    mesh->Verts[mesh->Num_Verts_Global+mesh->Lines[i]->Mark]->Coord[2]=   
    0.5*(mesh->Verts[mesh->Lines[i]->Verts[0]]->Coord[2]+mesh->Verts[mesh->Lines[i]->Verts[1]]->Coord[2]);
    //两条边
    mesh->Lines[mesh->Num_Lines_Global+mesh->Lines[i]->Mark]->Verts[0]=mesh->Num_Verts_Global
                                                                        +mesh->Lines[i]->Mark; 
    mesh->Lines[mesh->Num_Lines_Global+mesh->Lines[i]->Mark]->Verts[1]=mesh->Lines[i]->Verts[1];
    mesh->Lines[mesh->Num_Lines_Global+mesh->Lines[i]->Mark]->Mark=-1;
    mesh->Lines[mesh->Num_Lines_Global+mesh->Lines[i]->Mark]->Child[0]=-1;
    mesh->Lines[mesh->Num_Lines_Global+mesh->Lines[i]->Mark]->Child[1]=-1;
    mesh->Lines[mesh->Num_Lines_Global+mesh->Lines[i]->Mark]->Child[2]=-1;
    mesh->Lines[mesh->Num_Lines_Global+mesh->Lines[i]->Mark]->ID_Boundary=mesh->Lines[i]->ID_Boundary;
    mesh->Lines[mesh->Num_Lines_Global+mesh->Lines[i]->Mark]->Num_Volus_Owned=-1;
    mesh->Lines[mesh->Num_Lines_Global+mesh->Lines[i]->Mark]->Volus_Owned=NULL;
 
    mesh->Lines[i]->Verts[0]=mesh->Lines[i]->Verts[0];
    mesh->Lines[i]->Verts[1]=mesh->Num_Verts_Global+mesh->Lines[i]->Mark;
    mesh->Lines[i]->Child[0]=i;
    mesh->Lines[i]->Child[1]=mesh->Num_Lines_Global+mesh->Lines[i]->Mark;
    mesh->Lines[i]->Child[2]=mesh->Num_Verts_Global+mesh->Lines[i]->Mark;
   }
 }
  int Iter_volu=mesh->Num_Volus_Global;
  int idx_vert=mesh->Num_Verts_Global+nlines_mark;
  int idx_line=mesh->Num_Lines_Global+nlines_mark;
  int idx_face=mesh->Num_Faces_Global;
  int idx_volu=mesh->Num_Volus_Global;
  int v0,v1,v2,v3,l0,l1,l2,l3,l4,l5,f0,f1,f2,f3;
  for(i=0;i<3;++i)
 {//加密三次循环
    for(j=0;j<Iter_volu;++j)//Iter_volu是需要循环的单元个数，每一次循环后都要加上新产生的单元个数
    {//单元遍历
       volu=mesh->Volus[j];
       v0=volu->Verts[0];
       v1=volu->Verts[1];
       v2=volu->Verts[2];
       v3=volu->Verts[3];
       l0=volu->Lines[0];
       l1=volu->Lines[1];
       l2=volu->Lines[2];
       l3=volu->Lines[3];
       l4=volu->Lines[4];
       l5=volu->Lines[5];
       f0=volu->Faces[0];
       f1=volu->Faces[1];
       f2=volu->Faces[2];
       f3=volu->Faces[3];

      if(volu->Mark<0) continue;
      if(volu->Type==0)
      {//type0
          //首先是两个相邻面都需要加密的情况   
          if (  (mesh->Faces[volu->Faces[3]]->Child[0]<0&&mesh->Faces[volu->Faces[2]]->Child[0]<0)
	      ||(mesh->Faces[volu->Faces[3]]->Child[0]<0&&mesh->Faces[volu->Faces[2]]->Child[0]>=0
		  &&(mesh->Faces[volu->Faces[2]]->Child[2]==volu->Lines[0]||mesh->Faces[volu->Faces[2]]
	         ->Child[2]==volu->Lines[2]||mesh->Faces[volu->Faces[2]]->Child[2]==volu->Lines[4]))
	      ||(mesh->Faces[volu->Faces[2]]->Child[0]<0&&mesh->Faces[volu->Faces[3]]->Child[0]>=0
	         &&(mesh->Faces[volu->Faces[3]]->Child[2]==volu->Lines[0]||mesh->Faces[volu->Faces[3]]
		    ->Child[2]==volu->Lines[1]||mesh->Faces[volu->Faces[3]]->Child[2]==volu->Lines[3]) )
	      ||(mesh->Faces[volu->Faces[3]]->Child[0]>=0&&(mesh->Faces[volu->Faces[3]]->Child[2]==
		    volu->Lines[0]||mesh->Faces[volu->Faces[3]]
		    ->Child[2]==volu->Lines[1]||mesh->Faces[volu->Faces[3]]->Child[2]==volu->Lines[3]) 
		 &&mesh->Faces[volu->Faces[2]]->Child[0]>=0&&(mesh->Faces[volu->Faces[2]]->Child[2]==
		    volu->Lines[0]||mesh->Faces[volu->Faces[2]]
	         ->Child[2]==volu->Lines[2]||mesh->Faces[volu->Faces[2]]->Child[2]==volu->Lines[4])
		 )
	     )   
	  {//1
	   //更新面上两条新边的信息 
	   mesh->Lines[idx_line]->Verts[1]= mesh->Lines[l0]->Child[2];
           mesh->Lines[idx_line]->Mark=-1;
	   mesh->Lines[idx_line]->ID_Boundary=mesh->Faces[f3]->ID_Boundary;
	   mesh->Lines[idx_line]->Num_Volus_Owned=-1;
	   mesh->Lines[idx_line]->Volus_Owned=NULL;
           mesh->Lines[idx_line]->Child[0]=-1;
	   mesh->Lines[idx_line]->Child[1]=-1;
	   mesh->Lines[idx_line]->Child[2]=-1;
	   mesh->Lines[idx_line++]->Verts[0]=v2;

	   mesh->Lines[idx_line]->Verts[1]= mesh->Lines[l0]->Child[2];
           mesh->Lines[idx_line]->Mark=-1;
	   mesh->Lines[idx_line]->Num_Volus_Owned=-1;
	   mesh->Lines[idx_line]->Volus_Owned=NULL;
	   mesh->Lines[idx_line]->ID_Boundary=mesh->Faces[f2]->ID_Boundary;
	   mesh->Lines[idx_line]->Child[0]=-1;
	   mesh->Lines[idx_line]->Child[1]=-1;
	   mesh->Lines[idx_line]->Child[2]=-1;
           mesh->Lines[idx_line++]->Verts[0]=v3;
           //更新第一个面上的信息
           mesh->Faces[idx_face]->Verts[0]=v1;
           mesh->Faces[idx_face]->Verts[1]=v2;
           mesh->Faces[idx_face]->Verts[2]=mesh->Lines[l0]->Child[2];

	   mesh->Faces[idx_face]->Child[0]=-1;
           mesh->Faces[idx_face]->Child[1]=-1;
           mesh->Faces[idx_face]->Child[2]=-1;
           mesh->Faces[idx_face]->Old_Child[0]=-1;
           mesh->Faces[idx_face]->Old_Child[1]=-1;

	   mesh->Faces[idx_face]->ID_Boundary=mesh->Faces[f3]->ID_Boundary;

           mesh->Faces[idx_face]->Lines[0]=idx_line-2;
           mesh->Faces[idx_face]->Lines[1]=
	   v0==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
           mesh->Faces[idx_face++]->Lines[2]=l3;

	   mesh->Faces[f3]->Verts[0]=volu->Verts[0];
           mesh->Faces[f3]->Verts[1]=mesh->Lines[l0]->Child[2];
           mesh->Faces[f3]->Verts[2]=v2;
         
	   mesh->Faces[f3]->Lines[0]=idx_line-2;
           mesh->Faces[f3]->Lines[1]=l1;
           mesh->Faces[f3]->Lines[2]=
           v0==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];

	   mesh->Faces[f3]->Old_Child[0]=mesh->Faces[f3]->Child[1];
	   mesh->Faces[f3]->Old_Child[1]=mesh->Faces[f3]->Child[2];
	   mesh->Faces[f3]->Child[0]=f3;
           mesh->Faces[f3]->Child[1]=idx_face-1;
           mesh->Faces[f3]->Child[2]=idx_line-2;
	   //更新第二个面上的信息
           mesh->Faces[idx_face]->Verts[0]=v1;
           mesh->Faces[idx_face]->Verts[1]=v3;
           mesh->Faces[idx_face]->Verts[2]=mesh->Lines[l0]->Child[2];

	   mesh->Faces[idx_face]->Child[0]=-1;
           mesh->Faces[idx_face]->Child[1]=-1;
           mesh->Faces[idx_face]->Child[2]=-1;
      	   mesh->Faces[idx_face]->Old_Child[0]=-1;
           mesh->Faces[idx_face]->Old_Child[1]=-1;
           mesh->Faces[idx_face]->ID_Boundary=mesh->Faces[f2]->ID_Boundary;//?????????????????

           mesh->Faces[idx_face]->Lines[0]=idx_line-1;
           mesh->Faces[idx_face]->Lines[1]=
           v0==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
           mesh->Faces[idx_face++]->Lines[2]=l4;

	   mesh->Faces[f2]->Verts[0]=v0;
           mesh->Faces[f2]->Verts[1]=mesh->Lines[l0]->Child[2];
           mesh->Faces[f2]->Verts[2]=volu->Verts[3];
         
           mesh->Faces[f2]->Old_Child[0]=mesh->Faces[f2]->Child[1];
           mesh->Faces[f2]->Old_Child[1]=mesh->Faces[f2]->Child[2];
	   mesh->Faces[f2]->Child[0]=f2;
           mesh->Faces[f2]->Child[1]=idx_face-1;
           mesh->Faces[f2]->Child[2]=idx_line-1;

	   mesh->Faces[f2]->Lines[0]=idx_line-1;
           mesh->Faces[f2]->Lines[1]=l2;
           mesh->Faces[f2]->Lines[2]=
           v0==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];

           //更新中间的新面的信息
	   mesh->Faces[idx_face]->Verts[0]=mesh->Lines[l0]->Child[2];
           mesh->Faces[idx_face]->Verts[1]=v2;
           mesh->Faces[idx_face]->Verts[2]=v3;
         
	   mesh->Faces[idx_face]->Child[0]=-1;
           mesh->Faces[idx_face]->Child[1]=-1;
           mesh->Faces[idx_face]->Child[2]=-1;
           mesh->Faces[idx_face]->Old_Child[0]=-1;
           mesh->Faces[idx_face]->Old_Child[1]=-1;

	   mesh->Faces[idx_face]->ID_Boundary=0;

	   mesh->Faces[idx_face]->Lines[0]=l5;
           mesh->Faces[idx_face]->Lines[1]=idx_line-1;
           mesh->Faces[idx_face++]->Lines[2]=idx_line-2;

	      //更新单元信息
	   mesh->Volus[idx_volu]->Verts[0]=v1;
	   mesh->Volus[idx_volu]->Verts[1]=v3;
	   mesh->Volus[idx_volu]->Verts[2]=v2;
           mesh->Volus[idx_volu]->Verts[3]=mesh->Lines[l0]->Child[2];
           
	   mesh->Volus[idx_volu]->Lines[0]=l4;
	   mesh->Volus[idx_volu]->Lines[1]=l3;
	   mesh->Volus[idx_volu]->Lines[2]=
	   volu->Verts[0]==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];//yong volu->verts shi fou you wen ti????
 	   mesh->Volus[idx_volu]->Lines[3]=l5;
           mesh->Volus[idx_volu]->Lines[4]=idx_line-1;
 	   mesh->Volus[idx_volu]->Lines[5]=idx_line-2;

	   if(i==0)
	   {
            mesh->Volus[idx_volu]->Father=j;
            mesh->Volus[j]->Father=j;
	   }
           else
            mesh->Volus[idx_volu]->Father=mesh->Volus[j]->Father;
            mesh->Volus[idx_volu]->Ancestor=mesh->Volus[j]->Ancestor;

           mesh->Volus[idx_volu]->Type=1;

	   mesh->Volus[idx_volu]->Faces[0]=idx_face-1;
	   mesh->Volus[idx_volu]->Faces[1]=idx_face-3;
	   mesh->Volus[idx_volu]->Faces[2]=idx_face-2;
           mesh->Volus[idx_volu++]->Faces[3]=f0;

	   mesh->Volus[j]->Verts[0]=v0;
	   mesh->Volus[j]->Verts[1]=v2;
	   mesh->Volus[j]->Verts[2]=v3;
	   mesh->Volus[j]->Verts[3]=mesh->Lines[l0]->Child[2];


	   tmp_line=volu->Verts[0]==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];
           mesh->Volus[j]->Lines[0]=l1;
	   mesh->Volus[j]->Lines[1]=l2;
	   mesh->Volus[j]->Lines[2]=tmp_line;
	   mesh->Volus[j]->Lines[3]=l5;
           mesh->Volus[j]->Lines[4]=idx_line-2;
 	   mesh->Volus[j]->Lines[5]=idx_line-1;

	   tmp_face=volu->Faces[1];
	   mesh->Volus[j]->Faces[0]=idx_face-1;
	   mesh->Volus[j]->Faces[1]=f2;
	   mesh->Volus[j]->Faces[2]=f3;
	   mesh->Volus[j]->Faces[3]=tmp_face;

           mesh->Volus[j]->Type=1;

	   mesh->Volus[j]->Child[0]=j;
	   mesh->Volus[j]->Child[1]=idx_volu-1;

	  
	   if(MarkJudge(mesh,mesh->Volus[j])<0)
	   {
            mesh->Volus[j]->Mark=-1;
	   }
           if(MarkJudge(mesh,mesh->Volus[idx_volu-1])>=0) 
	   {
            mesh->Volus[idx_volu-1]->Mark=1;
	   }
	   else
	   mesh->Volus[idx_volu-1]->Mark=-1;

	  }//1
           //-----------------------------------------------------------------------------------
           //第二种情况：第一个面没加密，第二个面加密
	 else if(  (mesh->Faces[volu->Faces[3]]->Child[0]<0
		     &&mesh->Faces[volu->Faces[2]]->Child[0]>=0
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[0]
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[2]
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[4])
		  ||(mesh->Faces[volu->Faces[3]]->Child[0]>=0
		     &&(mesh->Faces[volu->Faces[3]]->Child[2]==
		    volu->Lines[0]||mesh->Faces[volu->Faces[3]]
		    ->Child[2]==volu->Lines[1]||mesh->Faces[volu->Faces[3]]->Child[2]==volu->Lines[3]) 
		     &&mesh->Faces[volu->Faces[2]]->Child[0]>=0
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[0]
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[2]
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[4]
		    )
	       )

	   //else if(mesh->Faces[volu->Faces[3]]->Child[0]<0&&mesh->Faces[volu->Faces[2]]->Child[0]>=0)
	   {//2
               mesh->Lines[idx_line]->Verts[1]=mesh->Lines[l0]->Child[2];
	       mesh->Lines[idx_line]->Child[0]=-1;
	       mesh->Lines[idx_line]->Child[1]=-1;
               mesh->Lines[idx_line]->Child[2]=-1;
	       mesh->Lines[idx_line]->Num_Volus_Owned=-1;
	       mesh->Lines[idx_line]->Volus_Owned=NULL;
               mesh->Lines[idx_line]->Mark=-1;
	       mesh->Lines[idx_line]->ID_Boundary=mesh->Faces[f3]->ID_Boundary;
	       mesh->Lines[idx_line++]->Verts[0]=v2;
	       	       
	       //更新第一个面上的新面
               mesh->Faces[f3]->Verts[0]=v0;
               mesh->Faces[f3]->Verts[1]=mesh->Lines[l0]->Child[2];
               mesh->Faces[f3]->Verts[2]=v2;
             
	       mesh->Faces[f3]->Old_Child[0]=mesh->Faces[f3]->Child[1];
               mesh->Faces[f3]->Old_Child[1]=mesh->Faces[f3]->Child[2];
               mesh->Faces[f3]->Child[0]=f3;
               mesh->Faces[f3]->Child[1]=idx_face;
               mesh->Faces[f3]->Child[2]=idx_line-1;

	       mesh->Faces[f3]->Lines[0]=idx_line-1;
               mesh->Faces[f3]->Lines[1]=l1;
               mesh->Faces[f3]->Lines[2]=
	       volu->Verts[0]==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];

               mesh->Faces[idx_face]->Verts[0]=volu->Verts[1];
               mesh->Faces[idx_face]->Verts[1]=volu->Verts[2];
               mesh->Faces[idx_face]->Verts[2]=mesh->Lines[l0]->Child[2];

               mesh->Faces[idx_face]->Child[0]=-1;
               mesh->Faces[idx_face]->Child[1]=-1;
               mesh->Faces[idx_face]->Child[2]=-1;

	       mesh->Faces[idx_face]->Old_Child[0]=-1;
               mesh->Faces[idx_face]->Old_Child[1]=-1;

               mesh->Faces[idx_face]->ID_Boundary=mesh->Faces[f3]->ID_Boundary;

               mesh->Faces[idx_face]->Lines[0]=idx_line-1;
               mesh->Faces[idx_face]->Lines[1]=
	       v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
               mesh->Faces[idx_face++]->Lines[2]=l3;
	       //更新中间新面的信息
               mesh->Faces[idx_face]->Verts[0]=mesh->Lines[l0]->Child[2];
               mesh->Faces[idx_face]->Verts[1]=v2;
               mesh->Faces[idx_face]->Verts[2]=v3;
              
               mesh->Faces[idx_face]->Child[0]=-1;
               mesh->Faces[idx_face]->Child[1]=-1;
               mesh->Faces[idx_face]->Child[2]=-1;

       	       mesh->Faces[idx_face]->Old_Child[0]=-1;
               mesh->Faces[idx_face]->Old_Child[1]=-1;

               mesh->Faces[idx_face]->ID_Boundary=0;

	       mesh->Faces[idx_face]->Lines[0]=l5;
               mesh->Faces[idx_face]->Lines[1]=mesh->Faces[f2]->Child[2];
               mesh->Faces[idx_face++]->Lines[2]=idx_line-1;
          
	       //更新两个单元的信息
               mesh->Volus[idx_volu]->Verts[0]=v1;
               mesh->Volus[idx_volu]->Verts[1]=v3;
               mesh->Volus[idx_volu]->Verts[2]=v2;
               mesh->Volus[idx_volu]->Verts[3]=mesh->Lines[l0]->Child[2];

               mesh->Volus[idx_volu]->Lines[0]=volu->Lines[4];
               mesh->Volus[idx_volu]->Lines[1]=volu->Lines[3];
               mesh->Volus[idx_volu]->Lines[2]=
	       volu->Verts[0]==mesh->Lines[volu->Lines[0]]->Verts[0]?
	       mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
               mesh->Volus[idx_volu]->Lines[3]=l5;
               /*
	       if(mesh->Faces[f2]->Old_Child[1]<0)
               mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f2]->Child[2];
	       else
               mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f2]->Old_Child[1];
	       */
	       mesh->Volus[idx_volu]->Lines[5]=idx_line-1;
               mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f2]->Child[2];

	       if(i==0)
	       {
                mesh->Volus[idx_volu]->Father=j;
                mesh->Volus[j]->Father=j;
	       }
	       else
	       mesh->Volus[idx_volu]->Father=mesh->Volus[j]->Father;
               mesh->Volus[idx_volu]->Ancestor=mesh->Volus[j]->Ancestor;
               mesh->Volus[idx_volu]->Type=1;

               mesh->Volus[idx_volu]->Faces[0]=idx_face-1;
               mesh->Volus[idx_volu]->Faces[1]=idx_face-2;
               /*
	       if(mesh->Faces[f2]->Old_Child[1]<0)
	       {
	       mesh->Volus[idx_volu]->Faces[2]=
               v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	       mesh->Faces[f2]->Child[1]:mesh->Faces[f2]->Child[0];
	       }
	       else
	      {
	       mesh->Volus[idx_volu]->Faces[2]=
               v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	       mesh->Faces[f2]->Old_Child[0]:mesh->Faces[f2]->Child[0];
	       }
	       */
               mesh->Volus[idx_volu]->Faces[2]=
               v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	       mesh->Faces[f2]->Child[1]:mesh->Faces[f2]->Child[0];

	       mesh->Volus[idx_volu++]->Faces[3]=volu->Faces[0];

	       mesh->Volus[j]->Verts[0]=v0;
               mesh->Volus[j]->Verts[1]=v2;
               mesh->Volus[j]->Verts[2]=v3;
               mesh->Volus[j]->Verts[3]=mesh->Lines[l0]->Child[2];

	       tmp_line= v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];
               mesh->Volus[j]->Lines[0]=l1;
               mesh->Volus[j]->Lines[1]=l2;
               mesh->Volus[j]->Lines[2]=tmp_line;
               mesh->Volus[j]->Lines[3]=l5;
               mesh->Volus[j]->Lines[4]=idx_line-1;
               /*
	       if(mesh->Faces[f2]->Old_Child[1]<0)//gaid
               mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f2]->Child[2];
	       else
               mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f2]->Old_Child[1];
               */
               mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Child[2];
	       tmp_face=volu->Faces[1];
	       mesh->Volus[j]->Faces[0]=idx_face-1;
	       /*
	       if(mesh->Faces[f2]->Old_Child[1]<0)
	       {
	       mesh->Volus[j]->Faces[1]=
	       v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	       mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Child[1];
               mesh->Volus[j]->Faces[2]=f3;
               mesh->Volus[j]->Faces[3]=tmp_face;
	       }
	       else
	       {
	       mesh->Volus[j]->Faces[1]=
               v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	       mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Old_Child[0];
               mesh->Volus[j]->Faces[2]=f3;
               mesh->Volus[j]->Faces[3]=tmp_face;
	       }
	       */
               mesh->Volus[j]->Faces[1]=
	       v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	       mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Child[1];
               mesh->Volus[j]->Faces[2]=f3;
               mesh->Volus[j]->Faces[3]=tmp_face;

	       mesh->Volus[j]->Type=1;

               mesh->Volus[j]->Child[0]=j;
	       mesh->Volus[j]->Child[1]=idx_volu-1;


	       //ReConstruct(mesh->Volus[idx_volu-1]);
               //ReConstruct(mesh->Volus[idx_volu-2]);
           if(MarkJudge(mesh,mesh->Volus[j])<0)
	   {
	      mesh->Volus[j]->Mark=-1; 
	   }
           if(MarkJudge(mesh,mesh->Volus[idx_volu-1])>=0) 
	   {
              mesh->Volus[idx_volu-1]->Mark=1;
	   }
           else
	   mesh->Volus[idx_volu-1]->Mark=-1;

        }//2
         //---------------------------------------------------------------------------------------------
	  else if ( (mesh->Faces[volu->Faces[3]]->Child[0]>=0
		     &&mesh->Faces[volu->Faces[3]]->Child[2]!=volu->Lines[0]
		     &&mesh->Faces[volu->Faces[3]]->Child[2]!=volu->Lines[1]
		     &&mesh->Faces[volu->Faces[3]]->Child[2]!=volu->Lines[3]
		     )
		  &&(mesh->Faces[volu->Faces[2]]->Child[0]>=0
		    &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[0]
		    &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[2]
		    &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[4]
		    )
		  ) 
	  //else if(mesh->Faces[volu->Faces[3]]->Child[0]>=0&&mesh->Faces[volu->Faces[2]]->Child[0]>=0)
	  //第三种情况：两个面都加密 
	   {//3
	   mesh->Faces[idx_face]->Verts[0]=mesh->Lines[l0]->Child[2];
           mesh->Faces[idx_face]->Verts[1]=v2;
           mesh->Faces[idx_face]->Verts[2]=v3;
         
	   mesh->Faces[idx_face]->Child[0]=-1;
           mesh->Faces[idx_face]->Child[1]=-1;
           mesh->Faces[idx_face]->Child[2]=-1;
 	   mesh->Faces[idx_face]->Old_Child[0]=-1;
           mesh->Faces[idx_face]->Old_Child[1]=-1;

           mesh->Faces[idx_face]->ID_Boundary=0;
	   

	   mesh->Faces[idx_face]->Lines[0]=volu->Lines[5];
           mesh->Faces[idx_face]->Lines[1]=mesh->Faces[f2]->Child[2];
           mesh->Faces[idx_face++]->Lines[2]=mesh->Faces[f3]->Child[2];

           //更新单元信息
	   mesh->Volus[idx_volu]->Verts[0]=v1;
	   mesh->Volus[idx_volu]->Verts[1]=v3;
	   mesh->Volus[idx_volu]->Verts[2]=v2;
           mesh->Volus[idx_volu]->Verts[3]=mesh->Lines[l0]->Child[2];

 	   mesh->Volus[idx_volu]->Lines[0]=l4;
           mesh->Volus[idx_volu]->Lines[1]=l3;
           mesh->Volus[idx_volu]->Lines[2]=
           v0==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
	   mesh->Volus[idx_volu]->Lines[3]=l5;
	   //mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[volu->Faces[2]]->Child[2];
	   //mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[volu->Faces[3]]->Child[2];

           /*
           if(mesh->Faces[volu->Faces[2]]->Old_Child[0]<0)    
	      mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f2]->Child[2];
	   else
              mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f2]->Old_Child[1];

           if(mesh->Faces[f3]->Old_Child[1]<0)//gaid     
	      mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f3]->Child[2];
	   else
              mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f3]->Old_Child[1];
           */
	    mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f2]->Child[2];
	    mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f3]->Child[2];

	   if(i==0)
	   {
            mesh->Volus[idx_volu]->Father=j;
            mesh->Volus[j]->Father=j;
	   }
           else
           mesh->Volus[idx_volu]->Father=mesh->Volus[j]->Father;
           mesh->Volus[idx_volu]->Ancestor=mesh->Volus[j]->Ancestor;

           mesh->Volus[idx_volu]->Type=1;

	   mesh->Volus[idx_volu]->Faces[0]=idx_face-1;
	   /*
	   if(mesh->Faces[f3]->Old_Child[0]<0)
	   mesh->Volus[idx_volu]->Faces[1]=
           volu->Verts[0]==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	   mesh->Faces[f3]->Child[1]:mesh->Faces[f3]->Child[0];
	   else
	   mesh->Volus[idx_volu]->Faces[1]=
           v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	   mesh->Faces[f3]->Old_Child[0]:mesh->Faces[f3]->Child[0];
           */
           mesh->Volus[idx_volu]->Faces[1]=
           volu->Verts[0]==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	   mesh->Faces[f3]->Child[1]:mesh->Faces[f3]->Child[0];

           /*
	   if(mesh->Faces[f2]->Old_Child[0]<0)
	   mesh->Volus[idx_volu]->Faces[2]=
           v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	   mesh->Faces[f2]->Child[1]:mesh->Faces[f2]->Child[0];
	   else
           mesh->Volus[idx_volu]->Faces[2]=
           v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	   mesh->Faces[f2]->Old_Child[0]:mesh->Faces[f2]->Child[0];
           */
           mesh->Volus[idx_volu]->Faces[2]=
           v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	   mesh->Faces[f2]->Child[1]:mesh->Faces[f2]->Child[0];

	   mesh->Volus[idx_volu++]->Faces[3]=f0;
           
	   mesh->Volus[j]->Verts[0]=v0;
	   mesh->Volus[j]->Verts[1]=v2;
	   mesh->Volus[j]->Verts[2]=v3;
           mesh->Volus[j]->Verts[3]=mesh->Lines[l0]->Child[2];
           tmp_line=volu->Verts[0]==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];
	   //mesh->Volus[j]->Lines[4]=mesh->Faces[volu->Faces[3]]->Child[2];
	   mesh->Volus[j]->Lines[3]=l5;
           //mesh->Volus[j]->Lines[5]=mesh->Faces[volu->Faces[2]]->Child[2];
	   mesh->Volus[j]->Lines[0]=l1;
	   mesh->Volus[j]->Lines[1]=l2;
           mesh->Volus[j]->Lines[2]=tmp_line;
	   /*   
	   if(mesh->Faces[volu->Faces[2]]->Old_Child[1]<0)//gaid    
	     mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Child[2];
   	   else
             mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Old_Child[1];
	   if(mesh->Faces[volu->Faces[3]]->Old_Child[1]<0)      
	     mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Child[2];
	   else
             mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Old_Child[1];
           */
	    mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Child[2];
	    mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Child[2];

	   tmp_face=f1;
	   mesh->Volus[j]->Faces[0]=idx_face-1;
	   /*
	   if(mesh->Faces[f2]->Old_Child[1]<0)
	   mesh->Volus[j]->Faces[1]=
           v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	   mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Child[1];
	   else
           mesh->Volus[j]->Faces[1]=
           volu->Verts[0]==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	   mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Old_Child[0];
           if(mesh->Faces[f3]->Child[1]<0)
	   mesh->Volus[j]->Faces[2]=
           v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	   mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Child[1];
	   else
           mesh->Volus[j]->Faces[2]=
           v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	   mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Old_Child[0];
           */
           mesh->Volus[j]->Faces[3]=tmp_face;
           mesh->Volus[j]->Faces[1]=
           v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	   mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Child[1];
           mesh->Volus[j]->Faces[2]=
           v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	   mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Child[1];

           mesh->Volus[j]->Type=1;

	   mesh->Volus[j]->Child[0]=j;
	   mesh->Volus[j]->Child[1]=idx_volu-1;

 	   if(MarkJudge(mesh,mesh->Volus[j])<0)
	   {
	    mesh->Volus[j]->Mark=-1;
	   }
           if(MarkJudge(mesh,mesh->Volus[idx_volu-1])>=0) 
	   {
            mesh->Volus[idx_volu-1]->Mark=1;
	   }
           else
	   mesh->Volus[idx_volu-1]->Mark=-1;

	   }//3
	   //-------------------------------------------------------------------------------------------
           //第四种情况：第一个面(012)加密，第二个面没加密 
	   else
	   {//4
               mesh->Lines[idx_line]->Verts[1]=mesh->Lines[l0]->Child[2];
	       mesh->Lines[idx_line]->Child[0]=-1;
	       mesh->Lines[idx_line]->Child[1]=-1;
               mesh->Lines[idx_line]->Child[2]=-1;
	       mesh->Lines[idx_line]->Num_Volus_Owned=-1;
	       mesh->Lines[idx_line]->Volus_Owned=NULL;
               mesh->Lines[idx_line]->Mark=-1;
	       mesh->Lines[idx_line]->ID_Boundary=mesh->Faces[f2]->ID_Boundary;
	       mesh->Lines[idx_line++]->Verts[0]=v3;
	       	       
	       //更新第二个面上的新面
               mesh->Faces[f2]->Verts[0]=v0;
               mesh->Faces[f2]->Verts[1]=mesh->Lines[l0]->Child[2];
               mesh->Faces[f2]->Verts[2]=v3;
             
               mesh->Faces[f2]->Old_Child[0]=mesh->Faces[f2]->Child[1];
               mesh->Faces[f2]->Old_Child[1]=mesh->Faces[f2]->Child[2];
               mesh->Faces[f2]->Child[0]=f2;
               mesh->Faces[f2]->Child[1]=idx_face;
               mesh->Faces[f2]->Child[2]=idx_line-1;

	       mesh->Faces[f2]->Lines[0]=idx_line-1;
               mesh->Faces[f2]->Lines[1]=l2;
               mesh->Faces[f2]->Lines[2]=
	       volu->Verts[0]==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];

               mesh->Faces[idx_face]->Verts[0]=v1;
               mesh->Faces[idx_face]->Verts[1]=v3;
               mesh->Faces[idx_face]->Verts[2]=mesh->Lines[l0]->Child[2];

               mesh->Faces[idx_face]->Child[0]=-1;
               mesh->Faces[idx_face]->Child[1]=-1;
               mesh->Faces[idx_face]->Child[2]=-1;

	       mesh->Faces[idx_face]->Old_Child[0]=-1;
               mesh->Faces[idx_face]->Old_Child[1]=-1;

               mesh->Faces[idx_face]->ID_Boundary=mesh->Faces[f2]->ID_Boundary;

               mesh->Faces[idx_face]->Lines[0]=idx_line-1;
               mesh->Faces[idx_face]->Lines[1]=
	       v0==mesh->Lines[volu->Lines[0]]->Verts[0]?
	       mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
               mesh->Faces[idx_face++]->Lines[2]=l4;

	       //更新中间新面的信息
               mesh->Faces[idx_face]->Verts[0]=mesh->Lines[l0]->Child[2];
               mesh->Faces[idx_face]->Verts[1]=v2;
               mesh->Faces[idx_face]->Verts[2]=v3;
              
               mesh->Faces[idx_face]->Child[0]=-1;
               mesh->Faces[idx_face]->Child[1]=-1;
               mesh->Faces[idx_face]->Child[2]=-1;

	       mesh->Faces[idx_face]->Old_Child[0]=-1;
               mesh->Faces[idx_face]->Old_Child[1]=-1;

               mesh->Faces[idx_face]->ID_Boundary=0;

	       mesh->Faces[idx_face]->Lines[0]=l5;
               mesh->Faces[idx_face]->Lines[1]=idx_line-1;
               mesh->Faces[idx_face++]->Lines[2]=mesh->Faces[f3]->Child[2];

	       //更新两个单元的信息
               mesh->Volus[idx_volu]->Verts[0]=v1;
               mesh->Volus[idx_volu]->Verts[1]=v3;
               mesh->Volus[idx_volu]->Verts[2]=v2;
               mesh->Volus[idx_volu]->Verts[3]=mesh->Lines[l0]->Child[2];

               mesh->Volus[idx_volu]->Lines[0]=l4;
               mesh->Volus[idx_volu]->Lines[1]=l3;
	       mesh->Volus[idx_volu]->Lines[2]=
	       v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
               
	       mesh->Volus[idx_volu]->Lines[3]=l5;
               mesh->Volus[idx_volu]->Lines[4]=idx_line-1;
               /*  
	       if(mesh->Faces[f3]->Old_Child[1]<0)//gaid
	       mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f3]->Child[2];
               else
               mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f3]->Old_Child[1];
               */
	       mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f3]->Child[2];
    	       if(i==0)
   	       {
                mesh->Volus[idx_volu]->Father=j;
                mesh->Volus[j]->Father=j;
	       }
               else
               mesh->Volus[idx_volu]->Father=mesh->Volus[j]->Father;
               mesh->Volus[idx_volu]->Ancestor=mesh->Volus[j]->Ancestor;

               mesh->Volus[idx_volu]->Type=1;

               mesh->Volus[idx_volu]->Faces[0]=idx_face-1;
               /*
	       if(mesh->Faces[f3]->Old_Child[1]<0)
               mesh->Volus[idx_volu]->Faces[1]=
               v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	       mesh->Faces[f3]->Child[1]:mesh->Faces[f3]->Child[0];
	       else
               mesh->Volus[idx_volu]->Faces[1]=
               v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	       mesh->Faces[f3]->Old_Child[0]:mesh->Faces[f3]->Child[0];
               */
               mesh->Volus[idx_volu]->Faces[1]=
               v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	       mesh->Faces[f3]->Child[1]:mesh->Faces[f3]->Child[0];

	       mesh->Volus[idx_volu]->Faces[2]=idx_face-2;
               mesh->Volus[idx_volu++]->Faces[3]=f0;

	       mesh->Volus[j]->Verts[0]=v0;
               mesh->Volus[j]->Verts[1]=v2;
               mesh->Volus[j]->Verts[2]=v3;
               mesh->Volus[j]->Verts[3]=mesh->Lines[l0]->Child[2];

	       tmp_line= v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];
               mesh->Volus[j]->Lines[0]=l1;
               mesh->Volus[j]->Lines[1]=l2;
               mesh->Volus[j]->Lines[2]=tmp_line;
               mesh->Volus[j]->Lines[3]=l5;
	       /*
	       if(mesh->Faces[f3]->Old_Child[0]<0)
               mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Child[2];
	       else
               mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Old_Child[1];
               */
               mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Child[2];
               mesh->Volus[j]->Lines[5]=idx_line-1;

	       tmp_face=volu->Faces[1];
               mesh->Volus[j]->Faces[0]=idx_face-1;
               mesh->Volus[j]->Faces[1]=volu->Faces[2];
	       /*
	       if( mesh->Faces[f3]->Old_Child[0]<0)
               mesh->Volus[j]->Faces[2]=
	       v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	       mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Child[1];
               else
               mesh->Volus[j]->Faces[2]=
	       v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	       mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Old_Child[0];
               */
	       mesh->Volus[j]->Faces[2]=
	       v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	       mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Child[1];

	       mesh->Volus[j]->Faces[3]=tmp_face;

               mesh->Volus[j]->Type=1;

	       mesh->Volus[j]->Child[0]=j;
	       mesh->Volus[j]->Child[1]=idx_volu-1;


	       //ReConstruct(mesh->Volus[idx_volu-1]);
               //ReConstruct(mesh->Volus[idx_volu-2]);
           
	      if(MarkJudge(mesh,mesh->Volus[j])<0)
	      {
	       mesh->Volus[j]->Mark=-1;
	      }
              if(MarkJudge(mesh,mesh->Volus[idx_volu-1])>=0) 
	      {
               mesh->Volus[idx_volu-1]->Mark=1;
	      }
              else
	       mesh->Volus[idx_volu-1]->Mark=-1;
            }//4
          }//type0
      
          else if(volu->Type==1)
          {//type1
	   //nprintf ( "%d,%d,%d,%d\n", volu->Faces[1],volu->Faces[2],volu->Faces[3],volu->Faces[4]);
           if (  (mesh->Faces[volu->Faces[3]]->Child[0]<0&&mesh->Faces[volu->Faces[2]]->Child[0]<0)
	      ||(mesh->Faces[volu->Faces[3]]->Child[0]<0&&mesh->Faces[volu->Faces[2]]->Child[0]>=0
		  &&(mesh->Faces[volu->Faces[2]]->Child[2]==volu->Lines[0]||mesh->Faces[volu->Faces[2]]
	         ->Child[2]==volu->Lines[2]||mesh->Faces[volu->Faces[2]]->Child[2]==volu->Lines[4]))
	      ||(mesh->Faces[volu->Faces[2]]->Child[0]<0&&mesh->Faces[volu->Faces[3]]->Child[0]>=0
	         &&(mesh->Faces[volu->Faces[3]]->Child[2]==volu->Lines[0]||mesh->Faces[volu->Faces[3]]
		    ->Child[2]==volu->Lines[1]||mesh->Faces[volu->Faces[3]]->Child[2]==volu->Lines[3]) )
	      ||(mesh->Faces[volu->Faces[3]]->Child[0]>=0&&(mesh->Faces[volu->Faces[3]]->Child[2]==
		    volu->Lines[0]||mesh->Faces[volu->Faces[3]]
		    ->Child[2]==volu->Lines[1]||mesh->Faces[volu->Faces[3]]->Child[2]==volu->Lines[3]) 
		 &&mesh->Faces[volu->Faces[2]]->Child[0]>=0&&(mesh->Faces[volu->Faces[2]]->Child[2]==
		    volu->Lines[0]||mesh->Faces[volu->Faces[2]]
	         ->Child[2]==volu->Lines[2]||mesh->Faces[volu->Faces[2]]->Child[2]==volu->Lines[4])
		 )
	     )   

          //if(mesh->Faces[volu->Faces[3]]->Child[0]<0
	  // &&mesh->Faces[volu->Faces[2]]->Child[0]<0)//第一种情况：两个面都没加密 
          {//1
	   //更新面上两条新边的信息 
	   mesh->Lines[idx_line]->Verts[1]= mesh->Lines[l0]->Child[2];
	   mesh->Lines[idx_line]->Num_Volus_Owned=-1;
	   mesh->Lines[idx_line]->Volus_Owned=NULL;
           mesh->Lines[idx_line]->Mark=-1;
	   mesh->Lines[idx_line]->ID_Boundary=mesh->Faces[f3]->ID_Boundary;
	   mesh->Lines[idx_line]->Child[0]=-1;
	   mesh->Lines[idx_line]->Child[1]=-1;
	   mesh->Lines[idx_line]->Child[2]=-1;
	   mesh->Lines[idx_line++]->Verts[0]=v2;

	   mesh->Lines[idx_line]->Verts[1]= mesh->Lines[l0]->Child[2];
	   mesh->Lines[idx_line]->Num_Volus_Owned=-1;
	   mesh->Lines[idx_line]->Volus_Owned=NULL;
           mesh->Lines[idx_line]->Mark=-1;
	   mesh->Lines[idx_line]->ID_Boundary=mesh->Faces[f2]->ID_Boundary;
	   mesh->Lines[idx_line]->Child[0]=-1;
	   mesh->Lines[idx_line]->Child[1]=-1;
	   mesh->Lines[idx_line]->Child[2]=-1;
           mesh->Lines[idx_line++]->Verts[0]=v3;
           //更新第一个面上的信息
           mesh->Faces[idx_face]->Verts[0]=v1;
           mesh->Faces[idx_face]->Verts[1]=v2;
           mesh->Faces[idx_face]->Verts[2]=mesh->Lines[l0]->Child[2];

	   mesh->Faces[idx_face]->Old_Child[0]=-1;
           mesh->Faces[idx_face]->Old_Child[1]=-1;

	   mesh->Faces[idx_face]->Child[0]=-1;
           mesh->Faces[idx_face]->Child[1]=-1;
           mesh->Faces[idx_face]->Child[2]=-1;
           mesh->Faces[idx_face]->ID_Boundary=mesh->Faces[f3]->ID_Boundary;

           mesh->Faces[idx_face]->Lines[0]=idx_line-2;
           mesh->Faces[idx_face]->Lines[1]=
	   volu->Verts[0]==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
           mesh->Faces[idx_face++]->Lines[2]=l3;

	   mesh->Faces[f3]->Verts[0]=v0;
           mesh->Faces[f3]->Verts[1]=mesh->Lines[l0]->Child[2];
           mesh->Faces[f3]->Verts[2]=v2;
        
           mesh->Faces[f3]->Old_Child[0]=mesh->Faces[f3]->Child[1];
	   mesh->Faces[f3]->Old_Child[1]=mesh->Faces[f3]->Child[2];
	   mesh->Faces[f3]->Child[0]=volu->Faces[3];
           mesh->Faces[f3]->Child[1]=idx_face-1;
           mesh->Faces[f3]->Child[2]=idx_line-2;

	   mesh->Faces[f3]->Lines[0]=idx_line-2;
           mesh->Faces[f3]->Lines[1]=l1;
           mesh->Faces[f3]->Lines[2]=
           volu->Verts[0]==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];

           //更新第二个面上的信息
           mesh->Faces[idx_face]->Verts[0]=v1;
           mesh->Faces[idx_face]->Verts[1]=v3;
           mesh->Faces[idx_face]->Verts[2]=mesh->Lines[l0]->Child[2];

           mesh->Faces[idx_face]->Old_Child[0]=-1;
           mesh->Faces[idx_face]->Old_Child[1]=-1;

	   mesh->Faces[idx_face]->Child[0]=-1;
           mesh->Faces[idx_face]->Child[1]=-1;
           mesh->Faces[idx_face]->Child[2]=-1;
           mesh->Faces[idx_face]->ID_Boundary=mesh->Faces[f2]->ID_Boundary;

           mesh->Faces[idx_face]->Lines[0]=idx_line-1;
           mesh->Faces[idx_face]->Lines[1]=
           volu->Verts[0]==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
           mesh->Faces[idx_face++]->Lines[2]=volu->Lines[4];

	   mesh->Faces[f2]->Verts[0]=v0;
           mesh->Faces[f2]->Verts[1]=mesh->Lines[l0]->Child[2];
           mesh->Faces[f2]->Verts[2]=v3;
        
           mesh->Faces[f2]->Old_Child[0]= mesh->Faces[f2]->Child[1];
           mesh->Faces[f2]->Old_Child[1]= mesh->Faces[f2]->Child[2];
	   mesh->Faces[f2]->Child[0]=f2;
           mesh->Faces[f2]->Child[1]=idx_face-1;
           mesh->Faces[f2]->Child[2]=idx_line-1;

	   mesh->Faces[f2]->Lines[0]=idx_line-1;
           mesh->Faces[f2]->Lines[1]=l2;
           mesh->Faces[f2]->Lines[2]=
           v0==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];

           //更新中间的新面的信息
	   mesh->Faces[idx_face]->Verts[0]=mesh->Lines[l0]->Child[2];
           mesh->Faces[idx_face]->Verts[1]=v2;
           mesh->Faces[idx_face]->Verts[2]=v3;
	   
	   mesh->Faces[idx_face]->Old_Child[0]=-1;
           mesh->Faces[idx_face]->Old_Child[1]=-1;

	   mesh->Faces[idx_face]->Child[0]=-1;
           mesh->Faces[idx_face]->Child[1]=-1;
           mesh->Faces[idx_face]->Child[2]=-1;
         
	   mesh->Faces[idx_face]->ID_Boundary=0;

	   mesh->Faces[idx_face]->Lines[0]=l5;
           mesh->Faces[idx_face]->Lines[1]=idx_line-1;
           mesh->Faces[idx_face++]->Lines[2]=idx_line-2;


	      //更新单元信息
	   mesh->Volus[idx_volu]->Verts[0]=v1;
	   mesh->Volus[idx_volu]->Verts[1]=v2;
	   mesh->Volus[idx_volu]->Verts[2]=v3;
           mesh->Volus[idx_volu]->Verts[3]=mesh->Lines[l0]->Child[2];
           
	   mesh->Volus[idx_volu]->Lines[0]=l3;
	   mesh->Volus[idx_volu]->Lines[1]=l4;
	   mesh->Volus[idx_volu]->Lines[2]=
	   v0==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
 	   mesh->Volus[idx_volu]->Lines[3]=l5;
           mesh->Volus[idx_volu]->Lines[4]=idx_line-2;
 	   mesh->Volus[idx_volu]->Lines[5]=idx_line-1;

	   if(i==0)
   	   {
            mesh->Volus[idx_volu]->Father=j;
            mesh->Volus[j]->Father=j;
           }
           else
           mesh->Volus[idx_volu]->Father=mesh->Volus[j]->Father;
           mesh->Volus[idx_volu]->Ancestor=mesh->Volus[j]->Ancestor;


           mesh->Volus[idx_volu]->Type=2;

	   mesh->Volus[idx_volu]->Faces[0]=idx_face-1;
	   mesh->Volus[idx_volu]->Faces[1]=idx_face-2;
	   mesh->Volus[idx_volu]->Faces[2]=idx_face-3;
           mesh->Volus[idx_volu++]->Faces[3]=volu->Faces[0];

	   mesh->Volus[j]->Verts[0]=v0;
	   mesh->Volus[j]->Verts[1]=v2;
	   mesh->Volus[j]->Verts[2]=v3;
	   mesh->Volus[j]->Verts[3]=mesh->Lines[l0]->Child[2];

	   tmp_line=volu->Verts[0]==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];
           mesh->Volus[j]->Lines[0]=l1;
	   mesh->Volus[j]->Lines[1]=l2;
	   mesh->Volus[j]->Lines[2]=tmp_line;
 	   mesh->Volus[j]->Lines[3]=l5;
           mesh->Volus[j]->Lines[4]=idx_line-2;
 	   mesh->Volus[j]->Lines[5]=idx_line-1;

	   tmp_face=volu->Faces[1];
	   mesh->Volus[j]->Faces[0]=idx_face-1;
	   mesh->Volus[j]->Faces[1]=f2;
	   mesh->Volus[j]->Faces[2]=f3;
	   mesh->Volus[j]->Faces[3]=tmp_face;

           mesh->Volus[j]->Type=2;

	   mesh->Volus[j]->Child[0]=j;
	   mesh->Volus[j]->Child[1]=idx_volu-1;

           //ReConstruct(mesh->Volus[idx_volu-1]);
           //ReConstruct(mesh->Volus[MarkSet[i]]);
           if(MarkJudge(mesh,mesh->Volus[j])<0)
	   {
            mesh->Volus[j]->Mark=-1;
	   }
           if(MarkJudge(mesh,mesh->Volus[idx_volu-1])>=0) 
	   {
            mesh->Volus[idx_volu-1]->Mark=1;
	   }
           else
	    mesh->Volus[idx_volu-1]->Mark=-1;
         }//1
           //-----------------------------------------------------------------------------------
           //第二种情况：第一个面没加密，第二个面加密 
	 else if(  (mesh->Faces[volu->Faces[3]]->Child[0]<0
		     &&mesh->Faces[volu->Faces[2]]->Child[0]>=0
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[0]
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[2]
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[4])
		  ||(mesh->Faces[volu->Faces[3]]->Child[0]>=0
		     &&(mesh->Faces[volu->Faces[3]]->Child[2]==
		    volu->Lines[0]||mesh->Faces[volu->Faces[3]]
		    ->Child[2]==volu->Lines[1]||mesh->Faces[volu->Faces[3]]->Child[2]==volu->Lines[3]) 
		     &&mesh->Faces[volu->Faces[2]]->Child[0]>=0
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[0]
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[2]
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[4]
		    )
	       )

	   {//2
	      //更新第一个面上的边	   
               mesh->Lines[idx_line]->Verts[1]=mesh->Lines[l0]->Child[2];
	       mesh->Lines[idx_line]->Child[0]=-1;
	       mesh->Lines[idx_line]->Child[1]=-1;
               mesh->Lines[idx_line]->Child[2]=-1;
	       mesh->Lines[idx_line]->Num_Volus_Owned=-1;
	       mesh->Lines[idx_line]->Volus_Owned=NULL;
               mesh->Lines[idx_line]->Mark=-1;
	       mesh->Lines[idx_line]->ID_Boundary=mesh->Faces[volu->Faces[3]]->ID_Boundary;
	       mesh->Lines[idx_line++]->Verts[0]=v2;
	       	       
	       //更新第一个面上的新面
               mesh->Faces[f3]->Verts[0]=v0;
               mesh->Faces[f3]->Verts[1]=mesh->Lines[l0]->Child[2];
               mesh->Faces[f3]->Verts[2]=v2;
            
	       mesh->Faces[f3]->Old_Child[0]=mesh->Faces[f3]->Child[1];
               mesh->Faces[f3]->Old_Child[1]=mesh->Faces[f3]->Child[2];
  
               mesh->Faces[f3]->Child[0]=f3;
               mesh->Faces[f3]->Child[1]=idx_face;
               mesh->Faces[f3]->Child[2]=idx_line-1;

	       mesh->Faces[f3]->Lines[0]=idx_line-1;
               mesh->Faces[f3]->Lines[1]=l1;
               mesh->Faces[f3]->Lines[2]=
	       v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];

               mesh->Faces[idx_face]->Verts[0]=v1;
               mesh->Faces[idx_face]->Verts[1]=v2;
               mesh->Faces[idx_face]->Verts[2]=mesh->Lines[l0]->Child[2];

               mesh->Faces[idx_face]->Child[0]=-1;
               mesh->Faces[idx_face]->Child[1]=-1;
               mesh->Faces[idx_face]->Child[2]=-1;

	       mesh->Faces[idx_face]->Old_Child[0]=-1;
               mesh->Faces[idx_face]->Old_Child[1]=-1;

               mesh->Faces[idx_face]->ID_Boundary=mesh->Faces[f3]->ID_Boundary;

               mesh->Faces[idx_face]->Lines[0]=idx_line-1;
               mesh->Faces[idx_face]->Lines[1]=
	       v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
               mesh->Faces[idx_face++]->Lines[2]=l3;

	       //更新中间新面的信息
               mesh->Faces[idx_face]->Verts[0]=mesh->Lines[l0]->Child[2];
               mesh->Faces[idx_face]->Verts[1]=v2;
               mesh->Faces[idx_face]->Verts[2]=v3;
              
               mesh->Faces[idx_face]->Child[0]=-1;
               mesh->Faces[idx_face]->Child[1]=-1;
               mesh->Faces[idx_face]->Child[2]=-1;

	       mesh->Faces[idx_face]->Old_Child[0]=-1;
               mesh->Faces[idx_face]->Old_Child[1]=-1;

               mesh->Faces[idx_face]->ID_Boundary=0;

	       mesh->Faces[idx_face]->Lines[0]=l5;
               mesh->Faces[idx_face]->Lines[1]=mesh->Faces[f2]->Child[2];
               mesh->Faces[idx_face++]->Lines[2]=idx_line-1;
          
	       //更新两个单元的信息
               mesh->Volus[idx_volu]->Verts[0]=v1;
               mesh->Volus[idx_volu]->Verts[1]=v2;
               mesh->Volus[idx_volu]->Verts[2]=v3;
               mesh->Volus[idx_volu]->Verts[3]=mesh->Lines[l0]->Child[2];

               mesh->Volus[idx_volu]->Lines[0]=l3;
               mesh->Volus[idx_volu]->Lines[1]=l4;
               mesh->Volus[idx_volu]->Lines[2]=
	       v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
               mesh->Volus[idx_volu]->Lines[3]=l5;
               /*  if(mesh->Faces[f2]->Old_Child[1]<0)
	       mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f2]->Child[2];
	       else
               mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f2]->Old_Child[1];
               */
	       mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f2]->Child[2];
	       mesh->Volus[idx_volu]->Lines[4]=idx_line-1;

               if(i==0)
   	       {
                mesh->Volus[idx_volu]->Father=j;
                mesh->Volus[j]->Father=j;
               }
               else
               mesh->Volus[idx_volu]->Father=mesh->Volus[j]->Father;
               mesh->Volus[idx_volu]->Ancestor=mesh->Volus[j]->Ancestor;


               mesh->Volus[idx_volu]->Type=2;

               mesh->Volus[idx_volu]->Faces[0]=idx_face-1;
               mesh->Volus[idx_volu]->Faces[2]=idx_face-2;
	       /*
	       if( mesh->Faces[volu->Faces[2]]->Old_Child[1]<0)
               mesh->Volus[idx_volu]->Faces[1]=
               v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	       mesh->Faces[f2]->Child[1]:mesh->Faces[f2]->Child[0];
	       else
	       mesh->Volus[idx_volu]->Faces[1]=
               v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	       mesh->Faces[f2]->Old_Child[0]:mesh->Faces[f2]->Child[0];
               */
               mesh->Volus[idx_volu]->Faces[1]=
               v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	       mesh->Faces[f2]->Child[1]:mesh->Faces[f2]->Child[0];

	       mesh->Volus[idx_volu++]->Faces[3]=f0;

	       mesh->Volus[j]->Verts[0]=v0;
               mesh->Volus[j]->Verts[1]=v2;
               mesh->Volus[j]->Verts[2]=v3;
               mesh->Volus[j]->Verts[3]=mesh->Lines[l0]->Child[2];

	       tmp_line=v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];
               mesh->Volus[j]->Lines[0]=l1;
               mesh->Volus[j]->Lines[1]=l2;
               mesh->Volus[j]->Lines[2]=tmp_line;
               mesh->Volus[j]->Lines[3]=l5;
               mesh->Volus[j]->Lines[4]=idx_line-1;
	       /* if(mesh->Faces[f2]->Old_Child[1]<0)
               mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Child[2];
	       else
               mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Old_Child[1];
               */
               mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Child[2];

	       tmp_face=f1;
               mesh->Volus[j]->Faces[0]=idx_face-1;
	       /*
	       if(mesh->Faces[f2]->Old_Child[1]<0)
               mesh->Volus[j]->Faces[1]=
	       v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	       mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Child[1];
	       else
	       mesh->Volus[j]->Faces[1]=
	       v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	       mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Old_Child[0];
               */
               mesh->Volus[j]->Faces[1]=
	       v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	       mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Child[1];

	       mesh->Volus[j]->Faces[2]=f3;
               mesh->Volus[j]->Faces[3]=tmp_face;

               mesh->Volus[j]->Type=2;


               mesh->Volus[j]->Child[0]=j;
	       mesh->Volus[j]->Child[1]=idx_volu-1;


	       //ReConstruct(mesh->Volus[idx_volu-1]);
               //ReConstruct(mesh->Volus[idx_volu-2]);
               if(MarkJudge(mesh,mesh->Volus[j])<0)
	       {
	        mesh->Volus[j]->Mark=-1; 
	       }
               if(MarkJudge(mesh,mesh->Volus[idx_volu-1])>=0) 
	       {
                mesh->Volus[idx_volu-1]->Mark=1;
	       }
               else
	        mesh->Volus[idx_volu-1]->Mark=-1;
           }//2
          //---------------------------------------------------------------------------------------------
	  else if ( (mesh->Faces[volu->Faces[3]]->Child[0]>=0
		     &&mesh->Faces[volu->Faces[3]]->Child[2]!=volu->Lines[0]
		     &&mesh->Faces[volu->Faces[3]]->Child[2]!=volu->Lines[1]
		     &&mesh->Faces[volu->Faces[3]]->Child[2]!=volu->Lines[3]
		     )
		  &&(mesh->Faces[volu->Faces[2]]->Child[0]>=0
		    &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[0]
		    &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[2]
		    &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[4]
		    )
		  ) 

	 // else if(mesh->Faces[volu->Faces[3]]->Child[0]>=0&&mesh->Faces[volu->Faces[2]]->Child[0]>=0)//第三种情况：两个面都加密 
	   {//3
           //更新中间的新面的信息
           mesh->Faces[idx_face]->Verts[0]=mesh->Lines[l0]->Child[2];
           mesh->Faces[idx_face]->Verts[1]=v2;
           mesh->Faces[idx_face]->Verts[2]=v3;
         
	   mesh->Faces[idx_face]->Old_Child[0]=-1;
           mesh->Faces[idx_face]->Old_Child[1]=-1;

	   mesh->Faces[idx_face]->Child[0]=-1;
           mesh->Faces[idx_face]->Child[1]=-1;
           mesh->Faces[idx_face]->Child[2]=-1;
 
           mesh->Faces[idx_face]->ID_Boundary=0;

	   mesh->Faces[idx_face]->Lines[0]=l5;
           mesh->Faces[idx_face]->Lines[1]=mesh->Faces[f2]->Child[2];
           mesh->Faces[idx_face++]->Lines[2]=mesh->Faces[f3]->Child[2];


           //更新单元信息
	   mesh->Volus[idx_volu]->Verts[0]=v1;
	   mesh->Volus[idx_volu]->Verts[1]=v2;
	   mesh->Volus[idx_volu]->Verts[2]=v3;
           mesh->Volus[idx_volu]->Verts[3]=mesh->Lines[l0]->Child[2];

 	   mesh->Volus[idx_volu]->Lines[0]=l3;
           mesh->Volus[idx_volu]->Lines[1]=l4;
           mesh->Volus[idx_volu]->Lines[2]=
           v0==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
	   mesh->Volus[idx_volu]->Lines[3]=l5;
	   /*
	   if(mesh->Faces[f2]->Old_Child[1]<0)
	   mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f2]->Child[2];
	   else
	   mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f2]->Old_Child[1];
           if(mesh->Faces[f3]->Old_Child[1]<0)
	   mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f3]->Child[2];
           else
           mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f3]->Old_Child[1];
          */
	   mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f3]->Child[2];
	   mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f2]->Child[2];

	   if(i==0)
   	       {
                mesh->Volus[idx_volu]->Father=j;
                mesh->Volus[j]->Father=j;
               }
           else
           mesh->Volus[idx_volu]->Father=mesh->Volus[j]->Father;
           mesh->Volus[idx_volu]->Ancestor=mesh->Volus[j]->Ancestor;


           mesh->Volus[idx_volu]->Type=2;

	   mesh->Volus[idx_volu]->Faces[0]=idx_face-1;
	   /*
	   if(mesh->Faces[f3]->Old_Child[0]<0)
	   mesh->Volus[idx_volu]->Faces[2]=
           v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	   mesh->Faces[f3]->Child[1]:mesh->Faces[f3]->Child[0];
	   else
           mesh->Volus[idx_volu]->Faces[2]=
           v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	   mesh->Faces[f3]->Old_Child[0]:mesh->Faces[f3]->Child[0];
	   if(mesh->Faces[f2]->Old_Child[0]<0)
	   mesh->Volus[idx_volu]->Faces[1]=
           v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	   mesh->Faces[f2]->Child[1]:mesh->Faces[f2]->Child[0];
	   else
	   mesh->Volus[idx_volu]->Faces[1]=
           v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	   mesh->Faces[f2]->Old_Child[0]:mesh->Faces[f2]->Child[0];
           */
           mesh->Volus[idx_volu]->Faces[2]=
           v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	   mesh->Faces[f3]->Child[1]:mesh->Faces[f3]->Child[0];
	   mesh->Volus[idx_volu]->Faces[1]=
           v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	   mesh->Faces[f2]->Child[1]:mesh->Faces[f2]->Child[0];
	   mesh->Volus[idx_volu++]->Faces[3]=f0;

           mesh->Volus[j]->Verts[0]=v0;
	   mesh->Volus[j]->Verts[1]=v2;
	   mesh->Volus[j]->Verts[2]=v3;
           mesh->Volus[j]->Verts[3]=mesh->Lines[l0]->Child[2];

	   tmp_line=v0==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];
	   mesh->Volus[j]->Lines[0]=l1;
	   mesh->Volus[j]->Lines[1]=l2;
           mesh->Volus[j]->Lines[2]=tmp_line;
	   mesh->Volus[j]->Lines[3]=l5;
           /*if(mesh->Faces[f3]->Old_Child[1]<0)
	   mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Child[2];
	   else
   	   mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Old_Child[1];
	   if(mesh->Faces[f2]->Old_Child[1]<0)
	   else
           mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Old_Child[1];
           */

	   mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Child[2];
           mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Child[2];

	   tmp_face=f1;
	   mesh->Volus[j]->Faces[0]=idx_face-1;
	   /*
	   if(mesh->Faces[f2]->Old_Child[1]<0)
	   mesh->Volus[j]->Faces[1]=
           v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	   mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Child[1];
	   else
	   mesh->Volus[j]->Faces[1]=
           v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	   mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Old_Child[0];
           if(mesh->Faces[f3]->Old_Child[1]<0)
	   mesh->Volus[j]->Faces[2]=
           v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	   mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Child[1];
	   else
	   mesh->Volus[j]->Faces[2]=
           v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	   mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Old_Child[0];
           */
	   mesh->Volus[j]->Faces[1]=
           v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	   mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Child[1];
	   mesh->Volus[j]->Faces[2]=
           v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	   mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Child[1];
	   mesh->Volus[j]->Faces[3]=tmp_face;

           mesh->Volus[j]->Type=2;

	   mesh->Volus[j]->Child[0]=j;
	   mesh->Volus[j]->Child[1]=idx_volu-1;


           //ReConstruct(mesh->Volus[idx_volu-1]);
           //ReConstruct(mesh->Volus[idx_volu-2]);
 	   if(MarkJudge(mesh,mesh->Volus[j])<0)
	   {
	    mesh->Volus[j]->Mark=-1;
	   }
           if(MarkJudge(mesh,mesh->Volus[idx_volu-1])>=0) 
	   {
            mesh->Volus[idx_volu-1]->Mark=1;
	   }
          else
	   mesh->Volus[idx_volu-1]->Mark=-1;


	   }//3
	   //-------------------------------------------------------------------------------------------
           //第四种情况：第一个面(012)加密，第二个面没加密 
	   else
	   {//4
	      //更新第二个面上的边	   
               mesh->Lines[idx_line]->Verts[1]=mesh->Lines[l0]->Child[2];
	       mesh->Lines[idx_line]->Child[0]=-1;
	       mesh->Lines[idx_line]->Child[1]=-1;
               mesh->Lines[idx_line]->Child[2]=-1;
	       mesh->Lines[idx_line]->Num_Volus_Owned=-1;
	       mesh->Lines[idx_line]->Volus_Owned=NULL;
               mesh->Lines[idx_line]->Mark=-1;
	       mesh->Lines[idx_line]->ID_Boundary=mesh->Faces[volu->Faces[2]]->ID_Boundary;
	       mesh->Lines[idx_line++]->Verts[0]=v3;
	       	       
	       //更新第二个面上的新面
               mesh->Faces[volu->Faces[2]]->Verts[0]=v0;
               mesh->Faces[volu->Faces[2]]->Verts[1]=mesh->Lines[l0]->Child[2];
               mesh->Faces[volu->Faces[2]]->Verts[2]=v3;
              
	       mesh->Faces[volu->Faces[2]]->Lines[0]=idx_line-1;
               mesh->Faces[volu->Faces[2]]->Lines[1]=l2;
               mesh->Faces[volu->Faces[2]]->Lines[2]=
	       volu->Verts[0]==mesh->Lines[volu->Lines[0]]->Verts[0]?
	       mesh->Lines[volu->Lines[0]]->Child[0]:mesh->Lines[l0]->Child[1];

               mesh->Faces[idx_face]->Verts[0]=v1;
               mesh->Faces[idx_face]->Verts[1]=v3;
               mesh->Faces[idx_face]->Verts[2]=mesh->Lines[l0]->Child[2];

	       mesh->Faces[idx_face]->Old_Child[0]=-1;
               mesh->Faces[idx_face]->Old_Child[1]=-1;

               mesh->Faces[idx_face]->Child[0]=-1;
               mesh->Faces[idx_face]->Child[1]=-1;
               mesh->Faces[idx_face]->Child[2]=-1;

               mesh->Faces[idx_face]->ID_Boundary=mesh->Faces[volu->Faces[2]]->ID_Boundary;

               mesh->Faces[idx_face]->Lines[0]=idx_line-1;
               mesh->Faces[idx_face]->Lines[1]=
	       v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
               mesh->Faces[idx_face++]->Lines[2]=l4;

	       mesh->Faces[f2]->Old_Child[0]=mesh->Faces[f2]->Child[1];
	       mesh->Faces[f2]->Old_Child[1]=mesh->Faces[f2]->Child[2];

	       mesh->Faces[f2]->Child[0]=f2;
	       mesh->Faces[f2]->Child[1]=idx_face-1;
	       mesh->Faces[f2]->Child[2]=idx_line-1;
	      
	       //更新中间新面的信息
               mesh->Faces[idx_face]->Verts[0]=mesh->Lines[l0]->Child[2];
               mesh->Faces[idx_face]->Verts[1]=v2;
               mesh->Faces[idx_face]->Verts[2]=v3;
              
	       mesh->Faces[idx_face]->Old_Child[0]=-1;
               mesh->Faces[idx_face]->Old_Child[1]=-1;

               mesh->Faces[idx_face]->Child[0]=-1;
               mesh->Faces[idx_face]->Child[1]=-1;
               mesh->Faces[idx_face]->Child[2]=-1;

               mesh->Faces[idx_face]->ID_Boundary=0;

	       mesh->Faces[idx_face]->Lines[0]=l5;
               mesh->Faces[idx_face]->Lines[1]=idx_line-1;
               mesh->Faces[idx_face++]->Lines[2]=mesh->Faces[f3]->Child[2];

	       //更新两个单元的信息
               mesh->Volus[idx_volu]->Verts[0]=v1;
               mesh->Volus[idx_volu]->Verts[1]=v2;
               mesh->Volus[idx_volu]->Verts[2]=v3;
               mesh->Volus[idx_volu]->Verts[3]=mesh->Lines[l0]->Child[2];

               mesh->Volus[idx_volu]->Lines[0]=l3;
               mesh->Volus[idx_volu]->Lines[1]=l4;
               mesh->Volus[idx_volu]->Lines[2]=
	       v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
               mesh->Volus[idx_volu]->Lines[3]=l5;
               mesh->Volus[idx_volu]->Lines[5]=idx_line-1;
	       /*  if(mesh->Faces[volu->Faces[3]]->Old_Child[1]<0)
               mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f3]->Child[2];
               else
               mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f3]->Old_Child[1];
               */
               mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f3]->Child[2];
	       if(i==0)
   	       {
                mesh->Volus[idx_volu]->Father=j;
                mesh->Volus[j]->Father=j;
               }
               else
               mesh->Volus[idx_volu]->Father=mesh->Volus[j]->Father;
               mesh->Volus[idx_volu]->Ancestor=mesh->Volus[j]->Ancestor;

               mesh->Volus[idx_volu]->Type=2;

               mesh->Volus[idx_volu]->Faces[0]=idx_face-1;
	       /*
	       if(mesh->Faces[f3]->Old_Child[0]<0)
               mesh->Volus[idx_volu]->Faces[2]=
               v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	       mesh->Faces[f3]->Child[1]:mesh->Faces[f3]->Child[0];
	       else
	       mesh->Volus[idx_volu]->Faces[2]=
               v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	       mesh->Faces[f3]->Old_Child[0]:mesh->Faces[f3]->Child[0];
               */
               mesh->Volus[idx_volu]->Faces[2]=
               v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	       mesh->Faces[f3]->Child[1]:mesh->Faces[f3]->Child[0];
               mesh->Volus[idx_volu]->Faces[1]=idx_face-2;
               mesh->Volus[idx_volu++]->Faces[3]=volu->Faces[0];

	       mesh->Volus[j]->Verts[0]=v0;
               mesh->Volus[j]->Verts[1]=v2;
               mesh->Volus[j]->Verts[2]=v3;
               mesh->Volus[j]->Verts[3]=mesh->Lines[l0]->Child[2];

	       tmp_line=v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];
               mesh->Volus[j]->Lines[0]=l1;
               mesh->Volus[j]->Lines[1]=l2;
               mesh->Volus[j]->Lines[2]=tmp_line;
               mesh->Volus[j]->Lines[3]=l5;
	       /*if(mesh->Faces[f3]->Old_Child[1]<0)
               mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Child[2];
               else
               mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Old_Child[1];
               */
               mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Child[2];
               mesh->Volus[j]->Lines[5]=idx_line-1;

	       tmp_face=volu->Faces[1];
               mesh->Volus[j]->Faces[0]=idx_face-1;
               mesh->Volus[j]->Faces[1]=f2;
	       /*
	       if(mesh->Faces[f3]->Old_Child[1]<0)
               mesh->Volus[j]->Faces[2]=
	       v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	       mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Child[1];
	       else
               mesh->Volus[j]->Faces[2]=
	       v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	       mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Old_Child[0];
               */
               mesh->Volus[j]->Faces[2]=
	       v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	       mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Child[1];
	       mesh->Volus[j]->Faces[3]=tmp_face;

               mesh->Volus[j]->Type=2;

	       mesh->Volus[j]->Child[0]=j;
	       mesh->Volus[j]->Child[1]=idx_volu-1;


	       //ReConstruct(mesh->Volus[idx_volu-1]);
               //ReConstruct(mesh->Volus[idx_volu-2]);
           
	       if(MarkJudge(mesh,mesh->Volus[j])<0)
	       {
	        mesh->Volus[j]->Mark=-1;
	       }
               if(MarkJudge(mesh,mesh->Volus[idx_volu-1])>=0) 
	       {
                mesh->Volus[idx_volu-1]->Mark=1;
	       }
               else
	        mesh->Volus[idx_volu-1]->Mark=-1;
            }//4
       }//type1
      else 
       {//type2
          if (  (mesh->Faces[volu->Faces[3]]->Child[0]<0&&mesh->Faces[volu->Faces[2]]->Child[0]<0)
	      ||(mesh->Faces[volu->Faces[3]]->Child[0]<0&&mesh->Faces[volu->Faces[2]]->Child[0]>=0
		  &&(mesh->Faces[volu->Faces[2]]->Child[2]==volu->Lines[0]||mesh->Faces[volu->Faces[2]]
	         ->Child[2]==volu->Lines[2]||mesh->Faces[volu->Faces[2]]->Child[2]==volu->Lines[4]))
	      ||(mesh->Faces[volu->Faces[2]]->Child[0]<0&&mesh->Faces[volu->Faces[3]]->Child[0]>=0
	         &&(mesh->Faces[volu->Faces[3]]->Child[2]==volu->Lines[0]||mesh->Faces[volu->Faces[3]]
		    ->Child[2]==volu->Lines[1]||mesh->Faces[volu->Faces[3]]->Child[2]==volu->Lines[3]) )
	      ||(mesh->Faces[volu->Faces[3]]->Child[0]>=0&&(mesh->Faces[volu->Faces[3]]->Child[2]==
		    volu->Lines[0]||mesh->Faces[volu->Faces[3]]
		    ->Child[2]==volu->Lines[1]||mesh->Faces[volu->Faces[3]]->Child[2]==volu->Lines[3]) 
		 &&mesh->Faces[volu->Faces[2]]->Child[0]>=0&&(mesh->Faces[volu->Faces[2]]->Child[2]==
		    volu->Lines[0]||mesh->Faces[volu->Faces[2]]
	         ->Child[2]==volu->Lines[2]||mesh->Faces[volu->Faces[2]]->Child[2]==volu->Lines[4])
		 )
	     )   

          {//1
	   //更新面上两条新边的信息 
	   mesh->Lines[idx_line]->Verts[1]= mesh->Lines[l0]->Child[2];
           mesh->Lines[idx_line]->Mark=-1;
	   mesh->Lines[idx_line]->Num_Volus_Owned=-1;
	   mesh->Lines[idx_line]->Volus_Owned=NULL;
           mesh->Lines[idx_line]->Child[0]=-1;
	   mesh->Lines[idx_line]->Child[1]=-1;
	   mesh->Lines[idx_line]->Child[2]=-1;
	   mesh->Lines[idx_line]->ID_Boundary=mesh->Faces[f3]->ID_Boundary;
	   mesh->Lines[idx_line++]->Verts[0]=v2;

	   mesh->Lines[idx_line]->Verts[1]= mesh->Lines[l0]->Child[2];
           mesh->Lines[idx_line]->Mark=-1;
	   mesh->Lines[idx_line]->Num_Volus_Owned=-1;
	   mesh->Lines[idx_line]->Volus_Owned=NULL;
	   mesh->Lines[idx_line]->Child[0]=-1;
	   mesh->Lines[idx_line]->Child[1]=-1;
	   mesh->Lines[idx_line]->Child[2]=-1;
           mesh->Lines[idx_line]->ID_Boundary=mesh->Faces[f2]->ID_Boundary;
	   mesh->Lines[idx_line++]->Verts[0]=v3;
           //更新第一个面上的信息
           mesh->Faces[idx_face]->Verts[0]=v1;
           mesh->Faces[idx_face]->Verts[1]=v2;
           mesh->Faces[idx_face]->Verts[2]=mesh->Lines[l0]->Child[2];

	   mesh->Faces[idx_face]->Old_Child[0]=-1;
           mesh->Faces[idx_face]->Old_Child[1]=-1;

	   mesh->Faces[idx_face]->Child[0]=-1;
           mesh->Faces[idx_face]->Child[1]=-1;
           mesh->Faces[idx_face]->Child[2]=-1;
           mesh->Faces[idx_face]->ID_Boundary=mesh->Faces[f3]->ID_Boundary;

           mesh->Faces[idx_face]->Lines[0]=idx_line-2;
           mesh->Faces[idx_face]->Lines[1]=
	   v0==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
           mesh->Faces[idx_face++]->Lines[2]=l3;

	   mesh->Faces[f3]->Verts[0]=v0;
           mesh->Faces[f3]->Verts[1]=mesh->Lines[l0]->Child[2];
           mesh->Faces[f3]->Verts[2]=v2;
         
	   mesh->Faces[f3]->Lines[0]=idx_line-2;
           mesh->Faces[f3]->Lines[1]=l1;
           mesh->Faces[f3]->Lines[2]=
           v0==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];
           
	   mesh->Faces[f3]->Old_Child[0]=mesh->Faces[f3]->Child[1];
	   mesh->Faces[f3]->Old_Child[1]=mesh->Faces[f3]->Child[2];
	   mesh->Faces[f3]->Child[0]=f3;
           mesh->Faces[f3]->Child[1]=idx_face-1;
           mesh->Faces[f3]->Child[2]=idx_line-2;

           //更新第二个面上的信息
           mesh->Faces[idx_face]->Verts[0]=v1;
           mesh->Faces[idx_face]->Verts[1]=v3;
           mesh->Faces[idx_face]->Verts[2]=mesh->Lines[l0]->Child[2];

           mesh->Faces[idx_face]->Old_Child[0]=-1;
           mesh->Faces[idx_face]->Old_Child[1]=-1;

	   mesh->Faces[idx_face]->Child[0]=-1;
           mesh->Faces[idx_face]->Child[1]=-1;
           mesh->Faces[idx_face]->Child[2]=-1;
           mesh->Faces[idx_face]->ID_Boundary=mesh->Faces[f2]->ID_Boundary;

           mesh->Faces[idx_face]->Lines[0]=idx_line-1;
           mesh->Faces[idx_face]->Lines[1]=
           v0==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
           mesh->Faces[idx_face++]->Lines[2]=l4;

	   mesh->Faces[f2]->Verts[0]=v0;
           mesh->Faces[f2]->Verts[1]=mesh->Lines[l0]->Child[2];
           mesh->Faces[f2]->Verts[2]=v3;
         

	   mesh->Faces[f2]->Lines[0]=idx_line-1;
           mesh->Faces[f2]->Lines[1]=l2;
           mesh->Faces[f2]->Lines[2]=
           v0==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];

	   mesh->Faces[f2]->Old_Child[0]=mesh->Faces[f2]->Child[1];
           mesh->Faces[f2]->Old_Child[1]=mesh->Faces[f2]->Child[2];

	   mesh->Faces[f2]->Child[0]=f2;
           mesh->Faces[f2]->Child[1]=idx_face-1;
           mesh->Faces[f2]->Child[2]=idx_line-1;


           //更新中间的新面的信息
	   mesh->Faces[idx_face]->Verts[0]=mesh->Lines[l0]->Child[2];
           mesh->Faces[idx_face]->Verts[1]=v2;
           mesh->Faces[idx_face]->Verts[2]=v3;
         
	   mesh->Faces[idx_face]->Old_Child[0]=-1;
           mesh->Faces[idx_face]->Old_Child[1]=-1;

	   mesh->Faces[idx_face]->Child[0]=-1;
           mesh->Faces[idx_face]->Child[1]=-1;
           mesh->Faces[idx_face]->Child[2]=-1;
         
	   mesh->Faces[idx_face]->ID_Boundary=0;

	   mesh->Faces[idx_face]->Lines[0]=l5;
           mesh->Faces[idx_face]->Lines[1]=idx_line-1;
           mesh->Faces[idx_face++]->Lines[2]=idx_line-2;


	      //更新单元信息
	   mesh->Volus[idx_volu]->Verts[0]=v1;
	   mesh->Volus[idx_volu]->Verts[1]=v2;
	   mesh->Volus[idx_volu]->Verts[2]=v3;
           mesh->Volus[idx_volu]->Verts[3]=mesh->Lines[l0]->Child[2];
           
	   mesh->Volus[idx_volu]->Lines[0]=l3;
	   mesh->Volus[idx_volu]->Lines[1]=l4;
	   mesh->Volus[idx_volu]->Lines[2]=
	   v0==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
 	   mesh->Volus[idx_volu]->Lines[3]=l5;
           mesh->Volus[idx_volu]->Lines[4]=idx_line-2;
 	   mesh->Volus[idx_volu]->Lines[5]=idx_line-1;

	   if(i==0)
   	       {
                mesh->Volus[idx_volu]->Father=j;
                mesh->Volus[j]->Father=j;
               }
           else
           mesh->Volus[idx_volu]->Father=mesh->Volus[j]->Father;
           mesh->Volus[idx_volu]->Ancestor=mesh->Volus[j]->Ancestor;


           mesh->Volus[idx_volu]->Type=0;

	   mesh->Volus[idx_volu]->Faces[0]=idx_face-1;
	   mesh->Volus[idx_volu]->Faces[1]=idx_face-2;
	   mesh->Volus[idx_volu]->Faces[2]=idx_face-3;
           mesh->Volus[idx_volu++]->Faces[3]=volu->Faces[0];

	   mesh->Volus[j]->Verts[0]=v0;
	   mesh->Volus[j]->Verts[1]=v2;
	   mesh->Volus[j]->Verts[2]=v3;
	   mesh->Volus[j]->Verts[3]=mesh->Lines[l0]->Child[2];

	   tmp_line=v0==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];
           mesh->Volus[j]->Lines[0]=l1;
	   mesh->Volus[j]->Lines[1]=l2;
	   mesh->Volus[j]->Lines[2]=tmp_line;
 	   mesh->Volus[j]->Lines[3]=l5;
           mesh->Volus[j]->Lines[4]=idx_line-2;
 	   mesh->Volus[j]->Lines[5]=idx_line-1;

	   tmp_face=f1;
	   mesh->Volus[j]->Faces[0]=idx_face-1;
	   mesh->Volus[j]->Faces[1]=f2;
	   mesh->Volus[j]->Faces[2]=f3;
	   mesh->Volus[j]->Faces[3]=tmp_face;

           mesh->Volus[j]->Type=0;

	   mesh->Volus[j]->Child[0]=j;
	   mesh->Volus[j]->Child[1]=idx_volu-1;

           //ReConstruct(mesh->Volus[idx_volu-1]);
           //ReConstruct(mesh->Volus[MarkSet[i]]);
           if(MarkJudge(mesh,mesh->Volus[j])<0)
	   {
            mesh->Volus[j]->Mark=-1;
	   }
           if(MarkJudge(mesh,mesh->Volus[idx_volu-1])>=0) 
	   {
            mesh->Volus[idx_volu-1]->Mark=1;
	   }
           else
	   mesh->Volus[idx_volu-1]->Mark=-1;
         }//1
           //-----------------------------------------------------------------------------------
           //第二种情况：第一个面没加密，第二个面加密 
	 else if(  (mesh->Faces[volu->Faces[3]]->Child[0]<0
		     &&mesh->Faces[volu->Faces[2]]->Child[0]>=0
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[0]
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[2]
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[4])
		  ||(mesh->Faces[volu->Faces[3]]->Child[0]>=0
		     &&(mesh->Faces[volu->Faces[3]]->Child[2]==
		    volu->Lines[0]||mesh->Faces[volu->Faces[3]]
		    ->Child[2]==volu->Lines[1]||mesh->Faces[volu->Faces[3]]->Child[2]==volu->Lines[3]) 
		     &&mesh->Faces[volu->Faces[2]]->Child[0]>=0
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[0]
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[2]
		     &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[4]
		    )
	       )
	   {//2
	      //更新第一个面上的边	   
	       mesh->Lines[idx_line]->Verts[1]=mesh->Lines[l0]->Child[2];
	       mesh->Lines[idx_line]->Child[0]=-1;
	       mesh->Lines[idx_line]->Child[1]=-1;
               mesh->Lines[idx_line]->Child[2]=-1;
	       mesh->Lines[idx_line]->Num_Volus_Owned=-1;
	       mesh->Lines[idx_line]->Volus_Owned=NULL;
               mesh->Lines[idx_line]->Mark=-1;
	       mesh->Lines[idx_line]->ID_Boundary=mesh->Faces[f3]->ID_Boundary;
	       mesh->Lines[idx_line++]->Verts[0]=v2;
	       	       
	       //更新第一个面上的新面
               mesh->Faces[f3]->Verts[0]=v0;
               mesh->Faces[f3]->Verts[1]=mesh->Lines[l0]->Child[2];
               mesh->Faces[f3]->Verts[2]=v2;
              
	       mesh->Faces[f3]->Lines[0]=idx_line-1;
               mesh->Faces[f3]->Lines[1]=l1;
               mesh->Faces[f3]->Lines[2]=
	       volu->Verts[0]==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];

               mesh->Faces[idx_face]->Verts[0]=v1;
               mesh->Faces[idx_face]->Verts[1]=v2;
               mesh->Faces[idx_face]->Verts[2]=mesh->Lines[l0]->Child[2];

	       mesh->Faces[idx_face]->Old_Child[0]=-1;
               mesh->Faces[idx_face]->Old_Child[1]=-1;
 
               mesh->Faces[idx_face]->Child[0]=-1;
               mesh->Faces[idx_face]->Child[1]=-1;
               mesh->Faces[idx_face]->Child[2]=-1;

               mesh->Faces[idx_face]->ID_Boundary=mesh->Faces[f3]->ID_Boundary;

               mesh->Faces[idx_face]->Lines[0]=idx_line-1;
               mesh->Faces[idx_face]->Lines[1]=
	       v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
               mesh->Faces[idx_face++]->Lines[2]=l3;

	       mesh->Faces[f3]->Old_Child[0]=mesh->Faces[f3]->Child[1];
	       mesh->Faces[f3]->Old_Child[1]=mesh->Faces[f3]->Child[2];
	       mesh->Faces[f3]->Child[0]=volu->Faces[3];
	       mesh->Faces[f3]->Child[1]=idx_face-1;
	       mesh->Faces[f3]->Child[2]=idx_line-1;
	      
	       //更新中间新面的信息
               mesh->Faces[idx_face]->Verts[0]=mesh->Lines[l0]->Child[2];
               mesh->Faces[idx_face]->Verts[1]=v2;
               mesh->Faces[idx_face]->Verts[2]=v3;
             
	       mesh->Faces[idx_face]->Old_Child[0]=-1;
               mesh->Faces[idx_face]->Old_Child[1]=-1;
 
               mesh->Faces[idx_face]->Child[0]=-1;
               mesh->Faces[idx_face]->Child[1]=-1;
               mesh->Faces[idx_face]->Child[2]=-1;

               mesh->Faces[idx_face]->ID_Boundary=0;

	       mesh->Faces[idx_face]->Lines[0]=l5;
               mesh->Faces[idx_face]->Lines[1]=mesh->Faces[f2]->Child[2];
               mesh->Faces[idx_face++]->Lines[2]=idx_line-1;
          
	       //更新两个单元的信息
               mesh->Volus[idx_volu]->Verts[0]=v1;
               mesh->Volus[idx_volu]->Verts[1]=v2;
               mesh->Volus[idx_volu]->Verts[2]=v3;
               mesh->Volus[idx_volu]->Verts[3]=mesh->Lines[l0]->Child[2];

               mesh->Volus[idx_volu]->Lines[0]=l3;
               mesh->Volus[idx_volu]->Lines[1]=l4;
               mesh->Volus[idx_volu]->Lines[2]=
	       v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
               mesh->Volus[idx_volu]->Lines[3]=l5;
	       /*if(mesh->Faces[f2]->Old_Child[1]<0)
               mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f2]->Child[2];
	       else
               mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f2]->Old_Child[1];
	       */
	       mesh->Volus[idx_volu]->Lines[4]=idx_line-1;
               mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f2]->Child[2];
	       if(i==0)
   	       {
                mesh->Volus[idx_volu]->Father=j;
                mesh->Volus[j]->Father=j;
               }
               else
               mesh->Volus[idx_volu]->Father=mesh->Volus[j]->Father;
               mesh->Volus[idx_volu]->Ancestor=mesh->Volus[j]->Ancestor;

               mesh->Volus[idx_volu]->Type=0;

               mesh->Volus[idx_volu]->Faces[0]=idx_face-1;
               mesh->Volus[idx_volu]->Faces[2]=idx_face-2;
	       /*
	       if(mesh->Faces[volu->Faces[2]]->Old_Child[0]<0)
               mesh->Volus[idx_volu]->Faces[1]=
               v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	       mesh->Faces[f2]->Child[1]:mesh->Faces[f2]->Child[0];
	       else
	        mesh->Volus[idx_volu]->Faces[1]=
               v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	       mesh->Faces[f2]->Old_Child[0]:mesh->Faces[f2]->Child[0];
               */
               mesh->Volus[idx_volu]->Faces[1]=
               v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	       mesh->Faces[f2]->Child[1]:mesh->Faces[f2]->Child[0];

	       mesh->Volus[idx_volu++]->Faces[3]=volu->Faces[0];

	       mesh->Volus[j]->Verts[0]=v0;
               mesh->Volus[j]->Verts[1]=v2;
               mesh->Volus[j]->Verts[2]=v3;
               mesh->Volus[j]->Verts[3]=mesh->Lines[l0]->Child[2];

	       tmp_line=v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];
               mesh->Volus[j]->Lines[0]=l1;
               mesh->Volus[j]->Lines[1]=l2;
               mesh->Volus[j]->Lines[2]=tmp_line;
               mesh->Volus[j]->Lines[3]=l5;
               mesh->Volus[j]->Lines[4]=idx_line-1;
	       /*if(mesh->Faces[f2]->Old_Child[1]<0)
               mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Child[2];
	       else
                mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Old_Child[1];
               */
               mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Child[2];
	       tmp_face=volu->Faces[1];
               mesh->Volus[j]->Faces[0]=idx_face-1;
	       /*
	       if(mesh->Faces[f2]->Old_Child[1]<0)
               mesh->Volus[j]->Faces[1]=
	       v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	       mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Child[1];
	       else
	       mesh->Volus[j]->Faces[1]=
	       v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	       mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Old_Child[0];
               */
	       mesh->Volus[j]->Faces[1]=
	       v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	       mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Child[1];

	       mesh->Volus[j]->Faces[2]=f3;
               mesh->Volus[j]->Faces[3]=tmp_face;

               mesh->Volus[j]->Type=0;

               mesh->Volus[j]->Child[0]=j;
	       mesh->Volus[j]->Child[1]=idx_volu-1;


	       //ReConstruct(mesh->Volus[idx_volu-1]);
               //ReConstruct(mesh->Volus[idx_volu-2]);
              if(MarkJudge(mesh,mesh->Volus[j])<0)
	      {
	       mesh->Volus[j]->Mark=-1; 
	      }
              if(MarkJudge(mesh,mesh->Volus[idx_volu-1])>=0) 
	      {
               mesh->Volus[idx_volu-1]->Mark=1;
	      }
              else
	       mesh->Volus[idx_volu-1]->Mark=-1;
         }//2
         //---------------------------------------------------------------------------------------------
           else if ( (mesh->Faces[volu->Faces[3]]->Child[0]>=0
		     &&mesh->Faces[volu->Faces[3]]->Child[2]!=volu->Lines[0]
		     &&mesh->Faces[volu->Faces[3]]->Child[2]!=volu->Lines[1]
		     &&mesh->Faces[volu->Faces[3]]->Child[2]!=volu->Lines[3]
		     )
		  &&(mesh->Faces[volu->Faces[2]]->Child[0]>=0
		    &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[0]
		    &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[2]
		    &&mesh->Faces[volu->Faces[2]]->Child[2]!=volu->Lines[4]
		    )
		  ) 
	   {//3
	    //更新中间的新面的信息
	   mesh->Faces[idx_face]->Verts[0]=mesh->Lines[l0]->Child[2];
           mesh->Faces[idx_face]->Verts[1]=v2;
           mesh->Faces[idx_face]->Verts[2]=v3;
        
	   mesh->Faces[idx_face]->Old_Child[0]=-1;
           mesh->Faces[idx_face]->Old_Child[1]=-1;
 
	   mesh->Faces[idx_face]->Child[0]=-1;
           mesh->Faces[idx_face]->Child[1]=-1;
           mesh->Faces[idx_face]->Child[2]=-1;
 
           mesh->Faces[idx_face]->ID_Boundary=0;

	   mesh->Faces[idx_face]->Lines[0]=l5;
           mesh->Faces[idx_face]->Lines[1]=mesh->Faces[f2]->Child[2];
           mesh->Faces[idx_face++]->Lines[2]=mesh->Faces[f3]->Child[2];

           //更新单元信息
           mesh->Volus[idx_volu]->Verts[0]=v1;
           mesh->Volus[idx_volu]->Verts[1]=v2;
           mesh->Volus[idx_volu]->Verts[2]=v3;
           mesh->Volus[idx_volu]->Verts[3]=mesh->Lines[l0]->Child[2];

 	   
	   mesh->Volus[idx_volu]->Lines[0]=l3;
           mesh->Volus[idx_volu]->Lines[1]=l4;
           mesh->Volus[idx_volu]->Lines[2]=
           v0==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
	   mesh->Volus[idx_volu]->Lines[3]=l5;
	   /*if(mesh->Faces[volu->Faces[2]]->Old_Child[1]<0)
	   mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f2]->Child[2];
	   else
           mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f2]->Old_Child[1];
           if(mesh->Faces[f3]->Old_Child[1]<0)
            mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f3]->Child[2];
	   else
	   mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f3]->Old_Child[1];
           */

	   mesh->Volus[idx_volu]->Lines[5]=mesh->Faces[f2]->Child[2];
           mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f3]->Child[2];

	   if(i==0)
   	       {
                mesh->Volus[idx_volu]->Father=j;
                mesh->Volus[j]->Father=j;
               }
           else
           mesh->Volus[idx_volu]->Father=mesh->Volus[j]->Father;
           mesh->Volus[idx_volu]->Ancestor=mesh->Volus[j]->Ancestor;

           mesh->Volus[idx_volu]->Type=0;

	   mesh->Volus[idx_volu]->Faces[0]=idx_face-1;
	   /*
	   if(mesh->Faces[f3]->Old_Child[0]<0)
	   mesh->Volus[idx_volu]->Faces[2]=
           v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	   mesh->Faces[f3]->Child[1]:mesh->Faces[f3]->Child[0];
	   else
	   mesh->Volus[idx_volu]->Faces[2]=
           v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	   mesh->Faces[f3]->Old_Child[0]:mesh->Faces[f3]->Child[0];
	   if(mesh->Faces[f2]->Old_Child[0]<0)
	   mesh->Volus[idx_volu]->Faces[1]=
           v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	   mesh->Faces[f2]->Child[1]:mesh->Faces[f2]->Child[0];
	   else
           mesh->Volus[idx_volu]->Faces[1]=
           v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	   mesh->Faces[f2]->Old_Child[0]:mesh->Faces[f2]->Child[0];
           */
           mesh->Volus[idx_volu]->Faces[2]=
           v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	   mesh->Faces[f3]->Child[1]:mesh->Faces[f3]->Child[0];
	   mesh->Volus[idx_volu]->Faces[1]=
           v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	   mesh->Faces[f2]->Child[1]:mesh->Faces[f2]->Child[0];

	   mesh->Volus[idx_volu++]->Faces[3]=f0;

           mesh->Volus[j]->Verts[0]=v0;
	   mesh->Volus[j]->Verts[1]=v2;
	   mesh->Volus[j]->Verts[2]=v3;
           mesh->Volus[j]->Verts[3]=mesh->Lines[l0]->Child[2];
	   //printf ( "test identity-=%d,%d,%d,%d\n", volu->Verts[0],volu->Verts[2],volu->Verts[3],mesh->Lines[volu->Lines[0]]->Child[2]);
	  
	   //printf ( "test identity-=%d,%d,%d,%d\n", mesh->Volus[j]->Verts[0],mesh->Volus[j]->Verts[1],mesh->Volus[j]->Verts[2],mesh->Volus[j]->Verts[3]);
	   tmp_line=v0==mesh->Lines[l0]->Verts[0]?
	   mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];
	   mesh->Volus[j]->Lines[0]=l1;
	   mesh->Volus[j]->Lines[1]=l2;
           mesh->Volus[j]->Lines[2]=tmp_line;
	   mesh->Volus[j]->Lines[3]=l5;
           mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Child[2];
           mesh->Volus[j]->Lines[5]=mesh->Faces[f2]->Child[2];

	   tmp_face=volu->Faces[1];
	   mesh->Volus[j]->Faces[0]=idx_face-1;
	   /*
	   if(mesh->Faces[f2]->Old_Child[1]<0)
	   mesh->Volus[j]->Faces[1]=
           v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	   mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Child[1];
	   else
           mesh->Volus[j]->Faces[1]=
           v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	   mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Old_Child[0];
           if(mesh->Faces[f3]->Old_Child[1]<0)
           mesh->Volus[j]->Faces[2]=
           v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	   mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Child[1];
	   else
	   mesh->Volus[j]->Faces[2]=
           v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	   mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Old_Child[0];
           */
           mesh->Volus[j]->Faces[1]=
           v0==mesh->Faces[mesh->Faces[f2]->Child[0]]->Verts[0]?
	   mesh->Faces[f2]->Child[0]:mesh->Faces[f2]->Child[1];
           mesh->Volus[j]->Faces[2]=
           v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	   mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Child[1];

	   mesh->Volus[j]->Faces[3]=tmp_face;

           mesh->Volus[j]->Type=0;

	   mesh->Volus[j]->Child[0]=j;
	   mesh->Volus[j]->Child[1]=idx_volu-1;


 	   if(MarkJudge(mesh,mesh->Volus[j])<0)
	   {
	    mesh->Volus[j]->Mark=-1;
	   }
	   if(MarkJudge(mesh,mesh->Volus[idx_volu-1])>=0) 
	   {
            mesh->Volus[idx_volu-1]->Mark=1;
	   }
          else
	   mesh->Volus[idx_volu-1]->Mark=-1;

	   }//3
	   //-------------------------------------------------------------------------------------------
           //第四种情况：第一个面(012)加密，第二个面没加密 
	   else
	   {//4
	      //更新第二个面上的边	   
	      mesh->Lines[idx_line]->Verts[1]=mesh->Lines[l0]->Child[2];
	       mesh->Lines[idx_line]->Child[0]=-1;
	       mesh->Lines[idx_line]->Child[1]=-1;
               mesh->Lines[idx_line]->Child[2]=-1;
	       mesh->Lines[idx_line]->Num_Volus_Owned=-1;
	       mesh->Lines[idx_line]->Volus_Owned=NULL;
               mesh->Lines[idx_line]->Mark=-1;
	       mesh->Lines[idx_line]->ID_Boundary=mesh->Faces[f2]->ID_Boundary;
	       mesh->Lines[idx_line++]->Verts[0]=v3;
	       	       
	       //更新第二个面上的新面
               mesh->Faces[f2]->Verts[0]=v0;
               mesh->Faces[f2]->Verts[1]=mesh->Lines[l0]->Child[2];
               mesh->Faces[f2]->Verts[2]=v3;
              
	       mesh->Faces[f2]->Lines[0]=idx_line-1;
               mesh->Faces[f2]->Lines[1]=volu->Lines[2];
               mesh->Faces[f2]->Lines[2]=
	       v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];

               mesh->Faces[idx_face]->Verts[0]=v1;
               mesh->Faces[idx_face]->Verts[1]=v3;
               mesh->Faces[idx_face]->Verts[2]=mesh->Lines[l0]->Child[2];

	       mesh->Faces[idx_face]->Old_Child[0]=-1;
               mesh->Faces[idx_face]->Old_Child[1]=-1;

               mesh->Faces[idx_face]->Child[0]=-1;
               mesh->Faces[idx_face]->Child[1]=-1;
               mesh->Faces[idx_face]->Child[2]=-1;

               mesh->Faces[idx_face]->ID_Boundary=mesh->Faces[volu->Faces[2]]->ID_Boundary;

               mesh->Faces[idx_face]->Lines[0]=idx_line-1;
               mesh->Faces[idx_face]->Lines[1]=
	       v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
               mesh->Faces[idx_face++]->Lines[2]=l4;

	       mesh->Faces[f2]->Old_Child[0]=mesh->Faces[f2]->Child[1];
	       mesh->Faces[f2]->Old_Child[1]=mesh->Faces[f2]->Child[2];

	       mesh->Faces[f2]->Child[0]=f2;
	       mesh->Faces[f2]->Child[1]=idx_face-1;
	       mesh->Faces[f2]->Child[2]=idx_line-1;
	      
	       //更新中间新面的信息
               mesh->Faces[idx_face]->Verts[0]=mesh->Lines[l0]->Child[2];
               mesh->Faces[idx_face]->Verts[1]=v2;
               mesh->Faces[idx_face]->Verts[2]=v3;
              
	       mesh->Faces[idx_face]->Old_Child[0]=-1;
               mesh->Faces[idx_face]->Old_Child[1]=-1;
               mesh->Faces[idx_face]->Child[0]=-1;
               mesh->Faces[idx_face]->Child[1]=-1;
               mesh->Faces[idx_face]->Child[2]=-1;

               mesh->Faces[idx_face]->ID_Boundary=0;

	       mesh->Faces[idx_face]->Lines[0]=l5;
               mesh->Faces[idx_face]->Lines[1]=idx_line-1;
               mesh->Faces[idx_face++]->Lines[2]=mesh->Faces[f3]->Child[2];

	       //更新两个单元的信息
               mesh->Volus[idx_volu]->Verts[0]=v1;
               mesh->Volus[idx_volu]->Verts[1]=v2;
               mesh->Volus[idx_volu]->Verts[2]=v3;
               mesh->Volus[idx_volu]->Verts[3]=mesh->Lines[l0]->Child[2];

               mesh->Volus[idx_volu]->Lines[0]=l3;
               mesh->Volus[idx_volu]->Lines[1]=l4;
               mesh->Volus[idx_volu]->Lines[2]=
	       v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[1]:mesh->Lines[l0]->Child[0];
               mesh->Volus[idx_volu]->Lines[3]=l5;
               mesh->Volus[idx_volu]->Lines[5]=idx_line-1;
	       /* if(mesh->Faces[f3]->Old_Child[1]<0)
               mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f3]->Child[2];
               else
	       mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f3]->Old_Child[1];
               */
               mesh->Volus[idx_volu]->Lines[4]=mesh->Faces[f3]->Child[2];
	       if(i==0)
   	       {
                mesh->Volus[idx_volu]->Father=j;
                mesh->Volus[j]->Father=j;
               }
               else
               mesh->Volus[idx_volu]->Father=mesh->Volus[j]->Father;
               mesh->Volus[idx_volu]->Ancestor=mesh->Volus[j]->Ancestor;

               mesh->Volus[idx_volu]->Type=0;

               mesh->Volus[idx_volu]->Faces[0]=idx_face-1;
	       /*
	       if(mesh->Faces[f3]->Old_Child[0]<0)
               mesh->Volus[idx_volu]->Faces[2]=
               v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	       mesh->Faces[f3]->Child[1]:mesh->Faces[f3]->Child[0];
               else
	       mesh->Volus[idx_volu]->Faces[2]=
               v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	       mesh->Faces[f3]->Old_Child[0]:mesh->Faces[f3]->Child[0];
               */
               mesh->Volus[idx_volu]->Faces[2]=
               v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	       mesh->Faces[f3]->Child[1]:mesh->Faces[f3]->Child[0];
	       mesh->Volus[idx_volu]->Faces[1]=idx_face-2;
               mesh->Volus[idx_volu++]->Faces[3]=f0;


	       mesh->Volus[j]->Verts[0]=v0;
               mesh->Volus[j]->Verts[1]=v2;
               mesh->Volus[j]->Verts[2]=v3;
               mesh->Volus[j]->Verts[3]=mesh->Lines[l0]->Child[2];

               tmp_line= v0==mesh->Lines[l0]->Verts[0]?
	       mesh->Lines[l0]->Child[0]:mesh->Lines[l0]->Child[1];
	       mesh->Volus[j]->Lines[0]=l1;
               mesh->Volus[j]->Lines[1]=l2;
               mesh->Volus[j]->Lines[2]=tmp_line;
               mesh->Volus[j]->Lines[3]=l5;
	       /*if(mesh->Faces[f3]->Old_Child[1]<0)
               mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Child[2];
	       else
               mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Old_Child[1];
               */
               mesh->Volus[j]->Lines[4]=mesh->Faces[f3]->Child[2];
               mesh->Volus[j]->Lines[5]=idx_line-1;

	       tmp_face=volu->Faces[1];
               mesh->Volus[j]->Faces[0]=idx_face-1;
               mesh->Volus[j]->Faces[1]=f2;
               /*  
	       if(mesh->Faces[f3]->Old_Child[1]<0)
               mesh->Volus[j]->Faces[2]=
	       v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	       mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Child[1];
	       else
	       mesh->Volus[j]->Faces[2]=
	       v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	       mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Old_Child[0];
               */
               mesh->Volus[j]->Faces[2]=
	       v0==mesh->Faces[mesh->Faces[f3]->Child[0]]->Verts[0]?
	       mesh->Faces[f3]->Child[0]:mesh->Faces[f3]->Child[1];

               mesh->Volus[j]->Faces[3]=tmp_face;

               mesh->Volus[j]->Type=0;


	       mesh->Volus[j]->Child[0]=j;
	       mesh->Volus[j]->Child[1]=idx_volu-1;


	       //ReConstruct(mesh->Volus[idx_volu-1]);
               //ReConstruct(mesh->Volus[idx_volu-2]);
           
	      if(MarkJudge(mesh,mesh->Volus[j])<0)
	      {
	       mesh->Volus[j]->Mark=-1;
	      }
              if(MarkJudge(mesh,mesh->Volus[idx_volu-1])>=0) 
	      {
               mesh->Volus[idx_volu-1]->Mark=1;
	      }
             else
	      mesh->Volus[idx_volu-1]->Mark=-1;
            }//4
         }//type2
      }//单元遍历
   Iter_volu=idx_volu;
  }//加密三次循环
  mesh->Num_Verts_Global=total_vert;
  mesh->Num_Lines_Global=total_line;
  mesh->Num_Faces_Global=total_face;
  mesh->Num_Volus_Global=total_volu;
}//MeshRefineBisection结束
 
   
int
MeshMarkElemNum(MESH *mesh, double *eta, double *theta)
{
  int i,j;
  int *ix=(int *)malloc(mesh->Num_Volus_Global);
  double curr_error=0.0,total_error=0.0;
  for(i=0;i<mesh->Num_Volus_Global;++i)
    ix[i]=i;
  QuickSort(eta,0,mesh->Num_Volus_Global-1,sizeof(double),ix,cmp);//按照eta的大小顺序调整ix的顺序
  for(i=0;i<mesh->Num_Volus_Global;++i)
    total_error+=eta[i];
  for(i=0;i<mesh->Num_Volus_Global;++i)
  {
   curr_error+=eta[ix[i]];
   if(curr_error>=*theta*total_error)
     break;
  }
  i++;
  free(ix);
  return i;
}

void
MeshMarkElem(MESH *mesh, int *M, double *eta, double *theta,int *num)
{
  int i,j;
  int *ix=(int *)malloc(mesh->Num_Volus_Global);
  double curr_error=0.0,total_error=0.0;
  for(i=0;i<mesh->Num_Volus_Global;++i)
   ix[i]=i;
  QuickSort(eta,0,mesh->Num_Volus_Global-1,sizeof(double),ix,cmp);//按照eta的大小顺序调整ix的顺序
  for(j=0;j<*num;++j)
   M[j]=ix[j];
  free(ix);
}

int* MeshMarkRefine(MESH *mesh, double *eta, double theta,int *num)
{
  int i,j;
  int *M;
  int *ix=(int *)malloc(mesh->Num_Volus_Global);
  double curr_error=0.0,total_error=0.0;
  for(i=0;i<mesh->Num_Volus_Global;++i)
   ix[i]=i;
  QuickSort(eta,0,mesh->Num_Volus_Global-1,sizeof(double),ix,cmp);//按照eta的大小顺序调整ix的顺序
  for(i=0;i<mesh->Num_Volus_Global;++i)
   total_error+=eta[i];
  for(i=0;i<mesh->Num_Volus_Global;++i)
  {
   curr_error+=eta[ix[i]];
   if(curr_error>=theta*total_error)
    break;
  }
  num[0]=i;
  M=(int *)malloc((num[0]));
  for(j=0;j<num[0];++j)
   M[j]=ix[j];
  free(ix);
  return M;
}
//进行times步二分加密
void MeshBisection(MESH *mesh,int times)
{
   int i;
   for(i=0;i<times;++i)
   {
      MeshRefineUniformBisectionNew1(mesh);
      MeshFullInfo4Line(mesh);
      MeshRenew(mesh);
   }
}
//进行times步一致加密
void MeshConsist(MESH *mesh,int times)
{
   int i;
   for(i=0;i<times;++i)
   {
      AdjustMesh(mesh);
      MeshRefineConsistent(mesh);
   }
}




//计算单元的体积
double 
GetVolume(MESH *Mesh,VOLU *volu)
{
  double volume;
  MESH *mesh;
  mesh=Mesh;
  double ax=mesh->Verts[volu->Verts[1]]->Coord[0]-mesh->Verts[volu->Verts[0]]->Coord[0];
  double ay=mesh->Verts[volu->Verts[1]]->Coord[1]-mesh->Verts[volu->Verts[0]]->Coord[1];
  double az=mesh->Verts[volu->Verts[1]]->Coord[2]-mesh->Verts[volu->Verts[0]]->Coord[2];
  double bx=mesh->Verts[volu->Verts[2]]->Coord[0]-mesh->Verts[volu->Verts[0]]->Coord[0];
  double by=mesh->Verts[volu->Verts[2]]->Coord[1]-mesh->Verts[volu->Verts[0]]->Coord[1];
  double bz=mesh->Verts[volu->Verts[2]]->Coord[2]-mesh->Verts[volu->Verts[0]]->Coord[2];
  double cx=mesh->Verts[volu->Verts[3]]->Coord[0]-mesh->Verts[volu->Verts[0]]->Coord[0];
  double cy=mesh->Verts[volu->Verts[3]]->Coord[1]-mesh->Verts[volu->Verts[0]]->Coord[1];
  double cz=mesh->Verts[volu->Verts[3]]->Coord[2]-mesh->Verts[volu->Verts[0]]->Coord[2];
  volume=cx*(ay*bz-az*by)+cy*(az*bx-ax*bz)+cz*(ax*by-ay*bx);
  volume= volume/6;
  return fabs(volume);
}

//求jacobi变换矩阵
void
GetJacobiMatrix(MESH *Mesh,VOLU *volu, double **B) //double b[3])
{
  MESH *mesh;
  mesh=Mesh;
  double x0=mesh->Verts[volu->Verts[0]]->Coord[0];
  double y0=mesh->Verts[volu->Verts[0]]->Coord[1];
  double z0=mesh->Verts[volu->Verts[0]]->Coord[2];
  double x1=mesh->Verts[volu->Verts[1]]->Coord[0];
  double y1=mesh->Verts[volu->Verts[1]]->Coord[1];
  double z1=mesh->Verts[volu->Verts[1]]->Coord[2];
  double x2=mesh->Verts[volu->Verts[2]]->Coord[0];
  double y2=mesh->Verts[volu->Verts[2]]->Coord[1];
  double z2=mesh->Verts[volu->Verts[2]]->Coord[2];
  double x3=mesh->Verts[volu->Verts[3]]->Coord[0];
  double y3=mesh->Verts[volu->Verts[3]]->Coord[1];
  double z3=mesh->Verts[volu->Verts[3]]->Coord[2];
  B[0][0]=x1-x0; 
  B[0][1]=x2-x0; 
  B[0][2]=x3-x0; 
  B[1][0]=y1-y0; 
  B[1][1]=y2-y0; 
  B[1][2]=y3-y0; 
  B[2][0]=z1-z0; 
  B[2][1]=z2-z0; 
  B[2][2]=z3-z0; 
}

//得到矩阵的逆
void
GetInverseMatrix(double **A,double **IA)
{
  double a11=A[0][0];
  double a12=A[0][1];
  double a13=A[0][2];
  double a21=A[1][0];
  double a22=A[1][1];
  double a23=A[1][2];
  double a31=A[2][0];
  double a32=A[2][1];
  double a33=A[2][2];
  double B=a11*a22*a33 - a11*a23*a32 - a12*a21*a33 + a12*a23*a31 + a13*a21*a32 - a13*a22*a31;
  IA[0][0]=  (a22*a33-a23*a32)/B;
  IA[0][1]= -(a12*a33-a13*a32)/B;
  IA[0][2]=  (a12*a23-a13*a22)/B;
  IA[1][0]= -(a21*a33-a23*a31)/B;
  IA[1][1]=  (a11*a33-a13*a31)/B;
  IA[1][2]= -(a11*a23-a13*a21)/B;
  IA[2][0]=  (a21*a32-a22*a31)/B;
  IA[2][1]= -(a11*a32-a12*a31)/B;
  IA[2][2]=  (a11*a22-a12*a21)/B;
}


/**从网格最终获得一个单元，用来进行有限元空间等的操作 */
void GetElement3D(MESH *mesh,int volu_num, ELEMENT3D *element3D)
{
  element3D->ID_Num = volu_num;
  element3D->Num_Verts = 4;
  int i;
  i = mesh->Volus[volu_num]->Verts[0];
  element3D->Vert_X[0] = mesh->Verts[i]->Coord[0];
  element3D->Vert_Y[0] = mesh->Verts[i]->Coord[1];
  element3D->Vert_Z[0] = mesh->Verts[i]->Coord[2];
  i = mesh->Volus[volu_num]->Verts[1];
  element3D->Vert_X[1] = mesh->Verts[i]->Coord[0];
  element3D->Vert_Y[1] = mesh->Verts[i]->Coord[1];
  element3D->Vert_Z[1] = mesh->Verts[i]->Coord[2];
  i = mesh->Volus[volu_num]->Verts[2];
  element3D->Vert_X[2] = mesh->Verts[i]->Coord[0];
  element3D->Vert_Y[2] = mesh->Verts[i]->Coord[1];  
  element3D->Vert_Z[2] = mesh->Verts[i]->Coord[2];
  i = mesh->Volus[volu_num]->Verts[3];
  element3D->Vert_X[3] = mesh->Verts[i]->Coord[0];
  element3D->Vert_Y[3] = mesh->Verts[i]->Coord[1];  
  element3D->Vert_Z[3] = mesh->Verts[i]->Coord[2];
  element3D->Volu = GetVolume(mesh,mesh->Volus[volu_num]);
  GetJacobiMatrix(mesh,mesh->Volus[volu_num],element3D->Jacobian); 
  GetInverseMatrix(element3D->Jacobian,element3D->Inv_Jacobian);
}


void InitialElem3D(MESH *mesh, ELEMENT3D *elem3D)
{ 
  int i,Num_Verts;
  Num_Verts = 4;
  elem3D->Vert_X = (double *)malloc(Num_Verts);
  elem3D->Vert_Y = (double *)malloc(Num_Verts);
  elem3D->Vert_Z = (double *)malloc(Num_Verts);

 
  elem3D->Inv_Jacobian = (double **)malloc(3);
  elem3D->Jacobian = (double **)malloc(3);
  for(i=0;i<3;i++)
   { 
    elem3D->Inv_Jacobian[i] = (double *)malloc(3);
   }
  for(i=0;i<3;i++)
   { 
    elem3D->Jacobian[i] = (double *)malloc(3);
   }
}

void FreeElem3D(ELEMENT3D *elem3D)
{
    free(elem3D->Vert_X);
    free(elem3D->Vert_Y);
    free(elem3D->Vert_Z);
    //printf("free2222\n");
    int i;
   
    for(i=0;i<3;i++)
    {
      //printf("%d\n",i);
      free(elem3D->Inv_Jacobian[i]);
      free(elem3D->Jacobian[i]);
      //printf("%d\n",i);
    }
    free(elem3D->Jacobian);
    free(elem3D->Inv_Jacobian);
    free(elem3D);
}



void GetElementCoord(ELEMENT3D *elem3D,double Quad_xi,double Quad_eta, double Quad_zeta,
                     double *Quad_x, double *Quad_y,double *Quad_z)
{
   // x-coordination
    Quad_x[0] = (1-Quad_xi-Quad_eta-Quad_zeta)*elem3D->Vert_X[0] + Quad_xi*elem3D->Vert_X[1] +
            Quad_eta*elem3D->Vert_X[2]+Quad_zeta*elem3D->Vert_X[3];
   //y-coordination
    Quad_y[0] = (1-Quad_xi-Quad_eta-Quad_zeta)*elem3D->Vert_Y[0] + Quad_xi*elem3D->Vert_Y[1] +
            Quad_eta*elem3D->Vert_Y[2]+Quad_zeta*elem3D->Vert_Y[3];
   //z-coordination
    Quad_z[0] = (1-Quad_xi-Quad_eta-Quad_zeta)*elem3D->Vert_Z[0] + Quad_xi*elem3D->Vert_Z[1] +
            Quad_eta*elem3D->Vert_Z[2]+Quad_zeta*elem3D->Vert_Z[3];
}

void FreeMesh(MESH *Mesh)
{
  int i;  
  for(i=0;i<Mesh->Num_Verts_Global;i++)  
    FreeVert(Mesh->Verts[i]);  
  free(Mesh->Verts);  
  for(i=0;i<Mesh->Num_Lines_Global;i++)
    FreeLine(Mesh->Lines[i]);
  free(Mesh->Lines);
  for(i=0;i<Mesh->Num_Faces_Global;i++)
    FreeFace(Mesh->Faces[i]);
  free(Mesh->Faces); 
  for(i=0;i<Mesh->Num_Volus_Global;i++)
  {
   FreeVolu(Mesh->Volus[i]);
  }
  free(Mesh->Volus);  
  free(Mesh);
}
void FreeVolu(VOLU *Volu)
{
  free(Volu);
}
void FreeFace(FACEE *Face)
{
   //free(Face->Volus_Owned);
   free(Face);
}

void FreeLine(LINE *Line)
{
   //free(Line->Faces_Owned);
   //free(Line->Volus_Owned);
   free(Line);
}

void FreeVert(VERT *Vert)
{
   //free(Vert->Lines_Owned);
   //free(Vert->Faces_Owned);
   //free(Vert->Volus_Owned);
   free(Vert);
}


double FACEE_AREA(MESH *mesh,FACEE *Face)
{
    double area;
   /* VERT *start, *end;
    double dx, dy;    
    int start_num, end_num;
    start_num = Line->Verts[0];
    end_num = Line->Verts[1];
    
    start = mesh->Verts[start_num];
    end = mesh->Verts[end_num];
    
    dx = end->Coord[0] - start->Coord[0];
    dy = end->Coord[1] - start->Coord[1];

    length = sqrt(dx*dx +dy*dy);
    
    return length;
*/
    return 1;
    }

//检查得到的网格是否协调
void CheckMesh(MESH *mesh)
{
   int i;
   int v1,v2,v3,v4,l1,l2,l3,l4,l5,l6,f1,f2,f3,f4;
   for(i=0;i<mesh->Num_Volus_Global;++i)
   {
     
     v1=mesh->Volus[i]->Verts[0];
     v2=mesh->Volus[i]->Verts[1];
     v3=mesh->Volus[i]->Verts[2];
     v4=mesh->Volus[i]->Verts[3];
     l1=mesh->Volus[i]->Lines[0];
     l2=mesh->Volus[i]->Lines[1];
     l3=mesh->Volus[i]->Lines[2];
     l4=mesh->Volus[i]->Lines[3];
     l5=mesh->Volus[i]->Lines[4];
     l6=mesh->Volus[i]->Lines[5];
     
     f1=mesh->Volus[i]->Faces[0];
     f2=mesh->Volus[i]->Faces[1];
     f3=mesh->Volus[i]->Faces[2];
     f4=mesh->Volus[i]->Faces[3];
     
     
     if((mesh->Lines[l1]->Verts[0]==v1&&mesh->Lines[l1]->Verts[1]==v2)||
          (mesh->Lines[l1]->Verts[1]==v1&&mesh->Lines[l1]->Verts[0]==v2))
     {
      
     }
     else
     {
       printf ("The first line of vould %d is wrong!\n",i);
     }
     
     if((mesh->Lines[l2]->Verts[0]==v1&&mesh->Lines[l2]->Verts[1]==v3)||
          (mesh->Lines[l2]->Verts[1]==v1&&mesh->Lines[l2]->Verts[0]==v3))
     {
      
     }
     else
     {
       printf ("The second line of vould %d is wrong!\n",i);
     }
     
     if((mesh->Lines[l3]->Verts[0]==v1&&mesh->Lines[l3]->Verts[1]==v4)||
          (mesh->Lines[l3]->Verts[1]==v1&&mesh->Lines[l3]->Verts[0]==v4))
     {
      
     }
     else
     {
       printf ("The third line of vould %d is wrong!\n",i);
     }
     
     if((mesh->Lines[l4]->Verts[0]==v2&&mesh->Lines[l4]->Verts[1]==v3)||
          (mesh->Lines[l4]->Verts[1]==v2&&mesh->Lines[l4]->Verts[0]==v3))
     {
      
     }
     else
     {
       printf ("The fourth line of vould %d is wrong!\n",i);
     }
     
          if((mesh->Lines[l5]->Verts[0]==v2&&mesh->Lines[l5]->Verts[1]==v4)||
          (mesh->Lines[l5]->Verts[1]==v2&&mesh->Lines[l5]->Verts[0]==v4))
     {
      
     }
     else
     {
       printf ("The fifth line of vould %d is wrong!\n",i);
     }
     
          if((mesh->Lines[l6]->Verts[0]==v3&&mesh->Lines[l6]->Verts[1]==v4)||
          (mesh->Lines[l6]->Verts[1]==v3&&mesh->Lines[l6]->Verts[0]==v4))
     {
      
     }
     else
     {
       printf ("The sixth line of vould %d is wrong!\n",i);
     }
     
     if(mesh->Faces[f1]->Verts[0]!=v1&&mesh->Faces[f1]->Verts[1]!=v1&&mesh->Faces[f1]->Verts[2]!=v1
       &&(mesh->Faces[f1]->Verts[0]==v2||mesh->Faces[f1]->Verts[0]==v3||mesh->Faces[f1]->Verts[0]==v4)
       &&(mesh->Faces[f1]->Verts[1]==v2||mesh->Faces[f1]->Verts[1]==v3||mesh->Faces[f1]->Verts[1]==v4)
       &&(mesh->Faces[f1]->Verts[2]==v2||mesh->Faces[f1]->Verts[2]==v3||mesh->Faces[f1]->Verts[2]==v4))
     {
       
     }
     else
      {
       printf ("The fisrt face of vould %d is wrong!\n",i);
     } 
     
     if(mesh->Faces[f2]->Verts[0]!=v2&&mesh->Faces[f2]->Verts[1]!=v2&&mesh->Faces[f2]->Verts[2]!=v2
       &&(mesh->Faces[f2]->Verts[0]==v1||mesh->Faces[f2]->Verts[0]==v3||mesh->Faces[f2]->Verts[0]==v4)
       &&(mesh->Faces[f2]->Verts[1]==v1||mesh->Faces[f2]->Verts[1]==v3||mesh->Faces[f2]->Verts[1]==v4)
       &&(mesh->Faces[f2]->Verts[2]==v1||mesh->Faces[f2]->Verts[2]==v3||mesh->Faces[f2]->Verts[2]==v4))
     {
       
     }
     else
      {
       printf ("The second face of vould %d is wrong!\n",i);
     } 
     
     if(mesh->Faces[f3]->Verts[0]!=v3&&mesh->Faces[f3]->Verts[1]!=v3&&mesh->Faces[f3]->Verts[2]!=v3
       &&(mesh->Faces[f3]->Verts[0]==v1||mesh->Faces[f3]->Verts[0]==v2||mesh->Faces[f3]->Verts[0]==v4)
       &&(mesh->Faces[f3]->Verts[1]==v1||mesh->Faces[f3]->Verts[1]==v2||mesh->Faces[f3]->Verts[1]==v4)
       &&(mesh->Faces[f3]->Verts[2]==v1||mesh->Faces[f3]->Verts[2]==v2||mesh->Faces[f3]->Verts[2]==v4))
     {
       
     }
     else
      {
       printf ("The third face of vould %d is wrong!\n",i);
     } 
     
     if(mesh->Faces[f4]->Verts[0]!=v4&&mesh->Faces[f4]->Verts[1]!=v4&&mesh->Faces[f4]->Verts[2]!=v4
       &&(mesh->Faces[f4]->Verts[0]==v1||mesh->Faces[f4]->Verts[0]==v2||mesh->Faces[f4]->Verts[0]==v3)
       &&(mesh->Faces[f4]->Verts[1]==v1||mesh->Faces[f4]->Verts[1]==v2||mesh->Faces[f4]->Verts[1]==v3)
       &&(mesh->Faces[f4]->Verts[2]==v1||mesh->Faces[f4]->Verts[2]==v2||mesh->Faces[f4]->Verts[2]==v3))
     {
       
     }
     else
      {
       printf ("The forth face of vould %d is wrong!\n",i);
     } 
     
     
     
      if((mesh->Faces[f1]->Lines[0]==l4||mesh->Faces[f1]->Lines[0]==l5||mesh->Faces[f1]->Lines[0]==l6)
       &&(mesh->Faces[f1]->Lines[1]==l4||mesh->Faces[f1]->Lines[1]==l5||mesh->Faces[f1]->Lines[1]==l6)
       &&(mesh->Faces[f1]->Lines[2]==l4||mesh->Faces[f1]->Lines[2]==l5||mesh->Faces[f1]->Lines[2]==l6))
     {
       
     }
     else
      {
       printf ("The fisrt line of face of vould %d is wrong!\n",i);
     } 
     
      if((mesh->Faces[f2]->Lines[0]==l2||mesh->Faces[f2]->Lines[0]==l3||mesh->Faces[f2]->Lines[0]==l6)
       &&(mesh->Faces[f2]->Lines[1]==l2||mesh->Faces[f2]->Lines[1]==l3||mesh->Faces[f2]->Lines[1]==l6)
       &&(mesh->Faces[f2]->Lines[2]==l2||mesh->Faces[f2]->Lines[2]==l3||mesh->Faces[f2]->Lines[2]==l6))
     {
       
     }
     else
      {
       printf ("The second line of face of vould %d is wrong!\n",i);
     } 
     if( (mesh->Faces[f3]->Lines[0]==l1||mesh->Faces[f3]->Lines[0]==l5||mesh->Faces[f3]->Lines[0]==l3)
       &&(mesh->Faces[f3]->Lines[1]==l1||mesh->Faces[f3]->Lines[1]==l5||mesh->Faces[f3]->Lines[1]==l3)
       &&(mesh->Faces[f3]->Lines[2]==l1||mesh->Faces[f3]->Lines[2]==l5||mesh->Faces[f3]->Lines[2]==l3))
     {
       
     }
     else
      {
       printf ("The third line of face of vould %d is wrong!\n",i);
     } 
     if((mesh->Faces[f4]->Lines[0]==l1||mesh->Faces[f4]->Lines[0]==l2||mesh->Faces[f4]->Lines[0]==l4)
       &&(mesh->Faces[f4]->Lines[1]==l1||mesh->Faces[f4]->Lines[1]==l2||mesh->Faces[f4]->Lines[1]==l4)
       &&(mesh->Faces[f4]->Lines[2]==l1||mesh->Faces[f4]->Lines[2]==l2||mesh->Faces[f4]->Lines[2]==l4))
     {
       
     }
     else
      {
       printf ("The fourth line of face of vould %d is wrong!\n",i);
     } 
     
    // printf("no problem-------------------------------------------------\n");

   
   }

}
               
/*  
void GetFaceInfo4Mesh(MESH *mesh)
{
  int *a,*b,idx_face,i,j,total;
  idx_face=0;
  VOLU *volu;
  int num_verts,num_volus;
  num_verts=mesh->Num_Verts_Global;
  num_volus=mesh->Num_Volus_Global;
  total=num_verts*num_verts*num_verts;
  a=malloc(total*sizeof(int));
  memset(a,-1,total*sizeof(int));
  for(i=0;i<num_volus;++i)
  {
     volu=mesh->Volus[i];
     for(j=0;j<4;++j)//四个面的循环
     {
      FaceInfo(volu->Verts[(j+1)%4],volu->Verts[(j+2)%4],volu->Verts[(j+3)%4],a,b,num_verts);
      if(*b<0)
       {
         volu->Faces[j]=idx_face;
	 mesh->Faces[idx_face]->Verts[0]=volu->Verts[(j+1)%4];
	 mesh->Faces[idx_face]->Verts[1]=volu->Verts[(j+2)%4];
	 mesh->Faces[idx_face]->Verts[2]=volu->Verts[(j+3)%4];
         switch(i)
	 {
	  case 0: 
	    mesh->Faces[idx_face]->Lines[0]=5;
	    mesh->Faces[idx_face]->Lines[1]=4;
	    mesh->Faces[idx_face]->Lines[2]=3;
	    break;
	 case 1: 
	    mesh->Faces[idx_face]->Lines[0]=2;
	    mesh->Faces[idx_face]->Lines[1]=1;
	    mesh->Faces[idx_face]->Lines[2]=5;
	    break;
         case 2: 
	    mesh->Faces[idx_face]->Lines[0]=0;
	    mesh->Faces[idx_face]->Lines[1]=4;
	    mesh->Faces[idx_face]->Lines[2]=2;
	    break;
         case 3: 
	    mesh->Faces[idx_face]->Lines[0]=3;
	    mesh->Faces[idx_face]->Lines[1]=1;
	    mesh->Faces[idx_face]->Lines[2]=0;
	    break;
	 }//switch结束
         *b=idx_face;   
         idx_face++;
       }
      else
         volu->Faces[j]=*b;
     
     }//面循环结束
  }//单元循环结束

  free(a);
}
     */ 

int FaceInfo(int a, int b, int c, int *array,int *point,int num)
{
  int a1,a2,a3;
  if(a<b&&b<c)
  {
    a1=a;
    a2=b;
    a3=c; 
  }
  else if(a<c&&c<b)
  {
    a1=a;
    a2=c;
    a3=b; 
  }
   else if(b<c&&c<a)
  {
    a1=b;
    a2=c;
    a3=a; 
  }
   else if(b<a&&a<c)
  {
    a1=b;
    a2=a;
    a3=c; 
  }
   else if(c<a&&a<b)
  {
    a1=c;
    a2=a;
    a3=b; 
  }
  else
  {
    a1=c;
    a2=b;
    a3=a; 
  }
point=array+a1*num+a2*num+a3;
}



void CopyMesh(MESH *mesh_dest, MESH *mesh_source)
{
  int num_verts,num_lines,num_faces,num_volus; 
  int New_verts,New_lines,New_faces,New_volus;
  VERT **Verts_new,**Verts_old;
  LINE **Lines_new,**Lines_old ;
  FACEE **Faces_new,**Faces_old;
  VOLU **Volus_new,**Volus_old;
  //VERT *vert;
  LINE *line, *new_line_1, *new_line_2, *new_line_3;
  FACEE *face, *new_face_1, *new_face_2;
  VOLU *volu;  
  int i,j,k, new_num, refine_ind;
  int elem_verts, elem_lines;  
  mesh_dest->Num_Verts_Global = mesh_source->Num_Verts_Global;
  mesh_dest->Num_Lines_Global = mesh_source->Num_Lines_Global;
  mesh_dest->Num_Faces_Global = mesh_source->Num_Faces_Global;
  mesh_dest->Num_Volus_Global = mesh_source->Num_Volus_Global;
  //复制所有的点，线，面的信息
  num_verts = mesh_dest->Num_Verts_Global;
  mesh_dest->Verts = (VERT **)malloc(num_verts);
  Verts_new = mesh_dest->Verts;
  Verts_old = mesh_source->Verts;
  //ShowMesh(mesh_source);
  for(i=0;i<num_verts;i++)
  {
    Verts_new[i] = (VERT *)malloc(1);
    CopyVert(Verts_new[i],Verts_old[i]);
  }
  
  num_lines = mesh_dest->Num_Lines_Global;
  mesh_dest->Lines = (LINE **)malloc(num_lines);
  Lines_new = mesh_dest->Lines;
  Lines_old = mesh_source->Lines;
  for(i=0;i<num_lines;i++)
  {
    Lines_new[i] = (LINE *)malloc(1);
    CopyLine(Lines_new[i],Lines_old[i]);
  }   

  num_faces = mesh_dest->Num_Faces_Global;
  mesh_dest->Faces = (FACEE **)malloc(num_faces);
  Faces_new = mesh_dest->Faces;
  Faces_old = mesh_source->Faces;
  for(i=0;i<num_faces;i++)
  {
    Faces_new[i] =(FACEE *)malloc(1);
    CopyFace(Faces_new[i],Faces_old[i]); 
  }  
  num_volus = mesh_dest->Num_Volus_Global;

  mesh_dest->Volus = (VOLU **)malloc(num_volus);
  Volus_new = mesh_dest->Volus;
  Volus_old = mesh_source->Volus;
  for(i=0;i<num_volus;i++)
  {
    Volus_new[i] = (VOLU *)malloc(1);
    CopyVolu(Volus_new[i],Volus_old[i]); 
  }  
}

void CopyVert(VERT *vert_dest, VERT *vert_source)
{
  vert_dest->Num = vert_source->Num;
  vert_dest->Coord[0] = vert_source->Coord[0];
  vert_dest->Coord[1] = vert_source->Coord[1];
  vert_dest->Coord[2] = vert_source->Coord[2];
}


void CopyLine(LINE *line_dest, LINE *line_source)
{
  line_dest->Num = line_source->Num;
  //memcpy(line_dest->Verts,line_source->Verts,2*sizeof(int));
  line_dest->Verts[0]=line_source->Verts[0];
  line_dest->Verts[1]=line_source->Verts[1];
}

void CopyFace(FACEE*face_dest, FACEE *face_source)
{
  face_dest->Num = face_source->Num;
  //memcpy(face_dest->Verts,face_source->Verts,3*sizeof(int));
  //memcpy(face_dest->Lines,face_source->Lines,3*sizeof(int));
  face_dest->Verts[0]=face_source->Verts[0];
  face_dest->Verts[1]=face_source->Verts[1];
  face_dest->Verts[2]=face_source->Verts[2];
  face_dest->Lines[0]=face_source->Lines[0];
  face_dest->Lines[1]=face_source->Lines[1];
  face_dest->Lines[2]=face_source->Lines[2];
  face_dest->ID_Boundary=face_source->ID_Boundary;
}

void CopyVolu(VOLU *volu_dest, VOLU *volu_source)
{
  volu_dest->Num = volu_source->Num;
  //memcpy(volu_dest->Verts,volu_source->Verts,4*sizeof(int));
  //memcpy(volu_dest->Lines,volu_source->Lines,6*sizeof(int));
  //memcpy(volu_dest->Faces,volu_source->Faces,4*sizeof(int));
  volu_dest->Verts[0]=volu_source->Verts[0];
  volu_dest->Verts[1]=volu_source->Verts[1];
  volu_dest->Verts[2]=volu_source->Verts[2];
  volu_dest->Verts[3]=volu_source->Verts[3];
  volu_dest->Lines[0]=volu_source->Lines[0];
  volu_dest->Lines[1]=volu_source->Lines[1];
  volu_dest->Lines[2]=volu_source->Lines[2];
  volu_dest->Lines[3]=volu_source->Lines[3];
  volu_dest->Lines[4]=volu_source->Lines[4];
  volu_dest->Lines[5]=volu_source->Lines[5];
  volu_dest->Faces[0]=volu_source->Faces[0];
  volu_dest->Faces[1]=volu_source->Faces[1];
  volu_dest->Faces[2]=volu_source->Faces[2];
  volu_dest->Faces[3]=volu_source->Faces[3];

  volu_dest->Father = volu_source->Father;
  volu_dest->Ancestor = volu_source->Ancestor; 
}

/*
void AddFefunct(double *fefunct,double *values_FE,int nrows,double *nrow)
{
 int i;
 for(i=0;i<nrows;++i)
 {
  fefunct[nrow[i]]=values_FE[i];  
 }
}
*/

