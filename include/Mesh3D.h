/*
 * =====================================================================================
 *
 *       Filename:  mesh.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2014年10月04日 15时14分08秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Xu Fei (Fly), xufei@lsec.cc.ac.cn
 *        Company:  
 *
 * =====================================================================================
 */
//点边面体的结构

#ifndef __MESH__
#define __MESH__

#define NVERTS_PER_FACEE 3
#define NLINES_PER_FACEE 3
#define NVERTS_PER_VOLU 4
#define NLINES_PER_VOLU 6
#define NFACEES_PER_VOLU 4
#define NVERTS_PER_LINE 2
//typedef unsigned int size_t;


//typedef double COORD[3]; 
typedef struct VERTS_ {
  int Num;
  double Coord[3];
  int Num_Lines_Owned; 
  int *Lines_Owned;
  int Num_Faces_Owned;
  int *Faces_Owned;
  int Num_Volus_Owned;
  int *Volus_Owned;
  
  //int ID_Boundary;
  
  int Count;
  int Mark;//yizhijiami  
} VERT;

typedef struct LINE_ {
  int Num;
  int Verts[2];
  int Num_Faces_Owned;
  int *Faces_Owned;
  int Num_Volus_Owned;
  int *Volus_Owned;
  int Count;
  int Mark;
  
  int ID_Boundary;
  
  int Child[3];  
} LINE;

typedef struct FACEE_ {
  int Num;
  int Verts[NVERTS_PER_FACEE];
  int Lines[NLINES_PER_FACEE];
  int Num_Volus_Owned;
  int *Volus_Owned;   
  
  int ID_Boundary;
  
  int Child[3];
  int Old_Child[2];
  
  int Count;
} FACEE;

typedef struct VOLU_{
  int Num;
  int Verts[NVERTS_PER_VOLU];
  int Lines[NLINES_PER_VOLU];
  int Faces[NFACEES_PER_VOLU];
  int Child [2];
  int Ancestor;
  int Father;
  
  int Sign;
  int Mark;
  int Type;
  int Trans;
} VOLU;


 //网格结构
typedef struct MESHES_ {
  int Num_Verts_Global;
  int Num_Lines_Global;
  int Num_Faces_Global;
  int Num_Volus_Global;

  VERT **Verts;
  LINE  **Lines;
  FACEE **Faces;
  VOLU **Volus;

  int get_full_info;
  int Mark;
} MESH;


typedef struct ELEMENT3D_ {
  int ID_Num;
  int Num_Verts;
  double *Vert_X;
  double *Vert_Y;  
  double *Vert_Z;
  double Volu;
  double **Inv_Jacobian;
  double **Jacobian;
}ELEMENT3D;


//sub programs 
void MeshRefineConsistent(MESH *mesh);

void MeshFullInfo(MESH *mesh);
void MeshFullInfo4Line(MESH *mesh);
void MeshFullInfo4Face(MESH *mesh);
void MeshFullInfo4Num(MESH *mesh);
//判断一个单元的边是否被Mark
int MarkJudge(MESH *mesh,VOLU *volu);
//释放网格所占有的内存
void MeshDestroy(MESH *mesh);
void ReadMesh(MESH *mesh,char *file);
void ReadMeshBis(MESH *mesh,char *file);
void ReadMeshLshape(MESH *mesh,char *file);
void ReadMeshConcave(MESH *mesh,char *file);
void WriteMesh(MESH *mesh, char *file);
void WriteMeshBis(MESH *mesh,char *file);
void WriteMesh4PHG(MESH *mesh, char *file);
void GetFacesVolusOwned(MESH *mesh);
void ShowVolu(MESH *mesh,VOLU *volu);
void ShowMesh(MESH *mesh);
void CheckMesh(MESH *mesh);
static int cmp(const double *a, const double *b);
//void QuickSort(double  *x, int low, int high, size_t size, int *ix, int (*cmp)(const double *a,const double *b));
//网格信息补全
void MeshFullInfoConsistent(MESH *mesh);
int FaceInfo(int a, int b, int c, int *array,int *point,int num);
//网格为下一次加密更新
void MeshRenew(MESH *mesh);
//线在单元里的局部编号
int GetLineNum(LINE *line,int line_num,VOLU *volu);
void MarkMesh(MESH *mesh, VOLU *volu, int volu_num,LINE *line, int line_num, int *nline_mark);
void MeshRefineBisection(MESH *mesh, int *M, int *num);
void MeshRefineUniformBisection(MESH *mesh);
void MeshRefineUniformBisection2(MESH *mesh);
void MeshRefineUniformBisectionNew(MESH *mesh);
void MeshRefineUniformBisectionNew1(MESH *mesh);

int MeshMarkElemNum(MESH *mesh, double *eta, double *theta);
int* MeshMarkRefine(MESH *mesh,double *eta,double theta,int *num);
void MeshMarkElem(MESH *mesh, int *M, double *eta, double *theta,int *num);
void MeshBisection(MESH *mesh,int times);
void MeshConsist(MESH *mesh,int times);
//计算单元的体积
double GetVolume(MESH *mesh,VOLU *volu);
//求jacobi变换矩阵
void GetJacobiMatrix(MESH *mesh,VOLU *volu, double **B);
void GetInverseMatrix(double **A,double **IA);
//单元操作
void GetElement3D(MESH *mesh,int volu_num, ELEMENT3D *element3D);
void InitialElem3D(MESH *mesh, ELEMENT3D *elem3D);
void InitialMesh(MESH *mesh);
void FreeElem3D(ELEMENT3D *elem3D);
void GetElementCoord(ELEMENT3D *elem3D,double Quad_xi,double Quad_eta, double Quad_zeta,
                     double *Quad_x, double *Quad_y,double *Quad_z);

//释放内存
void FreeMesh(MESH *Mesh);
void FreeVolu(VOLU *Volu);
void FreeFace(FACEE *Face);
void FreeLine(LINE *Line);
void FreeVert(VERT *Vert);
double FACEE_AREA(MESH *mesh,FACEE *Face);
//网格调整
void AdjustMesh(MESH *mesh);
void AdjustMeshL(MESH *mesh);
void AdjustMeshForRHF(MESH *mesh);
int GetLongestLine(double d1,double d2,double d3,double d4,double d5,double d6);

void CopyMesh(MESH* mesh_dest, MESH*mesh_source);
void CopyVert(VERT *vert_dest, VERT *vert_source);
void CopyLine(LINE *line_dest, LINE *line_source);
void CopyFace(FACEE*face_dest, FACEE *face_source);
void CopyVolu(VOLU *volu_dest, VOLU *volu_source);



void GetFace4Mesh(MESH *mesh);

#endif
