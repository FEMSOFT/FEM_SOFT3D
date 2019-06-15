// =======================================================================
// @(#)LinAlg.C        1.18 07/03/00
//
// Purpose:     basic routines for eigenvalue solving
//
// Author:      hhxie@lsec.cc.ac.cn
// =======================================================================

#include "EigenSolver.h"
#include "LinearAlg.h"
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "CG.h"
#include "DIRECTSOLVER.h""
#include "Matrix.h""

/**========================================================*/
/** Fortran subrutines for eigenvalues solving  */

//the subrutine for the symmetric matrix of the double precision
void dsaupd_(int *ido, //ç¨æ¥æç¤ºå¦ä½è¿è¡ééè®¯æ¥è¿è¡Arnoldiè¿­ä»£
	    const char *bmat, // è¡¨ç¤ºç¹å¾å¼é®é¢çç±»åï¼ å½Mä¸ºåä½ç©éµæ¯bmat='I'ï¼å¦å bmat=âGâ
	    int *n, // è¡¨ç¤ºæè¦è®¡ç®çç¹å¾å¼é®é¢çç»´æ°ï¼ç¹å¾åéçç»´æ°ï¼
	    const char *which, //æç¤ºæéè¦è®¡ç®çç¹å¾å¼ï¼ âLMâï¼âSMâ
	    int *nev, //è¡¨ç¤ºæéè¦è®¡ç®ç¹å¾å¼çä¸ªæ°
	    double *tol, //Arnoldiè¿­ä»£åæ­¢çè¯¯å·®é
	    double *resid, //æ®éåéï¼æ¯æ¬¡Arnoldiè¿­ä»£çåå§åé
	    int *ncv,     //Arnoldiè¿­ä»£è¿ç¨ä¸­ï¼Krylovå­ç©ºé´ä¸­çç»´æ° (æå¤å°ä¸ªåé)
	    double *v,    /* ç¨æ¥å­Krylovå­ç©ºé´çåéï¼å¤§å°æ¯ n*ncv */
	    int *ldv,  	  /* Krylovå­ç©ºé´ä¸­çåéVä¸­æ¯ä¸ªåéçé¿åº¦*/
	    int *iparam,  /* åæ°è¡¨ï¼éé¢è®°å½äºå®ä¹å·ä½ç®æ³çä¸äºåæ° */
	    int *ipntr,   /* WORKDåWORKLçèµ·å§ä½ç½®ï¼å¨å·ä½çç®æ³å®ç°ä¸­å¯ä»¥ä½ä¼å° */
	   /* Note (from ./PARPACK/SRC/MPI/pdsaupd.f):
	    *     -------------------------------------------------------------
	    *     IPNTR(1): pointer to the current operand vector X in WORKD.
	    *     IPNTR(2): pointer to the current result vector Y in WORKD.
	    *     IPNTR(3): pointer to the vector B * X in WORKD when used in
	    *               the shift-and-invert mode.
	    *     IPNTR(4): pointer to the next available location in WORKL
	    *               that is untouched by the program.
	    *     IPNTR(5): pointer to the NCV by 2 tridiagonal matrix T in WORKL.
	    *     IPNTR(6): pointer to the NCV RITZ values array in WORKL.
	    *     IPNTR(7): pointer to the Ritz estimates in array WORKL associated
	    *               with the Ritz values located in RITZ in WORKL.
	    *     IPNTR(11): pointer to the NP shifts in WORKL. See Remark 6 below.
	    *
	    *     Note: IPNTR(8:10) is only referenced by pdseupd . See Remark 2.
	    *     IPNTR(8): pointer to the NCV RITZ values of the original system.
	    *     IPNTR(9): pointer to the NCV corresponding error bounds.
	    *     IPNTR(10): pointer to the NCV by NCV matrix of eigenvectors
	    *                of the tridiagonal matrix T. Only referenced by
	    *                pdseupd  if RVEC = .TRUE. See Remarks. */
	    double *workd,  /* 3*nå¤§å°çä¸ä¸ªåéï¼æ¯æç®æ³å¨å¶ä¸­è¿è¡ work array of length 3*N */
	    double *workl,  /* work array of length LWORKL */
	    int *lworkl,  /* sizeof WORKL, at least NCV**2 + 8*NCV */
	    int *info    /* è¾å¥ä¿¡æ¯ Input:
			 *	0:  a randomly initial residual vector is used
			 *	!0: RESID contains the initial residual vector
			 * Output:
			 *	0: normal exit
			 *	1: maximum number of iterations taken.
			 *	... ... */
	    );
		    
void dseupd_(int *rvec,    /* FALSE: Ritz values only, TRUE: Ritz vectors */
	    char *All,     /* 'A': all, 'S': some (specified by SELECT[]) */
	    int *select,   /* logical array of dimension NCV */
	    double *d,     /* eigenvalues, array of dimension NEV (output) */
	    double *v1,    /* Ritz vectors, N by NEV array (may use V) */
	    int *ldv1,     /* leading dimension of Z */
	    double *sigma, /* represents the shift */
	    /* the following arguments are the same as for pdsaupd */
	    const char *bmat, 
	    int *n,
	    const char *which, 
	    int *nev,
	    double *tol, 
	    double *resid, 
	    int *ncv,
	    double *v,
	    int *ldv,
	    int *iparam, 
	    int *ipntr,
	    double *workd,
	    double *workl,
	    int *lworkl,
	    int *ierr
	    );		  
		  
//the subrutine for the non-symmetric matrix of the double precision	  		  
extern void dnaupd_(int *ido, const char *bmat, int *n, const char *which,
	      int *nev, double *tol, double *resid,
	      int *ncv, double *v, int *ldv,
	      int *iparam, int *ipntr, double *workd,
	      double *workl, int *lworkl, int *info);
	      
extern void dneupd_(int *rvec, char *All, int *select,
	      double *dr, double *di, double *v1, int *ldv1,
	      double *sigmar, double *sigmai, double *workev, const char *bmat,
	      int *n, const char *which, int *nev, double *tol,
	      double *resid, int *ncv,
	      double *v, int *ldv, int *iparam,
	      int *ipntr, double *workd,
	      double *workl, int *lworkl, int *ierr);			
		  

/** calculate the eigenvalues of the system using ARPACK routines*/
//å¨åç²¾åº¦ä¸æ±è§£å¯¹ç§°ç¹å¾å¼é®é¢ï¼æä»¬è¿éåè®¾å¯¹äºè¢«æ±è§£çååº¦ç©éµåè´¨éç©éµå·²ç»
// å¤çäºè¾¹çæ¡ä»¶ï¼å¨è¿éçè®¡ç®å°ä¸åå¤çè¾¹çæ¡ä»¶
void DSEigenSolver(MATRIX *A,  MATRIX *M, int nev, int indictor,double tau,
		double *Evals, double **Evecs, bool Change)
{
/* 
c     %--------------------------------------------------------%
c     | The work array WORKL is used in DSAUPD as         |
c     | workspace.  Its dimension LWORKL is set as         |
c     | illustrated below.  The parameter TOL determines   |
c     | the stopping criterion.  If TOL<=0, machine        |
c     | precision is used.  The variable IDO is used for  |  
c     | reverse communication and is initially set to 0.  |
c     | Setting INFO=0 indicates that a random vector is  |
c     | generated in DSAUPD to start the Arnoldi          |
c     | iteration.                                        | 
c     %-----------------------------------------------------% */
  int i,j,k,N_Eqn, ActiveBound, begin, end;
  int n = A->N_Rows; // ç¹å¾å¼è®¡ç®çç»´æ°
  int ido = 0; 
  const char *bmat; //è®°å½ç¹å¾å¼è®¡ç®çç±»åï¼âGâ æè'I'
  const char *which; //åªæ¯éè¦æ±è§£ç¹å¾å¼çèå´
  char All[4] = "All"; 
  double tol = 0.0;
  double *resid, *rhs;  //æ®éåé
  resid = malloc(sizeof(double)*n); 
  rhs = malloc(sizeof(double)*n);
  int ncv = 4*nev;   //Krylovå­ç©ºé´çåéçä¸ªæ°ï¼ä¸è¬ncv > 2*nev (nevè¡¨ç¤ºè¦æ±ç¹å¾å¼çä¸ªæ°)
  if (ncv>n) 
    ncv = n;
  double *v;  //ç¨æ¥å­Krylovå­ç©ºé´ä¸­çåé
  int ldv = n;
  //æ³¨æä¸é¢è¿ä¸ªåéå®ä¹ä»¥åå®çå¤§å°
  v = malloc(sizeof(double)*ldv*ncv);
  //ä¸é¢å®ä¹åæ°iparamï¼ç¨å°çå°±ç»åºå·ä½çå¼ï¼ä¸ç¨å°çå°±ä¸è¦ç®¡ï¼å ä¸ºåç®æ³æ§è¡è¿ç¨ä¸­ä¼èµå¼
  int *iparam;
  iparam = malloc(sizeof(int)*11); // new int[11];
  iparam[0] = 1;      //è¡¨ç¤ºç§»å¨ç­ç¥ï¼1ï¼è¡¨ç¤ºç¨äºshift
  iparam[2] = 3*n;   //æå¤§Arnoldiè¿­ä»£æ¬¡æ°
  iparam[6] = 3;     //åªæ¯è®¡ç®æ¨¡å¼
  int *ipntr;        //æç¤ºworkdåworklçä½ç½®åæ°ï¼ä¸Arnoldè¿­ä»£ä¹åçå­é®é¢çæ±è§£æå³
  ipntr = malloc(sizeof(int)*11); 
  double *workd; //æ³¨æworkdçå¤§å°
  workd = malloc(sizeof(double)*3*n); 
  double *workl; //worklæ¯ç¨æ¥å­å¨ä½QRåè§£çç©éµç
  workl = malloc(sizeof(double)*ncv*(ncv+8)); 
  int lworkl = ncv*(ncv+8);  //è®°å½worklçå¤§å°
  int info = 0;
  int rvec = 1;    //ï¼ç¨æ¥å¾å°ç¹å¾åéçæç¤ºå­ï¼ TRUEè¡¨ç¤ºè¦å¾å°ç¹å¾åéï¼FALSEï¼è¡¨ç¤ºåªéè¦å¾å°ç¹å¾å¼ï¼
  int *select;
  select = malloc(sizeof(int)*ncv); 
  double *d;  //ç¨æ¥å­å¨ç¹å¾å¼ ï¼QRè¿­ä»£çæ¶åä¹ç¨å®æ¥å­å¨ææå°ç©éµçç¹å¾å¼ï¼
  d = malloc(sizeof(double)*2*ncv); 
  double sigma; 
  sigma = tau;
  int ierr; 
  int *Row, *KCol;
  double *Values, *Valus_A;
/*
c     %---------------------------------------------------%
c     | This program uses exact shifts with respect to    |
c     | the current Hessenberg matrix (IPARAM(1) = 1).    |
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 3 specified in the      |
c     | documentation of DSAUPD is used (IPARAM(7) = 3).  |
c     | All these options may be changed by the user.     |
c     | For details, see the documentation in DSAUPD.     |
c     %---------------------------------------------------% */
  if((M) && (tau>0.0))
  {
   // MatrixAdd(A, M, -tau);  // è®¡ç®ç©éµ A=A-tau*M
  }
  if(indictor > 0)
    which = "LM";  //è®¡ç®æ¨¡æå¤§çç¹å¾å¼
  else 
    which = "SM";  //è®¡ç®æ¨¡æå°çç¹å¾å¼
  if(M==NULL)
    bmat = "I";   //ï¼å¦ææ²¡æè´¨éç©éµï¼è¡¨ç¤ºæ åç¹å¾å¼é®é¢
  else 
    bmat = "G";   //ï¼å¦ææè´¨éç©éµï¼è¡¨ç¤ºä¸è¬çç¹å¾å¼é®é¢
    
  //ä¸é¢è¿è¡ç¹å¾å¼æ±è§£è¿­ä»£
    /*
    c     %-----------------------------------------------%
    c     | M A I N   L O O P (Reverse communication) |
    c     %-----------------------------------------------%
    */
  do {
      /*
      c        %---------------------------------------------%
      c        | Repeatedly call the routine DSAUPD and take |
      c        | actions indicated by parameter IDO until    |
      c        | either convergence is indicated or maxitr   |
      c        | has been exceeded.                          |
      c        %---------------------------------------------% 
      */
      dsaupd_(&ido, bmat, &n, which, &nev, &tol, resid, 
	      &ncv, v, &ldv, iparam, ipntr, workd, workl,
	      &lworkl, &info); //ä¸æ¯æéçåéåé¢è¦å ä¸&ä»¥åå°å
      
      printf("Ipntr[1]=%d, Ipntr[2]=%d, Ipntr[3]=%d\n",ipntr[0],ipntr[1],ipntr[2]);
      
      switch(ido){
	/*
	c           %------------------------------------------------%
	c           | Perform  y <--- OP*x = inv[A-SIGMA*M]*M*x  |
	c           | to force the starting vector into the      |
	c           | range of OP.  The user should supply       |
	c           | his/her own matrix vector multiplication   |
	c           | routine and a linear system solver here.   |
	c           | The matrix vector multiplication routine   |
	c           | takes workd(ipntr(1)) as the input vector. |
	c           | The final result is returned to            |
	c           | workd(ipntr(2)).                           |
	c           %------------------------------------------------%
	*/
	case -1: 
	  //åè®¡ç® y = M*x
	  printf("ido=%d, Matrix-vector and solving\n",ido);
	  MatrixDotVec(M, workd+ipntr[0]-1, rhs);
	  //åè®¡ç® inv[A-SIGMA*M]y=w
	  CG(A, rhs, workd+ipntr[1]-1,1e-15,1000); 
	  break;
	case 1:
	  printf("ido=%d, solving\n",ido);
	  //æ workd(ipntr(3))æ·è´å° workd(ipntr(2))
	  //call dcopy ( n, workd(ipntr(3)), 1, workd(ipntr(2)), 1)
	  AssignVec(rhs,workd+ipntr[1]-1,n);
	  //è®¡ç® inv[A-SIGMA*M]y=w, ç»æå­å¨ workd(ipntr(2))
	  CG(A,rhs, workd+ipntr[1]-1,1e-15,1000); 
	  break;
	case 2:
	  //åä¸¾è¯ç¸ä¹ï¼è¾å¥workd(ipntr(1))ï¼è¾åºworkd(ipntr(2))
	  //workd(ipntr(2)) = M*workd(ipntr(1))
	  printf("ido=%d, Matrix multipy\n",ido);
	  //è®¡ç® y = M*x
	  MatrixDotVec(M, workd+ipntr[0]-1, workd+ipntr[1]-1);
	  break;	  
      }	
    } while ((ido==1)||(ido==-1)||(ido==2));

  if (info<0) {
         printf("Error with dsaupd, info = %d\n ",info);
         printf("Check documentation in dsaupd\n\n");
  } else {
    dseupd_(&rvec, All, select, d, v, &ldv, &sigma, bmat,
	    &n, which, &nev, &tol, resid, &ncv, v, &ldv,
	    iparam, ipntr, workd, workl, &lworkl, &ierr);    

    if (ierr!=0) {
      printf("Error with dseupd, info = %d\n",ierr); 
      printf("Check the documentation of dseupd.\n\n");
    } else if (info==1) {
       printf("Maximum number of iterations reached.\n\n");
    } else if (info==3) {
      printf("No shifts could be applied during implicit\n");
      printf("Arnoldi update, try increasing NCV.\n\n");
    }
    //printf("sigma= %d\n",sigma); 
    //printf("ldv= %d\n",ldv); 
    //printf("n= %d\n",n);
    int i, j;
    for (i=0; i<nev; i++)
    {
      Evals[i] = d[i];
      printf("The %d-th eigenvalue: %18.15f\n",i,Evals[i]);
    }
    
    if(Evecs)
    {
      Evecs = malloc(sizeof(double *)*nev); //new double *[nev];    
      for (i=0; i<nev; i++)
      {
	Evecs[i] = malloc(sizeof(double)*n); // new double [n];
	for (j=0; j<n; j++) Evecs[i][j] = v[i*n+j];
      }
    }
  }
  free(resid);
  free(v);
  free(iparam);
  free(ipntr);
  free(workd);
  free(workl);
  free(select);
  free(d);
}


/**   ====================================================*/

/** calculate the eigenvalues of the system using ARPACK routines for double non-symmetric matrix*/
//D åMè¡¨ç¤ºç¹å¾å¼é®é¢ä¸­çååº¦ç©éµåè´¨éç©éµ
//nevï¼ è¡¨ç¤ºè¦è®¡ç®ç¹å¾å¼çä¸ªæ°
//indicatorï¼ 
//tauï¼ åºç¨ shift-inverseè®¡ç®æ¨¡å¼
//Evals:å­å¨è®¡ç®åºçç¹å¾å¼
//Evecsï¼å­å¨è®¡ç®åºçç¹å¾åé
// Changeï¼ æ¯å¦éè¦å¯¹è´¨éç©éµè¿è¡æ¹åï¼ï¼
void DNEigenSolver(MATRIX *A,  MATRIX *M, int nev, int indictor, double tau, 
		   double *Evals, double **Evecs, bool Change)
{ 
  //OutPut("memory: " << setw(10) << GetMemory() << endl);
  int i,j,k,N_Eqn, ActiveBound, begin,end;
  int n = A->N_Rows;   //æè¦è®¡ç®çç©éµçç»´æ°
  int ido = 0;         //æç¤ºä¸ä¸æ­¥è¦åä»ä¹çåé
  const char *bmat;    //æ ç¤ºè¦åæ åç¹å¾å¼è®¡ç®æèæ¯å¹¿ä¹ç¹å¾å¼è®¡ç®
  const char *which;   //è¡¨ç¤ºè¦ç®åªé¨åçç¹å¾å¼ï¼æå¤§ï¼æèæå°ï¼æèç¦»æä¸ªæ°æè¿
  char All[4] = "All";
  double tol = 0.0;     //æ§å¶ç²¾åº¦
  double *resid;        //å­å¨åå§å¼çä¸ä¸ªåé
  resid = malloc(sizeof(double)*n); 
  int ncv = 4*nev;    
  if (ncv>n) ncv = n;
  double *v;
  int ldv = n;
  v = malloc(sizeof(double)*ldv*ncv); 
  int *iparam;
  iparam = malloc(sizeof(int)*11); 
  iparam[0] = 1;
  iparam[2] = 3*n;
  iparam[6] = 3;
  int *ipntr;
  ipntr = malloc(sizeof(int)*14); 
  double *workd;
  workd = malloc(sizeof(double)*3*n); 
  int lworkl = 3*ncv*ncv + 6*ncv;
  double *workl;
  workl = malloc(sizeof(double)*lworkl); 
  double *rwork;
  rwork = malloc(sizeof(double)*ncv); 
  int info = 0;
  int rvec = 1;      //Changed from above
  int *select;
  select = malloc(sizeof(int)*ncv); 
  double *dr, *di;
  dr = malloc(sizeof(double)*ncv); 
  di = malloc(sizeof(double)*ncv); 
  double sigmar, sigmai;
  sigmar = tau;
  sigmai = 0.0;
  double *workev;
  workev = malloc(sizeof(double)*3*ncv); 
  int ierr;
  double *rhs;
  rhs = malloc(sizeof(double)*n); 
  if(M==NULL)
  {
    bmat = "I";
    iparam[6] = 1;    
    if(indictor > 0)
      which = "LM";
    else 
      which = "SM";    
  }
  else
  {
    bmat = "G";
    iparam[6] = 3;
    if(indictor > 0)
    {
      printf("compute the smallest eigenvalue!\n");
      which = "SM";
    }
    else
    { 
      printf("compute the largest eigenvalue!\n");
      which = "LM"; 
    }
    if(tau!=0.0)
    {
      //è®¡ç® A =A-SIGMA*M
      MatrixAdd(A, M, -tau);
    }    
  }

  //=================================================================
  // start reverse communivation process
  //=================================================================
  do
  {
    dnaupd_(&ido, bmat, &n, which, &nev, &tol, resid, 
	    &ncv, v, &ldv, iparam, ipntr, workd, workl,
	    &lworkl,  &info);  
    //cout<<"iparam[6]= "<<iparam[6]<<endl;
    if(iparam[6] == 3)
    {
      //printf("ido= %d\n",ido); //<<endl;
      switch(ido)
      {	
	case -1:
	  /* Perform y <-- OP(x), where OP(x) is
	  inv[A-sigma*M]*M*x for shift-invert mode, or
	  inv[M]*A*x for regular-invert mode,
	  with x = WORKD(IPNTR(1)), y = WORKD(IPNTR(2))*/
	  MatrixDotVec(M, workd+ipntr[0]-1, rhs);	  
	   //DirectSolver(A, rhs, workd+ipntr[1]-1);
	  //MatVect(A, workd+ipntr[1]-1, rhs_temp);
	  //Daxpy(n,-1,rhs_temp,rhs);
	  //CG(A, rhs, workd+ipntr[1]-1,1e-8,100);
	  //AMGSOLVER(A, rhs, workd+ipntr[1]-1,1e-8,100);
	  DIRECTSOLVER(A, rhs, workd+ipntr[1]-1);
	  break;
	case 1:
	  /* Perform y <-- OP(x), where OP(x) is
	    inv[A-sigma*M]*M*x for shift-invert mode, or
	    inv[M]*A*x for regular-invert mode,
	    with x = WORKD(IPNTR(1)), y = WORKD(IPNTR(2))
	    Note: M*x is already stored at WORKD(IPNTR(3))
	    and need not be recomputed for shift-invert mode */    
	  //DirectSolver(A, workd+ipntr[2]-1, workd+ipntr[1]-1);
	  //MatVect(A,workd+ipntr[1]-1,rhs_temp);
	  //Daxpy(n,-1,rhs_temp,rhs);
	  //CG(A, workd+ipntr[2]-1, workd+ipntr[1]-1,1e-8,100);
	  //AMGSOLVER(A, workd+ipntr[2]-1, workd+ipntr[1]-1,1e-8,100);
	  DIRECTSOLVER(A, workd+ipntr[2]-1, workd+ipntr[1]-1);
	  break;
	case 2:
            MatrixDotVec(M, workd+ipntr[0]-1, workd+ipntr[1]-1);
	    break;
	default: 
	  //printf("Finish for Eigenvalue solver!\n");
	  //OutPut(endl);
	  break;
      }
    }
    else
      if ((ido==1)||(ido==-1)) 
  	MatrixDotVec(A, workd+ipntr[0]-1, workd+ipntr[1]-1);
      
  } while((ido==1)||(ido==-1)||(ido==2));
  //printf("Begin Compute The Eigenvector!\n");
  if (info<0)
  {
    printf("Error with dnaupd, info = %d\n",info);
    printf("Check documentation in znaupd\n\n");
  } 
  else
  {
    //è®¡ç®ç¹å¾åé
    dneupd_(&rvec, All, select, dr, di, v, &ldv, &sigmar,&sigmai, workev,
	    bmat, &n, which, &nev, &tol, resid, &ncv, v, &ldv,
	    iparam, ipntr, workd, workl, &lworkl, &ierr);
    //printf("info=%d\n",info);
    if (ierr!=0) 
    {
      printf("Error with zneupd, info = %d\n",ierr);
      printf("Check the documentation of zneupd.\n\n");
    } 
    else
      if(info==1) 
      {
	printf("Maximum number of iterations reached.\n\n");
      }
      else
	if (info==3)
	{
	  printf("No shifts could be applied during implicit\n");
	  printf("Arnoldi update, try increasing NCV.\n\n");
	}
    for (i=0; i<nev; i++)
    {
      Evals[i] = dr[i];
      //printf("the %d-th eigenvalue: %18.15f\n",i,Evals[i]);
    }
    if(Evecs)
    {
      //printf("Start save eigenfunction: ldv=%d, ncv=%d\n",ldv,ncv);
      for (i=0; i<nev; i++)
      {
	//cout<<"i= "<<i<<endl;
	//printf("i=%d\n",i);
	for (j=0; j<n; j++)
	{ 
	  //printf("j=%d\n",j);
	  //printf("v[%d]=%18.15f\n",i*n+j,v[i*n+j]);
	  Evecs[i][j] = v[i*n+j];  
	}
      }
    }
    //printf("End of save eigenfunction: \n");
    // Sort the eigenvalues and the corresponding eigenfunctions
    double temp;
    if(indictor>0)
    {
      //è¦æç¹å¾å¼åç¹å¾åéååæåº
      for (i=0; i<nev; i++)
      {
	for (j=i; j<nev; j++)
	{
	  if (Evals[j] > Evals[i])
	  {
	    temp = Evals[j];
	    Evals[j] = Evals[i];
	    Evals[i] = temp;
	    for (k=0; k<n; k++) 
	    {
	      temp = Evecs[i][k];
	      Evecs[i][k] = Evecs[j][k];
	      Evecs[j][k] = temp;	      
	    }	    
	  }	  
	}	
      }//end for nev
    }
    else
    {
      for (i=0; i<nev; i++)
      {
	for (j=i; j<nev; j++)
	{
	  if (Evals[j] < Evals[i])
	  {
	    temp = Evals[j];
	    Evals[j] = Evals[i];
	    Evals[i] = temp;
	    for (k=0; k<n; k++) 
	    {
	      temp = Evecs[i][k];
	      Evecs[i][k] = Evecs[j][k];
	      Evecs[j][k] = temp;	      
	    }	    
	  }	  
	}	
      }//end for nev
    }
  }//end for info
  //OutPut("memory: " << setw(10) << GetMemory() << endl);
  if(v)
   free(v);
    //OutPut("memory: " << setw(10) << GetMemory() << endl);
  if(iparam)
    free(iparam);  //OutPut("memory: " << setw(10) << GetMemory() << endl);
  if(ipntr)
    free(ipntr);
  //OutPut("memory: " << setw(10) << GetMemory() << endl);
  if(workd)
    free(workd);
  //OutPut("memory: " << setw(10) << GetMemory() << endl);
  if(workl)
    free(workl);
  //OutPut("memory: " << setw(10) << GetMemory() << endl);
  if(rwork)
    free(rwork);
  //OutPut("memory: " << setw(10) << GetMemory() << endl);
  if(select)
    free(select);
  //OutPut("memory: " << setw(10) << GetMemory() << endl);
  if(dr)
    free(dr);
  //OutPut("memory: " << setw(10) << GetMemory() << endl);
  if(di)
    free(di);
  //OutPut("memory: " << setw(10) << GetMemory() << endl);
  if(workev)
    free(workev);
  //OutPut("memory: " << setw(10) << GetMemory() << endl);
  if(rhs)
    free(rhs); 
  if(resid)
  {
    free(resid);
  }
}








/*
void DenseEigen(double **stiff_matrix,double **mass_matrix,  int dim, double *evals, double **evecs)
{

  int i,j;
  Engine *ep;
  if(!(ep = engOpen(NULL)))
  {
    printf("cannot open engine!");
    exit(1);
  }
  engSetVisible(ep,false);//å½æ°è°ç¨è¿ç¨ä¸­ä¸åºç°matlabè°ç¨ç»é¢
  double a[dim][dim],b[dim][dim];
  for(i=0;i<dim;++i)
    for(j=0;j<dim;++j)
    {
     a[i][j]=stiff_matrix[i][j];
     b[i][j]=mass_matrix[i][j];
    }
  double eigenvalues[dim][dim],eigenvector[dim][dim];
 
  mxArray *H = NULL,*T = NULL;
  H = mxCreateDoubleMatrix(dim,dim,mxREAL);
  T = mxCreateDoubleMatrix(dim,dim,mxREAL);
  memcpy((void *)mxGetPr(H),(void *)a, sizeof(double)*dim*dim);   
  memcpy((void *)mxGetPr(T),(void *)b, sizeof(double)*dim*dim);   
  engPutVariable(ep,"H",H);                                           
  engPutVariable(ep,"T",T);                                           

  //  char p[100]={1};
  //engOutputBuffer(ep,p,100);
  int ret1 = engEvalString(ep,"[V,D]=eig(H',T')");
  //printf("%s\n",p);
  mxArray *E=NULL,*F=NULL;
  E=engGetVariable(ep,"V");
  F=engGetVariable(ep,"D");
  memcpy((void *)eigenvector,(void *)mxGetPr(E), sizeof(double)*dim*dim);
  memcpy((void *)eigenvalues,(void *)mxGetPr(F), sizeof(double)*dim*dim);
  for(i=0;i<dim;++i)
     evals[i] = eigenvalues[i][i];
 
  for(i=0;i<dim;++i)
    for(j=0;j<dim;++j)
      evecs[i][j] = eigenvector[i][j];

  mxDestroyArray(F);
  mxDestroyArray(E);
  mxDestroyArray(H);
  mxDestroyArray(T);
  engClose(ep);

}

*/


