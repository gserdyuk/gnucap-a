/*
fortran-c calling for newton suite

covers BLAS, LAPACK and SLATEC routines mostly.
*/
#ifndef NEWTON_FC_HH
#define NEWTON_FC_HH

#include "inxnewt_api.hh"

// calling FORTAN from C
// some constants (it is impossible to use cinstant in calls)
static int ONE=1;
static int ZERO=0;
static double DONE=1.;
static double DZERO=0.;

extern "C" {

// BLAS
#define DCOPY dcopy_
#define DSCAL dscal_
#define DGEMV dgemv_
#define DAXPY daxpy_
#define DNRM2 dnrm2_

void DCOPY(int *N, double *x, int *i, double *y, int *j);
void DSCAL(int *N, double *a, double *x, int *i); 
void DGEMV(char *TRANS,int *N, int* N,double *alpha, double *A, int *N,
           double *X,int * incrx, double *beta, double *Y, int *incry);
void DAXPY(int *N, double *a, double *X, int *ix, double *Y, int* iy);
double DNRM2(int *N, double *x, int *ix);
           
//LAPACK

#define DGETRF dgetrf_
#define DGETRS dgetrs_
void DGETRF(int *M, int *N, double *A, int *LDA, int *IPIV, int *INFO);
void DGETRS(char* TRANS, /*int *TRANS_len, */ int *N, int *NRHS, double *A, int *LDA, int *iPiv, double *B, int *LDB, int *info);

//SLATEC

// OWN
#define MACHEPS macheps_
#define LINFNORMSC linfnormsc_
#define L2NORMSC l2normsc_

double MACHEPS(void);
double LINFNORMSC(int *N, double *a, double *s);
double L2NORMSC(int *N, double *a, double *s);

#define LINSRCDIV2 linsrcdiv2_
#define LINSRCDER  linsrcder_
#define LINSRC3P   linsrc3p_

void LINSRCDIV2(int *N, double *x, double *xp, double *dx, NEF_TYPE_INX NEF,
                void *p, double * fvp, double * fplus, double *sf, double *fc,
                int *retcode);

void LINSRCDER(int *N,double *x,double *fc,NEF_TYPE_INX NEF, double *grad,
               double *dx, double *sx, double *maxstep, double *steptol,
               int *retcode,double *xp, double *fplus, double *fvp,
               int *maxtaken, void *p, double *lambda,double *sf);
void LINSRC3P (int *N, double *x, double *fc,NEF_TYPE_INX NEF,double *dx,double *sx,
               double *maxstep,double *steptol, int *retcode, double *xp,
               double *fplus,double *fvp,int *maxtaken, void *p, double *lambda,double *sf);

#define FDJAC fdjac_
#define GRADJAC gradjac_
#define FDGRAD fdgrad_
#define FDJACDC fdjacdc_

void FDJAC  (int *N,double *x,double *fvc,NEF_TYPE_INX NEF,void *pF,double *sx,double *meps,double *jac);
void GRADJAC(int *N,double *jac,double *sf,double *fvc,double *grad);
void FDGRAD (int *N,double *sf,double *sx,NEF_TYPE_INX NEF,void* pF,double *grad);
void FDJACDC(int *N,double *x, double *f,NEF_TYPE_INX NEF,void* pp,double *sx, double *teta,double *jac, int *M);      

#define CLEAN clean_
void CLEAN(int *N, double *d);


#define DGMRES dgmres_
void DGMRES (int *N, double *fvc, double *sstep, int *nelt, int *ia, int *ja, double *a, int *isym,
             MATVECOP MATVEC, PRECONOP PRECON,            
            int *itol_lit, double *tol_lit, int *itmax_lit, int *iter_lit, double *err_lit, int *ierr_lit,
            int *iunit_lit,  double *sf, double *sx, double *rgwk, int *lrgw, int *igwk, int *ligw, 
	        double *rwork, int *iwork);

}

#endif
