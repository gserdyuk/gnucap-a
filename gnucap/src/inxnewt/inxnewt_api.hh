/*
Newton API 
function types
functions called or been called by inexact/exact/../ newton drivers
Gennady Serdyuk
gserdyuk@ridgetop-group.com
$Id: inxnewt_api.hh,v 1.9 2004/01/29 11:15:28 gserdyuk Exp $
*/

#ifndef NEWTON_API_HH
#define NEWTON_API_HH

#include "extradata.hh"   // class SolverExtraData;

#include "solver.hh"    // FVEC_PTR_TYPE
//typedef void (*FVEC_PTR_TYPE) (int *n, double* x, double* f, int *retcode, SolverExtraData *ed);
typedef FVEC_PTR_TYPE NEF_TYPE_INX;

//typedef void (*NEF_TYPE_INX)(int *N, double *x, double *fvc, int *retcode, SolverExtraData *P); 
// P may be passed through Fortran routines - so shoudl be represented as "integer"
// NEF_TYPE function may be called from C and Fortran - so all parameters are pointers to be conform to Fortran conventions

//---------------------------------------------------------------------------------

typedef void (*MATVECOP)(int *N, double *X, double *Y, 
                         int *nelt, int *ia, int *ja, double *a, int *isym);  // matrix-vector operation. may be called from fortran

void unitMV (int *N, double *X, double *Y, int *nelt, int *ia, int *ja, double *a, int *isym);
void multAX (int *N, double *X, double *Y, int *nelt, int *ia, int *ja, double *a, int *isym);
void dderF  (int *N, double *X, double *Y, int *nelt, int *ia, int *ja, double *a, int *isym);

// all those functions should be supplied

//---------------------------------------------------------------------------------
class MatVecProduct{                            // class to incapsulate matrix-vector product arguments
                                                // in simplest case may contain matrix . Or (for matrix-free) - contains parameters for function call 
    int N_;
public:                    
    MatVecProduct(int n){N_=n;};
    virtual ~MatVecProduct(){};
    virtual MATVECOP func(){return unitMV;};
    virtual int N(){return N_;};
    };
    
class MatVec_PlainJac: public MatVecProduct {      // simple jacobian matvec class with by-column (fortan-like) storage scheme
    double * jac_;                              // jacobian                        
public:
    MatVec_PlainJac(int N, double *j):MatVecProduct(N){ jac_=j;};
    virtual ~MatVec_PlainJac(){};
    virtual MATVECOP func(){ return multAX; };
    double *j(){ return jac_;};                 // for usage inside of MULT_AX()
    };

//---------------------------------------------------------------------------------        
class MatVec_DDer: public MatVecProduct {       // Directional Derivative matvec class - matrix free scheme
    NEF_TYPE_INX  func_;                        // function, 
    double eta_;                                // used precision, to calculate finite differences
    double *xc_;                                // current x-vector
    double *fc_;                                // furrent f-fector
    void *p_;
public:
    MatVec_DDer(int N, NEF_TYPE_INX f, double e, double * xc, double *fc, void* p):MatVecProduct(N){ 
                func_=f; eta_=e; xc_=xc; fc_=fc; p_=p;};
    virtual ~MatVec_DDer(){};
    virtual MATVECOP func() {return dderF; };
    NEF_TYPE_INX nef(){ return func_;};               // for usage inside of DDER()
    double *xc(){return xc_;};
    double *fc(){return fc_;};
    double eta(){return eta_;};
    void *p()   {return p_;};
    };
    
// note, that MatVec objects should be used with corespondent MATVECOP functions, which
// KNOW what matvec class is used. That is done because function is passed into FORTRAN 
// subroutine SEPARATELY from their data, so data and function should be agreed
// Indeed, MatVecData object pointer is passed inside MATVECOP as void pointer and 
// shoiuld be properly interpreted inside.

//---------------------------------------------------------------------------------
// same is valid for MSOLVEOP and PrecondData class. PrecondDGMRES class is data, which is passed 
// into (*PRECONOP) function as last argument (void *), so body PRECONOP and this class 
// should be agreed
typedef void (*PRECONOP)(int *N, double *R, double *Z, 
                         int *nelt, int *ia, int *ja, double *a, int *isym, 
                         double *rwork, int *iwork); // preconditioner function. may be called from fortran

void unitPR (int *N, double *R, double *Z, int *nelt, int *ia, int *ja,double *a, int *isym, double *rwork, int *iwork);
class PrecondDGMRES{   // class to incapsulate preconditioner arguments. same as 
                    // MatVecdata is virtual and contain soimethinf to pass into PRECON
                    // It is a root of all preconditioner classes. Contain not only methods are required 
                    // by DGMRES, byt others - preparation of preconditioner, initialization, etc.
protected:
    int N_;
public: 
    PrecondDGMRES(int n){N_=n;};
    virtual ~PrecondDGMRES(){};
    virtual PRECONOP func() {return &unitPR; };
    virtual int N() {return N_;};
    virtual int isize() { return 1;}; // should be argeed with func(); 1- default value
    virtual int rsize() { return 1;}; // same as isize
    virtual int create(void){return 1;};   // action, performed once for preconditioner - for example read it from file 
    virtual int update(double *x, double *f,NEF_TYPE_INX NEF,void *pp,double *sx,double teta){return 1;};    // action, being performed with preconditioner AT EACH STEP - for example update
    virtual int precondSide(void){return 0; };  // 0 - np reconditioning, >0- right, <0 - left.
    
    };

//---------------------------------------------------------------------------------
void LUSolveBlockPreconditioner (int *N, double *R, double *Z, int *nelt, int *ia, int *ja, double *a, int *isym, double *rwork, int *iwork);
class LUConstBlockPrecond: public PrecondDGMRES {  // LU factorised constant block (M = N/Nfreq) preconditioner
protected:
    int M_;         //< block matrix size
    int NF_;        //< number of frequencies
    double *bj_;    //< factorized block preconditioner pointer
    int *ipiv_;
    double *jac_;   //source of data to build preconditioner (generally - dc jacobian)
public:
    LUConstBlockPrecond(int M, int NF, double * j): PrecondDGMRES(2*M*NF)
            {M_=M; NF_=NF; bj_=new double[M_*M_]; ipiv_=new int[M_]; jac_=j;};  // store block size and matrix pointer
    virtual ~LUConstBlockPrecond(){delete[] bj_; delete[] ipiv_;};
    virtual double *bj(){return bj_;};
    virtual int *ipiv(){return ipiv_;};
    virtual PRECONOP func() { return &LUSolveBlockPreconditioner; };
    virtual int M() {return M_;};
    virtual int create(void); // redefine processing before iterations - factorize jacobian
    virtual int precondSide(void){return -1; };  // <0 - left.
    };
    
//---------------------------------------------------------------------------------
class LUVariableBlockPrecond: public LUConstBlockPrecond {  // LU factorised constant block (M = N/Nfreq) preconditioner
public:
    LUVariableBlockPrecond(int M, int NF): LUConstBlockPrecond(M,NF,0){};
    virtual ~LUVariableBlockPrecond(){};
    virtual int create(void){return 1;}; // do nothing 
    virtual int update(double *x, double *f,NEF_TYPE_INX NEF,void *pp,double* sx,double teta);     // update preconditioner 
    };
    
//---------------------------------------------------------------------------------
int  inxnslvr(int N, NEF_TYPE_INX NEF,double *fvc, double *x, double *typf,double *typx,void *P, void *PJ,
              double fvectol, double steptol, double mintol,double maxstep, int maxiter,
              int precisegrad, int &termcode, int globtype, int usejac, double *jac,
              PrecondDGMRES* precond, MatVecProduct* mv, int debuglevel,
              int iunit_lit, int maxl_lit_abs, double maxl_lit_rel, double maxl_lit_increase,
              int lit_restarts, int precond_update_n, int eta_update
              );  // inexact solver
          

#endif
