/*
 INexact newton solver project (INSol) Project 
 designed to debugging of different ways of globalizaton
 Gennady Serdyuk
 gsserdyuk@ridgetop-group.com

 Working Laguage - C++

 INterface
 INEXACT newton solver
 inxslv(int N,NEF,fvc,x,typf,typx,pf,pj,fvectol,steptol,mintol,maxstep,
          maxiter,precisegrad,termcode,globtype)

 int N          - system size 
 int (*NEF)(..) - subroutine to calculate vector-function
 double* fvc    - vector-function , length N
 double* x      - variables , size=N 
 double* typf, typx
                - vectors of typical values of F & X - used for scaling
                
 void *pf;      - pointer to pass extra parameters in NEF - 
 void *pj;      - pointer to pass extra parameters in NEJ (in addition to NEF)

 void (*NEJ)     - subroutine name for Jacobian
 void (*MULJX)   - subroutine name for multiplication Jacobian X Vector 
 void (*PRECON)  - subroutine fo preconitioning    

 double fvectol   - required absolute precision (l2 norm) - tolerance 
 double steptol   - minimal allowed absolute step (l2 norm) 
 double mintol    - parameter to detect saddle point
 double maxstep   - maxinum allowed step (to detect undefined function)
 int maxiter   - maxx number of iterations
 int precisegrad 
                  - yes or no. 1 - it is posssible to get precise gradient and 
                    perform gradient check
 int usejac       - use jacobian or not (matrix-free method)
 
 termcode - termination code
     = 0 - ok, continue iterations
     < 0 - error in input data:
      -1 - N<1
     > 0 - stop iterations
       1 - ||fvc|| < fvectol
       2 - ||dx||  < steptol (small step)
       3 - can't reduce ||fvc|| - globalization fails
       4 - iteration limit
       5 - 5 maxsteps have been made
       6 - local minimum detected

 globtype - type of globalization
       1 - linear search with division by 2
       2 - polinomial interpolation 2-nd or 3-rd order (requires gradient)
       3 - 3-point quadractic approximation linear search. see C.T.Kelley
           SIAM Frontiers..
debuglevel:
    0 - no msgs
    1 - main events- exit codes, etc
    2 - iterations
    3 - variables and errors at every iteration
    4 - print jacobian
    
*************************************************************************/
#include <iostream>
#include <math.h>
#include <assert.h>#include <malloc.h>

#include "inxnewt_fc.hh"
#include "inxnewt_api.hh"
#include "commondefs.hh"


using std::cout;
using std::cerr;

int incheck_inx  (int N, double macheps,double *typx,double *typf, double *SX, double *SF,
                   double fvectol, double steptol, double mintol, double maxstep, int itlimit);
int nestop0_inx  (int N,double *F, double *SF, double fvectol);
int nestop_inx   (int N,double *xc,double *xp, double *fvp, double fnorm,double *grad,double * sx,double *sf,
                   int retcode,double fvectol, double steptol, int itncount, int itnlimit, int &maxtaken,
                   double mintol,int &consecmax,int precisegrad);

void printv(int N, double *x){
    for (int i=0; i<N; i++)
        cout<<x[i]<<" ";
    cout<<endl;
    return;
    }
    
void printv(int N, int *x){
    for (int i=0; i<N; i++)
        cout<<x[i]<<" ";
    cout<<endl;
    return;
    }
    
void printvx(int N, double *x){
    for (int i=0; i<N; i++)
        if (x[i]==0)
            cout<<"_|";
        else 
            cout<<"x|";
    cout<<endl;
    return;
    }


int  inxnslvr(int N, NEF_TYPE_INX NEF,double *fvc, double *x, double *typf,double *typx,void *pf, void *pj,
              double fvectol, double steptol, double mintol,double maxstep, int maxiter,
              int precisegrad, int &termcode, int globtype, int usejac, double *jac,
              PrecondDGMRES *prc,MatVecProduct *mv, int debuglevel, 
              int iunit_lit, int maxl_lit_abs, double maxl_lit_rel, double maxl_increase, int lit_restarts,
              int precond_upd_n, int eta_update
              ) {
     

      double fc,fplus;
      double *xp   =dalloca(N);
      double *grad =dalloca(N);
      double *fvp  =dalloca(N);
// fc - weighted l2 norm of fvc
// xp - next point
// grad - vector gradient
// fvp - vector function (plus)
// fplus - weighted l2 norm of fvp

    double *sx  =dalloca(N); // sx - scale by x;  
    double *sf  =dalloca(N); // sf - scale by fv
    double *dx  =dalloca(N); // step ( corection vector)
    double *sstep=dalloca(N);// saved (previous) step - to call Krylov solver with nonzero initial approx

    double meps, lambda;
      
    int iter, consecmax=0, maxtaken=0, retcode;
            
    if (debuglevel >0) 
        cout<< " ::: INXsolvr :::"<<endl;
    meps=MACHEPS();
// check input data
    termcode=incheck_inx(N,meps,typx,typf,sx,sf,fvectol,steptol,
                         mintol,maxstep,maxiter);     
    if (debuglevel >=3) {
        cout<<"MEPS    ="<<meps<<endl;
        cout<<"globtype="<<globtype<<endl;
        cout<<"termcode="<<termcode<<endl;
        }
    
    if(termcode < 0) 
        return termcode; // bad input data
      
// start iterations  
    iter=0;
    int rc;   // return code;
    (*NEF)(&N,x,fvc,&rc,(SolverExtraData*)pf);
    fc=L2NORMSC(&N,fvc,sf);
    
    if (debuglevel >=3) {
        cout<<" NEF calculated\n";
        cout<<"   fc: "<<fc<<endl;      
        }
        
// check starting point    
    termcode=nestop0_inx(N,fvc,sf,fvectol);
            
// calculate  jacobian
      if (usejac == 1)
          FDJAC(&N,x,fvc,NEF,pf,sx,&meps,jac);
    
    if (debuglevel >=4) {
      cout<<"after jacobian calc "<<endl;
      cout<<"   x: "; printv(N,x);
      cout<<"   fvc: "; printv(N,fvc);
      if (usejac==1) {
        cout<<"   jacobian : \n";
        for (int j=0; j<N; j++) {
            cout<<" "<<j<<"\t";printvx(N,jac+j*N);
            }
        }
      }
      
// calculate gradient
    if (usejac == 1) 
        GRADJAC(&N,jac,sf,fvc,grad);
    else
        FDGRAD(&N,sf,sx,NEF,pf,grad);
    
    if (debuglevel >0)   
        cout<<"termcode= "<<termcode<<endl;

    CLEAN(&N,sstep);

//     === inner krylov precision control ====
    double eta;
    const double etamax=0.9999;
    const double gamma=0.9;
    const double ETA_DEFAULT=1.e-13; //sqrt(meps);  // ~ 1.e-9; value 1.e-14 caused a problem (see devel log for description, dd 26-01-2004)
    // int eta_update =0 set from prog interface
    switch (eta_update) {
	case 0: { // constant choise
		eta=ETA_DEFAULT; 
		break;
	}
	case 1: { // update eta. see Kelley, p. 105.
		eta=etamax;
		break;
	}
	default:{ // no update
	}
    }

//     === DGMRES PARAMETERS =========
// parameters, related to maxl will be moved in iteration loop , as maxl value may be changed during iterations
/*  igwk[1] = maxl - max krylow subsp size to find x-x0
    igwk[2]   = kmp - max krylow subspace size (ortho vectors): =maxl
    igwk[3] = jscale - scaling. 0 - no scaling, 1-x, 2-b, 3-both
    igwk[4] = jpre - preconditioning. 0- no, 
    igwk[5] = max no of restarts
            diagnostics:
    igwk[6] = mlwk - min len of rgwk
    igwk[7] = nms no calls of msolve
*/
    int itol_lit=0;             // lit - Linear ITerative
    int itmax_lit=0;            // see DGMRES - it is dummy
    double tol_lit, err_lit;
    int iter_lit,ierr_lit;
    int ligw=20; // default value for ligw - see GMRES description
    int igwk[ligw];

    int maxl=max(maxl_lit_abs, int(maxl_lit_rel*N));    // combined from absolute and relative values
    int lrgw=1+N*(maxl+6)+maxl*(maxl+3);                // size of wrk array = 1 + N*(MAXL+6) + MAXL*(MAXL+3)
    cout<<"***** maxl = "<<maxl<<endl;
    cout<<"***** lrgw = "<<lrgw<<endl;
//    double rgwk[lrgw];   - moved down - just before DGMRES

    int nelt=1;                                         // see DGMRES usage notes
    int *ia=ialloca(nelt);
    int *ja=ialloca(nelt);
    double a[nelt];
    int isym=1;    
    itol_lit    =0;
    // tol_lit   - defined later        
    igwk[1-1] = maxl;  
    igwk[2-1] = maxl;
    igwk[3-1] = 0; 
    igwk[4-1] = 0;  // preconditioner type; redefined later
    igwk[5-1] = lit_restarts;  // max number of restarts; -1 -> 0 restarts
    int isize=max(1,prc->isize()); // size of aux array, safeguarded; 
    int rsize=max(1,prc->rsize()); // size of aux array, safeguarded;
    int *iwork=ialloca(isize);
    double *rwork=dalloca(rsize);

    // control different matvec functions and preconditioners
    PRECONOP PRECON=prc->func();;
    igwk[4-1]=prc->precondSide();   // assuming <0 - using Left preconditioner; 0 - no preconditioner - see DGMRES
    ja[0]=int(prc);
    
    MATVECOP MATVEC=mv->func(); // MATVEC operation took from class - it is predefined (and related to matrix storage/presence 
    ia[0]=int(mv);

//     === END DGMRES PARS ======

    if (debuglevel >=3)
        cout<<" precondtype="<<prc->precondSide()<<"\n";  
        

    double *rgwk=new double[lrgw];
    
    while (termcode ==0) {
        iter++;
        tol_lit     =eta;
// starting iterations:
        if (debuglevel >=2) 
            cout<<"ITERATION   "<<iter<<endl<<" Err ="<<fc<<" "<<endl;
        if (debuglevel >=3) {
            cout<<"X ";printv(N,x);
            cout<<"F ";printv(N,fvc);
            }
//        print *,'   jac: ',jac
//        print *, 'fvc DGMRES=',fvc

/*
NOTES for GMRES usage.
Not to change much in GMRES, it is proposed to use NELT, IA, JA, A, RWORK and IWORK for teh following purposes:
NELT=1
A(1)    - is not used
IA(1)   - contains mv pointer 
JA(1)   - contains prc pointer
RWORK and IWORK should be alocated according to PRECONOP operation. default: [1] preconditioner and arrays are not used
*/        

        if (debuglevel >=3) {       
            cout<<"N="<<N<<endl;
            cout<<"lrgw before call="<<lrgw<<endl;    
            cout<<"IGWK(6) before call - MLWK (min size of RGWK[lrgw])="<<igwk[6-1]<<endl;    
            }
        if((iter-1)%precond_upd_n == 0){
            if (debuglevel >=3)
                cout<<"going to update precond @ iter"<<iter<<endl;;         
            prc->update(x,fvc,NEF,pf,sx,meps);   // update preconditioner, if necessary
            flush(cout);
            }
        int loop_again;
        
        do {
         // inner block - only in it scope "maxl"-dependent array is defined. size is dynamic in iteration loop
            loop_again=0;      // assume we will not loop unless if will be forced
            cout<<"lrgw="<<lrgw<<endl; //########################
            //double rgwk[lrgw]; - allocated in heap
            
            CLEAN(&N,sstep); // zero initial approximjation looks better - take less time to converge
/*                
            cout<<"BEFORE DGMRES-------------------------------------------------------\n";
            cout<<"fvc =   ";printv(N,fvc);
            cout<<"sstep = ";printv(N,sstep);
            cout<<"igwk  = ";printv(ligw,igwk);
            cout<<" N= "<<N<<" itol_lit="<<itol_lit<<" tol_lit (aka TOL)="<<tol_lit<<" iter_lit="<<iter_lit<<endl;
            cout<<"err_lit = "<<err_lit<<" ierr_lit= "<<ierr_lit<<endl;
//            cout<<"MODIFY ="<<modify<<endl;
            cout<<"--------------------------------------------------------------------\n";
/**/
            DGMRES (&N, fvc, sstep, &nelt, ia, ja, a, &isym, MATVEC, PRECON, 
                &itol_lit, &tol_lit, &itmax_lit, &iter_lit, &err_lit, &ierr_lit,
                &iunit_lit,  sf, sx, rgwk, &lrgw, igwk, &ligw, rwork, iwork);
/*
            cout<<"AFTER DGMRES-------------------------------------------------------\n";
            cout<<"fvc =   ";printv(N,fvc);
            cout<<"sstep = ";printv(N,sstep);
            cout<<"igwk  = ";printv(ligw,igwk);
            cout<<" N= "<<N<<" itol_lit="<<itol_lit<<" tol_lit="<<tol_lit<<" iter_lit="<<iter_lit<<endl;
            cout<<"err_lit = "<<err_lit<<" ierr_lit= "<<ierr_lit<<endl;
            cout<<"--------------------------------------------------------------------\n";
/**/
            if (debuglevel >=2)
                cout<<"KRYLOV: IERR, ERR, TOL="<<ierr_lit<<" "<<err_lit<<" "<<tol_lit<<endl;
            if (debuglevel >=3){                
                cout<<"lrgw after call="<<lrgw<<endl;    
                cout<<"IGWK(6) after call - MLWK (min size of RGWK[lrgw])="<<igwk[6-1]<<endl;    
                cout<<"IGWK(7) after call - no of MSOLVE calls="<<igwk[7-1]<<endl;
                }

            switch (ierr_lit) {
                case 1: {// reallocate rgwk and loop:
                        maxl=min(maxl+int(maxl_increase*N), N);                  // increase maxl by maxl_increase percents of N; no more than N
                        lrgw=1+N*(maxl+6)+maxl*(maxl+3);                    // size of wrk array = 1 + N*(MAXL+6) + MAXL*(MAXL+3)
                        igwk[1-1] = maxl;  
                        igwk[2-1] = maxl;
                        delete[] rgwk;
                        rgwk=new double[lrgw]; // allocate bigger size
                        if (debuglevel >=2){
                            cout<<" redefine MAXL - new one ="<<maxl;
                            cout<<"; restart current iteration"<<endl;
                            }
                        loop_again=1;  // restart DGMRES
                        if (maxl==N)    // #### error  - it should not happen
                            loop_again=0; 

                    break;   // 
                    }
                case 2:  break;
                case -1: break;
                case -2: break;
                case 0:  break ;// OK
                default:       ; //  should not occure - generate exception #### error
                }

        } while (loop_again == 1  );


        
        // handle failures:
        assert(ierr_lit==0);    //####
         
/*         
        if (debuglevel >=2)       
            cout<<"IERR, ERR, TOL="<<ierr_lit<<" "<<err_lit<<" "<<tol_lit<<endl;
        if (debuglevel >=3){
            cout<<"lrgw after call="<<lrgw<<endl;    
            cout<<"IGWK(6) after call - MLWK (min size of RGWK[lrgw])="<<igwk[6-1]<<endl;    
            cout<<"IGWK(7) after call - no of MSOLVE calls="<<igwk[7-1]<<endl;
            }
*/                
        DCOPY(&N,sstep,&ONE,dx,&ONE);
        double m_one=-1.;
        DSCAL(&N,&m_one,dx,&ONE);
        
        if (debuglevel >=3){
            cout<<"DX =";printv(N,dx);
            }
// check step (and find damping factor or so) 
        if (globtype == 1) 
            LINSRCDIV2(&N,x,xp,dx,NEF,&pf,fvp,&fplus,sf,&fc,&retcode);
        else if (globtype == 2) { 
            if (usejac != 1) 
                assert(0);
            LINSRCDER(&N,x,&fc,NEF,grad,dx,sx,&maxstep,&steptol,
                      &retcode,xp,&fplus,fvp,&maxtaken,
                      &pf,&lambda,sf);
            }
        else if (globtype == 3)
            LINSRC3P(&N,x,&fc,NEF,dx,sx,&maxstep,&steptol,
                     &retcode,xp,&fplus,fvp,&maxtaken,
                     &pf,&lambda,sf);

        if (debuglevel >=1)
            cout<<"FC, FPLUS    ="<<fc<<" "<<fplus<<endl;
        if (debuglevel >=3){
            cout<<"X\"updated\" =";printv(N,x);
            cout<<"XP updated =";printv(N,xp);
            cout<<"FVP= "; printv(N,fvp);
            }
	
	switch (eta_update) {
		case 0: { // constant choise
			eta=ETA_DEFAULT; 
            break;
		}
		case 1: { // update eta. see Kelley, p. 105.
        		double eta_a=gamma*(fplus/fc)*(fplus/fc); 
                if (debuglevel >=3)
                    cout<<"eta_a="<<eta_a<<endl;
			    // eta_b=min(etamax, eta_a);
        		double eta_c;
			    if((gamma*eta*eta) <= 0.1)
            		eta_c=min(etamax,eta_a);
        		else
            		eta_c=min(etamax,max(eta_a, gamma*eta*eta));
        		eta=min(etamax,max(eta_c,0.5*fvectol/fplus));
		        //eta=eta_a; // tmp1 #### - to force convergence on some linear cases
            if (debuglevel >=2)
                cout<<"resulting eta="<<eta<<endl;
			break;
            }
		default:{ // no update
			}
		}

// calculate  jacobian
        if (usejac == 1) 
            FDJAC(&N,xp,fvp,NEF,pf,sx,&meps,jac);
        
// check stopping condition
        termcode=nestop_inx(N,x,xp,fvp,fplus,grad,sx,sf,retcode,fvectol,
                steptol,iter,maxiter,maxtaken,mintol,
                consecmax,precisegrad);
     
        if(termcode == 1) { // finished OK. copy everything to whereever to return values        
            }
        
        DCOPY(&N,xp,&ONE,x,&ONE);
        DCOPY(&N,fvp,&ONE,fvc,&ONE);

        fc=fplus;
//###        
        if (debuglevel >=3){
            cout<<"at the end of iteration N "<<iter<<endl;
            cout<<" x= "; printv(N,x);
            cout<<" f= "; printv(N,fvc);
            }

        cout<<"-------------------------"<<endl<<"\t termcode="<<termcode<<endl;
        
    } // end WHILE
    delete[] rgwk;
    if (debuglevel >= 0)
        cout<<"exit. termcode =  "<<termcode<<endl;
      
    return termcode;
}

//**********************************************************
//      function incheck_inx(N,macheps,typx,typf,SX,SF,
//     +                 fvectol,steptol,mintol,maxstep,
//     +                 itlimit)
int incheck_inx  (int N, double macheps,double *typx,double *typf, double *SX, double *SF,
                   double fvectol, double steptol, double mintol, double maxstep, int itlimit){
// typx typf SX and SF are of size N
      
    int incheck_inx=0;
// check size of system            
    if(N<1)
        return -1;
        
// assign X scaling vector      
    for(int i=0; i<N; i++){
        SX[i]=1;
        if(typx[i] != 0)    
            SX[i]=1/typx[i];
        }

// assign F scaling vector      
    for (int i=0; i<N; i++){        
        SF[i]=1;
        if(typf[i] != 0)
            SF[i]=1/typf[i];
        }
      
// check fvectol,
    if(fvectol == 0) 
        fvectol=pow(macheps,(1./3.));

// check steptol
    if(steptol == 0) 
        steptol=pow(macheps,(2./3.));

// check mintol
    if(mintol == 0)  
        mintol= pow(macheps,(2./3.));

// check maxstep,
    if(maxstep == 0) 
        maxstep = 1000.;

// check itlimit
    if(itlimit == 0) 
        itlimit=100;
      
    return incheck_inx;
    }
      
//**********************************************************
//      function nestop0_inx(N,F,SF,fvectol)
int nestop0_inx(int N, double *F, double *SF, double fvectol){
    double ferr=0;
 
// max{SF(i)*|F(i)|} <=0.01*fvcetol
    for (int i=0; i<N; i++)
        ferr=max(ferr,SF[i]*abs(F[i]));

//      print *, 'nestop0_inx: ferr, fvectol',ferr, fvectol
    if(ferr <= 0.01*fvectol) 
        return 1;

    return 0;
    }
      
//**********************************************************
int nestop_inx(int N,double *xc,double *xp,double *fvp,double fnorm, double *grad,
               double *sx, double *sf,int retcode,double fvectol, double steptol,
               int itncount, int itnlimit, int &maxtaken, double mintol,
               int &consecmax, int precisegrad){
      
    int nestop_inx=0;
// nestop=1 - scaled |fp| is less then fvectol - solution is found
    double fpnorm=0;
    for (int i=0; i<N; i++)
        fpnorm=max(fpnorm,sf[i]*abs(fvp[i]));
    if(fpnorm<=fvectol){
        cout<<" err tolerance: err="<<fpnorm<<" tol="<<fvectol<<endl;
        return 1;
        }
//###    cout<<" &&&&&&&&&&&&&& err tolerance: err="<<fpnorm<<endl;
// nestop=2 - scaled |step| is less then steptol
    double stepnorm=0;
    for (int i=0; i<N; i++)
        stepnorm=max(stepnorm,abs(xp[i]-xc[i])/max(xp[i],1/sx[i]) );
    if(stepnorm <= steptol){
        cout<<" step tolerance: step="<<stepnorm<<" tol="<<steptol<<endl;
        return 2;
        }
        
// nestop=3 - can't reduce ||F|| - globalization fails 
// retcode is code from globalizator, if ==1 - glob. was failed
    if(retcode == 1)
        return 3;
        
// nestop=4 - iteration limit exceeded
    if(itncount>= itnlimit)
        return 4;
      
// nestop=5 - 5 maxsteps have been made - lower limit for F is in inf
    if(maxtaken == 1)
        consecmax=consecmax+1;
    if(consecmax == 5)
        return 5;
    else
        consecmax=0;
      
// nestop=6 - current point is local minimum for ||F||
    if(precisegrad==1) {
        double gradtol=0;
        for (int i=0; i<N; i++)
            gradtol=max(gradtol,abs(grad[i])*
                    max(xp[i],1./sx[i])/max(fnorm,double(N)/2) );
        if(gradtol<=mintol)
            return 6;
        }

      return 0;
    }
    
//------------ to call from fortran ---------------------------
/*
void inxnsf(int *N, int* NEF,double *fvc,double *x,double *typf,double *typx, void *pf, void* pj,
            double *fvectol,double *steptol,double *mintol,double *maxstep, int* maxiter,        
            int *precisegrad, int* termcode,int *globtype, int *usejac){
     
   // fortran adapter  
   inxnslvr(*N, (NEF_TYPE_INX) NEF,fvc, x, typf,typx,(void *)pf, (void *)pj,
              *fvectol, *steptol, *mintol,*maxstep, *maxiter,
              *precisegrad, *termcode, *globtype, *usejac);
    
    return;
    }
*/
