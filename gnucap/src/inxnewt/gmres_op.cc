/*
gserdyuk@ridhetop-group.com
$Id: gmres_op.cc,v 1.16 2004/03/17 23:59:21 gserdyuk Exp $
17-oct-2003
operations, required by GMRES.
*/

#include <iostream>
#include <malloc.h>

#include "inxnewt_api.hh"
#include "inxnewt_fc.hh"
#include "commondefs.hh"

void printv(int N, double *x);
// --- functions, which are used in MatVecData class -------------
// remember, that MatVecData structure is available via ia[0]
void unitMV (int *Np, double *X, double *Y, 
                         int *nelt, int *ia, int *ja, double *a, int *isym){
// represents multiplication onto UNIT matrix. it is just for example - not really used
// Y=[1]*X 
    MatVecProduct *mv=(MatVecProduct*)ia[0];             // that is for example
    DCOPY(Np,X,&ONE,Y,&ONE);                             // N is taken from parameters
    return;
}

void multAX (int *Np, double *X, double *Y, 
                         int *nelt, int *ia, int *ja, double *a, int *isym){
// represents multiplication onto J matrix. 
// Y=[1]*X 
    char trans='n';
    MatVec_PlainJac *mv=(MatVec_PlainJac*)ia[0];                  
    double *A=mv->j();
//cout<<" ===================== multAX ========\n";    
//cout<<"X: "; printv(*Np,X);
    DGEMV(&trans,Np,Np,&DONE,A,Np,X,&ONE,&DZERO,Y,&ONE);
//cout<<"Y: "; printv (*Np,Y);    
    return;                     
}


void printv(int n, double *x);
void dderF (int *Np, double *X, double *Y, 
                         int *nelt, int *ia, int *ja, double *a, int *isym){
// calculates dF/dx*p directly like dF/dp*|p|, p - required direction.                         
    MatVec_DDer *mv=(MatVec_DDer*)ia[0];     
// current point:
    double *xc=mv->xc();                // xc - current point; X - current step
    double *fc=mv->fc();                // function in current point
    double eta=mv->eta();
    NEF_TYPE_INX func=mv->nef();
    void *par=mv->p();
    double M1=-1.;
/*    
    double dx[*Np];
    double rd[*Np], ed[*Np];  // real difference, expected difference 
*/
    double *dx=dalloca(*Np);
    double *rd=dalloca(*Np);
    double *ed=dalloca(*Np);  // real difference, expected difference 

    double xnorm=DNRM2(Np,X,&ONE);
    double etadx=eta/xnorm;
    double xdeta=xnorm/eta;
    double etainv=1./eta;
    int retcode;
//cout<<" ===================== dder ==\n";
/* 
//cout<<" X: "; printv(*Np,X);
//cout<<" xc: "; printv(*Np,xc);
//cout<<" Y: "; printv(*Np,Y);
//cout<<" fc: "; printv(*Np,fc);
//cout<<"ia[0] %%%%%%%%%%%%%\t"<<ia[0]<<endl;
//cout<<"mv->eta() %%%%%%%%%\t"<<eta<<endl;
//cout<<"&par %%%%%%%%%%%%%%\t"<<int(par)<<endl;
//cout<<"ja[0] %%%%%%%%%%%%%\t"<<ja[0]<<endl;
cout<<"#"; flush(cout); // 
*/
    if (xnorm != 0) {
        DCOPY(Np,xc,&ONE,dx,&ONE);
//	cout<<" dx: xc copied: "; printv(*Np,xc);
//	cout<<"etadx = eta/||x||="<<etadx<<endl;
        DAXPY(Np,&etadx,X,&ONE,dx,&ONE);
//	cout<<" dx: eta/||X||*X added: "; printv(*Np,xc);
	

//        CLEAN(Np,ed);
//        CLEAN(Np,rd);
//        DAXPY(Np,&etadx,X,&ONE,ed,&ONE);   // expected difference = etadx*X
//        DCOPY(Np,dx,&ONE,rd,&ONE);         // real diff = dx - xc 
//        DAXPY(Np,&M1,xc,&ONE,rd,&ONE);     // 
//        cout<<" expected diff: "; printv(*Np,ed);
//        cout<<" real diff    : "; printv(*Np,rd);
//        cout<<" xc           : "; printv(*Np,xc);
                
        (*func)(Np,dx,Y,&retcode,(SolverExtraData*)par);
//	cout<<" f+ calculated:  "; printv(*Np,Y);	
        DAXPY(Np,&M1,fc,&ONE,Y,&ONE);
//	cout<<" f+ - fc  =   :  "; printv(*Np,Y);		
//        DSCAL(Np,&xdeta,Y,&ONE);
        DSCAL(Np,&xdeta,Y,&ONE);
//      cout<<"scaling 1/eta="<<xdeta<<endl;
//	cout<<" (f+-fc)/eta*|X|:"; printv(*Np,Y);		
        }
    else{
        for (int i=0; i<*Np; i++)
            Y[i]=0;
        } 
//
//cout<<" Y: "; printv(*Np,Y);
/********************************
            {  //            === via jacobian ===	    
//cout<<"jacobian:\n";
double dfdx[*Np**Np];
double sx[*Np];
double res[*Np], diff[*Np], reldiff[*Np];

for (int i=0; i<*Np; i++){
    sx[i]=1.;
    res[i]=0;
    diff[i]=0;
    }
double meps=eta*eta;
FDJAC(Np,xc,fc,func,par,sx,&meps,dfdx);
//for (int i=0; i<*Np; i++)
//    printv(*Np,&dfdx[i* *Np]);

for (int i=0; i<*Np; i++){
    double sum=0;
    for (int j=0; j<*Np; j++)
        sum+=dfdx[j+i**Np]*X[j];
    res[i]=sum;
    }            
            
//cout<<"Jacob*X(res)=";printv(*Np,res);
            
for (int i=0; i<*Np; i++){
    diff[i]=Y[i]-res[i];
    reldiff[i]=(Y[i]!=0?diff[i]/Y[i]:diff[i]);
    if (modify ==1   && (counter >8 && counter<12) )  {  // 
	Y[i]=res[i]; 
	}           
   }

//cout<<"***    DIFF *** :";printv(*Np,diff);
//cout<<"*** relDIFF *** :";printv(*Np,reldiff);
            }
*********************************************/

//cout<<" ----------------- dder  out -------\n";
//cout<<" X: "; printv(*Np,X);
//cout<<" xc: "; printv(*Np,xc);
//cout<<" Y: "; printv(*Np,Y);
//cout<<" fc: "; printv(*Np,fc);
//cout<<"ia[0] %%%%%%%%%%%%%"<<ia[0]<<endl;
//cout<<"mv->eta() %%%%%%%%%%%%%"<<eta<<endl;
//cout<<"&par %%%%%%%%%%%%%"<<int(par)<<endl;
//cout<<"ja[0] %%%%%%%%%%%%%"<<ja[0]<<endl;
    
    return;                 
}

// --- functions, which are used in PrecondDGMRES class -------------

void unitPR (int *Np, double *R, double *Z, int *nelt, int *ia, int *ja, 
             double *a, int *isym, double *rwork, int *iwork){
    PrecondDGMRES *prc=(PrecondDGMRES *)ja[0]; 
    DCOPY(Np,R,&ONE,Z,&ONE);
    return;
}


// --- functions, which are used in PrecondDGMRES class -------------

#include "vector_2d.hh"
void LUSolveBlockPreconditioner (int *N, double *R, double *Z, int *nelt, int *ia, int *ja, double *a, int *isym, double *rwork, int *iwork){
    LUConstBlockPrecond *prc=(LUConstBlockPrecond*) ja[0];
    int * ipiv=prc->ipiv();
    double * bj=prc->bj();
    int M=prc->M();                         // number of variables
    char trans='n';
    int info;
//    double y[*N];
    double *y=dalloca(*N);
    
    int nFreq=(*N)/(2*M);                   // number of frequencies
    vector_2d R2d(R,2*nFreq,*N);                // alternative addressing (freq,var), freq=1..2*nFreq
    vector_2d Z2d(Z,2*nFreq,*N);                // alternative addressing (freq,var)
    vector_2d y2d(y,M,*N);            // -"- (freq, var)
    
/*    cout<<"LUSolveBlockPreconditioner:\n"; 
/*    cout<<" \nR -=-=-=-=-=-=-=-=-=-===-\n";
      cout<<"N="<<*N<<" M="<<M<<" nFreq="<<nFreq<<endl;    printv(*N,R);
/**/
    
    for (int i=0; i<M; i++)              // for each frequency
        for (int j=0; j<nFreq; j++){
            y2d(i,2*j  )=R2d(2*j  ,i);
            y2d(i,2*j+1)=R2d(2*j+1,i);
            }            
/*    cout<<"R REORDERED: y -=-=-=-=-=-=-=-=-=-===-\n";    printv(*N,y); */

    int NRHS=2*nFreq;

    DGETRS(&trans, /*int *TRANS_len, */ &M, &NRHS, bj, &M, ipiv, y/*double *B*/, &M /*int *LDB*/, &info);
/*    cout<<"y solved:  -=-=-=-=-=-=-=-=-=-===-\n";    printv(*N,y);*/
    for (int i=0; i<M; i++)              // for each frequency
        for (int j=0; j<nFreq; j++){
            Z2d(2*j  ,i)= y2d(i,2*j);
            Z2d(2*j+1,i)= y2d(i,2*j+1);
            }
/*    cout<<"y reordered back: Z -=-=-=-=-=-=-=-=-=-===-\n";    printv(*N,Z); /**/
    return;
}

int LUConstBlockPrecond::create(){
        // read it
/*        ifstream in("JACOBIAN");
        int tmp1;
        double tmp2;
        in>>tmp1;
        for (int i=0; i<2*M_; i++) // skip last one 
            for( int j=0; j<2*M_; j++){
                in >> tmp2;
                if (i%2 ==0 && j%2 == 0)
                    bj_[M_*i/2+j/2]=tmp2;           // each first (bypass second)
                }
        in.close();        
/**/
//      cout<<"LUConstBlockPrecond starts ************* jac_="<<jac_<<endl;;

        for (int i=0; i<M_; i++) // skip last one 
            for( int j=0; j<M_; j++)
                bj_[M_*i+j]=jac_[2*M_*2*i+2*j];

/*      cout<<"is read ************* \n";
        for (int i=0; i<M_; i++){
            for (int j=0; j<M_; j++)
                cout<<bj_[j*M_+i]<<" ";
            cout<<endl;    
            }
/**/            
        int info;    
        DGETRF(&M_,&M_,bj_,&M_,ipiv_,&info);  // fortran
/*      cout<<"factorizED ************* \n";
        for (int i=0; i<M_; i++){
            for (int j=0; j<M_; j++)
                cout<<bj_[j*M_+i]<<" ";
            cout<<endl;    
            }        
/**/
        return 1;
}

int LUVariableBlockPrecond::update(double *x, double *f,NEF_TYPE_INX NEF,void *pp,double *sx,double teta){     // update preconditioner 
//        cout<<"var precond update \n";
        int NN=N();
        int MM=M();
        //cout<<"var precond update : calling FDJACDC, MM="<<MM<<endl;;
        
        FDJACDC(&NN,x, f,NEF,pp,sx, &teta,bj(),&MM);
        
/*        cout<<"is built  ************* \n";
        for (int i=0; i<M_; i++){
            for (int j=0; j<M_; j++)
                cout<<bj_[j*M_+i]<<" ";
            cout<<endl;    
        }
/**/        
        int info;
        //cout<<"var precond update : calling DGETRF\n";

        DGETRF(&M_,&M_,bj_,&M_,ipiv_,&info);      
/*        
         cout<<"factorizED ************* \n";
        for (int i=0; i<M_; i++){
            for (int j=0; j<M_; j++)
                cout<<bj_[j*M_+i]<<" ";
            cout<<endl;    
            }        
/**/        
        return 1;
}
