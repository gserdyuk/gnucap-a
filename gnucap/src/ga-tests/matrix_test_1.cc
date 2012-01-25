#include "m_matrix.h"

void print_matr(BSMATRIX<double>* a, std::string a){
  std::cout<<a<<"\n";
  a->printm();
}

void print_vect(int N, double* x, std::string a){
  std::cout<<a<<"\n";
  for (int ii = 0; ii<N; ii++) 
    std::cout<<"     x["<<ii<<"]="<<x[ii]<<"\n";
}
/


BOOST_AUTO_TEST_CASE(ap_test1){

BSMATRIX<double> _aa;
BSMATRIX<double> _lu;
double * X;
double * F;
int N = 3;
double gmin = 1.e-20;

// allocate _aa and _lu
_aa.init(N);

_aa.iwant(1,1);
_aa.iwant(2,2);
_aa.iwant(2,1);
_aa.iwant(1,2);

_aa.allocate();

_aa.dezero(gmin);
_aa.setminpivot(gmin);
 
_aa.m(1,1) = 1.e-12;
_aa.m(2,2) = -10;
_aa.m(1,2) = -1;
_aa.m(2,1) =  10;

//... and _lu
_lu.clone(_aa);
_lu.allocate();

// tbd
print_matr(_aa," _aa filled:");

BOOST_CHECK(1==1);

}
