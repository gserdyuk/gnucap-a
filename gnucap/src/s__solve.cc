/*$Id: s__solve.cc,v 26.133 2009/11/26 04:58:04 al Exp $ -*- C++ -*-
 * Copyright (C) 2001 Albert Davis
 * Author: Albert Davis <aldavis@gnu.org>
 *
 * This file is part of "Gnucap", the Gnu Circuit Analysis Package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 *------------------------------------------------------------------
 * solve one step of a transient or dc analysis
 */
//testing=script 2006.07.14
#include <vector>
#include <gsl/gsl_cblas.h>
#include "s__.h"
#include "e_cardlist.h"
#include "u_status.h"
#include "e_node.h"

/*--------------------------------------------------------------------------*/
//	bool	SIM::solve(int,int);
//	void	SIM::finish_building_evalq();
//	void	SIM::advance_time();
//	void	SIM::set_flags();
//	void	SIM::clear_arrays();
//	void	SIM::evaluate_models();
//	void	SIM::set_damp();
//	void	SIM::load_matrix();
//	void	SIM::solve_equations();
/*--------------------------------------------------------------------------*/
static bool converged = false;     // GS - refactor ? sttaic is a bad style
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
void print_rhs(CKT_BASE* s, std::string a){
  std::cout<<a<<"\n";
  std::cout<<"    == _sim->_i (rhs): \n";
  int N = s ->_sim->_total_nodes+1;
  for (int ii = 0; ii<N; ++ii) 
    std::cout<<"     _sim->_i["<<ii<<"]="<<s->_sim->_i[ii]<<"\n";
}
void print_sol(CKT_BASE* s, std::string a){
  std::cout<<a<<"\n";
  std::cout<<"    == _sim->_v0 (sol): \n";
  int N = s ->_sim->_total_nodes+1;
  for (int ii = 0; ii<N; ++ii) 
    std::cout<<"     _sim->_v0["<<ii<<"]="<<s->_sim->_v0[ii]<<"\n";
}
void print_matr(SIM* s, std::string a){
  std::cout<<a<<"\n";
  std::cout<<"    == _sim->_aa: \n";
  (s->_sim->_aa).printm();
}
void print_matr_lu(SIM* s, std::string a){
  std::cout<<a<<"\n";
  std::cout<<"    == _sim->_lu: \n";
  (s->_sim->_lu).printm();
}
void print_vect(int N, double* x, std::string a){
  std::cout<<a<<"\n";
  for (int ii = 0; ii<N; ii++) 
    std::cout<<"     x["<<ii<<"]="<<x[ii]<<"\n";
}
/*--------------------------------------------------------------------------*/
/*
 friend int conv_check_4stop(int* N, double* xc, double* xp, double* fplus, double* fplus_norm, 
                              double* grad, double* sx, double* sf,
                              int* retcode,
                              double* fvectol, double* steptol, int* itncount, 
                              int *itnlimit, 
                              bool *maxtaken, 
                              int* analjac,
                              int * analjac,
                              int*cheapF,
                              double* mintol)
{ 
   int termcode=0;   
   //
   //0 -  nothing happened
   //1 -  scaled norm of function is < of fvectol
   //2 -  scaled norm of step is smaller than steptol - probably solution
   //3 -  at last globalization step can not decrease |F| significantly. may be bad jacobian
   //4 -  iteration limit is exceeded
   //5 -  5 steps MAXSTEP is made - likley function asymptotivcally approaches to end value
   //6 -  xc probably is local minimum, which is not root, or mintol is too small
    
    if (*retcode == 1 ) {
        teremcode = 3;
        }
    else if (                          
                              
}
*/

extern "C"  {

void stop_(int* N, double* X, double *DX, double* F, double* FNOR, 
        double* GR, double* SX, double* SF,
        int* IRETCD,int* ITER, int* MAXTKN,
        int * KMAXDU, int* TERMCD, 
        double* EPSSOL, double* EPSDU, double* EPSMIN, double* MAXDU,
        int* LIMIT);        
}

/*------------------------------------------------------------------------- */
/*
J = dN/dX + Y
Fg= dN/dX * X - N - I

F = N + I + Y*Xc 

F = - (dN/dX * X - N - I) + ( dN/dX + Y ) * X = 
    - dN/dX * X  + N+ I + dN/dX * X + Y *X = 
    = Y * X + N + I

F = J* Xc - Fg

*/
void calc_error_func(SIM *s, double* F){
    int N = s->_sim->_total_nodes+1;
    
    std::fill_n(F, N, 0);
    
    std::cout<<" --- ################ calc error func \n";
    std::cout<<" F cleaned: \n";
    for (int k = 0; k<N; k++) 
      std::cout<<"     F["<<k<<"]="<<F[k]<<"\n";

    std::cout<<" V0 : \n";
    for (int k = 0; k<N; k++) 
      std::cout<<"     _v0["<<k<<"]="<<s->_sim->_v0[k]<<"\n";

        
    for (int j=0; j<N; j++){
        double sum=0;
        for (int k=0; k<N; k++)
            sum+= s->_sim->_aa.s(j,k) * s->_sim->_v0[k];
        F[j]+=sum;
        }

    std::cout<<" F =J*Xc: \n";
    for (int k = 0; k<N; k++) 
      std::cout<<"     F["<<k<<"]="<<F[k]<<"\n";

    for (int j=0; j<N; j++)
        F[j] -= s->_sim->_i[j];

    std::cout<<" F =J*Xc - Fg: \n";
    for (int k = 0; k<N; k++) 
      std::cout<<"     F["<<k<<"]="<<F[k]<<"\n";

    std::cout<<" === ################ END calc error func \n";
/*
    std::cout<<"\n";
    std::cout<<"    == F : \n";
    for (int k = 0; k<=N; k++) 
      std::cout<<"     F["<<k<<"]="<<F[k]<<"\n";
*/
        return;    
}

/*--------------------------------------------------------------------------*/

typedef void (*ErrorFunction_TYPE) (int* N, double* XN, double* FN, int* IFLAG, double* SF,
                        double* FNOR1, void* ADD_DATA);
                        

void calc_circuit (int* nn, double* x, double* f, int* iflag, double* sf,
                        double* fnor, void* add_data){
    std::cout<<" ##++++++++++++++++++++++++++++++++++ calc_circuit entering +\n";
    int  N=*nn;
    for (int i=0; i<N; i++)
        std::cout<<" x["<<i<<"]="<<x[i]<<"\n";
    std::cout<<"--\n";
    SIM  *simobject = (SIM*) add_data; 
                        
    std::cout<<"    ---------- copy X to _sim->_v0 \n";    
    cblas_dcopy(N, x , 1, simobject->_sim->_v0, 1);
    
    // TMP - clean _sim->_i and _aa
     //std::fill_n(simobject->_sim->_i,N, 0);
     simobject->clear_arrays();
         
    // simobject->finish_building_evalq();   ?
    std::cout<<"    ---------- evaluate_models \n";    
    simobject->evaluate_models();
    std::cout<<"    ---------- load_matrix \n";
    simobject->load_matrix();
    
    print_rhs(simobject, " matrix loaded");
    print_matr(simobject, " matrix loaded");
    
    
    std::cout<<"    ---------- calc_error_func \\n";
    calc_error_func(simobject, f);
    std::cout<<"    ---------- results f[] \n";
    for (int i=0; i<N; i++)
        std::cout<<" f["<<i<<"]="<<f[i]<<"\n";
    std::cout<<"--\n";   
    
    *iflag=0;

    for (int i=0; i<N; i++)         // scale
        f[i]*=sf[i];
        
    *fnor= cblas_ddot(N,f,1,f,1) / 2;  // fnor = F^2 / 2; fnor'=2F/2=F
    
    std::cout<<" ##================================== calc_circuit returning - fnor="<<*fnor<<"\n";
    
 }
 
extern "C"  {

void lsearch_(
        int* N, double* X, double* FNOR, 
        double* GR, double* Y, double* SX, double* SF,
        int* IRETCD,int* MAXTKN,
        double* XN, double* FN, 
        double* FNORN, 
        double* TLS, 
        ErrorFunction_TYPE FUNCT,
        double* EPSSOL, double* EPSDU, double* EPSMIN, double* MAXDU,
        int* LIMIT, 
        void* add_data);        
}


bool SIM::solve(OPT::ITL itl, TRACE trace)
{

  std::cout<<"==== SIM::solve/ linserach_solve \n";
  
  int N = this->_sim->_total_nodes+1;

  double *Xn = this->_sim->_v0;  // contains X, after solve_equations() contain Xn

  double *X;  // dc-tran voltage at the beginning of current iteration, 
  // GS - additional array
  X = new double[N];            
  std::fill_n(X,N, 0);
  
  double *f, *fn;   // dc-tran error function, previous and new. 
  fn = new double[N];              // F "new"
  std::fill_n(fn,N, 0);
  f = new double[N];               // F 
  std::fill_n(f,N, 0);

  double *sx, *sf;   // scaling vectors
  sx = new double[N];      // scale X: so far use 1 as scale. not tested with !=1  
  std::fill_n(sx,N, 1);
  sf = new double[N];      // scale F
  std::fill_n(sf,N, 1);
  
  double *gr;                           // gradient
  gr = new double[N];      
  std::fill_n(gr,N, 0);

  double *Y;                           // newton step
  Y = new double[N];      // step value = v0-v0_prev 
  std::fill_n(Y,N, 0);

  double lambda=1.;

// TBD
  double epssol = 1.e-10;       // OPT::abstol;
  double epsdu  = OPT::vntol;
  double maxdu  = 1.e12;
  double epsmin=1.e-17;         // 1.e-15;   - converged good, check for external solver
  int limit = 100;
  int retcode = 0;
  int maxtaken= 0;
  double fnor, fnor_n;    
  void *add_data=0;
  bool shall_stop;
  
  converged = false;                // GS old converged criterion - all devices are converged
  bool converged_iter=false;        // GS new converged ccritereion - system is converged
  int convergedcount = 0;
  
  _sim->reset_iteration_counter(iSTEP);
  advance_time();

  _sim->_damp = 1;   // not used now
  
  //GS  double epsim = macheps()     - compute machine eps
  //GS int termcode = ne_input_check(N, epsim, irnag, sf, sx, U, epssol, epsdu, epsmin, maxdu, limit)
  int termcode = 0;  
  //GS if  termcode <0 - stop
  
 
  do{
    std::cout<<"==== loop beginning,  iter="<<_sim->iteration_number()<<"\n";
    if (trace >= tITERATION) {
      print_results(static_cast<double>(-_sim->iteration_number()));
    }
    set_flags();            // GS check all conditions
    // clear_arrays(); - called in calc_circuit
    // clear local arrays TODO
    
    finish_building_evalq();
    
    // TODO - check what is this
    _sim->count_iterations(iPRINTSTEP);
    _sim->count_iterations(iSTEP);
    _sim->count_iterations(_sim->_mode);
    _sim->count_iterations(iTOTAL);
    
    //std::cout<<" === solve: evaluate_models \n";
    evaluate_models();
    //print_rhs(this," evaluate_models");
    //print_sol(this," evaluate_models");
    //print_matr(this," evaluate_models");

    // this poubt "converge==true" means that devices are converged


 /* TODO : GS - double check logic  
 
    converged       - 
    _sim->limiting  -
    convergedcount  -
       
    if (converged) {
      if (_sim->_limiting) {
	error(bDEBUG, "converged beyond limit, resetting limit\n");
	_sim->set_limit();
	convergedcount = 0;
      }else{
	++convergedcount;
      }
    }else{
      convergedcount = 0;
    }
    if (convergedcount <= OPT::itermin) {
      converged = false;
    }
 
 */      
    if (!converged_iter || !OPT::fbbypass || _sim->_damp < .99) {
      //std::cout<<" === solve: set_damp \n";

      //set_damp();
      _sim->_damp=1;

      //std::cout<<" === solve: load_matrix \n";
      load_matrix();

      print_rhs(this," load_matrix_print");
      print_sol(this," load_matrix_print");
      print_matr(this," load_matrix_print");

      calc_error_func(this, fn);      
      print_vect(N,fn,"*** FN= \n");
      
      fnor=cblas_ddot(N, fn, 1,fn,1)/2.;
      std::cout<<" *** fnor="<<fnor<<"\n";
      
      // calculate gradient  TODO
      //GRADIE(IRANG,N,DFDX,FN,SF,GR)
      this->_sim->_aa.gradient(fn,sf,gr);  
      print_vect(N,gr,"*** GR= \n");
    

      // keep current point      
      //cblas_daxpy(N, 1 , Xn, 1, X, 1);      // store current vector v0 in X 
      cblas_dcopy(N, Xn, 1, X, 1);      // store current vector v0 in X 
      
      // keep current func
      //cblas_daxpy(N, 1, fn, 1, f,1);         // store current vector f       
      cblas_dcopy(N, fn, 1, f,1);         // store current vector f       
      

      std::cout<<" === solve: solve_equations \n";
      solve_equations();

      print_sol(this," solution after linear solver");
      // now we have X_next
      
      // calculate  X next = X prev - Y; y = x prev - x next
      std::cout<<" *** calculating Y: \n";
      cblas_dcopy(N,  X,  1, Y, 1);
      print_vect(N,Y,"  X  = \n");
      print_vect(N,Xn," Xn = \n");
      cblas_daxpy(N, -1, Xn, 1, Y, 1);       // Y is full newton step now
                
      print_vect(N,Y,"*** X-Xn= Y= \n");

                        
      /* 
      N - defined           X - defined
      FNOR   - defined
      GR - gr - defiend     Y - defined
      SX = sx  defined      SF = sf - defined
      IRETCD,
      MAXTKN,
      XN = Xn - defiend     FN = fn - defined
      FNORN         - defined
      lambda,       - defined 
      FUNCT,
      EPSSOL        - abstol
      EPSDU         - vntol
      EPSMIN        - defined
      MAXDU         - defined
      LIMIT         - defined
      add_data
      */
      
      // linear search
      //LSEARCH(N,U,FNOR,GR,Y,SX,SF,IRETCD,MAXTKN,UN,FN,FNORN,TLS, FUNCT,EPSSOL,EPSDU,EPSMIN,MAXDU,LIMIT, add_data)
      
      lsearch_(&N,X,&fnor,gr,Y,sx,sf,
                &retcode,&maxtaken,
                Xn,fn,
                &fnor_n,&lambda, 
                &calc_circuit,
                &epssol,&epsdu,&epsmin,&maxdu,&limit, 
                (void*)this);
        
      std::cout<<" after lin search - lambda="<<lambda<<" retcode = "<<retcode<<"\n";  
      
      print_sol(this," solution after damping");
    
      // check final stop condition   
      int iterno=_sim->iteration_number();
      int maxtkn=0;
      int kmaxdu=5;
      stop_(&N, X, Y, fn, &fnor, 
            gr, sx, sf, &retcode, &iterno, &maxtkn,
            &kmaxdu, &termcode, 
            &epssol, &epsdu, &epsmin, &maxdu, &limit);        
      std::cout<<" ---- termocode ="<<termcode<<"\n";

      converged_iter  = (termcode == 1);
      shall_stop = converged_iter  || _sim->exceeds_iteration_limit(itl) || (termcode != 0)  ;
      std::cout<<"  converged_iter = "<<converged_iter<<" shall_stop ="<<shall_stop<<"\n";
      
      
    }
  }while (!shall_stop );

  return converged=converged_iter ;   // to pass outside "converged" for solve_with_homotopy
}
/*--------------------------------------------------------------------------*/
bool SIM::solve_with_homotopy(OPT::ITL itl, TRACE trace)
{
  solve(itl, trace);
  return converged;
  //GS - so far will not allow homotopy for this solver - untill understand when is it possible 
  //GS trace2("plain", ::status.iter[iSTEP], OPT::gmin);
  if (!converged && OPT::itl[OPT::SSTEP] > 0) {
    int save_itermin = OPT::itermin;
    OPT::itermin = 0;
    double save_gmin = OPT::gmin;
    OPT::gmin = 1;
    while (_sim->_iter[iPRINTSTEP] < OPT::itl[OPT::SSTEP] && OPT::gmin > save_gmin) {
      std::cout<<"while in solve_with_homotory\n";
      //CARD_LIST::card_list.precalc();
      _sim->set_inc_mode_no();
      solve(itl, trace);
      if (!converged) {
	trace2("fail", _sim->_iter[iSTEP], OPT::gmin);
	OPT::gmin *= 3.5;
      }else{
	trace2("success", _sim->_iter[iSTEP], OPT::gmin);
	OPT::gmin /= 4;
      }
    }
    OPT::itermin = save_itermin;
    OPT::gmin = save_gmin;
    //CARD_LIST::card_list.precalc();
    solve(itl, trace);
    if (!converged) {
      trace2("final fail", _sim->_iter[iSTEP], OPT::gmin);
    }else{
      trace2("final success", _sim->_iter[iSTEP], OPT::gmin);
    }
  }else{
  }
  return converged;
}
/*--------------------------------------------------------------------------*/
/* finish_building_evalq
 * This function scans the circuit to queue for evaluation anything 
 * that is relevant that the devices missed themselves.
 * For now, it scans the whole circuit, but this will change.
 * This is slow, but the result is still faster than not having a queue.
 * The plan is .. every node knows whether it needs eval or not, and
 * only those nodes needing eval will be scanned.
 * Its purpose is to catch nodes that wake up after being dormant
 */
void SIM::finish_building_evalq(void)
{
  ::status.queue.start();
  CARD_LIST::card_list.tr_queue_eval();
  ::status.queue.stop();
}
/*--------------------------------------------------------------------------*/
void SIM::advance_time(void)
{
  ::status.advance.start();
  static double last_iter_time;
  if (_sim->_time0 > 0) {
    if (_sim->_time0 > last_iter_time) {	/* moving forward */
      notstd::copy_n(_sim->_v0, _sim->_total_nodes+1, _sim->_vt1);
      CARD_LIST::card_list.tr_advance();
    }else{				/* moving backward */
      /* don't save voltages.  They're wrong! */
      /* instead, restore a clean start for iteration */
      notstd::copy_n(_sim->_vt1, _sim->_total_nodes+1, _sim->_v0);
      CARD_LIST::card_list.tr_regress();
    }
  }else{
    CARD_LIST::card_list.dc_advance();
  }
  last_iter_time = _sim->_time0;
  ::status.advance.stop();
}
/* last_iter_time is initially 0 by C definition.
 * On subsequent runs it will start with an arbitrary positive value.
 * _sim->_time0 starts at either 0 or the ending time of the last run.
 * In either case, (time0 > last_iter_time) is false on the first step.
 * This correctly results in "don't save voltages..."
 */
/*--------------------------------------------------------------------------*/
void SIM::set_flags()
{
  _sim->_limiting = false;
  _sim->_fulldamp = false;
  
  if (OPT::incmode == false) {
    _sim->set_inc_mode_no();
  }else if (_sim->inc_mode_is_bad()) {                            // TBD - check this
    _sim->set_inc_mode_no();
  }else if (_sim->is_iteration_number(OPT::itl[OPT::TRLOW])) {    // TBD - check this
    _sim->set_inc_mode_no();
  }else if (_sim->is_iteration_number(0)) {                       // TBD  - check this
    // leave it as is
  }else{
    _sim->set_inc_mode_yes();
  }

  _sim->_bypass_ok =                                                // TBD - check this 
    (is_step_rejected()  ||  _sim->_damp < OPT::dampmax*OPT::dampmax)
    ? false : bool(OPT::bypass);
    
}
/*--------------------------------------------------------------------------*/
void SIM::clear_arrays(void)
{
  if (!_sim->is_inc_mode()) {			/* Clear working array */
    _sim->_aa.zero();
    _sim->_aa.dezero(OPT::gmin);		/* gmin fudge */
    std::fill_n(_sim->_i, _sim->_aa.size()+1, 0);
  }
  _sim->_loadq.clear();
}
/*--------------------------------------------------------------------------*/
void SIM::evaluate_models()
{
  ::status.evaluate.start();
  if (OPT::bypass) {
    std::cout<<"    =* OPT::bypass == true\n";  
    converged = true;
    std::cout<<"    =* converged="<< converged <<"\n";
    swap(_sim->_evalq, _sim->_evalq_uc);
    while (!_sim->_evalq->empty()) {
    std::cout<<"    =* device queue1="<< _sim->_evalq->front()->long_label()<<" ->do_tr\n";
      converged &= _sim->_evalq->front()->do_tr();
    std::cout<<"    =* converged="<< converged <<"\n";
      _sim->_evalq->pop_front();
    }
  }else{
     std::cout<<"    =* OPT::bypass != true\n";  
    _sim->_evalq_uc->clear();
    converged = CARD_LIST::card_list.do_tr();
  }
  while (!_sim->_late_evalq.empty()) { //BUG// encapsulation violation
    std::cout<<"    =1 device queue2="<< _sim->_evalq->front()->long_label()<<"\n";
    converged &= _sim->_late_evalq.front()->do_tr_last();
    std::cout<<"    =1 converged="<< converged <<"\n";
    _sim->_late_evalq.pop_front();
  }
  ::status.evaluate.stop();
}
/*--------------------------------------------------------------------------*/
void SIM::set_damp()
{
  if (_sim->is_second_iteration() && !converged && OPT::dampstrategy&dsINIT) {
    _sim->_damp = OPT::dampmin;
  }else if (_sim->is_first_iteration()  ||  converged) {
    _sim->_damp = OPT::dampmax;
  }else if (_sim->_fulldamp) {
    _sim->_damp = OPT::dampmin;
  }else{
    _sim->_damp = OPT::dampmax;
  }
  
    // GS
  if (_sim->iteration_number() == 3 && false /**/) {
     std::cout<<" setting damp to 0.5 \n";
     _sim->_damp=0.5;
     }

  trace1("", _sim->_damp);
   std::cout<<" set_damp;  _sim->_damp="<< _sim->_damp<<"\n";
}
/*--------------------------------------------------------------------------*/
void SIM::load_matrix()
{
  std::cout<<"  ===# load_matrix entered \n";
  print_rhs(this, "point1");

  ::status.load.start();
  if (OPT::traceload && _sim->is_inc_mode()) {
    std::cout<<"    =# work with loadq \n";
    while (!_sim->_loadq.empty()) {
       std::cout<<"    =# device queue="<< _sim->_loadq.back()->long_label()<<" ->tr_load()\n";
      _sim->_loadq.back()->tr_load();
       print_rhs(this, "point2");
      _sim->_loadq.pop_back();
    }
  }else{
    std::cout<<"    =# work with ARD_LIST::card_list.tr_load() \n";
    _sim->_loadq.clear();
    CARD_LIST::card_list.tr_load();
  }
    print_rhs(this, "point \"load_matrix exit\"");

  ::status.load.stop();
}
/*--------------------------------------------------------------------------*/
void SIM::solve_equations()
{
  ::status.lud.start();
  print_matr_lu(this, "@@@ solve_equations, _lu before lu_decomp" );
  print_matr   (this, "@@@ solve_equations, _aa before lu_decomp" );
  
  _sim->_lu.lu_decomp(_sim->_aa, bool(OPT::lubypass && _sim->is_inc_mode()));

  print_matr_lu(this, "@@@ solve_equations, _lu AFTER lu_decomp" );
  print_matr   (this, "@@@ solve_equations, _aa AFTER lu_decomp" );
  ::status.lud.stop();

  ::status.back.start();
  print_rhs(this, "@@@ solve_equations, _sim->_i  before lu_fbsub");
  print_sol(this, "@@@ solve_equations, _sim->_v0 before lu_fbsub");
  
  _sim->_lu.fbsub(_sim->_v0, _sim->_i, _sim->_v0);

  print_rhs(this, "@@@ solve_equations, _sim->_i  AFTER lu_fbsub");
  print_sol(this, "@@@ solve_equations, _sim->_v0 AFTER lu_fbsub");

  ::status.back.stop();
  
  std::cout<<"    == in solve_equations: \n";
  print_rhs(this, "point3");
  print_sol(this, "point3");
   
  if (_sim->_nstat) {
    // mixed mode
    for (int ii = _sim->_lu.size(); ii >= 1; --ii) {
      _sim->_nstat[ii].set_a_iter();
    }
  }else{
    // pure analog
    untested();
  }
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
