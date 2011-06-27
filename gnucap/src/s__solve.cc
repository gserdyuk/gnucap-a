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
#include "e_cardlist.h"
#include "u_status.h"
#include "e_node.h"
#include "s__.h"
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
static bool converged = false;
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
void print_rhs(CKT_BASE* s, std::string a){
  std::cout<<a<<"\n";
  std::cout<<"    == _sim->_i (rhs): \n";
  for (int ii = 1; ii<=s->_sim->_lu.size(); ++ii) 
    std::cout<<"     _sim->_i["<<ii<<"]="<<s->_sim->_i[ii]<<"\n";
}
void print_sol(CKT_BASE* s, std::string a){
  std::cout<<a<<"\n";
  std::cout<<"    == _sim->_v0 (sol): \n";
  for (int ii = 1; ii<=s->_sim->_lu.size(); ++ii) 
    std::cout<<"     _sim->_v0["<<ii<<"]="<<s->_sim->_v0[ii]<<"\n";
}
void print_matr(SIM* s, std::string a){
  std::cout<<a<<"\n";
  std::cout<<"    == _sim->_aa: \n";
  (s->_sim->_aa).printm();
}
/*--------------------------------------------------------------------------*/
void error_func(SIM *s, F){
    int N = s->_sim->_total_nodes;
    
    std::fill_n(F, N+1, 0);

    for (int j=0; j<=N; j++)
        F[j]=-s->_sim->_i[j];
        
    for (int j=0; j<=N; j++){
        double sum=0;
        for (int k=0; k<=N; k++)
            sum+= s->_sim->_aa.s(j,k) * s->_sim->_v0[k];
        F[j]+=sum;
        }

/*
    std::cout<<"\n";
    std::cout<<"    == F : \n";
    for (int k = 0; k<=N; k++) 
      std::cout<<"     F["<<k<<"]="<<F[k]<<"\n";
*/
        return;    
}

/*--------------------------------------------------------------------------*/
bool SIM::solve(OPT::ITL itl, TRACE trace)
{

  double *v0_prev;  // dc-tran voltage at previous iteration, at the end of ityeration _v0_prev[..]=_sim->_v0[..]  
  double *f, *fn;   // dc-tran error function, previous and new. 

  // GS - additional array
  v0_prev = new double[_total_nodes+1];  
  std::fill_n(v0_prev,_total_nodes+1, 0);
  
  fn = new double[_total_nodes+1];  
  std::fill_n(fn,_total_nodes+1, 0);

  f = new double[_total_nodes+1];  
  std::fill_n(f,_total_nodes+1, 0);
  

  converged = false;
  int convergedcount = 0;
  
  _sim->reset_iteration_counter(iSTEP);
  advance_time();

  _sim->_damp = 1;   // not used now
  
  //GS  double epsim = macheps()     - compute machine eps
  //GS int termcode = ne_input_check(N, epsim, irnag, sf, sx, U, epssol, epsdu, epsmin, maxdu, limit)
  int tremcode = 0;
  
  int icode_return=termcode;
  
  //GS if  termcode <0 - stop
  
 
  do{
    std::cout<<"==== loop beginning,  iter="<<_sim->iteration_number()<<"\n";
    if (trace >= tITERATION) {
      print_results(static_cast<double>(-_sim->iteration_number()));
    }
    set_flags();            // GS check all conditions
    clear_arrays();
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
    if (!converged || !OPT::fbbypass || _sim->_damp < .99) {
      //std::cout<<" === solve: set_damp \n";

      //set_damp();
      _sim->damp=1;

      //std::cout<<" === solve: load_matrix \n";
      load_matrix();

      //print_rhs(this," load_matrix_print");
      //print_sol(this," load_matrix_print");
      //print_matr(this," load_matrix_print");

      error_func(this, fn);
      
      // calculate gradient  TODO

      // keep current point
      for (int i=0; i<_sim->_total_nodes+1; i++)
        v0_prev[i]=_sim->_v0[i];

      // keep current func
      for (int i=0; i<_sim->_total_nodes+1; i++)
        fn[i]=f[i];

      std::cout<<" === solve: solve_equations \n";
      solve_equations();

      print_sol(this," solution after linear solver");
      // now we have X_next
      // ...    
      
      // set voltage dump 
      double d1=1;
      if (_sim->iteration_number()==2 && false)   // fixing at the end of iteration 2 
        d1=0.5;
      
      for (int i=0; i<_sim->_total_nodes+1; i++)
        _sim->_v0[i] = _sim->_v0[i] *d1 + _sim->_v0_prev[i] * (1.-d1);

      print_sol(this," solution after damping");
    
    }
  }while (!converged && !_sim->exceeds_iteration_limit(itl));

  return converged;
}
/*--------------------------------------------------------------------------*/
bool SIM::solve_with_homotopy(OPT::ITL itl, TRACE trace)
{
  solve(itl, trace);
  //GS trace2("plain", ::status.iter[iSTEP], OPT::gmin);
  if (!converged && OPT::itl[OPT::SSTEP] > 0) {
    int save_itermin = OPT::itermin;
    OPT::itermin = 0;
    double save_gmin = OPT::gmin;
    OPT::gmin = 1;
    while (_sim->_iter[iPRINTSTEP] < OPT::itl[OPT::SSTEP] && OPT::gmin > save_gmin) {
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
  }else if (_sim->inc_mode_is_bad()) {
    _sim->set_inc_mode_no();
  }else if (_sim->is_iteration_number(OPT::itl[OPT::TRLOW])) {
    _sim->set_inc_mode_no();
  }else if (_sim->is_iteration_number(0)) {
    // leave it as is
  }else{
    _sim->set_inc_mode_yes();
  }

  _sim->_bypass_ok = 
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
  print_rhs(this);

  ::status.load.start();
  if (OPT::traceload && _sim->is_inc_mode()) {
    std::cout<<"    =# work with loadq \n";
    while (!_sim->_loadq.empty()) {
       std::cout<<"    =# device queue="<< _sim->_loadq.back()->long_label()<<" ->tr_load()\n";
      _sim->_loadq.back()->tr_load();
       print_rhs(this);
      _sim->_loadq.pop_back();
    }
  }else{
    std::cout<<"    =# work with all.tr_load \n";
    _sim->_loadq.clear();
    CARD_LIST::card_list.tr_load();
  }
  ::status.load.stop();
}
/*--------------------------------------------------------------------------*/
void SIM::solve_equations()
{
  ::status.lud.start();
  _sim->_lu.lu_decomp(_sim->_aa, bool(OPT::lubypass && _sim->is_inc_mode()));
  ::status.lud.stop();

  ::status.back.start();
  _sim->_lu.fbsub(_sim->_v0, _sim->_i, _sim->_v0);
  ::status.back.stop();
  
  std::cout<<"    == in solve_equations: \n";
  print_rhs(this);
  print_sol(this);
   
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
