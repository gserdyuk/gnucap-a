/*$Id$ -*- C++ -*-
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
 * behavioral modeling
 *  hspice/msim/whatever compatible "semiconductor resistor and capacitor"
 *
 * modified by G Serdyuk gena@analograils.com
 * 
 * Resistor description
 * --------------------
 * Rxxx n1 n2 <mname> <R=> resistance <<TC1=>val> <<TC2=>val>
 *            <SCALE=val> <M=val> <AC=val> <DTEMP=val> <L=val> <W=val> <C=val>
 *
 * .model mname R keywd=value <CRATIO=val>
 *
 *  keywd    UNITS   default   impl   descr
 *  BULK             gnd       n      default node for cap
 *  CAP      F       0         n      default capcaitance
 *  CAPSW    F/m     0         n      sidewall fringing cap
 *  COX      F/m^2   0         n      bottomwall cap
 *  DI               0         n      relative dielectric const
 *  DLR      m       0         n      Diff between drawn and actual length DLReff=DLR*SCALM
 *  DW       m       0         n      difference between drawn and actual width DWeff=DW*SCALM
 *  L        m       0         y      defualt length of wire; Lscaled=L*SRINK*SCALM
 *  LEVEL                      n      model selector (not used)
 *  RAC      Ohm               n      defualt AC redsistance (RACeff=Reff)
 *  RES      Ohm     0         n      default resistance
 *  RSH              0         y      sheet resistance
 *  SHRINK           1         n      shrink factor
 *  TC1C     1/deg   0         n      1st order tem coef for cap
 *  TC2C     1/deg^2 0         n      2nd order tem coef for cap
 *  TC1R     1/deg   0         y      1st order tem coef for res
 *  TC2R     1/deg^2 0         y      2nd order tem coef for res
 *  THICK    m       0         n      diel thinkness
 *  TREF     C       TNOM      y      reference temperature
 *  W        m       0         y      Default width Wscaled=W*shrink*scalm
 * 
 * 
 *  Weff = Wscaled - 2 * DWeff
 *  Leff = Lscaled - 2 * DLReff
 * if resisatnce is specififed:   Reff= R * scale (element) / M
 * if (weff * Leff * RSH ) >0:
 *  Reff = Leff * RSH * SCALE (element)/ ( M * Weff)
 * if (weff * Leff * RSH ) ==0:
 *  Reff = RES * SCALE (element)/  M 
 * if AC is specified in element:
 *  RACeff= AC*SCALE(element)/M
 * otehrwise, if AC is specified in model:
 *  RACeff=RAC*SCALE(eleemnt)/M
 * if neither:
 *  RACeff=Reff
 * if resistance is less then RESMIN, itis reset to RESMIN and warning is displayed
 *  RESMIN=1/ ( GMAX * 1000 * M )
 *
 *  SCALM is an option, shall be specified in .options scalm=val
 *  SCALM can be included in .model statement to override value for particular model 
 *  SCALE is also an option which shall be specified in .options SCALE=val
 *  can be included in device statement
 *  so - SCALM - for model, SCALE - for device
 * 
 * capacitor is not yet implemented this way - rather like in bmm_semi.cc
 *
 * Cxxx n1 n2 <mname> <C=> resistance <<TC1=>val> <<TC2=>val>
 *            <SCALE=val> <IC=val> <M=val> <DTEMP=val> <L=val> <W=val> 
 *
 * .model mname C keywd=value
 *
 *  device parameters
 *  keywd    UNITS   default   impl   descr
 *  CAP      F       0         n      default capcaitance
 *  CAPSW    F/m     0         n      sidewall fringing cap
 *  COX      F/m^2   0         n      bottomwall cap
 *  DEL      m       0         n      diff between drawn  and actual width/length DELeff=DEL*SCALM
 *  DI               0         n      relative dielectric const
 *  DLR      m       0         n      Diff between drawn and actual length DLReff=DLR*SCALM
 *  L        m       0         n      defualt length of cap
 *  SHRINK           1         n      shrink factor
 *  TC1C     1/deg   0         n      1st order tem coef for cap
 *  TC2C     1/deg^2 0         n      2nd order tem coef for cap
 *  THICK    m       0         n      diel thinkness
 *  TREF     C       TNOM      n      reference temperature
 *  W        m       0         n      Default width Wscaled=W*shrink*scalm
 *  scale            1         n      description - see ablove
 *
 * 
 * common parameters for the model
 *  par     R   C   both 
 *  BULK    +   
 *  CAP     +   +   +
 *  CAPSW   +   +   +
 *  COX     +   +   +
 *  DEL         +
 *  DI      +   +   +
 *  DLR     +   +   +
 *  DW      +
 *  L       +   +   +
 *  LEVEL   +
 *  RAC     +
 *  RES     +
 *  RSH     +
 *  SHRINK  +   +   +
 *  TC1C    +   +   +
 *  TC2C    +   +   +
 *  TC1R    +
 *  TC2R    +
 *  THICK   +   +   +
 *  TREF    +   +   +
 *  W       +   +   +
 *  scalm   +   +   +   description - see above 
 *
 * common parameters for the devices  - distribution in hierarchy
 *  COMMON_COMPONENT <-- EVAL_BM_BASE <-- EVAL_BM_ACTION_BASE <-- EVAL_BM_SEMI_BASE<-- EVAL_BM_SEMI_RESISTOR
 *                                                                                   - EVAL_BM_SEMI_CAPACITOR
 *
 *          Resistor        Capacitor
 * <R=>     resistance  
 * <C=>                     capacitance
 * <TC1=>   val             val             EVAL_BM_ACTION_BASE
 * <TC2=>   val             val             EVAL_BM_ACTION_BASE
 * SCALE    val             val             EVAL_BM_ACTION_BASE
 * IC                       val             EVAL_BM_ACTION_BASE
 * M        val             val             COMMON_COMPONENT
 * AC       val                             EVAL_BM_SEMI_RESISTOR
 * DTEMP    val             val             COMMON_COMPONENT
 * L        val             val             EVAL_BM_SEMI_BASE
 * W        val             val             EVAL_BM_SEMI_BASE
 * C        val                             EVAL_BM_SEMI_RESISTOR
 * scale    val             val             EVAL_BM_ACTION_BASE  -  implemented
 * 
 * COMMON_COMPONENT::parse_params_obsolete_callback()       tnom, dtemp, temp, m, mfactor
 * EVAL_BM_ACTION_BASE::parse_params_obsolete_callback()    bandwidth, delay, phase, ioffset, 
 *                                                          ooffset, scale, tc1, tc2, ic
 * EVAL_BM_SEMI_BASE::parse_params_obsolete_callback()      L, W,
 * EVAL_BM_SEMI_RESISTOR::...                               AC, C, [R],  [TC1], [TC2]
 * EVAL_BM_SEMI_CAPACITOR::...                              [C],  [TC1], [TC2]
 * [] - means parameter with optional name.
 *
 *
 * common parameters for the model  - distribution in hierarchy
 *  MODEL_CARD <-- MODEL_SEMI_BASE <-- MODEL_SEMI_RESISTPR
 *                                     MODEL_SEMI_CAPACITOR
 * 
 *  par     R   C   both 
 *  BULK    +   
 *  CAP     +   +   +                       MODEL_SEMI_BASE
 *  CAPSW   +   +   +                       MODEL_SEMI_BASE
 *  COX     +   +   +                       MODEL_SEMI_BASE
 *  DEL         +
 *  DI      +   +   +                       MODEL_SEMI_BASE
 *  DLR     +   +   +                       MODEL_SEMI_BASE
 *  DW      +                                   MODEL_SEMI_RESISTOR
 *  L       +   +   +                       MODEL_SEMI_BASE
 *  LEVEL   +                                   MODEL_SEMI_RESISTOR
 *  RAC     +                                   MODEL_SEMI_RESISTOR
 *  RES     +                                   MODEL_SEMI_RESISTOR
 *  RSH     +                                   MODEL_SEMI_RESISTOR
 *  SHRINK  +   +   +                       MODEL_SEMI_BASE
 *  TC1C    +   +   +                       MODEL_SEMI_BASE
 *  TC2C    +   +   +                       MODEL_SEMI_BASE
 *  TC1R    +                                   MODEL_SEMI_RESISTOR
 *  TC2R    +                                   MODEL_SEMI_RESISTOR
 *  THICK   +   +   +                       MODEL_SEMI_BASE
 *  TREF    +   +   +                   MODEL_CARD
 *  W       +   +   +                       MODEL_SEMI_BASE
 *  scalm   +   +   +                       MODEL_SEMI_BASE     - may be has to be moved upper  - but so far kept here
 *
 */
//
//testing
#include "u_lang.h"
#include "e_model.h" 
#include "bm.h"
/*
G. Serdyuk, 28-oct-2010:
Note, that classes EVAL_BM_SEMI_BASE, EVAL_BM_SEMI_RESISTOR, MODEL_BM_SEMI_BASE and
MODEL_BM_SEMI_RESISTOR in this file are differemy from classes with same names in bmm_semi.cc file
EVAL_BM_SEMI_CAPACITOR, MODEL_BM_SEMI_CAPACITOR coinside yet, but will be changed too.

With time those classed will be renames nbot to confuse them with classes in bmm_semi.cc

*/

/*--------------------------------------------------------------------------*/
class EVAL_BM_SEMI_BASE : public EVAL_BM_ACTION_BASE {
protected:
  PARAMETER<double> _length;
  PARAMETER<double> _width;
  double _value;
private:
  static double const _default_length;
  static double const _default_width;
  static double const _default_value;
protected:
  explicit EVAL_BM_SEMI_BASE(const EVAL_BM_SEMI_BASE& p);
           EVAL_BM_SEMI_BASE(const EVAL_BM_ACTION_BASE *a);   // constructor used to take parameters from EVAL_BM_MODEL 
           // (which is daugther from EVAM_BM_ACTION_BASE) 
           //  pointer originated from EVAL_BM_MODEL::expand() 
           // EVAL_BM_ACTIOIN_BASE* have choosen over EVAL_BM_MODEL - because this is basic class for TABLE SEMI... and MODEL itself
           // and would be good common denominator.
  explicit EVAL_BM_SEMI_BASE(int c=0);
  ~EVAL_BM_SEMI_BASE() {}
protected: // override virtual
  bool		operator==(const COMMON_COMPONENT&)const;
  COMMON_COMPONENT* clone()const = 0;
  void		print_common_obsolete_callback(OMSTREAM&, LANGUAGE*)const;

  void		precalc_first(const CARD_LIST*);
  void  	expand(const COMPONENT*);
  void		tr_eval(ELEMENT*)const;
  std::string	name()const	{untested();return modelname().c_str();}
  bool		ac_too()const		{untested();return false;}
  bool  	parse_params_obsolete_callback(CS&);
};
/*--------------------------------------------------------------------------*/
class EVAL_BM_SEMI_CAPACITOR : public EVAL_BM_SEMI_BASE {
private:
  explicit EVAL_BM_SEMI_CAPACITOR(const EVAL_BM_SEMI_CAPACITOR& p)
    :EVAL_BM_SEMI_BASE(p) {}
public:
  explicit EVAL_BM_SEMI_CAPACITOR(int c=0)
    :EVAL_BM_SEMI_BASE(c) {}
           EVAL_BM_SEMI_CAPACITOR(const EVAL_BM_ACTION_BASE *a)
    :EVAL_BM_SEMI_BASE(a) {}
  ~EVAL_BM_SEMI_CAPACITOR() {}
private: // override virtual
  bool		operator==(const COMMON_COMPONENT&)const;
  COMMON_COMPONENT* clone()const {return new EVAL_BM_SEMI_CAPACITOR(*this);}
  void  	expand(const COMPONENT*);
  void		precalc_last(const CARD_LIST*);
};
/*--------------------------------------------------------------------------*/
class EVAL_BM_SEMI_RESISTOR : public EVAL_BM_SEMI_BASE {
protected:
  PARAMETER<double> _resistance;
  PARAMETER<double> _capacitance;
  PARAMETER<double> _tc1;
  PARAMETER<double> _tc2;
  PARAMETER<double> _res_ac;
private:
  static double const _default_resistance;
  static double const _default_capacitance;
  static double const _default_tc1;
  static double const _default_tc2;
  static double const _default_res_ac;
private:
  explicit EVAL_BM_SEMI_RESISTOR(const EVAL_BM_SEMI_RESISTOR& p);  
public:
  explicit EVAL_BM_SEMI_RESISTOR(int c=0);
           EVAL_BM_SEMI_RESISTOR(const EVAL_BM_ACTION_BASE* a);
  ~EVAL_BM_SEMI_RESISTOR() {}
private: // override virtual
  bool		operator==(const COMMON_COMPONENT&)const;
  COMMON_COMPONENT* clone()const {return new EVAL_BM_SEMI_RESISTOR(*this);}  
  void		print_common_obsolete_callback(OMSTREAM&, LANGUAGE*)const;
  void  	expand(const COMPONENT*);
  void		precalc_last(const CARD_LIST*);
  bool  	parse_params_obsolete_callback(CS&);

};
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
class MODEL_SEMI_BASE : public MODEL_CARD {
public:
  PARAMETER<double> _narrow;
  PARAMETER<double> _defw;
  PARAMETER<double> _tc1;
  PARAMETER<double> _tc2;
  PARAMETER<double> _defl;
  
  PARAMETER<double> _cap;
  PARAMETER<double> _capsw;
  PARAMETER<double> _cox;
  PARAMETER<double> _di;
  PARAMETER<double> _dlr;
  PARAMETER<double> _shrink;
  PARAMETER<double> _thick;
  PARAMETER<double> _tref;

  PARAMETER<double> _scalm;

private:
  static double const _default_narrow;
  static double const _default_defw;
  static double const _default_tc1;
  static double const _default_tc2;
  static double const _default_defl;

  static double const _default_cap;
  static double const _default_capsw;
  static double const _default_cox;
  static double const _default_di;
  static double const _default_dlr;
  static double const _default_shrink;
  static double const _default_thick;
  static double const _default_tref;

  static double const _default_scalm;

protected:
  explicit MODEL_SEMI_BASE();
  explicit MODEL_SEMI_BASE(const MODEL_SEMI_BASE& p);
protected: // override virtual
  void  precalc_first();
  //void  precalc_last();
  //CARD* clone()const //MODEL_CARD/pure
  void		set_param_by_index(int, std::string&, int);
  bool		param_is_printable(int)const;
  std::string	param_name(int)const;
  std::string	param_name(int,int)const;
  std::string	param_value(int)const;
  int param_count()const {return (16 + MODEL_CARD::param_count());}
};
/*--------------------------------------------------------------------------*/
class MODEL_SEMI_CAPACITOR : public MODEL_SEMI_BASE {
public:
  PARAMETER<double> _cj;
  PARAMETER<double> _cjsw;
private:
  static double const _default_cj;
  static double const _default_cjsw;
private:
  explicit MODEL_SEMI_CAPACITOR(const MODEL_SEMI_CAPACITOR& p);
public:
  explicit MODEL_SEMI_CAPACITOR();
private: // override virtual
  std::string dev_type()const		{return "c";}
  void  precalc_first();
  //void  precalc_last();
  COMMON_COMPONENT* new_common()const {return new EVAL_BM_SEMI_CAPACITOR;}   // not used any more
  COMMON_COMPONENT* new_common(EVAL_BM_ACTION_BASE* a)const {return new EVAL_BM_SEMI_CAPACITOR;}
  CARD* clone()const		{return new MODEL_SEMI_CAPACITOR(*this);}
  void		set_param_by_index(int, std::string&, int);
  bool		param_is_printable(int)const;
  std::string	param_name(int)const;
  std::string	param_name(int,int)const;
  std::string	param_value(int)const;
  int param_count()const {return (2 + MODEL_SEMI_BASE::param_count());} 
};
/*--------------------------------------------------------------------------*/
class MODEL_SEMI_RESISTOR : public MODEL_SEMI_BASE {
public:
  PARAMETER<double> _rsh;
  PARAMETER<double> _res;
  PARAMETER<double> _rac;
  PARAMETER<int>    _level;
  PARAMETER<double> _dw;
  PARAMETER<double> _tc1r;
  PARAMETER<double> _tc2r;
private:
  static double const _default_rsh;
  static double const _default_res;
  static double const _default_rac;
  // static int const _default_level; // do not need it
  static double const _default_dw;
  static double const _default_tc1r;
  static double const _default_tc2r;
private:
  explicit MODEL_SEMI_RESISTOR(const MODEL_SEMI_RESISTOR& p);
public:
  explicit MODEL_SEMI_RESISTOR();
private: // override virtual
  std::string dev_type()const		{return "r";}
  void  precalc_first();
  //void  precalc_last();
  COMMON_COMPONENT* new_common()const {return new EVAL_BM_SEMI_RESISTOR;}   // not used any more
  COMMON_COMPONENT* new_common(EVAL_BM_ACTION_BASE* a)const {return new EVAL_BM_SEMI_RESISTOR (a);}
  CARD* clone()const		{return new MODEL_SEMI_RESISTOR(*this);}
  void		set_param_by_index(int, std::string&, int);
  bool		param_is_printable(int)const;
  std::string	param_name(int)const;
  std::string	param_name(int,int)const;
  std::string	param_value(int)const;
  int param_count()const {return (7 + MODEL_SEMI_BASE::param_count());}
};
/*--------------------------------------------------------------------------*/
double const EVAL_BM_SEMI_BASE::_default_length = NOT_INPUT;
double const EVAL_BM_SEMI_BASE::_default_width = NOT_INPUT;
double const EVAL_BM_SEMI_BASE::_default_value = NOT_INPUT;
/*--------------------------------------------------------------------------*/
static MODEL_SEMI_RESISTOR  p1;
static MODEL_SEMI_CAPACITOR p2;
static DISPATCHER<MODEL_CARD>::INSTALL
  d1(&model_dispatcher, "r|res", &p1),
  d2(&model_dispatcher, "c|cap", &p2);
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
EVAL_BM_SEMI_BASE::EVAL_BM_SEMI_BASE(int c)
  :EVAL_BM_ACTION_BASE(c),
   _length(_default_length),
   _width(_default_width),
   _value(_default_value)
{
}
/*--------------------------------------------------------------------------*/
EVAL_BM_SEMI_BASE::EVAL_BM_SEMI_BASE(const EVAL_BM_SEMI_BASE& p)
  :EVAL_BM_ACTION_BASE(p),
   _length(p._length),
   _width(p._width),
   _value(p._value)
{
}
/*--------------------------------------------------------------------------*/
EVAL_BM_SEMI_BASE::EVAL_BM_SEMI_BASE(const EVAL_BM_ACTION_BASE* a)
  :EVAL_BM_ACTION_BASE(*a),
   _length(_default_length),
   _width(_default_width),
   _value(_default_value)
{
}
/*--------------------------------------------------------------------------*/
bool EVAL_BM_SEMI_BASE::operator==(const COMMON_COMPONENT& x)const
{
  const EVAL_BM_SEMI_BASE* p = dynamic_cast<const EVAL_BM_SEMI_BASE*>(&x);
  bool rv = p
    && _length == p->_length
    && _width == p->_width
    && EVAL_BM_ACTION_BASE::operator==(x);
  if (rv) {
    untested();
  }else{
  }
  return rv;
}
/*--------------------------------------------------------------------------*/
void EVAL_BM_SEMI_BASE::print_common_obsolete_callback(OMSTREAM& o, LANGUAGE* lang)const
{
  assert(lang);
  o << modelname();
  print_pair(o, lang, "l", _length);
  print_pair(o, lang, "w", _width, _width.has_hard_value());
  EVAL_BM_ACTION_BASE::print_common_obsolete_callback(o, lang);
}
/*--------------------------------------------------------------------------*/
void EVAL_BM_SEMI_BASE::expand(const COMPONENT* d)
{
//    std::cout<<"EVAL_BM_SEMI_BASE::expand "<< modelname()<<" _temp_c="<<_temp_c<<"\n";
  EVAL_BM_ACTION_BASE::expand(d);
  attach_model(d);
}
/*--------------------------------------------------------------------------*/
void EVAL_BM_SEMI_BASE::precalc_first(const CARD_LIST* Scope)
{
  assert(Scope);
  EVAL_BM_ACTION_BASE::precalc_first(Scope);
  _length.e_val(_default_length, Scope);
  _width.e_val(_default_width, Scope);
}
/*--------------------------------------------------------------------------*/
void EVAL_BM_SEMI_BASE::tr_eval(ELEMENT* d)const
{
  tr_finish_tdv(d, _value);
}
/*--------------------------------------------------------------------------*/
bool EVAL_BM_SEMI_BASE::parse_params_obsolete_callback(CS& cmd)
{
  return ONE_OF
    || Get(cmd, "l",	&_length)
    || Get(cmd, "w",	&_width)
    || EVAL_BM_ACTION_BASE::parse_params_obsolete_callback(cmd)
    ;
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
bool EVAL_BM_SEMI_CAPACITOR::operator==(const COMMON_COMPONENT& x)const
{
  const EVAL_BM_SEMI_CAPACITOR*
    p = dynamic_cast<const EVAL_BM_SEMI_CAPACITOR*>(&x);
  bool rv = p
    && EVAL_BM_SEMI_BASE::operator==(x);
  if (rv) {
    untested();
  }else{
  }
  return rv;
}
/*--------------------------------------------------------------------------*/
void EVAL_BM_SEMI_CAPACITOR::expand(const COMPONENT* d)
{
  EVAL_BM_SEMI_BASE::expand(d);

  const MODEL_SEMI_CAPACITOR* m = dynamic_cast<const MODEL_SEMI_CAPACITOR*>(model());
  if (!m) {
    unreachable();
    throw Exception_Model_Type_Mismatch(d->long_label(), modelname(), "semi-capacitor (C)");
  }else{
  }
}
/*--------------------------------------------------------------------------*/
void EVAL_BM_SEMI_CAPACITOR::precalc_last(const CARD_LIST* Scope)
{
  assert(Scope);
  EVAL_BM_SEMI_BASE::precalc_last(Scope);

  const MODEL_SEMI_CAPACITOR* m = prechecked_cast<const MODEL_SEMI_CAPACITOR*>(model());

  double width = (_width == NOT_INPUT) ? m->_defw : _width;
  double eff_width = width - m->_narrow;
  //
  double eff_length = _length - m->_narrow;
  _value = m->_cj * eff_length * eff_width + 2. * m->_cjsw * (eff_length + eff_width);
  double tempdiff = (_temp_c - m->_tnom_c);
/*
  std::cout<<
  " CAP \n"<<
  " _temp_c   ="<<_temp_c<<"\n"<<
  " m->tnom_c ="<<m->_tnom_c<<"\n"<<
  " tempdiff  ="<<tempdiff<<"\n"<<
  "\n";
*/
  _value *= 1 + m->_tc1*tempdiff + m->_tc2*tempdiff*tempdiff;

  if (eff_width <= 0.) {untested();
    throw Exception_Precalc(modelname() + ": effective width is negative or zero\n");
  }else{
  }
  if (eff_length <= 0.) {untested();
    throw Exception_Precalc(modelname() + ": effective length is negative or zero\n");
  }else{
  }
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
double const EVAL_BM_SEMI_RESISTOR::_default_resistance = NOT_INPUT;
double const EVAL_BM_SEMI_RESISTOR::_default_tc1 = NOT_INPUT;
double const EVAL_BM_SEMI_RESISTOR::_default_tc2 = NOT_INPUT;
double const EVAL_BM_SEMI_RESISTOR::_default_capacitance = NOT_INPUT;
double const EVAL_BM_SEMI_RESISTOR::_default_res_ac = NOT_INPUT;
/*--------------------------------------------------------------------------*/
EVAL_BM_SEMI_RESISTOR::EVAL_BM_SEMI_RESISTOR(int c)
  :EVAL_BM_SEMI_BASE(c),
   _resistance(_default_resistance),
   _capacitance(_default_capacitance),
   _tc1(_default_tc1),
   _tc2(_default_tc2),
   _res_ac(_default_res_ac)
{
}
/*--------------------------------------------------------------------------*/
EVAL_BM_SEMI_RESISTOR::EVAL_BM_SEMI_RESISTOR(const EVAL_BM_SEMI_RESISTOR& p)
  :EVAL_BM_SEMI_BASE(p),
   _resistance(p._resistance),
   _capacitance(p._capacitance),
   _tc1(p._tc1),
   _tc2(p._tc2),
   _res_ac(p._res_ac)
{
}
/*--------------------------------------------------------------------------*/

EVAL_BM_SEMI_RESISTOR::EVAL_BM_SEMI_RESISTOR(const EVAL_BM_ACTION_BASE* a)
  :EVAL_BM_SEMI_BASE(a),  // SIC##
   _resistance(_default_resistance),
   _capacitance(_default_capacitance),
   _tc1(_default_tc1),
   _tc2(_default_tc2),
   _res_ac(_default_res_ac)
{
}

/*--------------------------------------------------------------------------*/
bool EVAL_BM_SEMI_RESISTOR::operator==(const COMMON_COMPONENT& x)const
{
  const EVAL_BM_SEMI_RESISTOR* p = dynamic_cast<const EVAL_BM_SEMI_RESISTOR*>(&x);
  bool rv = p
    && _resistance  == p->_resistance
    && _capacitance == p->_capacitance
    && _tc1         == p->_tc1
    && _tc2         == p->_tc2
    && _res_ac      == p->_res_ac
    && EVAL_BM_SEMI_BASE::operator==(x);
  if (rv) {
    untested();
  }else{
  }
  return rv;
}
/*--------------------------------------------------------------------------*/
void EVAL_BM_SEMI_RESISTOR::expand(const COMPONENT* d)
{
//   std::cout<<"EVAL_BM_SEMI_RESISTOR::expand "<< modelname()<<" _temp_c="<<_temp_c<<"\n";
  EVAL_BM_SEMI_BASE::expand(d);

  const MODEL_SEMI_RESISTOR* m = dynamic_cast<const MODEL_SEMI_RESISTOR*>(model());
  if (!m) {
    unreachable();
    throw Exception_Model_Type_Mismatch(d->long_label(), modelname(), "semi-resistor (R)");
  }else{
  }
}
/*--------------------------------------------------------------------------*/
void EVAL_BM_SEMI_RESISTOR::precalc_last(const CARD_LIST* Scope)
{

//  std::cout<<"EVAL_BM_SEMI_RESISTOR::precalc_last "<<modelname()<<"\n";
  assert(Scope);
  EVAL_BM_SEMI_BASE::precalc_last(Scope);
  _resistance.e_val(_default_resistance, Scope);
  _tc1.e_val(_default_tc1, Scope);
  _tc2.e_val(_default_tc2, Scope);
  _capacitance.e_val(_default_capacitance, Scope);
  _res_ac.e_val(_default_capacitance, Scope);

  const MODEL_SEMI_RESISTOR* m = prechecked_cast<const MODEL_SEMI_RESISTOR*>(model());

  double tc1= (_tc1 == NOT_INPUT ) ? m->_tc1: _tc1;
  double tc2= (_tc2 == NOT_INPUT ) ? m->_tc2: _tc2;

  double eff_width;
  double eff_length;   

  double resmin = 1/(OPT::gmax * _mfactor * 1000);

/*  
  std::cout<<" precalc_last: at input:"<<"\n"<<
    " _scale     ="<<_scale<<"\n"<<
    " _mfactor   ="<<_mfactor<<"\n"<<
    " _resistance="<<_resistance<<"\n"<<
    " _length    ="<<_length<<"\n"<<
    " _width     ="<<_width<<"\n"<<
    " m->_defl     ="<<m->_defl<<"\n"<<
    " m->_defw     ="<<m->_defw<<"\n"<<
    " m->_dlr      ="<<m->_dlr<<"\n"<<
    " m->_shrink   ="<<m->_shrink<<"\n"<<
    " m->_scalm    ="<<m->_scalm<<"\n"
    ;
*/ 
   
  _value = BIGBIG;   // if something will go wrong - very large resistance 
  
  bool scale_warning=false;
  bool mfactor_warning=false;
  if(_scale == 0.0){
    scale_warning=true;
    _scale=1.;
    }
  else {}
    
  if (_mfactor==0.0){
    mfactor_warning=true;
    _mfactor=1.;
    }  
  else{}
 
  bool rescalc_warning=false; 
  if (_resistance != NOT_INPUT) {               //resistance is specified
    _value=_resistance * _scale/ _mfactor;      // done - asssured that _scale !=0 , _m != 0 - 

/*  std::cout<<" precalc_last: branch 1"<<"\n"
    " _value="<<_value<<"\n";
*/    
    }
  else {
//  std::cout<<" precalc_last: branch 2 \n";

    double width    = (_width == NOT_INPUT)    ? m->_defw : _width;    
    width           = width * m->_shrink * m->_scalm;              // w_scaled
    eff_width       = width - m->_dw * m->_scalm;                  // w_eff = w_scaled-dw_eff, dw_eff=dw*scalm

/*    std::cout<<
    " _length    ="<<_length<<"\n"<<
    " m->_defl     ="<<m->_defl<<"\n";
*/
    double length   = (_length == NOT_INPUT)  ? m->_defl : _length;    
    length          = length * m->_shrink * m->_scalm;             // l_scaled
    eff_length      = length - m->_dlr * m->_scalm;               // l_eff = l_scaled-dlr_eff, dlr_eff=dlr*scalm
    
/*    std::cout<<" precalc_last: branch2"<<"\n"<<
    " eff_width ="<<eff_width<<"\n"<<
    " eff_length="<<eff_length<<"\n"<<
    " m->_rsh   ="<<m->_rsh<<"\n"<<
    " m->_shrink="<<m->_shrink<<"\n"<<
    " m->_scalm ="<<m->_scalm<<"\n"<<
    " m->_dw    ="<<m->_dw<<"\n"<<
    " m->_dlr   ="<<m->_dlr<<"\n";
*/

    if (eff_width * eff_length * m->_rsh > 0.) {
      _value = m->_rsh * eff_length / eff_width * _scale / _mfactor;
      }
    else if ( eff_width * eff_length * m->_rsh == 0.)  {
      _value = m->_res *  _scale / _mfactor;
      }
    else{
      rescalc_warning=true;
      unreachable(); // indeed - reachable. think about exception
      }
    
   }
   
 /*  todo - here - think which temperatures to use
    we have a lot of different tref (device and model) and tnom.
    
 */  
//  std::cout<<" _value="<<_value<<"\n";
  double tempdiff = (_temp_c - m->_tnom_c); 

/*  std::cout<<
  "Res \n"<<
  " _temp_c   ="<<_temp_c<<"\n"<<
  " m->tnom_c ="<<m->_tnom_c<<"\n"<<
  " tempdiff  ="<<tempdiff<<"\n"<<
  " tc1       ="<<tc1<<"\n"<<
  " tc2       ="<<tc2<<"\n"<<
  "\n";
*/
                                            
  _value *= 1 + tc1*tempdiff + tc2*tempdiff*tempdiff;
//  std::cout<<" tempdiff: _value="<<_value<<"\n";

  bool resmin_warning=false;
  if ( _value < resmin){
     resmin_warning=true;
    _value=resmin;
    }
  else{
    }
/*
  std::cout<<
  "  resmin_warning = " << resmin_warning  <<
  "  scale_warning  = " << scale_warning   <<
  "  mfactor_warning= " << mfactor_warning <<
  "  rescalc_warning= " << rescalc_warning <<
  "  _value="           << _value<<
  "\n";
*/  
  if (resmin_warning || scale_warning || mfactor_warning ||  rescalc_warning){
    std::string exception_msg="";
    if (resmin_warning)
         exception_msg += " at "+modelname() +": RS001 resistance is close to zero, resmin used\n";
    else {}
    
    if (scale_warning)
         exception_msg += " at "+modelname() +": RS002 scale is equal to zero, fixed to 1.0 \n";
    else{}
  
    if (mfactor_warning)
         exception_msg += " at "+modelname() +": RS003 mfactor is equal zero, fixed to 1.0\n";
    else{}
  
    if (rescalc_warning)
         exception_msg += " at "+modelname() +": RS004 resistance calculation failure, infinite (BIGBIG) is used\n";
    else{}
  
    throw Exception_Precalc (exception_msg);
    }
  else {}
}
/*--------------------------------------------------------------------------*/
bool EVAL_BM_SEMI_RESISTOR::parse_params_obsolete_callback(CS& cmd)
{    
  return ONE_OF
    || (GetToken(cmd, &_resistance) && GetToken(cmd, &_tc1) && GetToken(cmd, &_tc2))   // consumes (possible) first three num parameters
    || Get(cmd, "r",	&_resistance)
    || Get(cmd, "c",    &_capacitance)
    || EVAL_BM_SEMI_BASE::parse_params_obsolete_callback(cmd)
    ;
}
/*--------------------------------------------------------------------------*/
void EVAL_BM_SEMI_RESISTOR::print_common_obsolete_callback(OMSTREAM& o, LANGUAGE* lang)const
{
  assert(lang);
  EVAL_BM_SEMI_BASE::print_common_obsolete_callback(o, lang);
  //o << modelname();
  print_pair(o, lang, "r", _resistance);
  print_pair(o, lang, "c", _capacitance,_capacitance.has_hard_value());
  print_pair(o, lang, "tc1", _tc1,      _tc1.has_hard_value());
  print_pair(o, lang, "tc2", _tc2,      _tc2.has_hard_value());
  print_pair(o, lang, "rac", _res_ac,   _res_ac.has_hard_value());
  
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*
* common parameters for the model
 *  par     R   C   both             
 *  BULK    +                   -
 *  CAP     +   +   +           5 new
 *  CAPSW   +   +   +           6
 *  COX     +   +   +           7
 *  DEL         +               - 
 *  DI      +   +   +           8
 *  DLR     +   +   +           9
 *  DW      +                   -
 *  L       +   +   +           4
 *  LEVEL   +                   -
 *  RAC     +                   -
 *  RES     +                   -
 *  RSH     +                   -
 *  SHRINK  +   +   +           10
 *  TC1C    +   +   +           2 
 *  TC2C    +   +   +           3
 *  TC1R    +                   -
 *  TC2R    +                   -
 *  THICK   +   +   +           11
 *  TREF    +   +   +           12
 *  W       +   +   +            1
 *  scalm   +   +   +           15
 --------
 * narrow
* add names "defl" and "defw" in addition to "W" and "L" :
 * defl                         13
 * defw                         14

 *  par     R   C   both 
 *  BULK    +   
 *  CAP     +   +   +                       MODEL_SEMI_BASE
 *  CAPSW   +   +   +                       MODEL_SEMI_BASE
 *  COX     +   +   +                       MODEL_SEMI_BASE
 *  DEL         +
 *  DI      +   +   +                       MODEL_SEMI_BASE
 *  DLR     +   +   +                       MODEL_SEMI_BASE
 *  DW      +                                   MODEL_SEMI_RESISTOR
 *  L       +   +   +                       MODEL_SEMI_BASE
 *  LEVEL   +                                   MODEL_SEMI_RESISTOR
 *  RAC     +                                   MODEL_SEMI_RESISTOR
 *  RES     +                                   MODEL_SEMI_RESISTOR
 *  RSH     +                                   MODEL_SEMI_RESISTOR
 *  SHRINK  +   +   +                       MODEL_SEMI_BASE
 *  TC1C    +   +   +                       MODEL_SEMI_BASE
 *  TC2C    +   +   +                       MODEL_SEMI_BASE
 *  TC1R    +                                   MODEL_SEMI_RESISTOR
 *  TC2R    +                                   MODEL_SEMI_RESISTOR
 *  THICK   +   +   +                       MODEL_SEMI_BASE
 *  TREF    +   +   +                   MODEL_CARD
 *  W       +   +   +                       MODEL_SEMI_BASE
 *  scalm   +   +   +                       MODEL_SEMI_BASE

*/

double const MODEL_SEMI_BASE::_default_narrow = 0.;
double const MODEL_SEMI_BASE::_default_defw = 1e-6;
double const MODEL_SEMI_BASE::_default_tc1 = 0.;
double const MODEL_SEMI_BASE::_default_tc2 = 0.;
double const MODEL_SEMI_BASE::_default_defl = 0.;

double const MODEL_SEMI_BASE::_default_cap = 0.;  
double const MODEL_SEMI_BASE::_default_capsw = 0.;
double const MODEL_SEMI_BASE::_default_cox = 0.;
double const MODEL_SEMI_BASE::_default_di = 0.;
double const MODEL_SEMI_BASE::_default_dlr = 0.;
double const MODEL_SEMI_BASE::_default_shrink = 1.;
double const MODEL_SEMI_BASE::_default_thick = 0.;
// double const MODEL_SEMI_BASE::_default_tref = NOT_INPUT;  // done - no need it
// double const MODEL_SEMI_BASE::_default_scalm = NOT_INPUT; // done - no need it

/*--------------------------------------------------------------------------*/
MODEL_SEMI_BASE::MODEL_SEMI_BASE()
  :MODEL_CARD(NULL),
   _narrow(_default_narrow),
   _defw(_default_defw),
   _tc1(_default_tc1),
   _tc2(_default_tc2),
   _defl(_default_defl), 
   
   _cap(_default_cap),
   _capsw(_default_capsw),
   _cox(_default_cox),
   _di(_default_di),
   _dlr(_default_dlr),
   _shrink(_default_shrink),
   _thick(_default_thick),
   _tref(OPT::tnom_c),      // instead of _default_tref
   
   _scalm(OPT::scalm)       // instead of _default_scalm
{
}
/*--------------------------------------------------------------------------*/
MODEL_SEMI_BASE::MODEL_SEMI_BASE(const MODEL_SEMI_BASE& p)
  :MODEL_CARD(p),
   _narrow(p._narrow),
   _defw(p._defw),
   _tc1(p._tc1),
   _tc2(p._tc2),
   _defl(p._defl),

  _cap  (p._cap),
  _capsw(p._capsw),
  _cox  (p._cox),
  _di   (p._di),
  _dlr  (p._dlr),
  _shrink(p._shrink),
  _thick(p._thick),
  _tref (p._tref),

   _scalm(p._scalm)
{
}
/*--------------------------------------------------------------------------*/
void MODEL_SEMI_BASE::set_param_by_index(int i, std::string& value, int offset)
{
  switch (MODEL_SEMI_BASE::param_count() - 1 - i) {
  case 0: _narrow = value; break;
  case 1: _defw = value; break;
  case 2: _tc1  = value; break;
  case 3: _tc2  = value; break;
  case 4: _defl = value; break;

  case 5: _cap  = value; break;
  case 6: _capsw = value; break;
  case 7: _cox  = value; break;
  case 8: _di   = value; break;
  case 9: _dlr  = value; break;
  case 10: _shrink = value; break;
  case 11: _thick  = value; break;
  case 12: _tref   = value; break;

  case 13: _defw   = value; break;
  case 14: _defl   = value; break;

  case 15: _scalm = value; break;
    
  default: MODEL_CARD::set_param_by_index(i, value, offset); break;
  }
}
/*--------------------------------------------------------------------------*/
bool MODEL_SEMI_BASE::param_is_printable(int i)const
{
  switch (MODEL_SEMI_BASE::param_count() - 1 - i) {
  case 0: 
  case 1: 
  case 2: 
  case 3:
  case 4: 
  
  case 5: 
  case 6: 
  case 7: 
  case 8: 
  case 9: 
  case 10: 
  case 11: 
  case 12: 

  case 13: 
  case 14: 

  case 15: 
  return true;
  default: return MODEL_CARD::param_is_printable(i);
  }
}
/*--------------------------------------------------------------------------*/
std::string MODEL_SEMI_BASE::param_name(int i)const
{
  switch (MODEL_SEMI_BASE::param_count() - 1 - i) {
  case 0: return "narrow";
  case 1: return "w";
  case 2: return "tc1";
  case 3: return "tc2";
  case 4: return "l";

  case 5: return "cap";
  case 6: return "capsw";
  case 7: return "cox";
  case 8: return "di";
  case 9: return "dlr";
  case 10: return "shrink";
  case 11: return "thick";
  case 12: return "tref";

  case 13: return "defw";
  case 14: return "defl";

  case 15: return "scalm";
  default: return MODEL_CARD::param_name(i);
  }
}
/*--------------------------------------------------------------------------*/
std::string MODEL_SEMI_BASE::param_name(int i, int j)const
{
  if (j == 0) {untested();
    return param_name(i);
  }else if (i >= MODEL_CARD::param_count()) {
    return "";
  }else{
    return MODEL_CARD::param_name(i, j);
  }
}
/*--------------------------------------------------------------------------*/
std::string MODEL_SEMI_BASE::param_value(int i)const
{
  switch (MODEL_SEMI_BASE::param_count() - 1 - i) {
  case 0: return _narrow.string();
  case 1: return _defw.string();
  case 2: return _tc1.string();
  case 3: return _tc2.string();
  case 4: return _defl.string();

  case 5: return _cap.string();
  case 6: return _capsw.string();
  case 7: return _cox.string();
  case 8: return _di.string();
  case 9: return _dlr.string();
  case 10: return _shrink.string();
  case 11: return _thick.string();
  case 12: return _tref.string();

  case 13: return _defw.string();
  case 14: return _defl.string();

  case 15: return _scalm.string(); 
  default: return MODEL_CARD::param_value(i);
  }
}
/*--------------------------------------------------------------------------*/
void MODEL_SEMI_BASE::precalc_first()
{
  MODEL_CARD::precalc_first();

  const CARD_LIST* s = scope();
  assert(s);

  _narrow.e_val(_default_narrow, s);
  _defw.e_val(_default_defw, s);
  _tc1.e_val(_default_tc1, s);
  _tc2.e_val(_default_tc2, s);
  _defl.e_val(_default_defl, s);

  _cap.e_val(_default_cap,s);
  _capsw.e_val(_default_capsw,s);
  _cox.e_val(_default_cox,s);
  _di.e_val(_default_di,s);
  _dlr.e_val(_default_dlr,s);
  _shrink.e_val(_default_shrink,s);
  _thick.e_val(_default_thick, s);
  _tref.e_val(OPT::tnom_c,s);  // OPT::tnom_c is default for it
  _scalm.e_val(OPT::scalm,s);  // OPT::scalm is default for it
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
double const MODEL_SEMI_CAPACITOR::_default_cj = 0.;
double const MODEL_SEMI_CAPACITOR::_default_cjsw = 0.;
/*--------------------------------------------------------------------------*/
MODEL_SEMI_CAPACITOR::MODEL_SEMI_CAPACITOR()
  :MODEL_SEMI_BASE(),
   _cj(_default_cj),
   _cjsw(_default_cjsw)
{
}
/*--------------------------------------------------------------------------*/
MODEL_SEMI_CAPACITOR::MODEL_SEMI_CAPACITOR(const MODEL_SEMI_CAPACITOR& p)
  :MODEL_SEMI_BASE(p),
   _cj(p._cj),
   _cjsw(p._cjsw)
{
}
/*--------------------------------------------------------------------------*/
void MODEL_SEMI_CAPACITOR::set_param_by_index(int i, std::string& value, int offset)
{
  switch (MODEL_SEMI_CAPACITOR::param_count() - 1 - i) {
  case 0: _cj = value; break;
  case 1: _cjsw = value; break;
  default: MODEL_SEMI_BASE::set_param_by_index(i, value, offset); break;
  }
}
/*--------------------------------------------------------------------------*/
bool MODEL_SEMI_CAPACITOR::param_is_printable(int i)const
{
  switch (MODEL_SEMI_CAPACITOR::param_count() - 1 - i) {
  case 0: 
  case 1: return true;
  default: return MODEL_SEMI_BASE::param_is_printable(i);
  }
}
/*--------------------------------------------------------------------------*/
std::string MODEL_SEMI_CAPACITOR::param_name(int i)const
{
  switch (MODEL_SEMI_CAPACITOR::param_count() - 1 - i) {
  case 0: return "cj";
  case 1: return "cjsw";
  default: return MODEL_SEMI_BASE::param_name(i);
  }
}
/*--------------------------------------------------------------------------*/
std::string MODEL_SEMI_CAPACITOR::param_name(int i, int j)const
{
  if (j == 0) {
    return param_name(i);
  }else if (i >= MODEL_SEMI_BASE::param_count()) {
    return "";
  }else{
    return MODEL_SEMI_BASE::param_name(i, j);
  }
}
/*--------------------------------------------------------------------------*/
std::string MODEL_SEMI_CAPACITOR::param_value(int i)const
{
  switch (MODEL_SEMI_CAPACITOR::param_count() - 1 - i) {
  case 0: return _cj.string();
  case 1: return _cjsw.string();
  default: return MODEL_SEMI_BASE::param_value(i);
  }
}
/*--------------------------------------------------------------------------*/
void MODEL_SEMI_CAPACITOR::precalc_first()
{
  MODEL_SEMI_BASE::precalc_first();

  const CARD_LIST* s = scope();
  assert(s);

  _cj.e_val(_default_cj, s);
  _cjsw.e_val(_default_cjsw, s);
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
double const MODEL_SEMI_RESISTOR::_default_rsh  = 0;
double const MODEL_SEMI_RESISTOR::_default_dw   = 0;
double const MODEL_SEMI_RESISTOR::_default_rac  = 0;
double const MODEL_SEMI_RESISTOR::_default_res  = 0;
double const MODEL_SEMI_RESISTOR::_default_tc1r = 0;
double const MODEL_SEMI_RESISTOR::_default_tc2r = 0;
/*--------------------------------------------------------------------------*/
MODEL_SEMI_RESISTOR::MODEL_SEMI_RESISTOR()
  :MODEL_SEMI_BASE(),
   _rsh(_default_rsh),
   _dw(_default_dw),
   _rac(_default_rac),
   _res(_default_res),
   _tc1r(_default_tc1r),
   _tc2r(_default_tc2r),
   _level(0)            // not used
{
}
/*--------------------------------------------------------------------------*/
MODEL_SEMI_RESISTOR::MODEL_SEMI_RESISTOR(const MODEL_SEMI_RESISTOR& p)
  :MODEL_SEMI_BASE(p),  // GS - fix, it was :MODEL_SEMI_BASE(),
   _rsh(p._rsh),
   _dw(p._dw),
   _rac(p._rac),
   _res(p._res),
   _tc1r(p._tc1r),
   _tc2r(p._tc2r),
   _level(p._level)
{
}
/*--------------------------------------------------------------------------*/
void MODEL_SEMI_RESISTOR::set_param_by_index(int i, std::string& value, int offset)
{
  switch (MODEL_SEMI_RESISTOR::param_count() - 1 - i) {
  case 0: _rsh   = value; break;
  case 1: _dw    = value; break;
  case 2: _level = value; break;
  case 3: _rac   = value; break;
  case 4: _res   = value; break;
  case 5: _tc1r  = value; break;
  case 6: _tc2r  = value; break;

  default: MODEL_SEMI_BASE::set_param_by_index(i, value, offset); break;
  }
}
/*--------------------------------------------------------------------------*/
bool MODEL_SEMI_RESISTOR::param_is_printable(int i)const
{
  switch (MODEL_SEMI_RESISTOR::param_count() - 1 - i) {
  case 0: case 1: case 2: case 3: case 4: case 5: case 6: return true;
  default: return MODEL_SEMI_BASE::param_is_printable(i);
  }
}
/*--------------------------------------------------------------------------*/
std::string MODEL_SEMI_RESISTOR::param_name(int i)const
{
  switch (MODEL_SEMI_RESISTOR::param_count() - 1 - i) {
  case 0: return "rsh";
  case 1: return "dw";
  case 2: return "level";
  case 3: return "rac";
  case 4: return "res";
  case 5: return "tc1r";
  case 6: return "tc2r";
  default: return MODEL_SEMI_BASE::param_name(i);
  }
}
/*--------------------------------------------------------------------------*/
std::string MODEL_SEMI_RESISTOR::param_name(int i, int j)const
{
  if (j == 0) {
    return param_name(i);
  }else if (i >= MODEL_SEMI_BASE::param_count()) {
    return "";
  }else{
    return MODEL_SEMI_BASE::param_name(i, j);
  }
}
/*--------------------------------------------------------------------------*/
std::string MODEL_SEMI_RESISTOR::param_value(int i)const
{
  switch (MODEL_SEMI_RESISTOR::param_count() - 1 - i) {
  case 0: return _rsh.string();
  case 1: return _dw.string();
  case 2: return _level.string();
  case 3: return _rac.string();
  case 4: return _res.string();
  case 5: return _tc1r.string();
  case 6: return _tc2r.string();
  default: return MODEL_SEMI_BASE::param_value(i);
  }
}
/*--------------------------------------------------------------------------*/
void MODEL_SEMI_RESISTOR::precalc_first()
{
  MODEL_SEMI_BASE::precalc_first();

  const CARD_LIST* par_scope = scope();
  assert(par_scope);

  _rsh.e_val (_default_rsh,  par_scope);
  _dw.e_val  (_default_dw,   par_scope);
  _rac.e_val (_default_rac,  par_scope);
  _res.e_val (_default_res,  par_scope);
  _tc1r.e_val(_default_tc1r, par_scope);
  _tc2r.e_val(_default_tc2r, par_scope);
  
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
Done list:
1) Warnings in EVAL_SEMI_RESISTOR::precalc_last - concentrate them after main caclulation 
   and pass to outer scope as exceptions.
   See COMPONENT::precalc_last() to understand how to pass warnings from any COMMMON_COMPONENT to containing COMPONENT

2) watch for param_count() method - it should return correct value

3) understand process of transformation of model-containing common_component (bm_model) to proper class (EVAL_ACTION_BM_BASE hierarchy).
   a) usually any component (upd: some of them), which has model is created as bm_model
   b) during elaboration process, at expansion stage, component, which cotains bm_model (in form of common_component *) 
      calls bm_model->expand()
   c) in that moment member common_component::_model already has known type of the related ".model" statement. and when calls
      model()->new_common() which creates new  correspondent common_component
   
      for example. 
      parsing
      r1 a b model_r r=10
      .model model_r r ....
      
      a)
      COMPONENT:        = (r1 a b )
         |-common       = COMMON_COMPONENT*: BM_MODEL
               |- model = MODEL_SEMI_RESISTOR
               
      b) COMPONNET::expand()
           EVAL_BM_MODEL::expand
           
      c) at EVAL_BM_MODEL:
           common->model -> new_common()
           as common->model of proper model type (MODEL_SEMI_RESISTOR), 
           so called proper MODEL_SEMI_RESISTOR::new_common(), 
           which creates EVAL_BM_SEMI_RESISTOR for MODEL_SEMI_RESISTOR - 
           see bmm_semi(_hsp).cc for example
      ISSUE:
      EVAL_BM_MODEL already passed precalc_first and its temperature was set
      meanwhile EVAL_SEMI_RESISTOR - does not, so it has to be properly initialized.   
      
    What happends further:
    d) new common, created at step "c" (call it "c") is attached to EVAL_BM_MODLE::_func
    as I understand _func is special storate for it
    then we return from expand() and expect deflation of eval_bm_expand
    
    fixed - new constriuctors added, now "this" is passed to new_common and available parameters passed to 
    newly created common_component.      

4) understand how COMMON_COMPONENT passes data to COMPONENT - for example resistance
   Answer: _value is passed from COMMON_COMPONENT (namely - from EVAL_BM_ACTION_BASE) to containing device using metho
   EVAL_BM_ACTION_BASE::tr_finish_tdv
   it is called both in OP and AC:
   ac_eval -> tr_eval -> tr_finish_tdv

5) fix in bmm_semi.cc
   done
   
Todo list

1) Understand all temp coeficients - how, when and what to calculate

2) - param.2a param.2a-1, etc - make a w/o "'" calculated

3) write tests - step 1

4) complete capacitances -                  step 2

5) complete RAC - resistance for AC mode    step 2

*/

