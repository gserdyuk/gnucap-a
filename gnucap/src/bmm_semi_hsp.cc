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
 *  DW       m       0         n      
 *  L        m       0         y      defualt length of wire
 *  LEVEL                      n      model selector (n/a)
 *  RAC      Ohm               n      defualt AC redsistance
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
 *  also can be included in .model statement
 *
 * capacitor is not yet implemented this way - rather like in bmm_semi.cc
 *
 * Cxxx n1 n2 <mname> <C=> resistance <<TC1=>val> <<TC2=>val>
 *            <SCALE=val> <IC=val> <M=val> <DTEMP=val> <L=val> <W=val> 
 *
 * .model mname C keywd=value
 *
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
 *
 * common parameters for the devices
 *          Resistor        Capacitor
 * <R=>     resistance  
 * <C=>                     capacitance
 * <TC1=>   val             val             EVAL_BM_ACTION_BASE
 * <TC2=>   val             val             EVAL_BM_ACTION_BASE
 * SCALE    val             val             EVAL_BM_ACTION_BASE
 * IC                       val             EVAL_BM_ACTION_BASE
 * M        val             val             COMMON_COMPONENT
 * AC       val
 * DTEMP    val             val             COMMON_COMPONENT
 * L        val             val             EVAL_BM_SEMI_BASE
 * W        val             val             EVAL_BM_SEMI_BASE
 * C        val
 * 
 * COMMON_COMPONENT::parse_params_obsolete_callback()       tnom, dtemp, temp, m, mfactor
 * EVAL_BM_ACTION_BASE::parse_params_obsolete_callback()    bandwidth, delay, phase, ioffset, 
 *                                                          ooffset, scale, tc1, tc2, ic
 * EVAL_BM_SEMI_BASE::parse_params_obsolete_callback()      L, W,
 * EVAL_BM_SEMI_RESISTOR::...                               AC, C, [R],  [TC1], [TC2]
 * EVAL_BM_SEMI_CAPACITOR::...                              [C],  [TC1], [TC2]
 * [] - means parameter with optional name.
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
private:
  static double const _default_resistance;
  static double const _default_capacitance;
  static double const _default_tc1;
  static double const _default_tc2;
private:
  explicit EVAL_BM_SEMI_RESISTOR(const EVAL_BM_SEMI_RESISTOR& p);
  
public:
  explicit EVAL_BM_SEMI_RESISTOR(int c=0);
  
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
private:
  static double const _default_narrow;
  static double const _default_defw;
  static double const _default_tc1;
  static double const _default_tc2;
  static double const _default_defl;
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
  int param_count()const {return (4 + MODEL_CARD::param_count());}
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
  COMMON_COMPONENT* new_common()const {return new EVAL_BM_SEMI_CAPACITOR;}
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
private:
  static double const _default_rsh;
private:
  explicit MODEL_SEMI_RESISTOR(const MODEL_SEMI_RESISTOR& p);
public:
  explicit MODEL_SEMI_RESISTOR();
private: // override virtual
  std::string dev_type()const		{return "r";}
  void  precalc_first();
  //void  precalc_last();
  COMMON_COMPONENT* new_common()const {return new EVAL_BM_SEMI_RESISTOR;}
  CARD* clone()const		{return new MODEL_SEMI_RESISTOR(*this);}
  void		set_param_by_index(int, std::string&, int);
  bool		param_is_printable(int)const;
  std::string	param_name(int)const;
  std::string	param_name(int,int)const;
  std::string	param_value(int)const;
  int param_count()const {return (1 + MODEL_SEMI_BASE::param_count());}
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
/*--------------------------------------------------------------------------*/
EVAL_BM_SEMI_RESISTOR::EVAL_BM_SEMI_RESISTOR(int c)
  :EVAL_BM_SEMI_BASE(c),
   _resistance(_default_resistance),
   _capacitance(_default_capacitance),
   _tc1(_default_tc1),
   _tc2(_default_tc2)
{
}
/*--------------------------------------------------------------------------*/
EVAL_BM_SEMI_RESISTOR::EVAL_BM_SEMI_RESISTOR(const EVAL_BM_SEMI_RESISTOR& p)
  :EVAL_BM_SEMI_BASE(p),
   _resistance(p._resistance),
   _capacitance(p._capacitance),
   _tc1(p._tc1),
   _tc2(p._tc2)
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
  assert(Scope);
  EVAL_BM_SEMI_BASE::precalc_last(Scope);
  _resistance.e_val(_default_resistance, Scope);
  _tc1.e_val(_default_tc1, Scope);
  _tc2.e_val(_default_tc2, Scope);
  _capacitance.e_val(_default_capacitance, Scope);

  const MODEL_SEMI_RESISTOR* m = prechecked_cast<const MODEL_SEMI_RESISTOR*>(model());

  double tc1= (_tc1 == NOT_INPUT ) ? m->_tc1: _tc1;
  double tc2= (_tc2 == NOT_INPUT ) ? m->_tc2: _tc2;

  double eff_width;
  double eff_length;   

  _value = BIGBIG;   // if something will go wrong - very large resistance 

  if (_resistance != NOT_INPUT) {       //resistance is specified
    _value=_resistance;
    }
  else {
    double width    = (_width == NOT_INPUT)    ? m->_defw : _width;
    eff_width       = width - m->_narrow;
    double length   = (_length == NOT_INPUT)  ? m->_defl : _length;    
    eff_length      = _length - m->_narrow;      
    
    if (eff_width != 0.) {
      _value = m->_rsh * eff_length / eff_width;
      }
    else{
      _value = BIGBIG;
      }
    
    if (eff_width <= 0.) {untested();
        throw Exception_Precalc(modelname() + ": effective width is negative or zero\n");
    }else{
    }
    if (eff_length <= 0.) {
        throw Exception_Precalc(modelname() + ": effective length is negative or zero\n");
    }else{
    }
    if (m->_rsh == NOT_INPUT) {
        throw Exception_Precalc(modelname() + ": resistance is not set but needed \n");
    }else{
    }
      
   }
  
  double tempdiff = (_temp_c - m->_tnom_c);
  _value *= 1 + tc1*tempdiff + tc2*tempdiff*tempdiff;

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
  print_pair(o, lang, "c", _capacitance);
  print_pair(o, lang, "tc1", _tc1);
  print_pair(o, lang, "tc2", _tc2);
  
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
double const MODEL_SEMI_BASE::_default_narrow = 0.;
double const MODEL_SEMI_BASE::_default_defw = 1e-6;
double const MODEL_SEMI_BASE::_default_tc1 = 0.;
double const MODEL_SEMI_BASE::_default_tc2 = 0.;
double const MODEL_SEMI_BASE::_default_defl = 0.;
/*--------------------------------------------------------------------------*/
MODEL_SEMI_BASE::MODEL_SEMI_BASE()
  :MODEL_CARD(NULL),
   _narrow(_default_narrow),
   _defw(_default_defw),
   _tc1(_default_tc1),
   _tc2(_default_tc2),
   _defl(_default_defl)
{
}
/*--------------------------------------------------------------------------*/
MODEL_SEMI_BASE::MODEL_SEMI_BASE(const MODEL_SEMI_BASE& p)
  :MODEL_CARD(p),
   _narrow(p._narrow),
   _defw(p._defw),
   _tc1(p._tc1),
   _tc2(p._tc2),
   _defl(p._defl)  
{
}
/*--------------------------------------------------------------------------*/
void MODEL_SEMI_BASE::set_param_by_index(int i, std::string& value, int offset)
{
  switch (MODEL_SEMI_BASE::param_count() - 1 - i) {
  case 0: _narrow = value; break;
  case 1: _defw = value; break;
  case 2: _tc1 = value; break;
  case 3: _tc2 = value; break;
  case 4: _defl = value; break;
  
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
  case 4: return true;
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
double const MODEL_SEMI_RESISTOR::_default_rsh = NOT_INPUT;
/*--------------------------------------------------------------------------*/
MODEL_SEMI_RESISTOR::MODEL_SEMI_RESISTOR()
  :MODEL_SEMI_BASE(),
   _rsh(_default_rsh)
{
}
/*--------------------------------------------------------------------------*/
MODEL_SEMI_RESISTOR::MODEL_SEMI_RESISTOR(const MODEL_SEMI_RESISTOR& p)
  :MODEL_SEMI_BASE(p),  // GS - it was :MODEL_SEMI_BASE(),
   _rsh(p._rsh)
{
}
/*--------------------------------------------------------------------------*/
void MODEL_SEMI_RESISTOR::set_param_by_index(int i, std::string& value, int offset)
{
  switch (MODEL_SEMI_RESISTOR::param_count() - 1 - i) {
  case 0: _rsh = value; break;
  default: MODEL_SEMI_BASE::set_param_by_index(i, value, offset); break;
  }
}
/*--------------------------------------------------------------------------*/
bool MODEL_SEMI_RESISTOR::param_is_printable(int i)const
{
  switch (MODEL_SEMI_RESISTOR::param_count() - 1 - i) {
  case 0: return true;
  default: return MODEL_SEMI_BASE::param_is_printable(i);
  }
}
/*--------------------------------------------------------------------------*/
std::string MODEL_SEMI_RESISTOR::param_name(int i)const
{
  switch (MODEL_SEMI_RESISTOR::param_count() - 1 - i) {
  case 0: return "rsh";
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
  default: return MODEL_SEMI_BASE::param_value(i);
  }
}
/*--------------------------------------------------------------------------*/
void MODEL_SEMI_RESISTOR::precalc_first()
{
  MODEL_SEMI_BASE::precalc_first();

  const CARD_LIST* par_scope = scope();
  assert(par_scope);

  _rsh.e_val(_default_rsh, par_scope);
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
