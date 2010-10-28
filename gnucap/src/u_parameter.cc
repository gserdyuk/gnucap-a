/*$Id: u_parameter.cc,v 26.119 2009/09/09 13:27:53 al Exp $ -*- C++ -*-
 * Copyright (C) 2005 Albert Davis
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
 * A class for parameterized values
 * Used for .param statements
 * and passing arguments to models and subcircuits
 */
//testing=script,sparse 2006.07.14
#include "l_stlextra.h"
#include "u_parameter.h"
#include "u_lang.h"
#include "func_userdef.h"
/*
#include "u_function.h"
#include "ap.h"
#include "io_error.h"

//GS additon
class USERDEF_FUNCTION : public FUNCTION {
private:
    std::string _value;
    std::vector< std::string> _args;
public:
  USERDEF_FUNCTION(std::vector< std::string> arg, std::string val) { _args=arg; _value=val;}
  
  std::string eval(CS& fact_args, const CARD_LIST* Scope)const
  {
// evaluate position parameters in call using Scope - TBD

    std::vector<PARAMETER<double> >x(_args.size());
    PARAMETER<double> y;
    for (int i = 0; i<_args.size(); ++i)
        fact_args >> x[i]; 

    for (int i = 0; i<_args.size(); ++i)
        x[i].e_val(NOT_INPUT,Scope);

// create new environment - localScope - contans all position parameter only
    CARD_LIST localScope;
    PARAM_LIST* pl = localScope.params();
    
    for (int i=0; i<_args.size(); ++i) 
        pl->addParameter(_args[i],x[i]);

// create new command to evaluate
    CS::CS cm(CS::_STRING,_value);  // creates command with expression (_value is "1+3+a" - rhs of function) 
       
    Expression e(cm);             // Create expression
    //cm.check(bDanger, "expression syntax error");  // tbd   
    Expression r(e, &localScope); // (CARD_LIST*)NULL);  // reduced expression
    double res=r.eval();
    return to_string(res);
  }
} ;
*/
/*--------------------------------------------------------------------------*/

void PARAM_LIST::parse(CS& cmd)
{
  (cmd >> "real |integer "); // ignore type
  unsigned here = cmd.cursor();
  for (;;) {
    if (!(cmd.more() && (cmd.is_alpha() || cmd.match1('_')))) {
      break;
    }else{
    }
    std::string Name;
    PARAMETER<double> Value;
    // cmd >> Name >> '=' >> Value;
    cmd >> Name; 
    /* GS here split:
       if >>'=' - folow  current path
       if >>'(' - means .param name(....... - so read function
       in my unserstamndin, PARAM_LIST has to be split later onto 2 classes - old param_list and 
       func_list, or rather func_installer_list, which will manage function destruction when 
       function goes out of scope         
    */
    if (cmd.peek() == '=') {
      // parameter
      cmd >> '=' >> Value;;
      if (cmd.stuck(&here)) {untested();
        break;
      }else{
      }
      if (OPT::case_insensitive) {
        notstd::to_lower(&Name);
      }else{
      }
      _pl[Name] = Value;
    }
    else if (cmd.peek() == '('){
      // function
      std::vector<std::string> args;
      std::string a;
      
      cmd >> '(';
      
      if (OPT::case_insensitive) {
          notstd::to_lower(&Name);
      }else{
      }

      for(;;){       
        if (cmd.peek() == ')')  // todo: fix to handle absense of ")"
            break;
        cmd >> a;
        
        if (OPT::case_insensitive) {
            notstd::to_lower(&a);
        }else{
        }

        args.push_back(a);
        }
        
      cmd >> ')';
      cmd >> '=';
      cmd >> Value;

      // do not transform Value.string() to lowercase - it will be done in deep_search()
             
      // add here one function called; TODO - switch it out, use mstdout/OMSTREAM
      
      error(bDEBUG,"Userdef_function adding function: "+Name+"\n");

      USERDEF_FUNCTION *uf=new USERDEF_FUNCTION(args,Value.string());
      _fl[Name] = new DISPATCHER<FUNCTION>::INSTALL (&function_dispatcher, Name.c_str(), uf);
      _fal[Name]= args;  // tbd shall be string of arguments
    }
    else {
    /*issue*/
    }
  }
  cmd.check(bDANGER, "syntax error");

}
/*--------------------------------------------------------------------------*/
void PARAM_LIST::print(OMSTREAM& o, LANGUAGE* lang)const
{
  for (const_iterator i = _pl.begin(); i != _pl.end(); ++i) {
    if (i->second.has_hard_value()) {
      print_pair(o, lang, i->first, i->second);
    }else{
    }
  }
}
/*--------------------------------------------------------------------------*/

void PARAM_LIST::deep_print(OMSTREAM& o, LANGUAGE* lang, int i) const
{
  int j=i+1;
  print(o,lang);
  o<<"\n";
  const PARAM_LIST *pl=get_try_again();
  if(pl) pl->deep_print(o,lang,j);
  pl=get_upper_level_params();  
  if(pl) pl->deep_print(o,lang,j);
}

/*--------------------------------------------------------------------------*/
bool PARAM_LIST::is_printable(int i)const
{
  //BUG// ugly linear search
  int i_try = 0;
  for (const_iterator ii = _pl.begin(); ii != _pl.end(); ++ii) {
    if (i_try++ == i) {
      return ii->second.has_hard_value();
    }else{
    }
  }
  return false;
}
/*--------------------------------------------------------------------------*/
std::string PARAM_LIST::name(int i)const
{
  //BUG// ugly linear search
  int i_try = 0;
  for (const_iterator ii = _pl.begin(); ii != _pl.end(); ++ii) {
    if (i_try++ == i) {
      return ii->first;
    }else{
    }
  }
  return "";
}
/*--------------------------------------------------------------------------*/
std::string PARAM_LIST::value(int i)const
{
  //BUG// ugly linear search
  int i_try = 0;
  for (const_iterator ii = _pl.begin(); ii != _pl.end(); ++ii) {
    if (i_try++ == i) {
      return ii->second.string();
    }else{
    }
  }
  return "";
}
/*--------------------------------------------------------------------------*/
void PARAM_LIST::eval_copy(PARAM_LIST& p, const CARD_LIST* scope)
{
  assert(!_try_again);
  _try_again = p._try_again;
  assert(!_upper_level);
  _upper_level = p._upper_level;  

  for (iterator i = p._pl.begin(); i != p._pl.end(); ++i) {
    if (i->second.has_hard_value()) {
      if (_pl[i->first].has_hard_value()) {untested();
	_pl[i->first] = i->second.e_val(_pl[i->first], scope);
      }else{
	_pl[i->first] = i->second.e_val(NOT_INPUT, scope);
      }
    }else{
    }
  }
}
/*--------------------------------------------------------------------------*/
const PARAMETER<double>& PARAM_LIST::deep_lookup(std::string Name)const
{
// switching case is tricky.
// in case of regular parameter calculations everythiung works fine.
// but in case of function computation - list of parameters and name are handled in parse
// function body is handled due to this code
  if (OPT::case_insensitive) {
    notstd::to_lower(&Name);
  }else{
  }
  
  // normal branch - parhier=none or local
  if (OPT::parhier==parhNONE || parhLOCAL ) {
  PARAMETER<double> & rv = _pl[Name];
  if (rv.has_hard_value()) {
    // found a value, return it
    return rv;
    }
  
  if (_try_again) {
    // didn't find one, look in enclosing scope
    const PARAMETER<double> & rv2 = _try_again->deep_lookup(Name);
    if (rv2.has_hard_value() ) {
        return rv2;
        }
   }
   
  if (_upper_level) {
    // have to look at upper-level params
    const PARAMETER<double> & rv3 = get_upper_level_params()->deep_lookup(Name);
    if (rv3.has_hard_value() ) {
        return rv3;
        }
   }
  
  // no enclosing scope to look in
  // really didn't find it, give up
  // return garbage value (NOT_INPUT)
  return rv;
  }

  // specific branch - parhier=global - can be inefficeint  
  else if (OPT::parhier= parhGLOBAL) {
    unreachable();  
    // to be implemented parhier=global
  }
  
  else {
    unreachable();
  }
  
}
/*--------------------------------------------------------------------------*/
void PARAM_LIST::set(std::string Name, const std::string& Value)
{
  if (OPT::case_insensitive) {
    notstd::to_lower(&Name);
  }else{
  }
  _pl[Name] = Value;
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
bool Get(CS& cmd, const std::string& key, PARAMETER<bool>* val)
{
  if (cmd.umatch(key + ' ')) {
    if (cmd.skip1b('=')) {
      cmd >> *val;
    }else{
      *val = true;
    }
    return true;
  }else if (cmd.umatch("no" + key)) {
    *val = false;
    return true;
  }else{
    return false;
  }
}
/*--------------------------------------------------------------------------*/
bool Get(CS& cmd, const std::string& key, PARAMETER<int>* val)
{
  if (cmd.umatch(key + " {=}")) {
    *val = int(cmd.ctof());
    return true;
  }else{
    return false;
  }
}
/*--------------------------------------------------------------------------*/
INTERFACE bool GetToken(CS& cmd, PARAMETER<double>* val)
{
  if (cmd.is_float()) {                 // number
    *val = cmd.skipbl().ctof();
    return true;
  }
  else if (cmd.skipbl().match1("'")){   // expression 
    cmd>>*val;
    return true;
  }
  else{
    return false;
  }

}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
//GS used in lookup_solve to eliminate single quotes in expressions.
std::string s_char_subst(std::string s, char a, char b){
    for (int i=0; i<=s.size(); i++)
        if (s[i]==a)
            s[i]=b;
    return s;
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
