/*$Id: func_userdef.cc
 * Copyright (C) 2008 Albert Davis
 * Author: Gennadiy Serdyuk
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
 */
//#include "u_parameter.h"
//#include "u_lang.h"
#include "func_userdef.h"
//#include "u_parameter.h"
//#include "u_function.h"
//#include "ap.h"
#include "io_error.h"

/*--------------------------------------------------------------------------*/
//GS additon
/*--------------------------------------------------------------------------*/

std::string USERDEF_FUNCTION::eval(CS& fact_args, const CARD_LIST* Scope)const
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
    CS cm(CS::_STRING,_value);  // creates command with expression (_value is "1+3+a" - rhs of function) 
       
    Expression e(cm);             // Create expression
    //cm.check(bDanger, "expression syntax error");  // tbd   
    Expression r(e, &localScope); // (CARD_LIST*)NULL);  // reduced expression
    double res=r.eval();
    return to_string(res);
  }
/*--------------------------------------------------------------------------*/

