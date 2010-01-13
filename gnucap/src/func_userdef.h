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

#ifndef FUNC_USERDEF_H
#define FUNC_uSERDEF_H

//#include "u_parameter.h"
//#include "u_lang.h"

#include "u_parameter.h"
#include "u_function.h"
#include "ap.h"

/*--------------------------------------------------------------------------*/
//GS additon
/*--------------------------------------------------------------------------*/
class USERDEF_FUNCTION : public FUNCTION {
private:
    std::string _value;
    std::vector< std::string> _args;
public:
  USERDEF_FUNCTION(std::vector< std::string> arg, std::string val) { _args=arg; _value=val;}
  std::string eval(CS& fact_args, const CARD_LIST* Scope)const;
} ;
/*--------------------------------------------------------------------------*/

#endif
