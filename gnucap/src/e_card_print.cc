/*$Id: e_card.cc,v 26.96 2008/10/09 05:36:27 al Exp $ -*- C++ -*-
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
 * Base class for "cards" in the circuit description file
 */
//testing=script 2006.07.12
#include <stdio.h>
#include <string.h>
#include "e_card.h"
#include "e_cardlist.h"
#include "e_node.h"
using namespace std;
/*--------------------------------------------------------------------------*/
void CARD::print(std::string a){
cout<<"--CARD::print; \n";
cout<<"label="<<short_label()<<"\n";
//cout<<"ext nodes="<<ext_nodes()<<"\n";   - moved in e_compon
//cout<<"int nodes="<<int_nodes()<<"\n";
//cout<<"is_2port="<<is_2port()<<"\n";
//cout<<"is_source="<<is_source()<<"\n";
//cout<<"has_inode="<<has_inode()<<"\n";
cout<<"_net_nodes="<<_net_nodes<<"\n";
for (int i=0; i<_net_nodes; i++){
	cout<<" node_"<<i<<" = "<<_n[i].m_()<<"\n";		//????
	}
}
/*--------------------------------------------------------------------------*/
void CARD_LIST::print(std::string a){
cout<<"CARD_LIST::print; point: "+a+" \n";
for (iterator c=begin(); c!=end(); ++c){
	(*c)->print("");
	}
cout<<"CARD_LIST::print; end \n";
}

/*--------------------------------------------------------------------------*/
void CARD_LIST::print(std::string a, char *file, int line){
cout<<"CARD_LIST::print; point: "<<a<<" File: "<<file<<" line "<<line<<" \n";
for (iterator c=begin(); c!=end(); ++c){
	(*c)->print("");
}
cout<<"CARD_LIST::print; end \n";
}

/*--------------------------------------------------------------------------*/
