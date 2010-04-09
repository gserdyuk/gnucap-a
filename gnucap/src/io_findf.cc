/*$Id: io_findf.cc,v 26.81 2008/05/27 05:34:00 al Exp $ -*- C++ -*-
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
 * Modified by AD.  Sent to me by C-WARE
 * This file contains the routine to locate a file,
 *	using a path string for the directories to search.
 * Interface:
 *	findfile(filename, paths, mode)
 *	    filename is the name of the file to be searched for,
 *	    paths is the path to follow to find it.
 *	    mode is how you want to open the file
 *	returns full path name, if successful, else "".
 *
 * PATHSEP, ENDDIR are system dependent, defined in md.h
 */
//testing=script,sparse 2006.07.17
#include "l_lib.h"
/*--------------------------------------------------------------------------*/
std::string findfile(const std::string& filename, const std::string& path,
		     int mode)
{
#ifdef CHECK_LOCAL_FIRST
  if (OS::access_ok(filename, mode)) {
    untested(); 
    return filename;
  }else{
    untested(); 
  }
#endif
					// for each item in the path
  for (std::string::const_iterator
	 p_ptr=path.begin(); p_ptr!=path.end(); ++p_ptr) {
    // p_ptr changed internally                  ^^^^^ skip sep
    std::string target = "";
    while (*p_ptr != PATHSEP  &&  p_ptr != path.end()) { // copy 1 path item
      target += *p_ptr++;
    }
    if (!target.empty() &&  !strchr(ENDDIR,p_ptr[-1])) {
      target += *ENDDIR;		// append '/' if needed
    }else{untested();
    }
    
    target += filename;
    if (OS::access_ok(target, mode)) {untested();	// found it
      return target;
    }else if (p_ptr==path.end()) {	// ran out of path, didn't find it
      return "";
    }else{				// else try again
    }
  }
  return ""; // path doesn't exist - didn't go thru loop at all
}
/*--------------------------------------------------------------------------*/
// find files if it is not abs-path , looks in the list of locations
// if found - returns actuial full path
// if not found - return file name as it was initially passed

std::string findfile_paths(const std::string& filename, std::vector<std::string> paths, int mode){

// check original filename and return if it is asolute path  
  if (filename[0]=='/') 
    return filename;

//other - call findfile:  
  for (int i=0; i<paths.size(); i++) {
    if (paths[i]!="")  {
        std::string retname = findfile(filename, paths[i], mode);
        if (retname!="")
            return retname; // did found something - returning
    }
}
  
   return filename;  // did not find anything in both searches - returning initial filename
}
/*--------------------------------------------------------------------------*/
// expands $VARIABLE according to environment
std::string expand_filename(std::string filename){

  std::string delim(NAMEDELIM);  
  while(true){
    // look for $
    // look for nonalphanum after $
    //all what we found $.... (till nonalphanum) is variable name
    // look for such variable name and substitute
    std:size_t pos1=filename.find_first_of("$",0);
    if (pos1==std::string::npos) // this will happen finally, return is granted
      return filename;
  
    // get environment variable
    std::size_t pos2=filename.find_first_of(delim,pos1);
    std::string name=filename.substr(pos1+1,pos2-(pos1+1));  // if pos2==npos - good too  
    char * env=getenv(name.c_str());
    std::string env_subst;
    if (env==NULL)
      env_subst="";
    else
      env_subst=std::string(env);
      
    filename = filename.substr(0,pos1)+env_subst+filename.substr(pos2);
  }
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
