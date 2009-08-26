/*$Id: ap_construct.cc,v 26.110 2009/05/28 15:32:04 al Exp $ -*- C++ -*-
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
 * construction, copy, etc.
 */
//testing=script,sparse 2006.07.17
#include "u_opt.h"
#include "ap.h"

#if defined(HAVE_LIBREADLINE)
  #include <readline/readline.h>
  #include <readline/history.h>
#endif
/*--------------------------------------------------------------------------*/
static std::string getlines(FILE*);
OMSTREAM mout; // > file bitmap //BUG//encapsulation
OMSTREAM mlog; // log file bitmap
/*--------------------------------------------------------------------------*/
CS::CS(CS::STDIN)
  :_file(stdin),
   _name(),
   _cmd(),
   _cnt(0),
   _length(0),
   _begin_match(0),
   _end_match(0),
   _ms(msEXACT),
   _ok(true)
{
}
/*--------------------------------------------------------------------------*/
CS::CS(CS::INC_FILE, const std::string& name)
  :_file(fopen(name.c_str(), "r")),
   _name(name),
   _cmd(),
   _cnt(0),
   _length(0),
   _begin_match(0),
   _end_match(0),
   _ms(msEXACT),
   _ok(true)
{
  if (!_file) {itested();
    throw Exception_File_Open(name + ':' + strerror(errno));
  }else{
  }
}
/*--------------------------------------------------------------------------*/
CS::CS(CS::WHOLE_FILE, const std::string& name)
  :_file(NULL),
   _name(name),
   _cmd(),
   _cnt(0),
   _length(0),
   _begin_match(0),
   _end_match(0),
   _ms(msEXACT),
   _ok(true)
{
  int f = open(name.c_str(), O_RDONLY);
  if (f == EOF) {itested();
    throw Exception_File_Open(name + ':' + strerror(errno));
  }else{
  }
  _length = static_cast<unsigned>(lseek(f, off_t(0), SEEK_END));
  lseek(f, off_t(0), SEEK_SET);

  char* cmd = new char[_length+2];
  read(f, cmd, _length);
  cmd[_length++] = '\0';
  _cmd = cmd;

  close(f);
}
/*--------------------------------------------------------------------------*/
CS::CS(CS::STRING, const std::string& s)
  :_file(NULL),
   _name(),
   _cmd(s),
   _cnt(0),
   _length(static_cast<unsigned>(s.length())),
   _begin_match(0),
   _end_match(0),
   _ms(msEXACT),
   _ok(true)
{
}
/*--------------------------------------------------------------------------*/
#if 0
CS::CS(const CS& p)
  :_file(NULL),
   _name(p._name),
   _cmd(p._cmd),
   _cnt(p._cnt),
   _length(p._length),
   _begin_match(0),
   _end_match(0),
   _ms(p._ms),
   _ok(p._ok)
{untested();
}
#endif
/*--------------------------------------------------------------------------*/
CS& CS::operator=(const std::string& s)
{untested();
  assert(!_file);
  _cmd = s;
  _cnt = 0;
  _ok = true;
  _length = static_cast<unsigned>(s.length());
  return *this;
}
/*--------------------------------------------------------------------------*/
#if 0
CS& CS::operator=(const CS& p)
{untested();
  assert(&p != this);
  _name = p._name;
  _file = p._file;
  _cmd = p._cmd;
  _cnt = p._cnt;
  _ok = p._ok;
  _length = p._length;
  return *this;
}
#endif
/*--------------------------------------------------------------------------*/
CS& CS::get_line(const std::string& prompt)
{
  if (is_file()) {
    _cmd = getlines(_file);
    _cnt = 0;
    _length = static_cast<unsigned>(_cmd.length());
    _ok = true;
  }else{itested();
    assert(_file == stdin);
    char cmdbuf[BUFLEN];
    getcmd(prompt.c_str(), cmdbuf, BUFLEN);
    _cmd = cmdbuf;
    _cnt = 0;
    _length = static_cast<unsigned>(_cmd.length());
    _ok = true;
  }

  if (OPT::listing) {
    IO::mstdout << "\"" << fullstring() << "\"\n";
  }else{
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
/* getcmd: get a command.
 * if "fin" is stdin, display a prompt first.
 * Also, actually do logging, echo, etc.
 */
char *getcmd(const char *prompt, char *buffer, int buflen)
{
  assert(prompt);
  assert(buffer);
  if (isatty(fileno(stdin))) {
    // stdin is keyboard
#if defined(HAVE_LIBREADLINE)
    if (OPT::edit) {
      char* line_read = readline(prompt);
      if (!line_read) {itested();
	throw Exception_End_Of_Input("EOF on stdin");
      }else{
      }
      // readline gets a new buffer every time, so copy it to where we want it
      char* end_of_line = (char*)memccpy(buffer, line_read, 0, static_cast<size_t>(buflen-1));
      if (!end_of_line) {
	buffer[buflen-1] = '\0';
      }else{
	*end_of_line = '\0';
      }
      free(line_read);
      
      if (*buffer) {
	add_history(buffer);
      }else{
      }
    }else
#endif
      {
	IO::mstdout << prompt;	/* prompt & flush buffer */
	if (!fgets(buffer, buflen, stdin)) {untested();	/* read line */
	  throw Exception_End_Of_Input("EOF on stdin");
	}else{
	}
      }
    (IO::mstdout - mout) << '\r';	/* reset col counter */
    trim(buffer);
    (mlog + mout) << buffer << '\n';
    return buffer;
  }else{
    // stdin is file
    if (!fgets(buffer, buflen, stdin)) {itested();	/* read line */
      throw Exception_End_Of_Input("EOF on stdin");
    }else{
    }
    trim(buffer);
    (mlog + mout) << buffer << '\n';
    return buffer;
  }
}
/*--------------------------------------------------------------------------*/
static std::string getlines(FILE *fileptr)
{
  assert(fileptr);
  const int buffer_size = BIGBUFLEN;
  std::string s;

  bool more = true;  // get another line (extend)
  while (more) {
    char buffer[buffer_size+1];
    char* got = fgets(buffer, buffer_size, fileptr);
    if (!got) { // probably end of file
      more = false;
      if (s == "") {
	throw Exception_End_Of_Input("");
      }else{untested();
      }
    }else{
      trim(buffer);
      size_t count = strlen(buffer);
      if (buffer[count-1] == '\\') {
	itested();
	buffer[count-1] = '\0';
      }else{
	int c = fgetc(fileptr);
	if (c == '+') {
	  more = true;
	}else if (c == '\n') {
	  more = true;
	  ungetc(c,fileptr);
	}else{
	  more = false;
	  ungetc(c,fileptr);
	}
      }
      s += buffer;
      s += ' ';
    }
  }
  return s;
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
