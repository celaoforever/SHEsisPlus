/****************************************************************************************[Dimacs.h]
Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and
associated documentation files (the "Software"), to deal in the Software without
restriction,
including without limitation the rights to use, copy, modify, merge, publish,
distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom
the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
**************************************************************************************************/

#ifndef Minisat_Dimacs_h
#define Minisat_Dimacs_h

#include <stdio.h>
#include <boost/algorithm/string.hpp>
#include "ParseUtils.h"
#include "SolverTypes.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
namespace Minisat {

//=================================================================================================
// DIMACS Parser:
template <class Solver>
static void readClause(std::string& line, Solver& S, vec<Lit>& lits) {
  int parsed_lit, var;
  std::vector<std::string> strs;
  lits.clear();
  boost::split(strs, line, boost::is_any_of("\t "));
  if (!(strs[strs.size() - 1] == "0")) {
    exit(-1);
  }
  for (int i = 0; i < strs.size() - 1; i++) {
    parsed_lit = atoi(strs[i].c_str());
    if (parsed_lit == 0) break;
    var = abs(parsed_lit) - 1;
    while (var >= S.nVars()) S.newVar();
    lits.push((parsed_lit > 0) ? mkLit(var) : ~mkLit(var));
  }
}
template <class Solver>
static void parse_DIMACS_main(std::stringstream& in, Solver& S) {
  vec<Lit> lits;
  int vars = 0;
  int clauses = 0;
  int cnt = 0;
  std::string buf;
  while (std::getline(in, buf, '\n')) {
    if (buf == "") continue;
    cnt++;
    readClause(buf, S, lits);
    S.addClause_(lits);
  }
}
}

template <class Solver>
static void parse_DIMACS(std::stringstream& input_stream, Solver& S) {
  parse_DIMACS_main(input_stream, S);
}

#endif
