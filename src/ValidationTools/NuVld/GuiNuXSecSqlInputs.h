//_____________________________________________________________________________
/*!

\class    genie::nuvld::GuiNuXSecSqlInputs

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  August 25, 2004
*/
//_____________________________________________________________________________

#ifndef _V_XSEC_SQL_INPUTS_H_
#define _V_XSEC_SQL_INPUTS_H_

#include <string>

using std::string;

namespace genie {
namespace nuvld {

typedef struct 
{
  string _experiments;
  string _xsecs;
  string _nus;
  string _targets;

} GuiNuXSecSqlInputs_t ;

} // nuvld namespace
} // genie namespace

#endif
