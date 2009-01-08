//_____________________________________________________________________________
/*!

\class    genie::nuvld::e_diff_xsec_sql_inputs

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  August 25, 2004
*/
//_____________________________________________________________________________

#ifndef _E_DIFF_XSEC_SQL_INPUTS_H_
#define _E_DIFF_XSEC_SQL_INPUTS_H_

#include <string>

using std::string;

namespace genie {
namespace nuvld {

typedef struct
{
  string _experiments;
  string _targets;
  double _E_min;
  double _E_max;
  double _EP_min;
  double _EP_max;
  double _Theta_min;
  double _Theta_max;
  double _Q2_min;
  double _Q2_max;
  double _W2_min;
  double _W2_max;
  double _Nu_min;
  double _Nu_max;
  double _Gamma_min;
  double _Gamma_max;
  double _Epsilon_min;
  double _Epsilon_max;
  double _x_min;
  double _x_max;

} e_diff_xsec_sql_inputs_t;

} // nuvld namespace
} // genie namespace

#endif

