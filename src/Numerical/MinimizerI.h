//____________________________________________________________________________
/*!

\class    genie::MinimizerI

\brief    Numerical minimization/maximazation algorithm ABC

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  July 01, 2004

*/
//____________________________________________________________________________

#ifndef _MINIMIZER_I_H_
#define _MINIMIZER_I_H_

#include <vector>

#include "Algorithm/Algorithm.h"

using std::vector;

namespace genie {

class GSFunc;

class MinimizerI : public Algorithm
{
public:
  virtual ~MinimizerI();

  virtual void Minimize(
       GSFunc & gsfunc, vector<double> & x, bool min) const = 0;

protected:
  MinimizerI();
  MinimizerI(string name);
  MinimizerI(string name, string config);
};

}        // namespace
#endif   // _MINIMIZER_I_H_
