//____________________________________________________________________________
/*!

\class    genie::RESXSec

\brief    Computes the RES Cross Section.

          Is a concrete implementation of the XSecAlgorithmI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 04, 2004
 
*/
//____________________________________________________________________________

#ifndef _RES_XSEC_H_
#define _RES_XSEC_H_

#include "Base/XSecAlgorithmI.h"
#include "Utils/Range1.h"

namespace genie {

class IntegratorI;

class RESXSec : public XSecAlgorithmI {

public:

  RESXSec();
  RESXSec(const char * param_set);
  virtual ~RESXSec();

  //-- XSecAlgorithmI interface implementation

  double XSec (const Interaction * interaction) const;

private:

  const IntegratorI * Integrator (void) const;

  Range1D_t WRange  (const Interaction * interaction) const;
  Range1D_t Q2Range (const Interaction * interaction) const;

  int    NW     (void) const;
  int    NLogQ2 (void) const;
  double Wmin   (void) const;
  double Wmax   (void) const;
  double Q2min  (void) const;
  double Q2max  (void) const;
};

}       // genie namespace

#endif  // _RES_XSEC_H_
