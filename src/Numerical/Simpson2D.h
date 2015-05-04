//____________________________________________________________________________
/*!

\class    genie::Simpson2D

\brief    The 2-D extended Simpson rule (an open integration formula). The
          algorithm which is a direct extension of Simpson1D in 2-D, evaluates
          the numerical err and keeps improving its numerical estimate until it
          converges to the true value within some predefined margin of numerical
          accuracy.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  May 03, 2004

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _SIMPSON_2D_H_
#define _SIMPSON_2D_H_

#include "Numerical/IntegratorI.h"
#include "Numerical/GridSpacing.h"

namespace genie {

class GSFunc;
class FunctionMap;

class Simpson2D: public IntegratorI
{
public:
  Simpson2D();
  Simpson2D(string config);
  virtual ~Simpson2D();

  //! implement the IntegratorI interface
  double Integrate(GSFunc & gsfunc) const;

  //! override the Algorithm::Configure methods to load configuration
  //!  data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:

  //! actual Simpson integration rule
  double SimpsonRule(FunctionMap & func_map) const;

  //! load config data to private data members
  void LoadConfigData (void);

  //! actual config data
  unsigned int  fIMaxConv;   ///< max number of iterations before converging
  unsigned int  fNo;         ///< 2^No + 1 is the initial number of steps
  double        fMaxPcntErr; ///< max numerical error allowed (in %)
  GridSpacing_t fSpacing;    ///< grid points spacing rule
  bool          fFastDensityIncrease;
  bool          fForceFixedNBins;
  unsigned int  fNBinsD0;
  unsigned int  fNBinsD1;
};

}        // genie namespace
#endif   // _SIMPSON_2D_H_

