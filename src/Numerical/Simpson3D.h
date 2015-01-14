//____________________________________________________________________________
/*!

\class    genie::Simpson3D

\brief    The 3-D extended Simpson rule (an open integration formula). The
          algorithm which is a direct extension of Simpson2D in 3-D, evaluates
          the numerical err and keeps improving its numerical estimate until it
          converges to the true value within some predefined margin of numerical
          accuracy.

\author   Chris Marshall <marshall \at pas.rochester.edu>
          University of Rochester

\created  March 20, 2014

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _SIMPSON_3D_H_
#define _SIMPSON_3D_H_

#include "Numerical/IntegratorI.h"
#include "Numerical/GridSpacing.h"

namespace genie {

class GSFunc;
class FunctionMap;

class Simpson3D: public IntegratorI
{
public:
  Simpson3D();
  Simpson3D(string config);
  virtual ~Simpson3D();

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
  unsigned int  fNBinsD2;
};

}        // genie namespace
#endif   // _SIMPSON_3D_H_

