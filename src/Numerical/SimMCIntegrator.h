//____________________________________________________________________________
/*!

\class    genie::SimMCIntegrator

\brief    A simple MC integrator for scalar functions in a multi-dimensional
          parameter space. Note that the accuracy in simple MC integration
          increases only as the sqrt(N).

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  May 03, 2004

\cpright  Copyright (c) 2003-2011, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _SIM_MC_INTEGRATOR_H_
#define _SIM_MC_INTEGRATOR_H_

#include "Numerical/IntegratorI.h"

namespace genie {

class GSFunc;

class SimMCIntegrator: public IntegratorI
{
public:
  SimMCIntegrator();
  SimMCIntegrator(string config);
  virtual ~SimMCIntegrator();

  //! implement the IntegratorI interface
  double Integrate(GSFunc & gsfunc) const;

  //! override the Algorithm::Configure methods to load configuration
  //!  data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:

  //! load config data to private data members
  void LoadConfigData (void);

  //! actual config data
  unsigned int  fNMax;       ///< max allowed number of func evaluations
  unsigned int  fNMin;       ///< compute integral/error only every so many func evals
  double        fMaxPcntErr; ///< max numerical error allowed (in %)
  bool          fInLog;
};

}        // genie namespace
#endif   // _SIM_MC_INTEGRATOR_H_
