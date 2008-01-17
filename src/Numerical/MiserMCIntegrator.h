//____________________________________________________________________________
/*!

\class    genie::MiserMCIntegrator

\brief    The MISER adaptive MC integration algorithm based on recursive
          stratified sampling.

\ref      Numerical Recipes in C, Cambridge Univ. Press, 2002, page 316-328

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  February 23, 2006

\cpright  Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _MISER_MC_INTEGRATOR_H_
#define _MISER_MC_INTEGRATOR_H_

#include "Numerical/IntegratorI.h"

namespace genie {

class GSFunc;

class MiserMCIntegrator: public IntegratorI
{
public:
  MiserMCIntegrator();
  MiserMCIntegrator(string config);
  virtual ~MiserMCIntegrator();

  //! implement the IntegratorI interface
  double Integrate(GSFunc & gsfunc) const;

  //! override the Algorithm::Configure methods to load configuration
  //!  data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:

};

}        // genie namespace
#endif   // _MISER_MC_INTEGRATOR_H_
