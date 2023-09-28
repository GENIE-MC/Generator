//____________________________________________________________________________
/*!

\class    genie::PXSecOnElectron

\brief    nu/nubar + e- scattering differential cross section \n
          The cross section algorithm interfaces nu+e cross sections
          Final states are handled by sub modules 

          Is a concrete implementation of the XSecAlgorithmI interface. \n

\author   Brinden Carlson bcarlson1@ufl.edu
          University of Florida

\created  June 1, 2023

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _NU_ON_ELECTRON_PARTIAL_XSEC_H_
#define _NU_ON_ELECTRON_PARTIAL_XSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/NuclearState/ElectronVelocity.h"

namespace genie {

class IntegratorI;
class XSecIntegratorI;

class PXSecOnElectron : public XSecAlgorithmI {

public:
  PXSecOnElectron();
  PXSecOnElectron(string name,string config);
  virtual ~PXSecOnElectron();

  //-- XSecAlgorithmI interface implementation
  virtual double XSec    (const Interaction * i, KinePhaseSpace_t k) const override = 0;
  double Integral        (const Interaction * i) const final;
  bool   ValidProcess    (const Interaction * i) const override;
  bool   ValidKinematics (const Interaction * i) const override;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config) final;
  void Configure(string config) final;

private:

  
protected:
  virtual void LoadConfig (void);

  const XSecIntegratorI * fXSecIntegrator;
  const ElectronVelocity * fElectronVelocity;

  int fNIntegration; //Max number of integration samples
  double fErrTolerance; //Error tolerance acceptable before returning cross section average
};

}       // genie namespace
#endif  // _NU_ON_ELECTRON_PARTIAL_XSEC_H_

