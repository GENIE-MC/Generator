//____________________________________________________________________________
/*!

\class    genie::PXSecOnElectron

\brief    nu/nubar + e- scattering differential cross section \n
          The cross section algorithm interfaces nu+e cross sections
          Final states are handled by sub modules 

          Is a concrete implementation of the XSecAlgorithmI interface. \n

\ref      W.J.Marciano and Z.Parsa, Neutrino-electron scattering theory,
          J.Phys.G: Nucl.Part.Phys. 29 (2003) 2629-2645
\ref      Oleksandr Tomalak and Richard J. Hill, Theory of elastic neutrino-electron scattering,
          Phys. Rev. D 101, 033006 – Published 21 February 2020

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
  virtual double XSec            (const Interaction * i, KinePhaseSpace_t k) const = 0;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;
  bool   ValidKinematics (const Interaction * i) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  
protected:
  void LoadConfig (void);

  const XSecIntegratorI * fXSecIntegrator;
  const ElectronVelocity * fElectronVelocity;

  int fNIntegration; //Max number of integration samples
  double fErrTolerance; //Error tolerance acceptable before returning cross section average
};

}       // genie namespace
#endif  // _NU_ON_ELECTRON_PARTIAL_XSEC_H_

