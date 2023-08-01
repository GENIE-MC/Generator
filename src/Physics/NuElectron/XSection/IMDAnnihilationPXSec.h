//____________________________________________________________________________
/*!

\class    genie::IMDAnnihilationPXSec

\brief    nuebar + e- -> mu- + numubar [CC]
          scattering differential cross section \n

          Is a concrete implementation of the PXSecOnElectron interface. \n

\ref      W.J.Marciano and Z.Parsa, Neutrino-electron scattering theory,
          J.Phys.G: Nucl.Part.Phys. 29 (2003) 2629-2645

\author   Rosen Matev (r.matev@gmail.com)

          B. Carlson implemented changes to interface PXSecOnElectron

\created  October 3, 2011

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _IMD_ANNIHILATION_PXSEC_H_
#define _IMD_ANNIHILATION_PXSEC_H_

#include "Physics/NuElectron/XSection/PXSecOnElectron.h"

namespace genie {

class IntegratorI;
class XSecIntegratorI;

class IMDAnnihilationPXSec : public PXSecOnElectron {

public:
  IMDAnnihilationPXSec();
  IMDAnnihilationPXSec(string config);
  virtual ~IMDAnnihilationPXSec();

  //-- PXSecOnElectron interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const; 
  bool   ValidProcess    (const Interaction * i) const; 
  bool   ValidKinematics (const Interaction * i) const; 

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config); 
  void Configure(string config); 

private:
  void LoadConfig (void);

};

}       // genie namespace
#endif  // _IMD_ANNIHILATION_PXSEC_H_
