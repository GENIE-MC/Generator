//____________________________________________________________________________
/*!

\class    genie::Intranuke

\brief    Base class for INTRANUKE intranuclear hadron transport implementations

\created  September 15, 2026

\cpright  Copyright (c) 2003-2026, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org

*/
//____________________________________________________________________________

#ifndef _INTRANUKE_BASE_H_
#define _INTRANUKE_BASE_H_

#include <TGenPhaseSpace.h>

#include "Physics/NuclearState/NuclearModelI.h"

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/Conventions/GMode.h"

namespace genie {

class IntranukeBase : public EventRecordVisitorI {

public :

  virtual string GetINukeMode() const = 0;
  virtual string GetGenINukeMode() const = 0;

  // Setters used in reweighting
  void SetRemnA( int A ) = 0;
  void SetRemnZ( int Z ) = 0;

  double GetRemnA() const = 0;
  double GetRemnZ() const = 0;

  double GetR0() const = 0;
  double GetNR() const = 0;

  double GetDelRPion() const = 0;
  double GetDelRNucleon() const = 0;

  double GetNucRmvE() const = 0;
  double GetHadStep() const = 0;

  bool GetUseOset() const = 0;
  bool GetAltOset() const = 0;
  bool GetXsecNNCorr() const = 0;

};

}      // genie namespace

#endif // _INTRANUKE_BASE_H_
