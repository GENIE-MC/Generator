//____________________________________________________________________________
/*!

\class    genie::MKFFEM

\brief    Electromagnetic form factors for MK SPP model


\author   Igor Kakorin <kakorin@jinr.ru>, Joint Institute for Nuclear Research \n
          based on code of
          Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  Nov 12, 2019

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _MK_EM_FORM_FACTOR_MODEL_H_
#define _MK_EM_FORM_FACTOR_MODEL_H_

#include "Physics/QuasiElastic/XSection/QELFormFactorsModelI.h"
#include "Physics/QuasiElastic/XSection/ELFormFactors.h"

namespace genie {

class ELFormFactorsModelI;

class MKFFEM : public QELFormFactorsModelI {

public:
   MKFFEM();
   MKFFEM(string name);
  ~MKFFEM();

  double F1V     (const Interaction * interaction) const;
  double xiF2V   (const Interaction * interaction) const;
  double FA      (const Interaction * interaction) const;
  double Fp      (const Interaction * interaction) const;

  // Overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void LoadConfig (void);

  double tau    (const Interaction * interaction) const;

  double F1P    (const Interaction * interaction) const;  
  double F2P    (const Interaction * interaction) const;  
  double F1N    (const Interaction * interaction) const;  
  double F2N    (const Interaction * interaction) const;  
   
  const ELFormFactorsModelI   * fElFFModel;

  mutable ELFormFactors   fELFF;

};

}       // genie namespace

#endif  

