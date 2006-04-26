//____________________________________________________________________________
/*!

\class    genie::KNOHadronization2

\brief    An improved KNO hadronization model.

          An improved version of the KNO hadronization model reproduces (in
          addition) the observed XF-(bkw/fwd) multiplicity asymmetries &
          multiplicity correlations. \n

          Is a concrete implementation of the HadronizationModelI interface.\n

          -- Test & currently broken version --

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  November 25, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _KNO_HADRONIZATION_2_H_
#define _KNO_HADRONIZATION_2_H_

#include <map>

#include "Fragmentation/HadronizationModelI.h"
#include "Fragmentation/Multiplicity.h"

using std::map;

namespace genie {

class KNOHadronization2 : public HadronizationModelI {

public:

  KNOHadronization2();
  KNOHadronization2(string config);
  virtual ~KNOHadronization2();

  //-- define KNOHadronization2 interface

  void           Initialize   (void)                 const;
  TClonesArray * Hadronize    (const Interaction * ) const;

private:

  const HadronizationModelI * OriginalHadronizationModel(void) const;

  map<Multiplicity_t, double>
                ExpectedMultiplicities(const Interaction * interaction) const;

  TClonesArray * DecayExclusive(
          TLorentzVector & sp4, TClonesArray * particle_list, bool fwd) const;
                                         
  double MultiplicityParam  (char prm, Multiplicity_t mt, int v, int N) const;
  void   BoostToHCMS        (TClonesArray * pl) const;                                         
  bool   ComputeForwardness (TClonesArray * pl, const Interaction * in) const;  
  double PKick              (TClonesArray * pl, const Interaction * in) const;
  int    NParticles         (TClonesArray * pl, bool fwd) const;
  double TotalMass          (TClonesArray * pl, bool fwd) const;               
  void   Cooling            (TClonesArray * pl) const;
  void   Cooling2(TClonesArray * particle_list, double W) const;

};

}         // genie namespace

#endif    // _KNO_HADRONIZATION_2_H_

