//____________________________________________________________________________
/*!

\class    genie::AGCharmPythia8Hadro2023

\brief    Andreopoulos - Gallagher (AG) GENIE Charm Hadronization model.

          The model relies on empirical charm fragmentation and pT functions,
          as well as on experimentally-determined charm fractions, to produce
          the ID and 4-momentum of charmed hadron in charm production events.

          The remnant (non-charm) system is hadronised by a call to PYTHIA.

          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

          Hugh Gallagher <gallag@minos.phy.tufts.edu>
          Tufts University

\created  August 17, 2004

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _AGCHARM_PYTHIA8_HADRONIZATION_H_
#define _AGCHARM_PYTHIA8_HADRONIZATION_H_

#include "Physics/Hadronization/AGCharmPythiaBaseHadro2023.h"

#ifdef __GENIE_PYTHIA8_ENABLED__
#include "Pythia8/Pythia.h"
#endif

namespace genie {

class Spline;
class FragmentationFunctionI;

class AGCharmPythia8Hadro2023 : public AGCharmPythiaBaseHadro2023 {

public:
  AGCharmPythia8Hadro2023();
  AGCharmPythia8Hadro2023(string config);
  virtual ~AGCharmPythia8Hadro2023();


private:
  void           Initialize      (void)                                  const;
  bool           HadronizeRemnant(int qrkSyst1, int qrkSyst2, double WR, TLorentzVector p4R,
                                  unsigned int& rpos, TClonesArray * particle_list) const;


#ifdef __GENIE_PYTHIA8_ENABLED__
  mutable Pythia8::Pythia *             fPythia;      ///< remnant (non-charm) hadronizer
#endif
};

}         // genie namespace

#endif    // _AGCHARM_PYTHIA8_HADRONIZATION__H_
