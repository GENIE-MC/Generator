//____________________________________________________________________________
/*!

\class    genie::PythiaHadronization

\brief    Provides access to the PYTHIA hadronization models. \n
          Is a concrete implementation of the HadronizationModelI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  August 17, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _PYTHIA_HADRONIZATION_H_
#define _PYTHIA_HADRONIZATION_H_

#include <TPythia6.h>

#include "Physics/Hadronization/HadronizationModelBase.h"

namespace genie {

class DecayModelI;
class PythiaHadronization : public HadronizationModelBase {

public:
  PythiaHadronization();
  PythiaHadronization(string config);
  virtual ~PythiaHadronization();

  //-- implement the HadronizationModelI interface
  void           Initialize       (void)                                  const;
  TClonesArray * Hadronize        (const Interaction*)                    const;
  double         Weight           (void)                                  const;
  PDGCodeList *  SelectParticles  (const Interaction*)                    const;
  TH1D *         MultiplicityProb (const Interaction*, Option_t* opt="")  const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void LoadConfig     (void);
  bool AssertValidity (const Interaction * i) const;
/*
  void SwitchDecays   (int pdgc, bool on_off) const;
  void HandleDecays   (TClonesArray * plist) const;
*/
  mutable TPythia6 * fPythia;   ///< PYTHIA6 wrapper class

  const DecayModelI * fDecayer;

  //-- configuration parameters
  //   Note: additional configuration parameters common to all hadronizers 
  //   (Wcut,Rijk,...) are declared one layer down in the inheritance tree
  double fSSBarSuppression;   ///< ssbar suppression
  double fGaussianPt2;        ///< gaussian pt2 distribution width
  double fNonGaussianPt2Tail; ///< non gaussian pt2 tail parameterization
  double fRemainingECutoff;   ///< remaining E cutoff for stopping fragmentation
};

}         // genie namespace

#endif    // _PYTHIA_HADRONIZATION__H_

