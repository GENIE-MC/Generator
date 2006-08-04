//____________________________________________________________________________
/*!

\class    genie::PythiaHadronization

\brief    Provides access to the PYTHIA hadronization models. \n
          Is a concrete implementation of the HadronizationModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  August 17, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _PYTHIA_HADRONIZATION_H_
#define _PYTHIA_HADRONIZATION_H_

#include <TPythia6.h>

#include "Fragmentation/HadronizationModelBase.h"

namespace genie {

class PythiaHadronization : public HadronizationModelBase {

public:

  PythiaHadronization();
  PythiaHadronization(string config);
  virtual ~PythiaHadronization();

  //! implement the HadronizationModelI interface
  void           Initialize       (void)                                  const;
  TClonesArray * Hadronize        (const Interaction*)                    const;
  double         Weight           (void)                                  const;
  PDGCodeList *  SelectParticles  (const Interaction*)                    const;
  TH1D *         MultiplicityProb (const Interaction*, Option_t* opt="")  const;

  //! overload the Algorithm::Configure() methods to load private data
  //! members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void LoadConfig     (void);
  bool AssertValidity (const Interaction * i) const;
  void SyncSeeds      (void) const;

  mutable TPythia6 * fPythia;   ///< PYTHIA6 wrapper class
  mutable long int   fCurrSeed; ///< always keep PYTHIA's & GENIE's seeds in sync

  //! configuration parameters
  //! Note: additional configuration parameters common to all hadronizers 
  //! (Wcut,Rijk,...) are declared one layer down in the inheritance tree
  double fSSBarSuppression;   ///< ssbar suppression
  double fGaussianPt2;        ///< gaussian pt2 distribution width
  double fNonGaussianPt2Tail; ///< non gaussian pt2 tail parameterization
  double fRemainingECutoff;   ///< remaining E cutoff for stopping fragmentation
};

}         // genie namespace

#endif    // _PYTHIA_HADRONIZATION__H_

