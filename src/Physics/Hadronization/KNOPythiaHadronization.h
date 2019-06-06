//____________________________________________________________________________
/*!

\class    genie::KNOPythiaHadronization

\brief    A 'composite' hadronization model using a KNO-based hadronization
          model at low W and PYTHIA/JETSET at higher W.
          Contains no new hadronization code but merely a configurable KNO to
          PYTHIA transition scheme.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  June 08, 2006

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _KNO_PYTHIA_HADRONIZATION_H_
#define _KNO_PYTHIA_HADRONIZATION_H_

#include "Physics/Hadronization/HadronizationModelI.h"

namespace genie {

class KNOPythiaHadronization : public HadronizationModelI {

public:

  KNOPythiaHadronization();
  KNOPythiaHadronization(string config);
  virtual ~KNOPythiaHadronization();

  //-- implement the HadronizationModelI interface
  void           Initialize       (void)                                 const;
  TClonesArray * Hadronize        (const Interaction* )                  const;
  double         Weight           (void)                                 const;
  PDGCodeList *  SelectParticles  (const Interaction*)                   const;
  TH1D *         MultiplicityProb (const Interaction*, Option_t* opt="") const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

 //-- private methods & mutable parameters

  void LoadConfig (void);
  const HadronizationModelI * SelectHadronizer(const Interaction *) const;

  mutable double fWeight; ///< weight for generated event

  //-- configuration

  const HadronizationModelI * fKNOHadronizer;    ///< KNO Hadronizer
  const HadronizationModelI * fPythiaHadronizer; ///< PYTHIA Hadronizer

  int    fMethod;       ///< KNO -> PYTHIA transition method
  double fWminTrWindow; ///< min W in transition region (pure KNO    < Wmin)
  double fWmaxTrWindow; ///< max W in transition region (pure PYTHIA > Wmax)
};

}         // genie namespace

#endif    // _KNO_PYTHIA_HADRONIZATION_H_

