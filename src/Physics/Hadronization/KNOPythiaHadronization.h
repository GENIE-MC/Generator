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

#include "Framework/Interaction/Interaction.h"
#include "Framework/EventGen/EventRecordVisitorI.h"

namespace genie {

class KNOPythiaHadronization : protected EventRecordVisitorI {

public:

  KNOPythiaHadronization();
  KNOPythiaHadronization(string config);
  virtual ~KNOPythiaHadronization();

  // Implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event) const;

  // Overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  virtual void Configure(const Registry & config);
  virtual void Configure(string config);

private:

  void LoadConfig (void);

  const EventRecordVisitorI * SelectHadronizer(const Interaction *) const;

  //-- configuration

  const EventRecordVisitorI * fKNOHadronizer;    ///< KNO Hadronizer
  const EventRecordVisitorI * fPythiaHadronizer; ///< PYTHIA Hadronizer

  int    fMethod;       ///< KNO -> PYTHIA transition method
  double fWminTrWindow; ///< min W in transition region (pure KNO    < Wmin)
  double fWmaxTrWindow; ///< max W in transition region (pure PYTHIA > Wmax)
};

}         // genie namespace

#endif    // _KNO_PYTHIA_HADRONIZATION_H_
