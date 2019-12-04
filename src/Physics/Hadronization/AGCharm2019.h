//____________________________________________________________________________
/*!

\class    genie::AGCharm2019

\brief    Andreopoulos - Gallagher (AG) GENIE Charm Hadronization model.

          The model relies on empirical charm fragmentation and pT functions,
          as well as on experimentally-determined charm fractions, to produce
          the ID and 4-momentum of charmed hadron in charm production events.

          The remnant (non-charm) system is hadronised by a call to PYTHIA.

          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

          Hugh Gallagher <gallag@minos.phy.tufts.edu>
          Tufts University

\created  August 17, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _CHARM_HADRONIZATION_H_
#define _CHARM_HADRONIZATION_H_

#include <TGenPhaseSpace.h>

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/Interaction/Interaction.h"

class TPythia6;
class TF1;

namespace genie {

class Spline;
class FragmentationFunctionI;

class AGCharm2019 : public EventRecordVisitorI {

public:
  AGCharm2019();
  AGCharm2019(string config);
  virtual ~AGCharm2019();

  // Implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event) const;

  // Overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  //
  void Configure(const Registry & config);
  void Configure(string config);


private:

  void           LoadConfig          (void);
  void           Initialize          (void)                                    const ;
  TClonesArray * Hadronize           (const Interaction* )                     const ;
  int            GenerateCharmHadron (int nupdg, double EvLab)                 const ;

  double         Weight              (void)                                    const ;

  mutable TGenPhaseSpace fPhaseSpaceGenerator; ///< a phase space generator

  // Configuration parameters
  //
  bool                           fCharmOnly;   ///< don't hadronize non-charm blob
  TF1 *                          fCharmPT2pdf; ///< charm hadron pT^2 pdf
  const FragmentationFunctionI * fFragmFunc;   ///< charm hadron fragmentation func

  double fFracMaxEnergy ;                      ///< Maximum energy available for the Meson fractions 

  Spline *                       fD0FracSpl;   ///< nu charm fraction vs Ev: D0
  Spline *                       fDpFracSpl;   ///< nu charm fraction vs Ev: D+
  Spline *                       fDsFracSpl;   ///< nu charm fraction vs Ev: Ds+
  double                         fD0BarFrac;   ///< nubar \bar{D0} charm fraction
  double                         fDmFrac;      ///< nubar D- charm fraction
  mutable TPythia6 *             fPythia;      ///< remnant (non-charm) hadronizer
};

}         // genie namespace

#endif    // _CHARM_HADRONIZATION__H_

