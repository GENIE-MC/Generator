//____________________________________________________________________________
/*!

\class    genie::AGCharmPythiaBaseHadro2023

\brief    Andreopoulos - Gallagher (AG) GENIE Charm Hadronization model.

          The model relies on empirical charm fragmentation and pT functions,
          as well as on experimentally-determined charm fractions, to produce
          the ID and 4-momentum of charmed hadron in charm production events.

          The remnant (non-charm) system is hadronised by a call to PYTHIA.

          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

          Hugh Gallagher <gallag@minos.phy.tufts.edu>
          Tufts University

\created  August 17, 2004

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _CHARM_HADRONIZATION_H_
#define _CHARM_HADRONIZATION_H_

#include <TGenPhaseSpace.h>

#include "Framework/Conventions/GBuild.h"
#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/Interaction/Interaction.h"

class TF1;

namespace genie {

class Spline;
class FragmentationFunctionI;

class AGCharmPythiaBaseHadro2023 : public EventRecordVisitorI {

public:

  // Implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event) const;

protected:

  AGCharmPythiaBaseHadro2023();
  AGCharmPythiaBaseHadro2023(string config);
  AGCharmPythiaBaseHadro2023(string name, string config);
  virtual ~AGCharmPythiaBaseHadro2023();

protected:
  void           LoadConfig          (void);
  void           Initialize          (void)                                    const ;

  virtual bool HadronizeRemnant(int qrkSyst1, int qrkSyst2, double WR, TLorentzVector p4R,
                                unsigned int& rpos, TClonesArray * particle_list)      const = 0;
  // needs to append to TClonesArray of GHepParticle starting at rpos

  // Overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  //
  void Configure(const Registry & config);
  void Configure(string config);

private:

  TClonesArray * Hadronize           (const Interaction* )                     const ;
  // returns TClonesArray of GHepParticle

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
};

}         // genie namespace

#endif    // _CHARM_HADRONIZATION__H_
