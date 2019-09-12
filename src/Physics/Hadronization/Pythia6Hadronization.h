//____________________________________________________________________________
/*!

\class    genie::Pythia6Hadronization

\brief    Provides access to the PYTHIA hadronization models. \n
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  August 17, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _PYTHIA6_HADRONIZATION_H_
#define _PYTHIA6_HADRONIZATION_H_

#include <TPythia6.h>

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/ParticleData/PDGCodeList.h" 

namespace genie {

class GHepParticle;

class Pythia6Hadronization : protected EventRecordVisitorI {

public:
  Pythia6Hadronization();
  Pythia6Hadronization(string config);
  virtual ~Pythia6Hadronization();

  // Implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void           Initialize         (void)                    const;
  TClonesArray * Hadronize          (const Interaction*)      const;
  double         Weight             (void)                    const;
  PDGCodeList *  SelectParticles    (const Interaction*)      const;
  TH1D *         MultiplicityProb   (const Interaction*)      const;

  bool           AssertValidity     (const Interaction * i)   const;
  double         MaxMult            (const Interaction * i)   const;
  TH1D *         CreateMultProbHist (double maxmult)          const;
  double         Wmin               (void)                    const;

  void LoadConfig                   (void);

  mutable TPythia6 * fPythia;   ///< PYTHIA6 wrapper class

  //-- configuration parameters
  //   Note: additional configuration parameters common to all hadronizers
  //   (Wcut,Rijk,...) are declared one layer down in the inheritance tree
  double fSSBarSuppression;   ///< ssbar suppression
  double fGaussianPt2;        ///< gaussian pt2 distribution width
  double fNonGaussianPt2Tail; ///< non gaussian pt2 tail parameterization
  double fRemainingECutoff;   ///< remaining E cutoff for stopping fragmentation
  double fDiQuarkSuppression; ///< di-quark suppression parameter
  double fLightVMesonSuppression; ///< Light vector meon suppression
  double fSVMesonSuppression ; ///< Strange vector meson suppression
  double fLunda; ///< Lund a parameter
  double fLundb; ///< Lund b parameter
  double fLundaDiq; ///< adjustment of Lund a for di-quark
};

}         // genie namespace

#endif    // _PYTHIA6_HADRONIZATION_H_
