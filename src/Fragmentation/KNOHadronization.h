//____________________________________________________________________________
/*!

\class    genie::KNOHadronization

\brief    The KNO hadronization model.
          This hadronization scheme is similar to the one originally used
          in NeuGEN by G.Barr, G.F.Pearce, H.Gallagher et al. \n
          Improvements were made by C.Andreopoulos, H.Gallagher, T.Yang. \n

          Both the 'historical' version and the new versions are supported.
          See the algorithms configuration for details.

          Is a concrete implementation of the HadronizationModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  August 17, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _KNO_HADRONIZATION_H_
#define _KNO_HADRONIZATION_H_

#include <vector>

#include <TGenPhaseSpace.h>

#include "Fragmentation/HadronizationModelI.h"

using std::vector;

class TF1;

namespace genie {

class MultiplicityProbModelI;
class DecayModelI;

class KNOHadronization : public HadronizationModelI {

public:

  KNOHadronization();
  KNOHadronization(string config);
  virtual ~KNOHadronization();

  //! implement the HadronizationModelI interface
  void           Initialize   (void)                 const;
  TClonesArray * Hadronize    (const Interaction * ) const;
  double         Weight       (void)                 const;

  //! overload the Algorithm::Configure() methods to load private data
  //! members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void          LoadConfig            (void);
  vector<int> * GenerateFSHadronCodes (int mult, int maxQ, double W) const;
  int           GenerateBaryonPdgCode (int mult, int maxQ)           const;
  int           HadronShowerCharge    (const Interaction * proc)     const;
  void          HandleDecays          (TClonesArray * particle_list) const;
  double        ReWeightPt2           (const vector<int> & pdgcv)    const;

  TClonesArray* DecayMethod1    (double W, const vector<int> & pdgv) const;
  TClonesArray* DecayMethod2    (double W, const vector<int> & pdgv) const;
  TClonesArray* DecayBackToBack (double W, const vector<int> & pdgv) const;

  void PhaseSpaceDecay(
             TClonesArray & pl, TLorentzVector & pd, 
	                     const vector<int> & pdgv, int offset=0) const;

  mutable TGenPhaseSpace fPhaseSpaceGenerator;
  mutable double         fWeight;

  //! configuration parameters

  const MultiplicityProbModelI * fMultProbModel;
  const DecayModelI *            fDecayer;

  double fPpi0;                ///< pi0 pi0 production probability
  double fPpic;                ///< pi+ pi- production probability
  double fPKc;                 ///< K+  K- production probability
  double fPK0;                 ///< K0  K0bar production probability
  bool   fUseBaryonXfPt2Param; ///< Generate baryon xF,pT2 from experimental parameterization?
  bool   fReWeightDecays;      ///< Reweight phase space decays to reproduce exp. pT^2 distributions?
  bool   fForceDecays;         ///< force decays of unstable hadrons produced?
  bool   fForceMinMult;        ///< force minimum multiplicity if (at low W) generated less?
  bool   fGenerateWeighted;    ///< generate weighted events?
  TF1 *  fBaryonXFpdf;         ///< baryon xF PDF
  TF1 *  fBaryonPT2pdf;        ///< baryon pT^2 PDF
};

}         // genie namespace

#endif    // _KNO_HADRONIZATION_H_

