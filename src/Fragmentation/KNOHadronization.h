//____________________________________________________________________________
/*!

\class    genie::KNOHadronization

\brief    The KNO hadronization model.
          This hadronization scheme originates from the NEUGEN's KNO model
          (G.Barr, G.F.Pearce, H.Gallagher et al.) \n
          Improvements were made by C.Andreopoulos, H.Gallagher, T.Yang. \n

          Both the 'historical' version and the new versions are supported.
          See the algorithms configuration for details.

          Is a concrete implementation of the HadronizationModelI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

          Hugh Gallagher <gallag@minos.phy.tufts.edu>
          Tufts University

          Tinjun Yang <tjyang@stanford.edu>
          Stanford University

\created  August 17, 2004

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _KNO_HADRONIZATION_H_
#define _KNO_HADRONIZATION_H_

#include <TGenPhaseSpace.h>

#include "Fragmentation/HadronizationModelBase.h"

class TF1;

namespace genie {

class DecayModelI;
//class Spline;

class KNOHadronization : public HadronizationModelBase {

public:
  KNOHadronization();
  KNOHadronization(string config);
  virtual ~KNOHadronization();

  // implement the HadronizationModelI interface
  void           Initialize       (void)                                    const;
  TClonesArray * Hadronize        (const Interaction* )                     const;
  double         Weight           (void)                                    const;
  PDGCodeList *  SelectParticles  (const Interaction*)                      const;
  TH1D *         MultiplicityProb (const Interaction*, Option_t* opt = "")  const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  // private methods & mutable parameters

  void          LoadConfig            (void);
  bool          AssertValidity        (const Interaction * i)        const;
  PDGCodeList * GenerateHadronCodes   (int mult, int maxQ, double W) const;
  int           GenerateBaryonPdgCode (int mult, int maxQ, double W) const;
  int           HadronShowerCharge    (const Interaction * )         const;
  double        KNO                   (int nu, int nuc, double z)    const;
  double        AverageChMult         (int nu, int nuc, double W)    const;
  void          HandleDecays          (TClonesArray * particle_list) const;
  double        ReWeightPt2           (const PDGCodeList & pdgcv)    const;

  TClonesArray* DecayMethod1    (double W, const PDGCodeList & pdgv, bool reweight_decays) const;
  TClonesArray* DecayMethod2    (double W, const PDGCodeList & pdgv, bool reweight_decays) const;
  TClonesArray* DecayBackToBack (double W, const PDGCodeList & pdgv) const;

  bool PhaseSpaceDecay(
         TClonesArray & pl, TLorentzVector & pd, 
	   const PDGCodeList & pdgv, int offset=0, bool reweight=false) const;

  mutable TGenPhaseSpace fPhaseSpaceGenerator; ///< a phase space generator
  mutable double         fWeight;              ///< weight for generated event

  // Configuration parameters
  // Note: additional configuration parameters common to all hadronizers
  // (Wcut,Rijk,...) are declared one layer down in the inheritance tree

  const DecayModelI * fDecayer;  ///< decay algorithm
  bool     fForceNeuGenLimit;    ///< force upper hadronic multiplicity to NeuGEN limit
//bool     fUseLegacyKNOSpline;  ///< use legacy spline instead of Levy
  bool     fUseIsotropic2BDecays;///< force isotropic, non-reweighted 2-body decays for consistency with neugen/daikon
  bool     fUseBaryonXfPt2Param; ///< Generate baryon xF,pT2 from experimental parameterization?
  bool     fReWeightDecays;      ///< Reweight phase space decays?
  bool     fForceDecays;         ///< force decays of unstable hadrons produced?
  bool     fForceMinMult;        ///< force minimum multiplicity if (at low W) generated less?
  bool     fGenerateWeighted;    ///< generate weighted events?
  double   fPhSpRwA;             ///< parameter for phase space decay reweighting
  double   fPpi0;                ///< {pi0 pi0  } production probability
  double   fPpic;                ///< {pi+ pi-  } production probability
  double   fPKc;                 ///< {K+  K-   } production probability
  double   fPK0;                 ///< {K0  K0bar} production probability
  double   fPpi0eta;             ///< {Pi0 eta} production probability
  double   fPeta;                ///< {eta eta} production probability
  double   fAvp;                 ///< offset in average charged hadron multiplicity = f(W) relation for vp
  double   fAvn;                 ///< offset in average charged hadron multiplicity = f(W) relation for vn
  double   fAvbp;                ///< offset in average charged hadron multiplicity = f(W) relation for vbp
  double   fAvbn;                ///< offset in average charged hadron multiplicity = f(W) relation for vbn
  double   fBvp;                 ///< slope  in average charged hadron multiplicity = f(W) relation for vp
  double   fBvn;                 ///< slope  in average charged hadron multiplicity = f(W) relation for vn
  double   fBvbp;                ///< slope  in average charged hadron multiplicity = f(W) relation for vbp
  double   fBvbn;                ///< slope  in average charged hadron multiplicity = f(W) relation for vbn
  double   fAhyperon;            ///< parameter controlling strange baryon production probability via associated production (P=a+b*lnW^2)
  double   fBhyperon;            ///< see above
  double   fCvp;                 ///< Levy function parameter for vp
  double   fCvn;                 ///< Levy function parameter for vn
  double   fCvbp;                ///< Levy function parameter for vbp
  double   fCvbn;                ///< Levy function parameter for vbn
  TF1 *    fBaryonXFpdf;         ///< baryon xF PDF
  TF1 *    fBaryonPT2pdf;        ///< baryon pT^2 PDF
//Spline * fKNO;                 ///< legacy KNO distribution (superseded by the Levy func)
};

}         // genie namespace

#endif    // _KNO_HADRONIZATION_H_

