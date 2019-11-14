//____________________________________________________________________________
/*!

\class    genie::KNOHadronization

\brief    A KNO-based hadronization model.

          Is a concrete implementation of the EventRecordVisitorI interface.

\author   The main authors of this model are:

          - Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
            University of Liverpool & STFC Rutherford Appleton Lab

          - Hugh Gallagher <gallag@minos.phy.tufts.edu>
            Tufts University

          - Tinjun Yang <tjyang@stanford.edu>
            Stanford University

          This is an improved version of the legacy neugen3 KNO-based model.
          Giles Barr, Geoff Pearce, and Hugh Gallagher were listed as authors
          of the original neugen3 model.

          Strange baryon production was implemented by Keith Hofmann and
          Hugh Gallagher (Tufts)

          Production of etas was added by Ji Liu (W&M)

          Changes required to implement the GENIE Boosted Dark Matter module
          were installed by Josh Berger (Univ. of Wisconsin)

\created  August 17, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _KNO_HADRONIZATION_H_
#define _KNO_HADRONIZATION_H_

#include <TGenPhaseSpace.h>

#include "Physics/Decay/Decayer.h"
#include "Framework/EventGen/EventRecordVisitorI.h"


class TF1;

namespace genie {

class Decayer;
//class Spline;

class KNOHadronization : public EventRecordVisitorI {

public:

  KNOHadronization();
  KNOHadronization(string config);
  virtual ~KNOHadronization();

  // Implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event) const;

  // Overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  virtual void Configure(const Registry & config);
  virtual void Configure(string config);

  friend class KNOTunedQPMDISPXSec ;

private:

  void           LoadConfig            (void);
  void           Initialize            (void)                                        const;
  TClonesArray * Hadronize             (const Interaction* )                         const;
  double         Weight                (void)                                        const;
  PDGCodeList *  SelectParticles       (const Interaction*)                          const;
  TH1D *         MultiplicityProb      (const Interaction*, Option_t* opt = "")      const;
  bool           AssertValidity        (const Interaction * i)                       const;
  PDGCodeList *  GenerateHadronCodes   (int mult, int maxQ, double W)                const;
  int            GenerateBaryonPdgCode (int mult, int maxQ, double W)                const;
  int            HadronShowerCharge    (const Interaction * )                        const;
  double         KNO                   (int nu, int nuc, double z)                   const;
  double         AverageChMult         (int nu, int nuc, double W)                   const;
  void           HandleDecays          (TClonesArray * particle_list)                const;
  double         ReWeightPt2           (const PDGCodeList & pdgcv)                   const;
  double         MaxMult               (const Interaction * i)                       const;
  TH1D *         CreateMultProbHist    (double maxmult)                              const;
  void           ApplyRijk             (const Interaction * i, bool norm, TH1D * mp) const;
  double         Wmin                  (void)                                        const;

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

  // nuegen parameters
  double   fWcut;      ///< Rijk applied for W<Wcut (see DIS/RES join scheme)
  double   fRvpCCm2;   ///< Rijk: vp,  CC, multiplicity = 2
  double   fRvpCCm3;   ///< Rijk: vp,  CC, multiplicity = 3
  double   fRvpNCm2;   ///< Rijk: vp,  NC, multiplicity = 2
  double   fRvpNCm3;   ///< Rijk: vp,  NC, multiplicity = 3
  double   fRvnCCm2;   ///< Rijk: vn,  CC, multiplicity = 2
  double   fRvnCCm3;   ///< Rijk: vn,  CC, multiplicity = 3
  double   fRvnNCm2;   ///< Rijk: vn,  NC, multiplicity = 2
  double   fRvnNCm3;   ///< Rijk: vn,  NC, multiplicity = 3
  double   fRvbpCCm2;  ///< Rijk: vbp, CC, multiplicity = 2
  double   fRvbpCCm3;  ///< Rijk: vbp, CC, multiplicity = 3
  double   fRvbpNCm2;  ///< Rijk: vbp, NC, multiplicity = 2
  double   fRvbpNCm3;  ///< Rijk: vbp, NC, multiplicity = 3
  double   fRvbnCCm2;  ///< Rijk: vbn, CC, multiplicity = 2
  double   fRvbnCCm3;  ///< Rijk: vbn, CC, multiplicity = 3
  double   fRvbnNCm2;  ///< Rijk: vbn, NC, multiplicity = 2
  double   fRvbnNCm3;  ///< Rijk: vbn, NC, multiplicity = 3

};

}         // genie namespace

#endif    // _KNO_HADRONIZATION_H_
