//____________________________________________________________________________
/*!

\class    genie::GHepSummaryBuilder

\brief    An object that knows how to look at the GHEP event record and extract
          summary information.
          The summary builder is used for building the plain (PR) ntuple which
          can also be generated for GHEP records converted from other generator
          outputs.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  July 10, 2005

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _GHEP_SUMMARY_BUILDER_H_
#define _GHEP_SUMMARY_BUILDER_H_

class TLorentzVector;

#include "Interaction/Interaction.h"
#include "Interaction/InteractionType.h"
#include "Interaction/ScatteringType.h"

namespace genie {

class GHepRecord;

class GHepSummaryBuilder {

public :

  GHepSummaryBuilder();
  ~GHepSummaryBuilder();

  //-- analyze the event record to extract summary information
  void AnalyzeEventRecord(const GHepRecord & evrec);

  //-- pdgc codes of various particles participating in the event
  int ProbePdgC     (void) const { return fProbePdgC; }
  int FslPdgC       (void) const { return fFslPdgC;   }
  int TgtPdgC       (void) const { return fTgtPdgC;   }
  int HitNuclPdgC   (void) const { return fNuclPdgC;  }
  int HitQuarkPdgC  (void) const { return fIQrkPdgC;  }
  int OutQuarkPdgC  (void) const { return fFQrkPdgC;  }
  int ResPdgC       (void) const { return fResPdgC;   }
  int CharmHadPdgC  (void) const { return fChHadPdgC; }

  //-- nuclear target Z,A,N
  int NuclTgtZ      (void) const { return fTgtZ;      }
  int NuclTgtA      (void) const { return fTgtA;      }
  int NuclTgtN      (void) const { return fTgtN;      }

  //-- kinematical variables
  double KineNu     (void) const { return fNu;        }
  double KineX      (void) const { return fX;         }
  double KineY      (void) const { return fY;         }
  double FragmZ     (void) const { return fZ;         }
  double KineQ2     (void) const { return fQ2;        }
  double KineW      (void) const { return fW;         }
  double EmFrac     (void) const { return fEmFrac;    }

  //-- scattering and interaction type
  ScatteringType_t  ScatType (void) const { return fScatType;  }
  InteractionType_t ProcType (void) const { return fProcType;  }

  //-- number of final state particles
  unsigned int NProton  (void) const { return fNProton;  }
  unsigned int NNeutron (void) const { return fNNeutron; }
  unsigned int NPi0     (void) const { return fNPi0;     }
  unsigned int NPiPlus  (void) const { return fNPiPlus;  }
  unsigned int NPiMinus (void) const { return fNPiMinus; }
  unsigned int NK0      (void) const { return fNK0;      }
  unsigned int NKPlus   (void) const { return fNKPlus;   }
  unsigned int NKMinus  (void) const { return fNKMinus;  }

  //-- vertices and particle 4-momenta
  const TLorentzVector & Vtx       (void) const { return *fVtx;      }
  const TLorentzVector & q4p       (void) const { return *fq4p;      }
  const TLorentzVector & Probe4P   (void) const { return *fProbe4P;  }
  const TLorentzVector & HitNucl4P (void) const { return *fNucl4P;   }
  const TLorentzVector & Fsl4P     (void) const { return *fFsl4P;    }
  const TLorentzVector & HadShw4P  (void) const { return *fHadShw4P; }

private:

  void Init();    ///< initialize data members
  void CleanUp(); ///< cleans up results from analyzing the previous GHEP record

  int                fProbePdgC;///< Probe (v,e...) PDG code
  int                fFslPdgC;  ///< Final state primary lepton PDG code
  int                fTgtPdgC;  ///< Nuclear target PDG code
  int                fNuclPdgC; ///< Struck nucleon PDG code
  int                fIQrkPdgC; ///< Struck quark PDG code
  int                fFQrkPdgC; ///< Outgoing quark PDG code
  int                fResPdgC;  ///< Resonance PDG code
  int                fChHadPdgC;///< Charm hadron PDG code
  int                fTgtZ;     ///< Nuclear target: Z
  int                fTgtA;     ///< Nuclear target: A
  int                fTgtN;     ///< Nuclear target: N
  ScatteringType_t   fScatType; ///< Scattering type (QEL, DIS, RES, COH, ...)
  InteractionType_t  fProcType; ///< Process type (E/M, CC, NC, ...)
  double             fNu;       ///< v
  double             fX;        ///< Bjorken x
  double             fY;        ///< Inelasticity
  double             fZ;        ///< Fragmentation z
  double             fQ2;       ///< Momentum transfer (>0 = -q2)
  double             fW;        ///< Hadronic invariant mass
  double             fEmFrac;   ///< E/M fraction
  TLorentzVector *   fVtx;      ///< Interaction vertex (x,y,z,t)
  TLorentzVector *   fProbe4P;  ///< Probe 4-P (px,py,pz,E) in LAB
  TLorentzVector *   fNucl4P;   ///< Struck nucleon 4-P (px,py,py,E) in LAB
  TLorentzVector *   fFsl4P;    ///< Final state primary lepton 4-P (px,py,py,E) in LAB
  TLorentzVector *   f2Fsl4P;   ///<
  TLorentzVector *   fHadShw4P; ///< Hadronic system 4-P (px,py,py,E) in LAB
  TLorentzVector *   fq4p;      ///<
  unsigned int       fNProton;  ///< Number of generated protons
  unsigned int       fNNeutron; ///< Number of generated neutrons
  unsigned int       fNPi0;     ///< Number of generated pi^0
  unsigned int       fNPiPlus;  ///< Number of generated pi^+
  unsigned int       fNPiMinus; ///< Number of generated pi^-
  unsigned int       fNK0;      ///< Number of generated K^0
  unsigned int       fNKPlus;   ///< Number of generated K^+
  unsigned int       fNKMinus;  ///< Number of generated K^-
};

}      // genie namespace

#endif // _GHEP_SUMMARY_BUILDER_H_
