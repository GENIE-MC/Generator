//____________________________________________________________________________
/*!

\class   genie::NtpGHepAnalyzer

\brief   An object that knows how to look at the GHEP event record and extract
         summary information for the generated event so as to fill the output
         ntuple's NtpMCSummary.

         This information also exists at the Interaction summary attached to
         the event record, but the attached Interaction is NULL for the GHEP
         records not generated within the GENIE framework but extracted &
         translated from external generators.
         Being able to analyze the GHEP record allows us to properly generate
         the GENIE ntuple for other generators too facilitating cross-generator
         comparisons.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created July 10, 2005

*/
//____________________________________________________________________________

#ifndef _NTP_GHEP_ANALYZER_H_
#define _NTP_GHEP_ANALYZER_H_

class TLorentzVector;

#include "Interaction/InteractionType.h"
#include "Interaction/ScatteringType.h"

namespace genie {

class EventRecord;

class NtpGHepAnalyzer {

public :

  NtpGHepAnalyzer();
  ~NtpGHepAnalyzer();

  void AnalyzeEventRecord(const EventRecord & evrec);

  int   ProbePdgCode     (void) const { return fProbePdgC; }
  int   FslPdgCode       (void) const { return fFslPdgC;   }
  int   TgtPdgCode       (void) const { return fTgtPdgC;   }
  int   HitNuclPdgCode   (void) const { return fNuclPdgC;  }
  int   HitQuarkPdgCode  (void) const { return fIQrkPdgC;  }
  int   OutgQuarkPdgCode (void) const { return fFQrkPdgC;  }
  int   ResPdgCode       (void) const { return fResPdgC;   }
  int   CharmHadPdgCode  (void) const { return fChHadPdgC; }
  int   NuclTgtZ         (void) const { return fTgtZ;      }
  int   NuclTgtA         (void) const { return fTgtA;      }
  int   NuclTgtN         (void) const { return fTgtN;      }
  int   ScatType         (void) const { return fScatType;  }
  int   ProcType         (void) const { return fProcType;  }
  float KineNu           (void) const { return fNu;        }
  float KineX            (void) const { return fX;         }
  float KineY            (void) const { return fY;         }
  float FragmZ           (void) const { return fZ;         }
  float KineQ2           (void) const { return fQ2;        }
  float KineW            (void) const { return fW;         }
  float XSec             (void) const { return fXSec;      }
  float dXSec            (void) const { return fdXSec;     }
  float EmFrac           (void) const { return fEmFrac;    }

  unsigned int NumProton  (void) const { return fNProton;  }
  unsigned int NumNeutron (void) const { return fNNeutron; }
  unsigned int NumPi0     (void) const { return fNPi0;     }
  unsigned int NumPiPlus  (void) const { return fNPiPlus;  }
  unsigned int NumPiMinus (void) const { return fNPiMinus; }
  unsigned int NumK0      (void) const { return fNK0;      }
  unsigned int NumKPlus   (void) const { return fNKPlus;   }
  unsigned int NumKMinus  (void) const { return fNKMinus;  }

  const TLorentzVector & Vtx       (void) const { return *fVtx;      }
  const TLorentzVector & q4p       (void) const { return *fq4p;      }
  const TLorentzVector & Probe4P   (void) const { return *fProbe4P;  }
  const TLorentzVector & HitNucl4P (void) const { return *fNucl4P;   }
  const TLorentzVector & Fsl4P     (void) const { return *fFsl4P;    }
  const TLorentzVector & HadShw4P  (void) const { return *fHadShw4P; }

private:

  void Init();    ///< initialize data members
  void CleanUp(); ///< cleans up results from analyzing the previous GHEP record

  unsigned int NEntries(int pdgc) const;

  const EventRecord * fEventRec; ///< latest attached EventRecord

  int               fProbePdgC;///< Probe (v,e...) PDG code
  int               fFslPdgC;  ///< Final state primary lepton PDG code
  int               fTgtPdgC;  ///< Nuclear target PDG code
  int               fNuclPdgC; ///< Struck nucleon PDG code
  int               fIQrkPdgC; ///< Struck quark PDG code
  int               fFQrkPdgC; ///< Outgoing quark PDG code
  int               fResPdgC;  ///< Resonance PDG code
  int               fChHadPdgC;///< Charm hadron PDG code
  int               fTgtZ;     ///< Nuclear target: Z
  int               fTgtA;     ///< Nuclear target: A
  int               fTgtN;     ///< Nuclear target: N
  ScatteringType_t  fScatType; ///< Scattering type (QEL, DIS, RES, COH, ...)
  InteractionType_t fProcType; ///< Process type (E/M, CC, NC, ...)
  float             fNu;       ///< v
  float             fX;        ///< Bjorken x
  float             fY;        ///< Inelasticity
  float             fZ;        ///< Fragmentation z
  float             fQ2;       ///< Momentum transfer (>0 = -q2)
  float             fW;        ///< Hadronic invariant mass
  float             fXSec;     ///< Cross section for this interaction sig=f(E)
  float             fdXSec;    ///< Cross section for the selected kinematics
  float             fEmFrac;   ///< E/M fraction
  TLorentzVector *  fVtx;      ///< Interaction vertex (x,y,z,t)
  TLorentzVector *  fProbe4P;  ///< Probe 4-P (px,py,pz,E) in LAB
  TLorentzVector *  fNucl4P;   ///< Struck nucleon 4-P (px,py,py,E) in LAB
  TLorentzVector *  fFsl4P;    ///< Final state primary lepton 4-P (px,py,py,E) in LAB
  TLorentzVector *  f2Fsl4P;   ///<
  TLorentzVector *  fHadShw4P; ///< Hadronic system 4-P (px,py,py,E) in LAB
  TLorentzVector *  fq4p;      ///<
  unsigned int      fNProton;  ///< Number of generated protons
  unsigned int      fNNeutron; ///< Number of generated neutrons
  unsigned int      fNPi0;     ///< Number of generated pi^0
  unsigned int      fNPiPlus;  ///< Number of generated pi^+
  unsigned int      fNPiMinus; ///< Number of generated pi^-
  unsigned int      fNK0;      ///< Number of generated K^0
  unsigned int      fNKPlus;   ///< Number of generated K^+
  unsigned int      fNKMinus;  ///< Number of generated K^-
};

}      // genie namespace

#endif // _NTP_GHEP_ANALYZER_H_
