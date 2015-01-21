//____________________________________________________________________________
/*!

\class   MerenyiNuEvent

\brief   Neutrino event description for the Merenyi test

\author  Pauli Kehayias (Tufts Univ)

\created Jan 13, 2009

\cpright  Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _MERENYI_NEUTRINO_EVENT_H_
#define _MERENYI_NEUTRINO_EVENT_H_

#include <vector>

#include <TParticle.h>

using namespace std;

namespace genie       {
namespace vld_merenyi { 
   
class MerenyiNuEvent
{
public:
  MerenyiNuEvent();
  MerenyiNuEvent(int numHadrons);
  MerenyiNuEvent(int numHadrons, int interactionType, double nuEnergy, double muPx, double muPy, double muPz, double muEnergy);

  void zeroAll();
    
  vector<TParticle> getHadList(void) {return hadList;}  ///<r eturns hadList
  void   addParticle(double particleID, double px, double py, double pz, double E); ///< add a particle to hadList
  double getLepAngCos(); ///< calculate lepton angle cosine
  
  int    getNumHadrons (void) { return numHad;      }
  int    getIntType    (void) { return intType;     }
  double getNuE        (void) { return nuE;         }
  double getLepPx      (void) { return lepton.Px(); }
  double getLepPy      (void) { return lepton.Py(); }
  double getLepPz      (void) { return lepton.Pz(); }
  int    getNumGam     (void) { return numGam;      }
  int    getNumPi0     (void) { return numPi0;      }
  int    getNumPiP     (void) { return numPiP;      }
  int    getNumPiM     (void) { return numPiM;      }
  int    getNumN       (void) { return numN;        }
  int    getNumP       (void) { return numP;        }
  int    getNumStr0    (void) { return numStr0;     }
  int    getNumStrP    (void) { return numStrP;     }
  int    getNumStrM    (void) { return numStrM;     }
  
  void   setNumHadrons (int numHadrons)                            { numHadrons = numHad;               }
  void   setNuE        (double nuEnergy)                           { nuE = nuEnergy;                    }
  void   setIntType    (int type)                                  { intType = type;                    }
  void   setLepP       (double px, double py, double pz, double E) { lepton.SetMomentum(px, py, pz, E); }

  void   addGam        (void) { numGam++;  }
  void   addPi0        (void) { numPi0++;  }
  void   addPiP        (void) { numPiP++;  }
  void   addPiM        (void) { numPiM++;  }
  void   addN          (void) { numN++;    }
  void   addP          (void) { numP++;    }
  void   addStr0       (void) { numStr0++; }
  void   addStrP       (void) { numStrP++; }
  void   addStrM       (void) { numStrM++; }

protected:
  vector<TParticle> hadList;  ///< contains particle IDs, momenta, energy for shower hadrons
  int               numHad;   ///< number of hadrons in the shower
  TParticle         lepton;   ///< resulting lepton
  int               intType;  ///< interaction type (1 = true QE, 2 = resonance production, 3 = deep inelastic scattering)
  double            nuE;      ///< neutrino energy
  int               numGam;   ///< number of gammas
  int               numPi0;   ///< number of pi0
  int               numPiP;   ///< number of pi+
  int               numPiM;   ///< number of pi-
  int               numN;     ///< number of neutrons
  int               numP;     ///< number of protons
  int               numStr0;  ///< number of strange 0,
  int               numStrP;  ///< number of strange +
  int               numStrM;  ///< number of strange -
};

} // vld_merenyi
} // genie

#endif
