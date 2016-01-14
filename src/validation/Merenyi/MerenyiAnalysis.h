//____________________________________________________________________________
/*!

\class   MerenyiAnalysis

\brief   Merenyi analysis base class

\author  Pauli Kehayias (Tufts Univ)

\created Jan 13, 2009

\cpright  Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _MERENYI_ANALYSIS_H_
#define _MERENYI_ANALYSIS_H_

#include <fstream>

#include "validation/Merenyi/MerenyiNuEvent.h"

using namespace std;

namespace genie       {
namespace vld_merenyi {

class MerenyiAnalysis {

public:
  MerenyiAnalysis ();
  virtual ~MerenyiAnalysis () { }

  void setFile(string file) { myFile = file.c_str(); }
  void readData();

  // display populations
  //
  void dispPops      (void);
  void dispNuEHist   (void);
  void dispMuAngHist (void);
  void dispPiMomHist (void);

  // write populations to file
  //
  void printPops      (void);
  void printNuEHist   (void);
  void printMuAngHist (void);
  void printPiMomHist (void);

  // returns a vector of the populations
  //
  vector <double> getPops      (void);
  vector <double> getNuEHist   (void);
  vector <double> getMuAngHist (void);
  vector <double> getPiMomHist (void);

  // returns i'th element of the neutrino energy histogram
  //
  double getNuEi   (int i);                                            
  double getMuAngi (int i);
  double getPiMomi (int i);

  int getNumEvents   (void) { return numEvents;               }
  int getNumCCEvents (void) { return (numEvents - numTossed); }
  int getNumChPions  (void) { return numPions;                } ///< return the number of charged pions

protected:

  ifstream in;
  int      numEvents;
  int      numTossed;

  void findPiMom(TParticle pion);
  virtual void checkCandP(MerenyiNuEvent& /*myEvent*/) {}  //check for candidate proton
  virtual void classifyHadrons(MerenyiNuEvent& /*myEvent*/) {}

private:

  string       myFile;
  vector <int> pEventTypes;    ///< 6 types of events with proton targets
  vector <int> nEventTypes;    ///< 8 types of events with neutron targets
  int          ill1EventTypes; ///< charge illegal event, type 1 (last line of Merenyi et al. Table 1)
  int          ill2EventTypes; ///< charge illegal event, type 2 (not on the table)

  void clearAll (void);

  void readNuEnergy (MerenyiNuEvent& myEvent);
  void readMuMom    (MerenyiNuEvent& myEvent);
  void processEvent (MerenyiNuEvent  myEvent);
  void getReaction  (MerenyiNuEvent  myEvent);
  void pTarget      (MerenyiNuEvent  myEvent);  ///< classifies proton target events
  void nTarget      (MerenyiNuEvent  myEvent);  ///< classifies neutron target events
  void illTarget1   (MerenyiNuEvent  myEvent);  ///< classifies first type of charge illegal event
  void illTarget2   (void);              ///< classifies second type of charge illegal event

  const double binSizeE;        ///< for neutrino energies
  const double binSizeAng;      ///< for muon angles
  const double binSizeMom;      ///< for pi+- momenta
  vector <int> nuECounts;
  vector <int> muAngCounts;
  vector <int> piMomCounts;
  int numPions;
};


} // vld_merenyi
} // genie

#endif


