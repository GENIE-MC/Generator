//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightFGM

\brief    Reweighting the Fermi Gas nuclear model

\author   Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Apr 26, 2010

\cpright  Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_FGM_H_
#define _G_REWEIGHT_FGM_H_

#include <map>

#include "ReWeight/GReWeightI.h"

using std::map;

using namespace genie::rew;
using namespace genie;

class TH1D;
class TNtupleD;
class TFile;

namespace genie {

class NuclearModelI;

namespace rew   {

 class GReWeightFGM : public GReWeightI 
 {
 public:
   GReWeightFGM();
  ~GReWeightFGM();

   // implement the GReWeightI interface
   bool   IsHandled      (GSyst_t syst);
   void   SetSystematic  (GSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     (const EventRecord & event);
   double CalcChisq      (void);

 private:

   void Init(void);

   double RewCCQEPauliSupViaKF   (const EventRecord & event);
   double RewCCQEMomDistroFGtoSF (const EventRecord & event);

   double fKFTwkDial;
   double fMomDistroTwkDial;

   const NuclearModelI * fFG;
   const NuclearModelI * fSF;

   map<int, TH1D *> fMapFGn;
   map<int, TH1D *> fMapFGp;
   map<int, TH1D *> fMapSFn;
   map<int, TH1D *> fMapSFp;

   TFile *    fTestFile;
   TNtupleD * fTestNtp;
 };

} // rew
} // genie

#endif

