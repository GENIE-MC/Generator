//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightINuke

\brief    Reweighting GENIE INTRANUKE/hA hadron transport model.

	  The reweighting code considers two sets of physics changes:
          * Change in the hadron mean free path \lambda, i.e. change in  
            the total rescattering probability P_{rescat}.
          * Changes in probabilty for rescattering mode X, given a fixed
            total rescattering probability P(X | \lambda).
            X = {elastic, inelastic, charge exchange, pion production,
                 absorption + multi-nucleon knockout}.        

          Physics changes are considered separately for pions and nucleons.
          Unitarity is explicitly conserved.

\author   Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Sep 10, 2009

\cpright  Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_INUKE_H_
#define _G_REWEIGHT_INUKE_H_

//#define _G_REWEIGHT_INUKE_DEBUG_NTP_

#include "ReWeight/GReWeightI.h"
#include "ReWeight/GReWeightINukeParams.h"

using namespace genie::rew;
using namespace genie;

class TFile;
class TNtuple;
class TLorentzVector;

namespace genie {
namespace rew   {

 class GReWeightINuke : public GReWeightI 
 {
 public:
   GReWeightINuke();
  ~GReWeightINuke();

   // implement the GReWeightI interface
   bool   IsHandled      (GSyst_t syst);
   void   SetSystematic  (GSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     (const EventRecord & event);
   double CalcChisq      (void);

 private:

   GReWeightINukeParams fINukeRwParams;
   TFile *              fTestFile;
   TNtuple *            fTestNtp;
 };

} // rew
} // genie

#endif

