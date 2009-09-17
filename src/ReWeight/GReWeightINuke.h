//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightINuke

\brief    Reweighting GENIE INTRANUKE/hA hadron transport model.

\author   Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Sep 10, 2009

\cpright  Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
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

class TLorentzVector;
class TNtuple;

namespace genie {
namespace rew   {

 class GReWeightINuke : public GReWeightI 
 {
 public:
   GReWeightINuke();
  ~GReWeightINuke();

   // implement the GReWeightI interface
   void   SetSystematic  (GSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     (const EventRecord & event);
   double CalcChisq      (void);

 private:

   GReWeightINukeParams fINukeRwParams;
   TNtuple *            fTestNtp;
 };

} // rew
} // genie

#endif

