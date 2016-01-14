//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightDISNuclMod

\brief    Reweighting the DIS nuclear modification model

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

#ifndef _G_REWEIGHT_DISNUCLMOD_H_
#define _G_REWEIGHT_DISNUCLMOD_H_

#include "ReWeight/GReWeightI.h"

using namespace genie::rew;
using namespace genie;

namespace genie {
namespace rew   {

 class GReWeightDISNuclMod : public GReWeightI 
 {
 public:
   GReWeightDISNuclMod();
  ~GReWeightDISNuclMod();

   // implement the GReWeightI interface
   bool   IsHandled      (GSyst_t syst);
   void   SetSystematic  (GSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     (const EventRecord & event);
   double CalcChisq      (void);

 private:
   void Init(void);

   double fNuclModTwkDial;
 };

} // rew
} // genie

#endif

