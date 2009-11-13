//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightNuXSec

\brief    Reweighting GENIE neutrino cross sections

\author   Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Aug 1, 2009

\cpright  Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_NEUTRINO_CROSS_SECTION_H_
#define _G_REWEIGHT_NEUTRINO_CROSS_SECTION_H_

#include <map>

#include "Algorithm/AlgFactory.h"
#include "Algorithm/AlgConfigPool.h"
#include "Registry/Registry.h"
#include "Interaction/Interaction.h"
#include "ReWeight/GReWeightNuXSecParams.h"
#include "ReWeight/GReWeightNuXSecHelper.h"
#include "ReWeight/GReWeightI.h"

using std::map;

namespace genie {
namespace rew   {

 class GReWeightNuXSec : public GReWeightI 
 {
 public:
   GReWeightNuXSec();
  ~GReWeightNuXSec();

   // implement the GReWeightI interface
   bool   IsHandled      (GSyst_t syst);
   void   SetSystematic  (GSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     (const EventRecord & event);
   double CalcChisq      (void);
      
 private:

   void   Init (void);

   double WghtStandardXSec  (const EventRecord & event);
   double WghtNonResBkgXSec (const EventRecord & event);

   GReWeightNuXSecHelper fXSecRwHelper;      ///<
   GReWeightNuXSecParams fXSecRwParams;      ///<
 };

} // rew   namespace
} // genie namespace

#endif

