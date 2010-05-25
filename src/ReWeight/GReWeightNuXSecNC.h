//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightNuXSecNC

\brief    

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

\created  May 25, 2010

\cpright  Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_NU_XSEC_NC_H_
#define _G_REWEIGHT_NU_XSEC_NC_H_

#include "ReWeight/GReWeightI.h"

namespace genie {
namespace rew   {

 class GReWeightNuXSecNC : public GReWeightI 
 {
 public:
   GReWeightNuXSecNC();
  ~GReWeightNuXSecNC();

   // implement the GReWeightI interface
   bool   IsHandled      (GSyst_t syst);
   void   SetSystematic  (GSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     (const EventRecord & event);
   double CalcChisq      (void);

   // various config options
   void RewNue       (bool   tf)  { fRewNue     = tf; }
   void RewNuebar    (bool   tf)  { fRewNuebar  = tf; }
   void RewNumu      (bool   tf)  { fRewNumu    = tf; }
   void RewNumubar   (bool   tf)  { fRewNumubar = tf; }

 private:

   void   Init (void);

   bool   fRewNue;         ///< reweight nu_e?
   bool   fRewNuebar;      ///< reweight nu_e_bar?
   bool   fRewNumu;        ///< reweight nu_mu?
   bool   fRewNumubar;     ///< reweight nu_mu_bar?

   double fNCTwkDial;      ///<
 };

} // rew   namespace
} // genie namespace

#endif

