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

#include "Algorithm/AlgFactory.h"
#include "Algorithm/AlgConfigPool.h"
#include "Registry/Registry.h"
#include "Interaction/Interaction.h"
#include "ReWeight/GReWeightNuXSecHelper.h"
#include "ReWeight/GSystErrors.h"
#include "ReWeight/GReWeightI.h"

namespace genie {
namespace rew   {

 class GReWeightNuXSec : public GReWeightI 
 {
 public:
   GReWeightNuXSec();
  ~GReWeightNuXSec();

   void   SetSystematic  (GSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     (const EventRecord & event);
   double CalcChisq      (void);
      
 private:

   void Init (void);

   GSystSet *            fSystSet;           ///<
   AlgFactory *          fAlgFactory;        ///<
   AlgConfigPool *       fConfigPool;        ///<
   Registry *            fUserPhysicsConfig; ///<
   GReWeightNuXSecHelper fXSecRwHelper;      ///<

   bool   fTweaked;

   double fMaQEL;                       ///<
   double fMvQEL;                       ///<
   double fMaRES;                       ///<
   double fMvRES;                       ///<
   double fRvpCC1pi;                    ///<
   double fRvpCC2pi;                    ///<
   double fRvpNC1pi;                    ///<
   double fRvpNC2pi;                    ///<
   double fRvnCC1pi;                    ///<
   double fRvnCC2pi;                    ///<
   double fRvnNC1pi;                    ///<
   double fRvnNC2pi;                    ///<
   double fRvbarpCC1pi;                 ///<
   double fRvbarpCC2pi;                 ///<
   double fRvbarpNC1pi;                 ///<
   double fRvbarpNC2pi;                 ///<
   double fRvbarnCC1pi;                 ///<
   double fRvbarnCC2pi;                 ///<
   double fRvbarnNC1pi;                 ///<
   double fRvbarnNC2pi;                 ///<

   double fMaQEL_def;
   double fMvQEL_def;
   double fMaRES_def;
   double fMvRES_def;
   double fRvpCC1pi_def;    
   double fRvpCC2pi_def;   
   double fRvpNC1pi_def;   
   double fRvpNC2pi_def;   
   double fRvnCC1pi_def;   
   double fRvnCC2pi_def;   
   double fRvnNC1pi_def;   
   double fRvnNC2pi_def;   
   double fRvbarpCC1pi_def;
   double fRvbarpCC2pi_def;
   double fRvbarpNC1pi_def;
   double fRvbarpNC2pi_def;
   double fRvbarnCC1pi_def;
   double fRvbarnCC2pi_def;
   double fRvbarnNC1pi_def;
   double fRvbarnNC2pi_def;

   double fNuXSec_MaQELTwkDial;
   double fNuXSec_MvQELTwkDial;
   double fNuXSec_MaRESTwkDial;
   double fNuXSec_MvRESTwkDial;
   double fNuXSec_RvpCC1piTwkDial;    
   double fNuXSec_RvpCC2piTwkDial;   
   double fNuXSec_RvpNC1piTwkDial;   
   double fNuXSec_RvpNC2piTwkDial;   
   double fNuXSec_RvnCC1piTwkDial;   
   double fNuXSec_RvnCC2piTwkDial;   
   double fNuXSec_RvnNC1piTwkDial;   
   double fNuXSec_RvnNC2piTwkDial;   
   double fNuXSec_RvbarpCC1piTwkDial;
   double fNuXSec_RvbarpCC2piTwkDial;
   double fNuXSec_RvbarpNC1piTwkDial;
   double fNuXSec_RvbarpNC2piTwkDial;
   double fNuXSec_RvbarnCC1piTwkDial;
   double fNuXSec_RvbarnCC2piTwkDial;
   double fNuXSec_RvbarnNC1piTwkDial;
   double fNuXSec_RvbarnNC2piTwkDial;

   bool fMaQEL_included;
   bool fMvQEL_included;
   bool fMaRES_included;
   bool fMvRES_included;
   bool fRvpCC1pi_included;    
   bool fRvpCC2pi_included;   
   bool fRvpNC1pi_included;   
   bool fRvpNC2pi_included;   
   bool fRvnCC1pi_included;   
   bool fRvnCC2pi_included;   
   bool fRvnNC1pi_included;   
   bool fRvnNC2pi_included;   
   bool fRvbarpCC1pi_included;
   bool fRvbarpCC2pi_included;
   bool fRvbarpNC1pi_included;
   bool fRvbarpNC2pi_included;
   bool fRvbarnCC1pi_included;
   bool fRvbarnCC2pi_included;
   bool fRvbarnNC1pi_included;
   bool fRvbarnNC2pi_included;     
 };

} // rew   namespace
} // genie namespace

#endif

