//____________________________________________________________________________ 
/*
 Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author:  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Oct 20, 2009 - CA
   Was first added in v2.5.1. Included method to marginalize N-1 "hadron fate" 
   nuisance parameter. The code can be used to estimate the effect of tweaking 
   one of the absorption, charge exchange, inelastic or pion production fate 
   fraction *irrespective* of the other fate fractions (within the unitarity constraint).

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Conventions/Units.h"
#include "GHEP/GHepParticle.h"
#include "HadronTransport/INukeHadroData.h"
#include "HadronTransport/INukeHadroFates.h"
#include "HadronTransport/INukeUtils.h"
#include "Messenger/Messenger.h"
#include "Numerical/Spline.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGCodes.h"
#include "ReWeight/GReWeightMarginalize.h"
#include "Utils/GSLUtils.h"

using namespace genie;
using namespace genie::rew;

//____________________________________________________________________________
double genie::rew::margin::MarginalizeFates(
  double fixed_dial, genie::rew::GSyst_t fixed_fate, 
  const genie::EventRecord * event, int n)
{
  bool is_pion_fate    = 
         genie::rew::GSyst::IsINukePionFateSystematic(fixed_fate); 
  bool is_nucleon_fate = 
         genie::rew::GSyst::IsINukeNuclFateSystematic(fixed_fate);

  // Peek fwd to avoid wasting time
  // Return 1 if margninalizing pion rescattering fates but there are
  // no primary pions in this event (similarly for nucleons)
  TIter event_iter(event);
  GHepParticle * p = 0;
  int npi  = 0;
  int nnuc = 0;
  while((p=dynamic_cast<GHepParticle *>(event_iter.Next())))
  {
    if(p->Status() != kIStHadronInTheNucleus) continue;
    if(pdg::IsPion   (p->Pdg())) { npi++;  }
    if(pdg::IsNucleon(p->Pdg())) { nnuc++; }
  }

  bool skip = (is_pion_fate    && npi  == 0) || 
              (is_nucleon_fate && nnuc == 0);
  if(skip) {
    return 1.;
  }

  genie::rew::GReWeight rw;
  genie::rew::GSystSet & systematics = rw.Systematics();

  genie::rew::GSyst_t cushion_term = kSystNull;

  if(is_pion_fate) {
     systematics.Include ( genie::rew::kSystINuke_CExTwk_pi    );
     systematics.Include ( genie::rew::kSystINuke_InelTwk_pi   );
     systematics.Include ( genie::rew::kSystINuke_AbsTwk_pi    );
     systematics.Include ( genie::rew::kSystINuke_PiProdTwk_pi );
     cushion_term = genie::rew::kSystINuke_ElTwk_pi;
  } 
  else
  if(is_nucleon_fate) {
     systematics.Include ( genie::rew::kSystINuke_CExTwk_N    );
     systematics.Include ( genie::rew::kSystINuke_InelTwk_N   );
     systematics.Include ( genie::rew::kSystINuke_AbsTwk_N    );
     systematics.Include ( genie::rew::kSystINuke_PiProdTwk_N );
     cushion_term = genie::rew::kSystINuke_ElTwk_N;
  }

  RandomGen * rnd = RandomGen::Instance();

  systematics.SetCurValue (fixed_fate, fixed_dial);

  double tot_wgt = 0;
  int itry=0;

  while(itry < n) {
    itry++;
    LOG("ReW", pNOTICE) 
        << "Marginalization progress : " << itry << "/" << n;
    if(is_pion_fate) {
      for(int i=0; i<5; i++) {
        genie::rew::GSyst_t syst = genie::rew::GSyst::NextPionFateSystematic(i);
        if(syst == fixed_fate  ) continue;
        if(syst == cushion_term) continue;
        systematics.SetCurValue(syst, rnd->RndGen().Gaus(0,1));
      }
    } 
    else
    if(is_nucleon_fate) {
      for(int i=0; i<5; i++) {
        genie::rew::GSyst_t syst = genie::rew::GSyst::NextNuclFateSystematic(i);
        if(syst == fixed_fate  ) continue;
        if(syst == cushion_term) continue;
        systematics.SetCurValue(syst, rnd->RndGen().Gaus(0,1));
      }
    }

    // Reconfigure GENIE to use the new settings and calculate a weight
    rw.Reconfigure();
    double w = rw.CalcWeight(*event);
    tot_wgt += w;

    LOG("ReW", pNOTICE) 
        << "Current weight = " << w << " (sum = " << tot_wgt << ")";

  }//n

  double avg_wgt = tot_wgt/n;
  return avg_wgt;
}
//____________________________________________________________________________
