//____________________________________________________________________________ 
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author:  Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 11, 2009 - CA
   Was first added in v2.5.1. Adapted from the T2K-specific version of the
   GENIE reweighting tool.

*/
//____________________________________________________________________________

#include <TMath.h>
//#include <Math/AllIntegrationTypes.h>
//#include <Math/IntegratorMultiDim.h>

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
#ifdef __GENIE_GSL_ENABLED__
/*
  ROOT::Math::IBaseFunctionMultiDim * func = 
    new genie::rew::margin::INukeFateIntgrd(fixed_dial, fixed_fate, event);
  ROOT::Math::IntegrationMultiDim::Type ig_type = 
    genie::utils::gsl::IntegrationNDimTypeFromString("vegas");

  ROOT::Math::IntegratorMultiDim ig(ig_type);
  ig.SetFunction(*func);

  double min[3] = { -1, -1, -1 };
  double max[3] = { +1, +1, +1 };

  double wght = ig.Integral(min, max);
  return wght; 
*/
#endif

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

  double min = -1;
  double max = +1;

  systematics.SetCurValue (fixed_fate, fixed_dial);

  double w_wght_sum = 0;
  int itry=0;

  while(itry < n) {
    itry++;

    if(is_pion_fate) {
      for(int i=0; i<5; i++) {
        genie::rew::GSyst_t syst = genie::rew::GSyst::NextPionFateSystematic(i);
        if(syst == fixed_fate  ) continue;
        if(syst == cushion_term) continue;
        double dial = min + (max-min) * rnd->RndGen().Rndm();
        systematics.SetCurValue(syst, dial);
      }
    } 
    else
    if(is_nucleon_fate) {
      for(int i=0; i<5; i++) {
        genie::rew::GSyst_t syst = genie::rew::GSyst::NextNuclFateSystematic(i);
        if(syst == fixed_fate  ) continue;
        if(syst == cushion_term) continue;
        double dial = min + (max-min) * rnd->RndGen().Rndm();
        systematics.SetCurValue(syst, dial);
      }
    }

    rw.Reconfigure();

    double wght = rw.CalcWeight(*event);

    double prob = 1.;

    if(is_pion_fate) {
      for(int i=0; i<5; i++) {
        genie::rew::GSyst_t syst = genie::rew::GSyst::NextPionFateSystematic(i);
        prob *= (TMath::Gaus( systematics.CurValue(syst), 0., 1., true ));
      }
    }
    else
    if(is_nucleon_fate) {
      for(int i=0; i<5; i++) {
        genie::rew::GSyst_t syst = genie::rew::GSyst::NextNuclFateSystematic(i);
        prob *= (TMath::Gaus( systematics.CurValue(syst), 0., 1., true ));
      }
    }

    w_wght_sum += (wght * prob);
  }//n

  double w_wght_avg = w_wght_sum/n;

  return w_wght_avg;
}
//____________________________________________________________________________
/*
genie::rew::margin::INukeFateIntgrd::INukeFateIntgrd(
  double fixed_twk_dial, genie::rew::GSyst_t fixed_syst, 
  const genie::EventRecord * event) :
ROOT::Math::IBaseFunctionMultiDim()
{
  fFixedTwKDial  = fixed_twk_dial;
  fFixedFateSyst = fixed_syst;
  fEvent         = event;

  fRW = new genie::rew::GReWeight;
  genie::rew::GSystSet & systematics = fRW->Systematics();

  if(genie::rew::GSyst::IsINukePionFateSystematic(fFixedFateSyst)) {
     systematics.Include ( kSystINuke_CExTwk_pi    );
     systematics.Include ( kSystINuke_InelTwk_pi   );
     systematics.Include ( kSystINuke_AbsTwk_pi    );
     systematics.Include ( kSystINuke_PiProdTwk_pi );
  } 
  else
  if(genie::rew::GSyst::IsINukeNuclFateSystematic(fFixedFateSyst)) {
     systematics.Include ( kSystINuke_CExTwk_N    );
     systematics.Include ( kSystINuke_InelTwk_N   );
     systematics.Include ( kSystINuke_AbsTwk_N    );
     systematics.Include ( kSystINuke_PiProdTwk_N );
  }
}
unsigned int genie::rew::margin::INukeFateIntgrd::NDim(void) const
{
  return 3; // 5 fates - 1 cushion term - 1 fixed term at each integration 
}
double genie::rew::margin::INukeFateIntgrd::DoEval(const double * xin) const
{
  genie::rew::GSystSet & systematics = fRW->Systematics();

  systematics.SetCurValue (fFixedFateSyst, fFixedTwKDial);

  if(GSyst::IsINukePionFateSystematic(fFixedFateSyst)) {
     systematics.SetCurValue ( genie::rew::kSystINuke_CExTwk_pi    , xin[0]);
     systematics.SetCurValue ( genie::rew::kSystINuke_InelTwk_pi   , xin[1]);
     systematics.SetCurValue ( genie::rew::kSystINuke_PiProdTwk_pi , xin[2]);
  }

  fRW->Reconfigure();

  const genie::EventRecord & event = *fEvent;

  double wght = fRW->CalcWeight(event);
  double prob = 1.;
  {
    prob *= (TMath::Gaus( systematics.CurValue(genie::rew::kSystINuke_CExTwk_pi   ), 0., 1., true ));
    prob *= (TMath::Gaus( systematics.CurValue(genie::rew::kSystINuke_InelTwk_pi  ), 0., 1., true ));
    prob *= (TMath::Gaus( systematics.CurValue(genie::rew::kSystINuke_PiProdTwk_pi), 0., 1., true ));
    prob *= (TMath::Gaus( systematics.CurValue(genie::rew::kSystINuke_ElTwk_pi    ), 0., 1., true ));
  }
  double w_wght = wght * prob;
  return w_wght;
}
ROOT::Math::IBaseFunctionMultiDim * 
   genie::rew::margin::INukeFateIntgrd::Clone() const
{
  return 
    new genie::rew::margin::INukeFateIntgrd(
          fFixedTwKDial, fFixedFateSyst, fEvent);
}
//____________________________________________________________________________
*/
