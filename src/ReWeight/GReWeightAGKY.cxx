//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 10, 2009 - CA
   Skeleton first included in v2.5.1.

*/
//____________________________________________________________________________

#include <TLorentzVector.h>

#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "ReWeight/GReWeightAGKY.h"
#include "ReWeight/GReWeightUtils.h"

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GReWeightAGKY::GReWeightAGKY() :
GReWeightI()
{

}
//_______________________________________________________________________________________
GReWeightAGKY::~GReWeightAGKY()
{

}
//_______________________________________________________________________________________
void GReWeightAGKY::SetSystematic(GSyst_t /*syst*/, double /*val*/)
{

}
//_______________________________________________________________________________________
void GReWeightAGKY::Reset(void)
{

}
//_______________________________________________________________________________________
void GReWeightAGKY::Reconfigure(void)
{

}
//_______________________________________________________________________________________
double GReWeightAGKY::CalcWeight(const EventRecord & event) 
{  
  // Skip events not handled by the AGKY hadronization model
  if(! utils::rew::HadronizedByAGKY(event)) return 1.0;

  // Skip high-W events hadronized by AGKY/PYTHIA
  if(utils::rew::HadronizedByAGKYPythia(event)) return 1.0;

  //
  // Reweight events hadronized by AGKY/KNO
  //

  double event_weight = 1.0;

  // Get the hadronic system 4-momentum and boost velocity
  TLorentzVector p4hadsyst = utils::rew::Hadronic4pLAB(event);
  TVector3 bv = -1 * p4hadsyst.BoostVector();
  double W = p4hadsyst.M();

  // Loop over stdhep entries and only calculate weights for particles.
  GHepParticle * p = 0;
  TIter event_iter(&event);
  while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {

	int imom = p->FirstMother();

        if(imom==-1) continue;
        if(event.Particle(imom)->Pdg() != kPdgHadronicSyst) continue;

        int pdgc = p->Pdg();

        TLorentzVector p4lab (p->Px(), p->Py(), p->Pz(), p->E());
        TLorentzVector p4hcm(p4lab);
        p4hcm.Boost(bv);

        double pT2 = p4hcm.Px()*p4hcm.Px() + p4hcm.Py()*p4hcm.Py();
        double xF  = p4hcm.Pz() / (W/2.);

        event_weight *= utils::rew::AGKYWeight(pdgc, xF, pT2);
  }

  return event_weight;
}
//_______________________________________________________________________________________
double GReWeightAGKY::CalcChisq(void)
{
  return 0.;
}
//_______________________________________________________________________________________
