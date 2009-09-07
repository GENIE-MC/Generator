//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 25, 2008 - CA
   This event generation modules was first added in version 2.3.1 as part of
   the new event generation thread handling amonaly-mediated single gamma
   interactions. 

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Conventions/Constants.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"
#include "VHE/GlashowResonanceGenerator.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
GlashowResonanceGenerator::GlashowResonanceGenerator() :
EventRecordVisitorI("genie::GlashowResonanceGenerator")
{

}
//___________________________________________________________________________
GlashowResonanceGenerator::GlashowResonanceGenerator(string config) :
EventRecordVisitorI("genie::GlashowResonanceGenerator", config)
{

}
//___________________________________________________________________________
GlashowResonanceGenerator::~GlashowResonanceGenerator()
{

}
//___________________________________________________________________________
void GlashowResonanceGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
  this -> PickHitElectron (evrec);
  this -> AddResonance    (evrec);
}
//___________________________________________________________________________
void GlashowResonanceGenerator::PickHitElectron(GHepRecord * evrec) const
{
//
//

}
//___________________________________________________________________________
void GlashowResonanceGenerator::AddResonance(GHepRecord * evrec) const
{
//
//
/*
  TLorentzVector p4_nu (* evrec -> Probe()       -> P4());
  TLorentzVector p4_e  (* evrec -> HitElectron() -> P4());

  double etot = p4_nu.E() + p4_e.E();

  TVector3 cm_boost = -1 * (p4_nu.Vect() + p4_e.Vect()) / etot;

  p4_nu.Boost(cm_boost);
  p4_e. Boost(cm_boost);

  TLorentzVector p4_W = p4_nu + p4_e;

  p4

  evrec->
*/
}
//___________________________________________________________________________
