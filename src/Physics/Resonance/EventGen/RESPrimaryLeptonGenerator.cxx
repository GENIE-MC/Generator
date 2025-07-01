//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool 
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/GHEP/GHepRecord.h"
#include "Framework/Interaction/Interaction.h"
#include "Physics/Resonance/EventGen/RESPrimaryLeptonGenerator.h"

using namespace genie;

//___________________________________________________________________________
RESPrimaryLeptonGenerator::RESPrimaryLeptonGenerator() :
PrimaryLeptonGenerator("genie::RESPrimaryLeptonGenerator")
{

}
//___________________________________________________________________________
RESPrimaryLeptonGenerator::RESPrimaryLeptonGenerator(string config) :
PrimaryLeptonGenerator("genie::RESPrimaryLeptonGenerator", config)
{

}
//___________________________________________________________________________
RESPrimaryLeptonGenerator::~RESPrimaryLeptonGenerator()
{

}
//___________________________________________________________________________
void RESPrimaryLeptonGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// This method generates the final state primary lepton in RES events

  // no modification is required to the std implementation
  PrimaryLeptonGenerator::ProcessEventRecord(evrec);
}
//___________________________________________________________________________
