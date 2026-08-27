//____________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org


 Author: Zach Orr
         Colorado State University

         Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Physics/BoostedDarkMatter/EventGen/DMRESOutgoingDarkGenerator.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Messenger/Messenger.h"

using namespace genie;

//___________________________________________________________________________
DMRESOutgoingDarkGenerator::DMRESOutgoingDarkGenerator() :
OutgoingDarkGenerator("genie::DMRESOutgoingDarkGenerator")
{

}
//___________________________________________________________________________
DMRESOutgoingDarkGenerator::DMRESOutgoingDarkGenerator(string config) :
OutgoingDarkGenerator("genie::DMRESOutgoingDarkGenerator", config)
{

}
//___________________________________________________________________________
DMRESOutgoingDarkGenerator::~DMRESOutgoingDarkGenerator()
{

}
//___________________________________________________________________________
void DMRESOutgoingDarkGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// This method generates the final state primary lepton in RES events

  // no modification is required to the std implementation
  OutgoingDarkGenerator::ProcessEventRecord(evrec);
}
//___________________________________________________________________________
