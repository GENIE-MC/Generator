//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - October 03, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Mar 03, 2009 - CA
   Moved into the new RES package from its previous location (EVGModules)

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

