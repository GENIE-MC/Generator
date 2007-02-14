//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - October 03, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>

#include "EVGModules/RESPrimaryLeptonGenerator.h"
#include "GHEP/GHepRecord.h"
#include "Interaction/Interaction.h"

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

