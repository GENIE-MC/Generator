//____________________________________________________________________________
/*
 Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - Sep 22, 2008

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 22, 2008 - CA
   This event generation modules was first added in version 2.5.1 as part of
   the new event generation thread handling MEC interactions. 

*/
//____________________________________________________________________________

#include <TMath.h>

#include "MEC/MECGenerator.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"

using namespace genie;

//___________________________________________________________________________
MECGenerator::MECGenerator() :
EventRecordVisitorI("genie::MECGenerator")
{

}
//___________________________________________________________________________
MECGenerator::MECGenerator(string config) :
EventRecordVisitorI("genie::MECGenerator", config)
{

}
//___________________________________________________________________________
MECGenerator::~MECGenerator()
{

}
//___________________________________________________________________________
void MECGenerator::ProcessEventRecord(GHepRecord * /*evrec*/) const
{

}
//___________________________________________________________________________
