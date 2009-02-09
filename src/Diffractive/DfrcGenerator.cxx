//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - October 03, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Conventions/Constants.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "Messenger/Messenger.h"
#include "Diffractive/DfrcGenerator.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
DfrcGenerator::DfrcGenerator() :
EventRecordVisitorI("genie::DfrcGenerator")
{

}
//___________________________________________________________________________
DfrcGenerator::DfrcGenerator(string config) :
EventRecordVisitorI("genie::DfrcGenerator", config)
{

}
//___________________________________________________________________________
DfrcGenerator::~DfrcGenerator()
{

}
//___________________________________________________________________________
void DfrcGenerator::ProcessEventRecord(GHepRecord * /*evrec*/) const
{

}
//___________________________________________________________________________
