//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Messenger/Messenger.h"
#include "PartonModel/DISStructureFuncModelEM.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using namespace genie;

//____________________________________________________________________________
DISStructureFuncModelEM::DISStructureFuncModelEM() :
DISStructureFuncModel("genie::DISStructureFuncModelEM")
{

}
//____________________________________________________________________________
DISStructureFuncModelEM::DISStructureFuncModelEM(string config):
DISStructureFuncModel("genie::DISStructureFuncModelEM", config)
{

}
//____________________________________________________________________________
DISStructureFuncModelEM::~DISStructureFuncModelEM()
{

}
//____________________________________________________________________________
