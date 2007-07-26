//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - October 03, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "BodekYang/BYStructureFuncModelCC.h"
#include "Conventions/Constants.h"
#include "Messenger/Messenger.h"
#include "PDF/PDF.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
BYStructureFuncModelCC::BYStructureFuncModelCC() :
BYStructureFuncModel("genie::BYStructureFuncModelCC")
{

}
//____________________________________________________________________________
BYStructureFuncModelCC::BYStructureFuncModelCC(string config):
BYStructureFuncModel("genie::BYStructureFuncModelCC", config)
{

}
//____________________________________________________________________________
BYStructureFuncModelCC::~BYStructureFuncModelCC()
{

}
//____________________________________________________________________________
