//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - May 02, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Messenger/Messenger.h"
#include "PartonModel/DISStructureFuncModelCC.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using namespace genie;

//____________________________________________________________________________
DISStructureFuncModelCC::DISStructureFuncModelCC() :
DISStructureFuncModel("genie::DISStructureFuncModelCC")
{

}
//____________________________________________________________________________
DISStructureFuncModelCC::DISStructureFuncModelCC(string config):
DISStructureFuncModel("genie::DISStructureFuncModelCC", config)
{

}
//____________________________________________________________________________
DISStructureFuncModelCC::~DISStructureFuncModelCC()
{

}
//____________________________________________________________________________
