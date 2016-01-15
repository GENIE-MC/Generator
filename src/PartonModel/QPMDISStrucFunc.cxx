//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Oct 09, 2009 - CA
   Renamed to QPMDISStrucFunc following code re-organization

*/
//____________________________________________________________________________

#include "Messenger/Messenger.h"
#include "PartonModel/QPMDISStrucFunc.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using namespace genie;

//____________________________________________________________________________
QPMDISStrucFunc::QPMDISStrucFunc() :
QPMDISStrucFuncBase("genie::QPMDISStrucFunc")
{

}
//____________________________________________________________________________
QPMDISStrucFunc::QPMDISStrucFunc(string config):
QPMDISStrucFuncBase("genie::QPMDISStrucFunc", config)
{

}
//____________________________________________________________________________
QPMDISStrucFunc::~QPMDISStrucFunc()
{

}
//____________________________________________________________________________
