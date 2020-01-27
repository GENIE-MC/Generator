//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
*/
//____________________________________________________________________________

#include "Framework/Messenger/Messenger.h"
#include "Physics/DeepInelastic/XSection/QPMDISStrucFunc.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"

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
