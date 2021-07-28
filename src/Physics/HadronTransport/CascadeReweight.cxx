//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 
 Julia Tena-Vidal <j.tena-vidal \at liverpool.ac.uk>
 University of Liverpool 

*/
//____________________________________________________________________________

#include <cstdlib>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/Constants.h"
#include "Physics/HadronTransport/HadronTransporter.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"

using namespace genie;
using namespace genie::utils;
using namespace genie::constants;

//___________________________________________________________________________
CascadeReweight::CascadeReweight() :
EventRecordVisitorI("genie::CascadeReweight")
{

}
//___________________________________________________________________________
CascadeReweight::CascadeReweight(string config) :
EventRecordVisitorI("genie::CascadeReweight", config)
{

}
//___________________________________________________________________________
CascadeReweight::~CascadeReweight()
{

}
//___________________________________________________________________________
void CascadeReweight::ProcessEventRecord(GHepRecord * evrec) const
{

}
//___________________________________________________________________________
void CascadeReweight::GetEventWeight (const GHepRecord & ev) const{

}
//___________________________________________________________________________
void CascadeReweight::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void CascadeReweight::Configure(string param_set)
{
  Algorithm::Configure(param_set);

  Registry * algos = AlgConfigPool::Instance() -> GlobalParameterList() ;
  Registry r( "CascadeReweight_specific", false ) ;

  Algorithm::Configure(r) ;

  this->LoadConfig();
}
//___________________________________________________________________________
void CascadeReweight::LoadConfig(void)
{

}
//___________________________________________________________________________


