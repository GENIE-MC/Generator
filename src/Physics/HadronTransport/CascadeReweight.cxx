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
  // Get Associated weight
  // Set weight 

}
//___________________________________________________________________________
void CascadeReweight::GetEventWeight (const GHepRecord & ev) const{
  // Get faith
  // Return weight given the faith
  // List of faiths

  /*
  //__________________________________________________________________________
  static string AsString(INukeFateHN_t fate) {
     switch (fate) {
      case kIHNFtUndefined : return "** Undefined HN-mode fate **"; break;
      case kIHNFtCEx       : return "HN-mode / cex";    break;
      case kIHNFtElas      : return "HN-mode / elas";   break;
      case kIHNFtInelas    : return "HN-mode / inelas"; break;
      case kIHNFtAbs       : return "HN-mode / abs";    break;
      case kIHNFtCmp   : return "HN-mode / compound"; break;
      case kIHNFtNoInteraction : return "HN-mode / no interaction"; break;
      default              : break; 
     }
     return "** Undefined HN-mode fate **"; 
  }
  //__________________________________________________________________________
  static string AsString(INukeFateHA_t fate) {
     switch (fate) {
      case kIHAFtUndefined : return "** Undefined HA-mode fate **"; break;
      case kIHAFtNoInteraction : return "HA-mode / no interaction"; break;
      case kIHAFtCEx       : return "HA-mode / cex";            break;
      //      case kIHAFtElas      : return "HA-mode / elas";           break;
      case kIHAFtInelas    : return "HA-mode / inelas";         break;
      case kIHAFtAbs       : return "HA-mode / abs";            break;
      case kIHAFtKo        : return "HA-mode / knock-out";      break;
      case kIHAFtCmp       : return "HA-mode / compound";       break;
      case kIHAFtPiProd    : return "HA-mode / pi-production" ; break;
      case kIHAFtInclPip   : return "HA-mode / pi-prod incl pi+";   break;
      case kIHAFtInclPim   : return "HA-mode / pi-prod incl pi-";   break;
      case kIHAFtInclPi0   : return "HA-mode / pi-prod incl pi0";   break;
      case kIHAFtDCEx      : return "HA-mode / dcex";           break;
      default              : break;
     }
     return "** Undefined HA-mode fate **"; 
  }
   */
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

  // Read xml configuration
  // Here we need to store the scaling factors associated to each weight

}
//___________________________________________________________________________


