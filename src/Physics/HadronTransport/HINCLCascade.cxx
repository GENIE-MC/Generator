#include <cstdlib>
#include <sstream>

#include <iostream>
#include <iomanip>
#include <string>

#include <cassert>
#include <cstdlib>
#include <map>
#include <TMath.h>
#include "G4INCLConfig.hh"
#include "G4INCLVersion.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLUnorderedVector.hh"
#include "G4INCLParticle.hh"
#include "G4INCLEventInfo.hh"
#include "G4INCLStore.hh"
#include "G4INCLIAvatar.hh"
#include "G4INCLCascade.hh"
#include "G4INCLConfigEnums.hh"
#include "DatafilePaths.hh" 
#include "G4INCLParticle.hh"


#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Physics/HadronTransport/INukeException.h"
#include "Physics/HadronTransport/Intranuke2018.h"
#include "Physics/HadronTransport/HAIntranuke2018.h"
#include "Physics/HadronTransport/INukeHadroData2018.h"
#include "Physics/HadronTransport/INukeUtils2018.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Physics/NuclearState/NuclearModelI.h"
#include "Physics/NuclearState/NuclearModelMap.h"
#include "G4INCLCascade.hh"
#include "G4INCLIAvatar.hh"
#include "G4INCLIPropagationModel.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLStandardPropagationModel.hh"
#include "G4INCLRandom.hh"
#include "G4INCLRanecu.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLKinematicsUtils.hh"

#include "HINCLCascade.h"
#include "INCLCascade.h"

// signal handler (for Linux and GCC)

using namespace genie;
using namespace genie::utils;

using namespace genie::constants;
using namespace genie::controls;
using namespace G4INCL;
using std::ostringstream;

HINCLCascade::HINCLCascade() :
INCLCascade("genie::HINCLCascade")
{

}
//___________________________________________________________________________
HINCLCascade::HINCLCascade(string config) :
INCLCascade("genie::HINCLCascade",config)
{

}
//___________________________________________________________________________
HINCLCascade::~HINCLCascade()
{

}
void HINCLCascade::ProcessEventRecord(GHepRecord * evrec) const{
  LOG("HINCLCascade", pNOTICE)
     << "************ Running HINCLCascade MODE INTRANUKE ************";

 INCLCascade::ProcessEventRecord(evrec);
//INCLCascade::INCLcascade();
  LOG("HINCLCascade", pINFO) << "Done with this event";
  //this->TransportHadrons(evrec);
}

void HINCLCascade::LoadConfig(void)
{
  LOG("HINCLCAscade", pINFO) << "Settings for INCL++ mode: " ;
 }