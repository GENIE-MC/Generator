/*

coppying heavily from MECPXSec.cxx



 */

// includes
#include <TMath.h>

#include "Conventions/Constants.h"
#include "Conventions/Controls.h"
#include "Conventions/GBuild.h"
#include "Conventions/Units.h"
#include "GHEP/GHepParticle.h"
#include "Messenger/Messenger.h"
#include "MECTensor/MECTensorPXSec.h"
#include "MECTensor/MECLoadHadTensor.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"
#include "Utils/KineUtils.h"


using namespace genie;
using namespace genie::constants;
using namespace genie::controls;

//____________________________________________________________
MECTensorPXSec::MECTensorPXSec() :
  XSecAlgorithmI("genie::MECTensorPXSec")
{

}

//_____________________________________________________________
MECTensorPXSec::MECTensorPXSec(string config) :
  XSecAlgorithmI("genie::MECPSXec",config)
{

}

//_____________________________________________________________
MECTensorPXSec::~MECTensorPXSec()
{

}

//_____________________________________________________________
double MECTensorPXSec::XSec(
		     const Interaction * interaction, KinePhaseSpace_t kps) const
{
  // respond to a request for a double-differential cross-section.
  // computes Tmu and CosTheta
  // and loads the appropriate hadron tensor (if not already loaded).

  // initial state
  int    tgtpdg = interaction->InitState().Tgt().Pdg();
  int    nupdg  = interaction->InitState().ProbePdg();
  double Enu    = interaction->InitState().ProbeE(kRfLab);
  TLorentzVector * v4Nu = interaction->InitState().GetProbeP4(kRfLab);

  // final state
  const Kinematics kinematics = interaction->Kine();
  TLorentzVector v4lep = kinematics.FSLeptonP4();
  double Mlep = interaction->FSPrimLepton()->Mass();

  // kinematics 
  double Tmu = v4lep.E() - Mlep;
  double CosTheta = cos(v4lep.Theta() - v4Nu->Theta());
  TVector3 lep3 = v4lep.Vect();
  TVector3 nu3 = v4Nu->Vect();  // why are these not both pointers?
  CosTheta = (lep3 * nu3) / (lep3.Mag() * nu3.Mag());
  //doubl CosTheta = dotproduct / magnitudes
  
  // get hadron tensor
  MECLoadHadTensor * hadtensor = MECLoadHadTensor::Instance();

  //double xsec = hadtensor->XSecFullAll(tgtpdg, nupdg, Enu, Tmu, CosTheta);

  double xsec = 0.0;
  double FourXSec[4];

  hadtensor->XSecFour(tgtpdg, nupdg, Enu, Tmu, CosTheta, FourXSec, 0);
  xsec = FourXSec[0];

  return xsec;
}

//_____________________________________________________________
double MECTensorPXSec::Integral(const Interaction * interaction) const {
  // splinemakers and denominators rejoice!
  // return integrated xsec at the given neutrino energy
  // this should return zero for bogus nuclei, so test for it.

  // initial state
  int    tgtpdg = interaction->InitState().Tgt().Pdg();
  int    nupdg  = interaction->InitState().ProbePdg();
  double Enu = interaction->InitState().ProbeE(kRfHitNucRest);

  // get hadron tensor
  MECLoadHadTensor * hadtensor = MECLoadHadTensor::Instance();

  double xSecAtE = hadtensor->TotalXsecAtE(tgtpdg, nupdg, Enu);
  if (!xSecAtE) return 0;
  else return xSecAtE;  

}


//_____________________________________________________________
bool MECTensorPXSec::ValidProcess(const Interaction * interaction) const {
  if (interaction->TestBit(kISkipProcessChk)) return true;

  const ProcessInfo & proc_info = interaction->ProcInfo();
  if (!proc_info.IsMEC()) {
      return false;
  }
  return true;
}

//_____________________________________________________________
void MECTensorPXSec::Configure(const Registry & config){
  Algorithm::Configure(config);
  this->LoadConfig();
}
//_____________________________________________________________
void MECTensorPXSec::LoadConfig(void){
  //  fTensorModel = 0;
  //fPDDComponent = 0;
  //dynamic_cast<const XSecAlgorithmI *> (this->SubAlg("
}
