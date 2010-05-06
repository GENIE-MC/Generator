//____________________________________________________________________________
/*
 Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Apr 27, 2010 - CA
   First included in v2.7.1.

*/
//____________________________________________________________________________

#include "Conventions/Controls.h"
#include "EVGCore/EventRecord.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "ReWeight/GReWeightFGM.h"

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GReWeightFGM::GReWeightFGM() :
GReWeightI()
{

}
//_______________________________________________________________________________________
GReWeightFGM::~GReWeightFGM()
{

}
//_______________________________________________________________________________________
bool GReWeightFGM::IsHandled(GSyst_t syst)
{
  switch(syst) {
    case ( kSystNucl_CCQEPauliSupViaKF   ) : 
    case ( kSystNucl_CCQEMomDistroFGtoSF ) : 
        return true;
        break;
    default:
        return false;
        break;
  }
  return false;
}
//_______________________________________________________________________________________
void GReWeightFGM::SetSystematic(GSyst_t syst, double val)
{
  switch(syst) {
    case ( kSystNucl_CCQEPauliSupViaKF   ) : 
        fKFTwkDial = val;
        break;
    case ( kSystNucl_CCQEMomDistroFGtoSF ) : 
        fMomDistroTwkDial = val;
        break;
    default:
        return;
        break;
  }
}
//_______________________________________________________________________________________
void GReWeightFGM::Reset(void)
{
  fKFTwkDial        = 1.;
  fMomDistroTwkDial = 1.;
}
//_______________________________________________________________________________________
void GReWeightFGM::Reconfigure(void)
{
  
}
//_______________________________________________________________________________________
double GReWeightFGM::CalcWeight(const EventRecord & event) 
{ 
  double wght =
    this->RewCCQEPauliSupViaKF   (event) *
    this->RewCCQEMomDistroFGtoSF (event);
 
  return wght;
}
//_______________________________________________________________________________________
double GReWeightFGM::CalcChisq(void)
{
  return 0.;
}
//_______________________________________________________________________________________
double GReWeightFGM::RewCCQEPauliSupViaKF(const EventRecord & event) 
{
  bool kF_tweaked = (TMath::Abs(fKFTwkDial) < controls::kASmallNum);
  if(!kF_tweaked) return 1.;

  bool is_qe = event.Summary()->ProcInfo().IsQuasiElastic();
  bool is_cc = event.Summary()->ProcInfo().IsWeakCC();
  if(!is_qe || !is_cc) return 1.;

  double wght = 1.;

//  double kF_def = 0;/////
//  double kF_err = 0;/////
//  double kF_twk = kF_def * (1 + fKFTwkDial * kF_err);

  //
  // ...
  // ...
  // ...
  //


  return wght;
}
//_______________________________________________________________________________________
double GReWeightFGM::RewCCQEMomDistroFGtoSF(const EventRecord & event) 
{
  bool momdistro_tweaked = (TMath::Abs(fMomDistroTwkDial) < controls::kASmallNum);
  if(!momdistro_tweaked) return 1.;

  bool is_qe = event.Summary()->ProcInfo().IsQuasiElastic();
  bool is_cc = event.Summary()->ProcInfo().IsWeakCC();
  if(!is_qe || !is_cc) return 1.;

  double wght = 1.;

  //
  // ...
  // ...
  // ...
  //


  return wght;
}
//_______________________________________________________________________________________
