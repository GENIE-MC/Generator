//_____________________________________________________________________________
/*!

\class    genie::nuvld::NeuGenFitParamsDialog

\brief    A GUI dialog for selecting NeuGEN physics parameters to be fitted

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#include <sstream>

#include "Messenger/Messenger.h"
#include "NuVldGUI/NeuGenFitParamsDialog.h"
#include "NuVldGUI/SysLogSingleton.h"
#include "Utils/StringUtils.h"

using std::ostringstream;

using namespace genie::string_utils;
using namespace genie::nuvld;

ClassImp(NeuGenFitParamsDialog)

//______________________________________________________________________________
bool NeuGenFitParamsDialog::fHaveNoHistory = true;
//______________________________________________________________________________
NeuGenFitParamsDialog::NeuGenFitParamsDialog(
             const TGWindow * p,const TGWindow * main,
                   UInt_t w, UInt_t h, UInt_t options, NeuGenFitParams * ngfp) :
fNGFP(ngfp)
{
  fMain = new TGTransientFrame(p, main, w, h, options);
  fMain->Connect("CloseWindow()",
                  "genie::nuvld::NeuGenFitParamsDialog", this, "CloseWindow()");
  
  fFitParamsGrpFrm = new TGGroupFrame(
                    fMain, "Select NeuGEN fit param and range", kVerticalFrame);

  fFitParamsGrpFrm->SetTitlePos(TGGroupFrame::kLeft);
  fFitParamsGrpFrm->SetLayoutManager(
                                new TGMatrixLayout(fFitParamsGrpFrm, 0, 4, 4));

  fFitParamLb   = new TGLabel( fFitParamsGrpFrm,
                       new TGString(" Select param: ") );
  fFitSpacerLb  = new TGLabel( fFitParamsGrpFrm, new TGString(" "));

  fMinFitParamLb = new TGLabel(
                          fFitParamsGrpFrm, new TGString(" min: "));
  fMaxFitParamLb = new TGLabel(
                          fFitParamsGrpFrm, new TGString(" max: "));
  fStepFitParamLb = new TGLabel(
                         fFitParamsGrpFrm, new TGString(" step: "));

  for(int iparam = 0; iparam < kNNGFitParams; iparam++) {

     // add a check-button for selecting parameters to be be fitted

     fFitParamChkB[iparam] = new TGCheckButton(fFitParamsGrpFrm,
                          fNGFP->ParamAsString(iparam).c_str(), 200 + iparam);

     // add option to set a range for the parameters to be fitted

     fFitParamMinNmEV[iparam] = new TGNumberEntry(
                         fFitParamsGrpFrm, 0, 8, 3, TGNumberFormat::kNESReal);
     fFitParamMaxNmEV[iparam] = new TGNumberEntry(
                         fFitParamsGrpFrm, 0, 8, 3, TGNumberFormat::kNESReal);
     fFitParamStepNmEV[iparam] = new TGNumberEntry(
                         fFitParamsGrpFrm, 0, 8, 3, TGNumberFormat::kNESReal);
  }

  fBtnLt     = new TGLayoutHints(kLHintsBottom | kLHintsCenterX, 2, 2, 2, 2);

  fBtnCmpFrm = new TGCompositeFrame(fMain, 3, 3, kHorizontalFrame);

  fOkBtn     = new TGTextButton(fBtnCmpFrm, "&Ok",     1);
  fCancelBtn = new TGTextButton(fBtnCmpFrm, "&Cancel", 2);
  fResetBtn  = new TGTextButton(fBtnCmpFrm, "&Reset",  3);

  fOkBtn->Connect("Clicked()",
                         "genie::nuvld::NeuGenFitParamsDialog", this, "Ok()");
  fCancelBtn->Connect("Clicked()",
                     "genie::nuvld::NeuGenFitParamsDialog", this, "Cancel()");
  fResetBtn->Connect("Clicked()",
                      "genie::nuvld::NeuGenFitParamsDialog", this, "Reset()");
                     
  fFitParamsGrpFrm -> AddFrame( fFitSpacerLb );
  fFitParamsGrpFrm -> AddFrame( fFitSpacerLb );
  fFitParamsGrpFrm -> AddFrame( fFitSpacerLb );
  fFitParamsGrpFrm -> AddFrame( fFitSpacerLb );

  fFitParamsGrpFrm -> AddFrame( fFitParamLb     );
  fFitParamsGrpFrm -> AddFrame( fMinFitParamLb  );
  fFitParamsGrpFrm -> AddFrame( fMaxFitParamLb  );
  fFitParamsGrpFrm -> AddFrame( fStepFitParamLb );

  fFitParamsGrpFrm -> AddFrame( fFitSpacerLb );
  fFitParamsGrpFrm -> AddFrame( fFitSpacerLb );
  fFitParamsGrpFrm -> AddFrame( fFitSpacerLb );
  fFitParamsGrpFrm -> AddFrame( fFitSpacerLb );
  
  for(int iparam = 0; iparam < kNNGFitParams; iparam++) {

      fFitParamsGrpFrm -> AddFrame( fFitParamChkB[iparam]     );
      fFitParamsGrpFrm -> AddFrame( fFitParamMinNmEV[iparam]  );
      fFitParamsGrpFrm -> AddFrame( fFitParamMaxNmEV[iparam]  );
      fFitParamsGrpFrm -> AddFrame( fFitParamStepNmEV[iparam] );
  }
  fFitParamsGrpFrm -> AddFrame( fFitSpacerLb );
                     
  fBtnCmpFrm->AddFrame (fOkBtn);
  fBtnCmpFrm->AddFrame (fCancelBtn);
  fBtnCmpFrm->AddFrame (fResetBtn);

  fMain->AddFrame (fFitParamsGrpFrm);
  fMain->AddFrame (fBtnCmpFrm, fBtnLt);

  fMain->MapSubwindows();
  fMain->Resize();

  this->PositionRelativeToParent(main); 

  fMain->SetWindowName("NeuGEN Fit Parameters");


  if(fHaveNoHistory) this->Reset();
  else               this->LoadLastEntries();

  fHaveNoHistory = false;
  
  fMain->MapWindow();
}
//______________________________________________________________________________
NeuGenFitParamsDialog::~NeuGenFitParamsDialog()
{
// the order is important here

  LOG("NuVld", pDEBUG) << "Deleting the NeuGenFitParamsDialog...";

  delete fOkBtn;
  delete fCancelBtn;
  delete fResetBtn;

  //delete [] fFitParamChkB;
  //delete [] fFitParamMinNmEV;
  //delete [] fFitParamMaxNmEV;
  //delete [] fFitParamStepNmEV;
   
  for(int iparam = 0; iparam < kNNGFitParams; iparam++)
  {
     delete fFitParamChkB[iparam];
     delete fFitParamMinNmEV[iparam];
     delete fFitParamMaxNmEV[iparam];    
     delete fFitParamStepNmEV[iparam];    
  }
  
  delete fFitParamLb;
  delete fFitSpacerLb;
  delete fMinFitParamLb;
  delete fMaxFitParamLb;  
  delete fStepFitParamLb;  
  delete fBtnLt;
  delete fBtnCmpFrm;
  delete fFitParamsGrpFrm;
  //delete fMain;

  LOG("NuVld", pINFO) << "...Done";   
}
//______________________________________________________________________________
void NeuGenFitParamsDialog::PositionRelativeToParent(const TGWindow * main)
{
// position relative to the parent's window

  int ax, ay;
  Window_t wdum;

  gVirtualX->TranslateCoordinates(main->GetId(), fMain->GetParent()->GetId(),
             (Int_t)(((TGFrame *) main)->GetWidth() - fMain->GetWidth()) >> 1,
             (Int_t)(((TGFrame *) main)->GetHeight() - fMain->GetHeight()) >> 1,
             ax, ay, wdum);
  fMain->Move(ax, ay);
}
//______________________________________________________________________________
void NeuGenFitParamsDialog::Reset(void)
{
  for(int i = 0; i < kNNGFitParams; i++) {

    fFitParamChkB     [i] -> SetState(kButtonUp);
    fFitParamMinNmEV  [i] -> SetNumber(0);
    fFitParamMaxNmEV  [i] -> SetNumber(0);  
    fFitParamStepNmEV [i] -> SetNumber(0);  
  }
}
//______________________________________________________________________________
void NeuGenFitParamsDialog::Ok(void)
{
  this->SetFitParameters();

  this->Report();
    
  fMain->SendCloseMessage();
}
//______________________________________________________________________________
void NeuGenFitParamsDialog::SetFitParameters(void)
{
  for (int ip = 0; ip < kNNGFitParams; ip++)
  {
    fNGFP->fIsFitted[ip] = (fFitParamChkB     [ip] -> GetState() == kButtonDown);
    fNGFP->fRangeMin[ip] =  fFitParamMinNmEV  [ip] -> GetNumber();
    fNGFP->fRangeMax[ip] =  fFitParamMaxNmEV  [ip] -> GetNumber();
    fNGFP->fStep[ip]     =  fFitParamStepNmEV [ip] -> GetNumber();
  }
}
//______________________________________________________________________________
void NeuGenFitParamsDialog::LoadLastEntries(void)
{
  this->Reset();
  
  for (int ip = 0; ip < kNNGFitParams; ip++)
  {
    if(fNGFP->fIsFitted[ip]) 
                   fFitParamChkB [ip] -> SetState(kButtonDown);
    
    fFitParamMinNmEV  [ip] -> SetNumber( fNGFP->fRangeMin[ip] );
    fFitParamMaxNmEV  [ip] -> SetNumber( fNGFP->fRangeMax[ip] );

    fFitParamStepNmEV [ip] -> SetNumber( fNGFP->fStep[ip]     );
  }
}
//______________________________________________________________________________
void NeuGenFitParamsDialog::Report(void) const
{
  SysLogSingleton * syslog = SysLogSingleton::Instance();

  syslog->StatusBar()->SetText( "Setting NeuGEN fit parameters",   0 );
  syslog->Log()->AddLine("Setting NeuGEN fit parameters");

  for (int ip = 0; ip < kNNGFitParams; ip++)
  {
     ostringstream entry;

     entry << "IP: " << ip << ", PARAM: " << fNGFP->ParamAsString(ip);

     if(fNGFP->fIsFitted[ip]) {
           entry << " *** FITTED *** "
                   << " range: (" << fNGFP->fRangeMin[ip]
                               << ", " << fNGFP->fRangeMax[ip] 
                                       << "), step: " << fNGFP->fStep[ip];
    }
    syslog->Log()->AddLine( entry.str().c_str() );
  }
}
//______________________________________________________________________________

