//_____________________________________________________________________________
/*!

\class    genie::nuvld::NeuGenConfigDialog

\brief    A GUI dialog for reading NeuGEN physics configuration parameters

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#include "NuVldGUI/NeuGenConfigDialog.h"
#include "NuVldGUI/NeuGenCards.h"
#include "NuVldGUI/SysLogSingleton.h"
#include "Utils/StringUtils.h"

using namespace genie::utils::str;
using namespace genie::nuvld;

ClassImp(NeuGenConfigDialog)

//______________________________________________________________________________
NeuGenConfigDialog::NeuGenConfigDialog(const TGWindow * p,
                      const TGWindow * main, UInt_t w, UInt_t h, UInt_t options)
{
  fMain = new TGTransientFrame(p, main, w, h, options);
  fMain->Connect("CloseWindow()",
                    "genie::nuvld::NeuGenConfigDialog", this, "CloseWindow()");

  // Build sub-Frames 
                    
  this->BuildAxialMassFrame();     //! "Axial Mass" Frame
  this->BuildQelParamsFrame();     //! Etc QEL Params (FA(Q2=0), Eta) Frame
  this->BuildResParamsFrame();     //! Etc RES Params (Omega, Z) Frame
  this->BuildCohParamsFrame();     //! Etc COH Params (Nucl scale, pi Re/Im) Frame
  this->BuildPdfSetFrame();        //! "Parton Density Functions" Frame
  this->BuildKnoParamsFrame();     //!  KNO Hadronization model params Frame
  this->BuildDisResTuningFrame();  //! "DIS/Resonance Tuning" Frame
  this->BuildButtonsFrame();       //! "Buttons" Frame

  // Add Frames to Main Frame

  fMain -> AddFrame( fMaGrpFrm,     fMaFrmLt     );
  fMain -> AddFrame( fQelPrmGrpFrm, fQelPrmFrmLt );
  fMain -> AddFrame( fResPrmGrpFrm, fResPrmFrmLt );
  fMain -> AddFrame( fCohPrmGrpFrm, fCohPrmFrmLt );
  fMain -> AddFrame( fPdfGrpFrm,    fPdfFrmLt    );
  fMain -> AddFrame( fKnoGrpFrm,    fKnoFrmLt    );
  fMain -> AddFrame( fDisResGrpFrm, fDisResFrmLt );
  fMain -> AddFrame( fBtnHFrm,      fBtnFrmLt    );

  fMain->MapSubwindows();
  fMain->Resize();

  this->RestoreDefaults(); // load 'default' values

  this->PositionRelativeToParent(main); // position relative to the parent's window

  fMain->SetWindowName("NeuGEN Config Dialog");

  fMain->MapWindow();
}
//______________________________________________________________________________
NeuGenConfigDialog::~NeuGenConfigDialog()
{
   delete fMaQelNum;
   delete fMaResNum;
   delete fMaCohNum;
   delete fFa0QelNum;
   delete fEtaQelNum;
   delete fOmegaResNum;
   delete fZResNum;
   delete fR0CohNum;
   delete fREICohNum;
   delete fPdfGroupNum;
   delete fPdfSetNum;   
   delete fKnoAvpNum;
   delete fKnoAvnNum;
   delete fKnoAvbpNum;
   delete fKnoAvbnNum;
   delete fKnoBNum;
   delete fDisResM2vpNum;
   delete fDisResM2vnNum;
   delete fDisResM2vbpNum;
   delete fDisResM2vbnNum;
   delete fDisResM3vpNum;
   delete fDisResM3vnNum;
   delete fDisResM3vbpNum;
   delete fDisResM3vbnNum;
   delete fMaQelLbl;
   delete fMaResLbl;
   delete fMaCohLbl;
   delete fFa0QelLbl;
   delete fEtaQelLbl;
   delete fOmegaResLbl;
   delete fZResLbl;
   delete fR0CohLbl;
   delete fREICohLbl;
   delete fPdfGroupLbl;
   delete fPdfSetLbl;
   delete fKnoISvpLbl;
   delete fKnoISvnLbl;
   delete fKnoISvbpLbl;
   delete fKnoISvbnLbl;
   delete fKnoBScaleLbl;
   delete fMultLbl;
   delete fMultEq2Lbl;
   delete fMultEq3Lbl;
   delete fISvpLbl;
   delete fISvnLbl;
   delete fISvbpLbl;
   delete fISvbnLbl;
   delete fUseInputsBtn;
   delete fCancelBtn;
   delete fRestoreDefaultBtn;
   delete fUseDefaultsBtn;
   delete fMaFrmLt;
   delete fQelPrmFrmLt;   
   delete fResPrmFrmLt;
   delete fCohPrmFrmLt;
   delete fPdfFrmLt;
   delete fKnoFrmLt;
   delete fDisResFrmLt;
   delete fBtnFrmLt;
   delete fMaGrpFrm;
   delete fQelPrmGrpFrm;   
   delete fResPrmGrpFrm;
   delete fCohPrmGrpFrm;
   delete fPdfGrpFrm;
   delete fKnoGrpFrm;
   delete fDisResGrpFrm;
   delete fBtnHFrm;
   delete fMain;
}
//______________________________________________________________________________
void NeuGenConfigDialog::BuildAxialMassFrame(void)
{
//! widgets for importing Axial Mass values

  fMaGrpFrm = new TGGroupFrame(fMain, "Axial Vector Mass Ma", kVerticalFrame);
  
  fMaGrpFrm->SetTitlePos(TGGroupFrame::kRight); // right aligned
  fMaGrpFrm->SetLayoutManager(new TGMatrixLayout(fMaGrpFrm, 0, 2, 3));

  fMaQelLbl  = new TGLabel(fMaGrpFrm, new TGString( "Quasi-Elastic Scattering (GeV):"));
  fMaResLbl  = new TGLabel(fMaGrpFrm, new TGString( "Resonance Scattering (GeV)    :"));
  fMaCohLbl  = new TGLabel(fMaGrpFrm, new TGString( "Coherent Scattering (GeV)     :"));

  fMaQelNum = new TGNumberEntry(fMaGrpFrm, 0, 11, 5, TGNumberFormat::kNESReal);
  fMaResNum = new TGNumberEntry(fMaGrpFrm, 0, 11, 5, TGNumberFormat::kNESReal);
  fMaCohNum = new TGNumberEntry(fMaGrpFrm, 0, 11, 5, TGNumberFormat::kNESReal);

  fMaGrpFrm -> AddFrame( fMaQelLbl );
  fMaGrpFrm -> AddFrame( fMaQelNum );
  fMaGrpFrm -> AddFrame( fMaResLbl );
  fMaGrpFrm -> AddFrame( fMaResNum );
  fMaGrpFrm -> AddFrame( fMaCohLbl );
  fMaGrpFrm -> AddFrame( fMaCohNum );

  fMaFrmLt = new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2);
}
//______________________________________________________________________________
void NeuGenConfigDialog::BuildQelParamsFrame(void)
{
//! widgets for importing etc params for QEL models

  fQelPrmGrpFrm = new TGGroupFrame(fMain,
                                 "QEL Scattering model params", kVerticalFrame);

  fQelPrmGrpFrm->SetTitlePos(TGGroupFrame::kRight); // right aligned
  fQelPrmGrpFrm->SetLayoutManager(new TGMatrixLayout(fQelPrmGrpFrm, 0, 2, 2));

  fFa0QelLbl  = new TGLabel(fQelPrmGrpFrm, new TGString( "Axial Form Factor FA (Q^2=0): "));
  fEtaQelLbl  = new TGLabel(fQelPrmGrpFrm, new TGString( "Elastic Scattering Parameter: "));

  fFa0QelNum = new TGNumberEntry(fQelPrmGrpFrm, 0, 11, 5, TGNumberFormat::kNESReal);
  fEtaQelNum = new TGNumberEntry(fQelPrmGrpFrm, 0, 11, 5, TGNumberFormat::kNESReal);

  fQelPrmGrpFrm -> AddFrame( fFa0QelLbl );
  fQelPrmGrpFrm -> AddFrame( fFa0QelNum );
  fQelPrmGrpFrm -> AddFrame( fEtaQelLbl );
  fQelPrmGrpFrm -> AddFrame( fEtaQelNum );

  fQelPrmFrmLt = new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2);
}
//______________________________________________________________________________
void NeuGenConfigDialog::BuildResParamsFrame(void)
{
//! widgets for importing etc params for QEL models

  fResPrmGrpFrm = new TGGroupFrame(fMain,
                                 "RES Scattering model params", kVerticalFrame);

  fResPrmGrpFrm->SetTitlePos(TGGroupFrame::kRight); // right aligned
  fResPrmGrpFrm->SetLayoutManager(new TGMatrixLayout(fResPrmGrpFrm, 0, 2, 2));

  fOmegaResLbl = new TGLabel(fResPrmGrpFrm, new TGString( "Rein-Seghal parameter Omega: "));
  fZResLbl     = new TGLabel(fResPrmGrpFrm, new TGString( "Rein-Seghal parameter Z   :"));

  fOmegaResNum = new TGNumberEntry(fResPrmGrpFrm, 0, 11, 5, TGNumberFormat::kNESReal);
  fZResNum     = new TGNumberEntry(fResPrmGrpFrm, 0, 11, 5, TGNumberFormat::kNESReal);

  fResPrmGrpFrm -> AddFrame( fOmegaResLbl );
  fResPrmGrpFrm -> AddFrame( fOmegaResNum );
  fResPrmGrpFrm -> AddFrame( fZResLbl );
  fResPrmGrpFrm -> AddFrame( fZResNum );

  fResPrmFrmLt = new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2);
}
//______________________________________________________________________________
void NeuGenConfigDialog::BuildCohParamsFrame(void)
{
//! widgets for importing etc params for QEL models

  fCohPrmGrpFrm = new TGGroupFrame(fMain,
                                 "COH Scattering model params", kVerticalFrame);

  fCohPrmGrpFrm->SetTitlePos(TGGroupFrame::kRight); // right aligned
  fCohPrmGrpFrm->SetLayoutManager(new TGMatrixLayout(fCohPrmGrpFrm, 0, 2, 2));

  fR0CohLbl  = new TGLabel(fCohPrmGrpFrm, new TGString( "Nuclear size scale param R0:    "));
  fREICohLbl = new TGLabel(fCohPrmGrpFrm, new TGString( "Pion scattering amplitude Re/Im:"));

  fR0CohNum  = new TGNumberEntry(fCohPrmGrpFrm, 0, 11, 5, TGNumberFormat::kNESReal);
  fREICohNum = new TGNumberEntry(fCohPrmGrpFrm, 0, 11, 5, TGNumberFormat::kNESReal);

  fCohPrmGrpFrm -> AddFrame( fR0CohLbl  );
  fCohPrmGrpFrm -> AddFrame( fR0CohNum  );
  fCohPrmGrpFrm -> AddFrame( fREICohLbl );
  fCohPrmGrpFrm -> AddFrame( fREICohNum );

  fCohPrmFrmLt = new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2);
}
//______________________________________________________________________________
void NeuGenConfigDialog::BuildPdfSetFrame(void)
{
//! widgets for importing PDF Group/Set

  fPdfGrpFrm = new TGGroupFrame(fMain,"Parton Density Functions", kVerticalFrame);
  
  fPdfGrpFrm->SetTitlePos(TGGroupFrame::kRight); // right aligned
  fPdfGrpFrm->SetLayoutManager(new TGMatrixLayout(fPdfGrpFrm, 0, 2, 2));

  fPdfGroupLbl  = new TGLabel(fPdfGrpFrm, new TGString( "  PDGLIB Group code  :") );
  fPdfSetLbl    = new TGLabel(fPdfGrpFrm, new TGString( "  PDFLIB Set code    :") );

  fPdfGroupNum = new TGNumberEntry(fPdfGrpFrm, 0, 12, 0, TGNumberFormat::kNESInteger);
  fPdfSetNum   = new TGNumberEntry(fPdfGrpFrm, 0, 12, 0, TGNumberFormat::kNESInteger);

  fPdfGrpFrm -> AddFrame( fPdfGroupLbl );
  fPdfGrpFrm -> AddFrame( fPdfGroupNum );
  fPdfGrpFrm -> AddFrame( fPdfSetLbl   );
  fPdfGrpFrm -> AddFrame( fPdfSetNum   );

  fPdfFrmLt = new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2);
}
//______________________________________________________________________________
void NeuGenConfigDialog::BuildKnoParamsFrame(void)
{
//! widgets for tuning DIS/Resonances

  fKnoGrpFrm = new TGGroupFrame(fMain,
                              "KNO Hadronization model params", kVerticalFrame);

  fKnoGrpFrm->SetTitlePos(TGGroupFrame::kRight); // right aligned
  fKnoGrpFrm->SetLayoutManager(new TGMatrixLayout(fKnoGrpFrm, 0, 5, 2));

  fKnoISvpLbl    = new TGLabel(fKnoGrpFrm, new TGString( " A(v-p) ")    );
  fKnoISvnLbl    = new TGLabel(fKnoGrpFrm, new TGString( " A(v-n) ")    );
  fKnoISvbpLbl   = new TGLabel(fKnoGrpFrm, new TGString( " A(vbar-p) ") );
  fKnoISvbnLbl   = new TGLabel(fKnoGrpFrm, new TGString( " A(vbar-n) ") );
  fKnoBScaleLbl  = new TGLabel(fKnoGrpFrm, new TGString( " B ") );
  
  fKnoAvpNum  = new TGNumberEntry(fKnoGrpFrm, 0, 6, 3, TGNumberFormat::kNESReal);
  fKnoAvnNum  = new TGNumberEntry(fKnoGrpFrm, 0, 6, 3, TGNumberFormat::kNESReal);
  fKnoAvbpNum = new TGNumberEntry(fKnoGrpFrm, 0, 6, 3, TGNumberFormat::kNESReal);
  fKnoAvbnNum = new TGNumberEntry(fKnoGrpFrm, 0, 6, 3, TGNumberFormat::kNESReal);
  fKnoBNum    = new TGNumberEntry(fKnoGrpFrm, 0, 6, 3, TGNumberFormat::kNESReal);

  fKnoGrpFrm -> AddFrame( fKnoISvpLbl   );
  fKnoGrpFrm -> AddFrame( fKnoISvnLbl   );
  fKnoGrpFrm -> AddFrame( fKnoISvbpLbl  );
  fKnoGrpFrm -> AddFrame( fKnoISvbnLbl  );
  fKnoGrpFrm -> AddFrame( fKnoBScaleLbl );
  fKnoGrpFrm -> AddFrame( fKnoAvpNum    );
  fKnoGrpFrm -> AddFrame( fKnoAvnNum    );
  fKnoGrpFrm -> AddFrame( fKnoAvbpNum   );
  fKnoGrpFrm -> AddFrame( fKnoAvbnNum   );
  fKnoGrpFrm -> AddFrame( fKnoBNum      );

  fKnoFrmLt = new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2);
}
//______________________________________________________________________________
void NeuGenConfigDialog::BuildDisResTuningFrame(void)
{
//! widgets for tuning DIS/Resonances

  fDisResGrpFrm = new TGGroupFrame(fMain, "DIS / Renonance tuning", kVerticalFrame);
  
  fDisResGrpFrm->SetTitlePos(TGGroupFrame::kRight); // right aligned
  fDisResGrpFrm->SetLayoutManager(new TGMatrixLayout(fDisResGrpFrm, 0, 5, 3));

  fMultLbl    = new TGLabel(fDisResGrpFrm, new TGString( " multiplicity "));
  fMultEq2Lbl = new TGLabel(fDisResGrpFrm, new TGString( " 2 ")      );
  fMultEq3Lbl = new TGLabel(fDisResGrpFrm, new TGString( " 3 ")      );
  fISvpLbl    = new TGLabel(fDisResGrpFrm, new TGString( " v-p ")    );
  fISvnLbl    = new TGLabel(fDisResGrpFrm, new TGString( " v-n ")    );
  fISvbpLbl   = new TGLabel(fDisResGrpFrm, new TGString( " vbar-p ") );
  fISvbnLbl   = new TGLabel(fDisResGrpFrm, new TGString( " vbar-n ") );

  fDisResM2vpNum  = new TGNumberEntry(fDisResGrpFrm, 0, 6, 1, TGNumberFormat::kNESReal);
  fDisResM2vnNum  = new TGNumberEntry(fDisResGrpFrm, 0, 6, 1, TGNumberFormat::kNESReal);
  fDisResM2vbpNum = new TGNumberEntry(fDisResGrpFrm, 0, 6, 1, TGNumberFormat::kNESReal);
  fDisResM2vbnNum = new TGNumberEntry(fDisResGrpFrm, 0, 6, 1, TGNumberFormat::kNESReal);
  fDisResM3vpNum  = new TGNumberEntry(fDisResGrpFrm, 0, 6, 1, TGNumberFormat::kNESReal);
  fDisResM3vnNum  = new TGNumberEntry(fDisResGrpFrm, 0, 6, 1, TGNumberFormat::kNESReal);
  fDisResM3vbpNum = new TGNumberEntry(fDisResGrpFrm, 0, 6, 1, TGNumberFormat::kNESReal);
  fDisResM3vbnNum = new TGNumberEntry(fDisResGrpFrm, 0, 6, 1, TGNumberFormat::kNESReal);

  fDisResGrpFrm -> AddFrame( fMultLbl         );
  fDisResGrpFrm -> AddFrame( fISvpLbl         );
  fDisResGrpFrm -> AddFrame( fISvnLbl         );
  fDisResGrpFrm -> AddFrame( fISvbpLbl        );
  fDisResGrpFrm -> AddFrame( fISvbnLbl        );
  fDisResGrpFrm -> AddFrame( fMultEq2Lbl      );
  fDisResGrpFrm -> AddFrame( fDisResM2vpNum   );
  fDisResGrpFrm -> AddFrame( fDisResM2vnNum   );
  fDisResGrpFrm -> AddFrame( fDisResM2vbpNum  );
  fDisResGrpFrm -> AddFrame( fDisResM2vbnNum  );
  fDisResGrpFrm -> AddFrame( fMultEq3Lbl      );
  fDisResGrpFrm -> AddFrame( fDisResM3vpNum   );
  fDisResGrpFrm -> AddFrame( fDisResM3vnNum   );
  fDisResGrpFrm -> AddFrame( fDisResM3vbpNum  );
  fDisResGrpFrm -> AddFrame( fDisResM3vbnNum  );

  fDisResFrmLt = new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2);
}
//______________________________________________________________________________
void NeuGenConfigDialog::BuildButtonsFrame(void)
{
//! NeuGEN Config dialog - control buttons

  fBtnHFrm = new TGHorizontalFrame(fMain, 10, 10);

  fUseInputsBtn      = new TGTextButton(fBtnHFrm, "Use &Input Values", 1);
  fCancelBtn         = new TGTextButton(fBtnHFrm, "&Cancel",           2);
  fRestoreDefaultBtn = new TGTextButton(fBtnHFrm, "&Restore Defaults", 3);
  fUseDefaultsBtn    = new TGTextButton(fBtnHFrm, "Use &Defaults",     4);

  fUseInputsBtn      -> Connect("Clicked()", "genie::nuvld::NeuGenConfigDialog", this, "UseNewValues()");
  fCancelBtn         -> Connect("Clicked()", "genie::nuvld::NeuGenConfigDialog", this, "Cancel()");
  fRestoreDefaultBtn -> Connect("Clicked()", "genie::nuvld::NeuGenConfigDialog", this, "RestoreDefaults()");
  fUseDefaultsBtn    -> Connect("Clicked()", "genie::nuvld::NeuGenConfigDialog", this, "UseDefaults()");

  fBtnHFrm->AddFrame( fUseInputsBtn       );
  fBtnHFrm->AddFrame( fUseDefaultsBtn     );
  fBtnHFrm->AddFrame( fRestoreDefaultBtn );
  fBtnHFrm->AddFrame( fCancelBtn           );

  fBtnFrmLt = new TGLayoutHints(kLHintsBottom | kLHintsCenterX,  2, 2, 2, 2);
}
//______________________________________________________________________________
void NeuGenConfigDialog::RestoreDefaults(void)
{
  NeuGenConfig * default_conf = new NeuGenConfig("std_neugen");

  fMaQelNum       -> SetNumber( default_conf -> MaQel()          );
  fMaResNum       -> SetNumber( default_conf -> MaRes()          );
  fMaCohNum       -> SetNumber( default_conf -> MaCoh()          );
  
  fFa0QelNum      -> SetNumber( default_conf -> Fa0Qel()         );
  fEtaQelNum      -> SetNumber( default_conf -> EtaQel()         );

  fOmegaResNum    -> SetNumber( default_conf -> OmegaRes()       );
  fZResNum        -> SetNumber( default_conf -> ZRes()           );

  fR0CohNum       -> SetNumber( default_conf -> R0Coh()          );
  fREICohNum      -> SetNumber( default_conf -> REICoh()         );

  fPdfGroupNum    -> SetNumber( default_conf -> PdfGroup()       );
  fPdfSetNum      -> SetNumber( default_conf -> PdfSet()         );

  fKnoAvpNum      -> SetNumber( default_conf -> KnoA(e_vp)       );
  fKnoAvnNum      -> SetNumber( default_conf -> KnoA(e_vn)       );
  fKnoAvbpNum     -> SetNumber( default_conf -> KnoA(e_vbp)      );
  fKnoAvbnNum     -> SetNumber( default_conf -> KnoA(e_vbn)      );
  fKnoBNum        -> SetNumber( default_conf -> KnoB()           );

  fDisResM2vpNum  -> SetNumber( default_conf -> DisRes(2, e_vp)  );
  fDisResM2vnNum  -> SetNumber( default_conf -> DisRes(2, e_vn)  );
  fDisResM2vbpNum -> SetNumber( default_conf -> DisRes(2, e_vbp) );
  fDisResM2vbnNum -> SetNumber( default_conf -> DisRes(2, e_vbn) );
  fDisResM3vpNum  -> SetNumber( default_conf -> DisRes(3, e_vp)  );
  fDisResM3vnNum  -> SetNumber( default_conf -> DisRes(3, e_vn)  );
  fDisResM3vbpNum -> SetNumber( default_conf -> DisRes(3, e_vbp) );
  fDisResM3vbnNum -> SetNumber( default_conf -> DisRes(3, e_vbn) );

  delete default_conf;
}
//______________________________________________________________________________
void NeuGenConfigDialog::UseNewValues(void)
{
  NeuGenCards * cards = NeuGenCards::Instance();

  cards -> CurrConfig() -> SetMaQel    ( this->ReadMaQel()     );
  cards -> CurrConfig() -> SetMaRes    ( this->ReadMaRes()     );
  cards -> CurrConfig() -> SetMaCoh    ( this->ReadMaCoh()     );
  
  cards -> CurrConfig() -> SetFa0Qel    ( this->ReadQelFa0()   );
  cards -> CurrConfig() -> SetEtaQel    ( this->ReadQelEta()   );

  cards -> CurrConfig() -> SetOmegaRes  ( this->ReadResOmega() );
  cards -> CurrConfig() -> SetZRes      ( this->ReadResZ()     );

  cards -> CurrConfig() -> SetR0Coh     ( this->ReadCohR0()    );
  cards -> CurrConfig() -> SetREICoh    ( this->ReadCohREI()   );

  cards -> CurrConfig() -> SetPdfGroup ( this->ReadPdfGroup()  );
  cards -> CurrConfig() -> SetPdfSet   ( this->ReadPdfSet()    );
  
  cards -> CurrConfig() -> SetKnoA ( e_vp,  this->ReadKnoAvp()  );
  cards -> CurrConfig() -> SetKnoA ( e_vn,  this->ReadKnoAvn()  );
  cards -> CurrConfig() -> SetKnoA ( e_vbp, this->ReadKnoAvbp() );
  cards -> CurrConfig() -> SetKnoA ( e_vbn, this->ReadKnoAvbn() );
  cards -> CurrConfig() -> SetKnoB ( this->ReadKnoB()           );

  cards -> CurrConfig() -> SetDisRes   ( 2, e_vp,   this->ReadDR2vp()  );
  cards -> CurrConfig() -> SetDisRes   ( 2, e_vn,   this->ReadDR2vn()  );
  cards -> CurrConfig() -> SetDisRes   ( 2, e_vbp,  this->ReadDR2vbp() );
  cards -> CurrConfig() -> SetDisRes   ( 2, e_vbn,  this->ReadDR2vbn() );
  cards -> CurrConfig() -> SetDisRes   ( 3, e_vp,   this->ReadDR3vp()  );
  cards -> CurrConfig() -> SetDisRes   ( 3, e_vn,   this->ReadDR3vn()  );
  cards -> CurrConfig() -> SetDisRes   ( 3, e_vbp,  this->ReadDR3vbp() );
  cards -> CurrConfig() -> SetDisRes   ( 3, e_vbn,  this->ReadDR3vbn() );

  this->Report();
  
  fMain->SendCloseMessage();
}
//______________________________________________________________________________
void NeuGenConfigDialog::UseDefaults(void)
{
  NeuGenCards * cards = NeuGenCards::Instance();

  cards->SetDefaultConfig();

  this->Report();
  
  fMain->SendCloseMessage();
}
//______________________________________________________________________________
void NeuGenConfigDialog::PositionRelativeToParent(const TGWindow * main)
{
//! position relative to the parent's window

  int ax, ay;
  Window_t wdum;

  gVirtualX->TranslateCoordinates(main->GetId(), fMain->GetParent()->GetId(),
             (Int_t)(((TGFrame *) main)->GetWidth() - fMain->GetWidth()) >> 1,
             (Int_t)(((TGFrame *) main)->GetHeight() - fMain->GetHeight()) >> 1,
             ax, ay, wdum);

  fMain->Move(ax, ay);
}
//______________________________________________________________________________
void NeuGenConfigDialog::Report(void) const
{
  SysLogSingleton * syslog = SysLogSingleton::Instance();

  syslog->StatusBar()->SetText( "Configuring NeuGEN",   0 );
  syslog->Log()->AddLine("Configuring NeuGEN");

  syslog->Log()->AddLine( Concat("Axial Mass (QES)... ", this->ReadMaQel()   ) );
  syslog->Log()->AddLine( Concat("Axial Mass (RES)... ", this->ReadMaRes()   ) );
  syslog->Log()->AddLine( Concat("Axial Mass (COH)... ", this->ReadMaCoh()   ) );
  syslog->Log()->AddLine( Concat("QEL FA(Q^2 = 0).... ", this->ReadQelFa0()  ) );
  syslog->Log()->AddLine( Concat("QEL Eta............ ", this->ReadQelEta()  ) );
  syslog->Log()->AddLine( Concat("RES R-S Omega...... ", this->ReadResOmega()) );
  syslog->Log()->AddLine( Concat("RES R-S Z...........", this->ReadResZ()    ) ); 
  syslog->Log()->AddLine( Concat("COH Nucl R0 scale.. ", this->ReadCohR0()   ) );
  syslog->Log()->AddLine( Concat("COH pi Re/Im Ampl...", this->ReadCohREI()  ) );
  syslog->Log()->AddLine( Concat("PDF Group...........", this->ReadPdfGroup()) );
  syslog->Log()->AddLine( Concat("PDF Set.............", this->ReadPdfSet()  ) );
  syslog->Log()->AddLine( Concat("KNO Param A /v -p...", this->ReadKnoAvp()  ) );
  syslog->Log()->AddLine( Concat("KNO Param A /v -n...", this->ReadKnoAvn()  ) );
  syslog->Log()->AddLine( Concat("KNO Param A /vb-p...", this->ReadKnoAvbp() ) );
  syslog->Log()->AddLine( Concat("KNO Param A /vb-n...", this->ReadKnoAvbn() ) );
  syslog->Log()->AddLine( Concat("KNO Param B.........", this->ReadKnoB()    ) );
  syslog->Log()->AddLine( Concat("DIS/RES - m=2/v -p..", this->ReadDR2vp()   ) );
  syslog->Log()->AddLine( Concat("DIS/RES - m=2/v -n..", this->ReadDR2vn()   ) );
  syslog->Log()->AddLine( Concat("DIS/RES - m=2/vb-p..", this->ReadDR2vbp()  ) );
  syslog->Log()->AddLine( Concat("DIS/RES - m=2/vb-n..", this->ReadDR2vbn()  ) );
  syslog->Log()->AddLine( Concat("DIS/RES - m=3/v -p..", this->ReadDR3vp()   ) );
  syslog->Log()->AddLine( Concat("DIS/RES - m=3/v -n..", this->ReadDR3vn()   ) );
  syslog->Log()->AddLine( Concat("DIS/RES - m=3/vb-p..", this->ReadDR3vbp()  ) );
  syslog->Log()->AddLine( Concat("DIS/RES - m=3/vb-n..", this->ReadDR3vbn()  ) );  
}
//______________________________________________________________________________
int NeuGenConfigDialog::ReadPdfGroup(void) const
{
  return (int) fPdfGroupNum->GetNumber();
}  
//______________________________________________________________________________
int NeuGenConfigDialog::ReadPdfSet(void) const
{
  return (int) fPdfSetNum->GetNumber();
}
//______________________________________________________________________________
float NeuGenConfigDialog::ReadMaQel(void) const
{
  return (float) fMaQelNum->GetNumber();
}
//______________________________________________________________________________
float NeuGenConfigDialog::ReadMaRes(void) const  
{
  return (float) fMaResNum->GetNumber();
}
//______________________________________________________________________________
float NeuGenConfigDialog::ReadMaCoh(void) const
{
  return (float) fMaCohNum->GetNumber();
}
//______________________________________________________________________________
float NeuGenConfigDialog::ReadQelFa0(void) const
{
  return (float) fFa0QelNum->GetNumber();
}
//______________________________________________________________________________
float NeuGenConfigDialog::ReadQelEta(void) const
{
  return (float) fEtaQelNum->GetNumber();
}
//______________________________________________________________________________
float NeuGenConfigDialog::ReadResOmega(void) const
{
  return (float) fOmegaResNum->GetNumber();
}
//______________________________________________________________________________
float NeuGenConfigDialog::ReadResZ(void) const
{
  return (float) fZResNum->GetNumber();
}
//______________________________________________________________________________
float NeuGenConfigDialog::ReadCohR0(void) const
{
  return (float) fR0CohNum->GetNumber();
}
//______________________________________________________________________________
float NeuGenConfigDialog::ReadCohREI(void) const
{
  return (float) fREICohNum->GetNumber();
}
//______________________________________________________________________________
float NeuGenConfigDialog::ReadKnoAvp(void) const
{
  return (float) fKnoAvpNum->GetNumber();
}
//______________________________________________________________________________
float NeuGenConfigDialog::ReadKnoAvn(void) const
{
  return (float) fKnoAvnNum->GetNumber();
}
//______________________________________________________________________________
float NeuGenConfigDialog::ReadKnoAvbp(void) const
{
  return (float) fKnoAvbpNum->GetNumber();
}
//______________________________________________________________________________
float NeuGenConfigDialog::ReadKnoAvbn(void) const
{
  return (float) fKnoAvbnNum->GetNumber();
}
//______________________________________________________________________________
float NeuGenConfigDialog::ReadKnoB(void) const
{
  return (float) fKnoBNum->GetNumber();
}
//______________________________________________________________________________
float NeuGenConfigDialog::ReadDR2vp(void) const
{
  return (float) fDisResM2vpNum->GetNumber();
}
//______________________________________________________________________________
float NeuGenConfigDialog::ReadDR2vn(void) const 
{
  return (float) fDisResM2vnNum->GetNumber();
}
//______________________________________________________________________________
float NeuGenConfigDialog::ReadDR2vbp(void) const
{
  return (float) fDisResM2vbpNum->GetNumber();
}
//______________________________________________________________________________
float NeuGenConfigDialog::ReadDR2vbn(void) const
{
  return (float) fDisResM2vbnNum->GetNumber();
}
//______________________________________________________________________________
float NeuGenConfigDialog::ReadDR3vp(void) const
{
  return (float) fDisResM3vpNum->GetNumber();
}
//______________________________________________________________________________
float NeuGenConfigDialog::ReadDR3vn(void) const
{
  return (float) fDisResM3vnNum->GetNumber();
}
//______________________________________________________________________________
float NeuGenConfigDialog::ReadDR3vbp(void) const
{
  return (float) fDisResM3vbpNum->GetNumber();
}
//______________________________________________________________________________
float NeuGenConfigDialog::ReadDR3vbn(void) const
{
  return (float) fDisResM3vbnNum->GetNumber();
}
//______________________________________________________________________________

