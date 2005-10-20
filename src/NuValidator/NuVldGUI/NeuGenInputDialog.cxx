//_____________________________________________________________________________
/*!

\class    genie::nuvld::NeuGenInputDialog

\brief    A GUI dialog for specifying what information to extract from a
          neutrino generator.

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#include <iostream>

#include <TGraph.h>
#include <TCanvas.h>
#include <TGFrame.h>
#include <TGLabel.h>
#include <TGButton.h>
#include <TGComboBox.h>
#include <TGListBox.h>
#include <TGNumberEntry.h>

#include "NuVldGUI/NeuGenInputDialog.h"
#include "NuVldGUI/SysLogSingleton.h"
#include "NuVldGUI/BrowserSingleton.h"
#include "NuVldGUI/NeuGenCards.h"
#include "NuVldGUI/NeuGenConstants.h"
#include "Utils/StringUtils.h"
#include "Utils/GUIUtils.h"

using std::cout;
using std::endl;

using namespace genie;
using namespace genie::nuvld;
using namespace genie::utils::str;
using namespace genie::nuvld::facades;

ClassImp(NeuGenInputDialog)

//______________________________________________________________________________
bool NeuGenInputDialog::fHaveNoHistory = true;
//______________________________________________________________________________
NeuGenInputDialog::NeuGenInputDialog(const TGWindow * p,
                      const TGWindow * main, UInt_t w, UInt_t h, UInt_t options)
{
  fMain = new TGTransientFrame(p, main, w, h, options);
  fMain->Connect("CloseWindow()",
                      "genie::nuvld::NeuGenInputDialog", this, "CloseWindow()");

  // Build frames

  fPlotTypeGrpf    = this->BuildXSecTypeFrame();
  fTotXSecConfGrpf = this->BuildConfigTotalXSecFrame();
  fDifXSecConfGrpf = this->BuildConfigDiffXSecFrame();
  fSFConfGrpf      = this->BuildConfigSFFrame();
  fPlotVarGrpf     = this->BuildPlotVarFrame();

  fInteractionGrpf = this->BuildInteractionFrame();

  fCutsSumsHFrm   = new TGHorizontalFrame(fMain, 10, 10);

  fCutsGrpf = this->BuildCutsFrame();
  fSumsGrpf = this->BuildSumsFrame();

  fCutsSumsHFrm -> AddFrame( fCutsGrpf );
  fCutsSumsHFrm -> AddFrame( fSumsGrpf );

  fBtnHFrm = this->BuildButtonsFrame();

  // Build frame layouts

  fPlotTypeFrmLt    = new TGLayoutHints(kLHintsTop | kLHintsCenterX,  2, 2, 2, 2);
  fTotXSecConfFrmLt = new TGLayoutHints(kLHintsTop | kLHintsCenterX,  2, 2, 2, 2);
  fDifXSecConfFrmLt = new TGLayoutHints(kLHintsTop | kLHintsCenterX,  2, 2, 2, 2);
  fSFConfFrmLt      = new TGLayoutHints(kLHintsTop | kLHintsCenterX,  2, 2, 2, 2);
  fPlotVarFrmLt     = new TGLayoutHints(kLHintsTop | kLHintsCenterX,  2, 2, 2, 2);
  fInteractionFrmLt = new TGLayoutHints(kLHintsTop | kLHintsCenterX,  2, 2, 2, 2);
  fCutsSumsFrmLt    = new TGLayoutHints(kLHintsTop | kLHintsCenterX,  2, 2, 2, 2);
  fBtnFrmLt         = new TGLayoutHints(kLHintsTop | kLHintsCenterX,  2, 2, 2, 2);
  fBtnFrmLt         = new TGLayoutHints(kLHintsTop | kLHintsCenterX,  2, 2, 2, 2);

  // Add Frames to Main Frame

  fMain -> AddFrame( fPlotTypeGrpf,    fPlotTypeFrmLt    );
  fMain -> AddFrame( fTotXSecConfGrpf, fTotXSecConfFrmLt );
  fMain -> AddFrame( fDifXSecConfGrpf, fDifXSecConfFrmLt );
  fMain -> AddFrame( fSFConfGrpf,      fSFConfFrmLt      );
  fMain -> AddFrame( fPlotVarGrpf,     fPlotVarFrmLt     );
  fMain -> AddFrame( fInteractionGrpf, fInteractionFrmLt );
  fMain -> AddFrame( fCutsSumsHFrm,    fCutsSumsFrmLt    );
  fMain -> AddFrame( fBtnHFrm,         fBtnFrmLt         );

  fMain->MapSubwindows();
  fMain->Resize();

  this->PositionRelativeToParent(main); // position relative to the parent's window

  if (fHaveNoHistory) this->Defaults();        // load some default values, or
  else                this->LoadLastEntries(); // last known set of user inputs

  fHaveNoHistory = false;

  fMain->SetWindowName("NeuGEN Inputs Dialog");

  fMain->MapWindow();
}
//______________________________________________________________________________
NeuGenInputDialog::~NeuGenInputDialog()
{
//-- note: the order is significant

   delete fPlotTypeLb;
   delete fPlotTypeCbx;
   delete fNpLb;
   delete fNPointsNmE;
   delete fFluxLb;
   delete fFluxCbx;
   delete fPlotTypeHFrm;

   delete fMinEnergyNmE;
   delete fMaxEnergyNmE;
   delete fMinEnergyLb;
   delete fMaxEnergyLb;
   delete fTotXSecSpacer;
   delete fTotXSecConfHFrm;

   delete fEnergyNmE;
   delete fEnergyLb;
   delete fPlotRangeLb;
   delete fPlotRangeCbx;

   delete fSFFixVarNmE;
   delete fSFFixVarLb;
   delete fSFLb;
   delete fSFTypeCbx;
   delete fSFRawDisCkb;

   delete fPlotVarLb;
   delete fMinVarLb;
   delete fMaxVarLb;
   delete fPlotVarCbx;
   delete fMinVarNmE;
   delete fMaxVarNmE;

   delete fProbeCbx;
   delete fWkCurrCbx;
   delete fInitStateCbx;
   delete fProbeLb;
   delete fWkCurrLb;
   delete fInitStateLb;
   delete fInteractionSpacer;
   delete fFinStateLb;
   delete fALb;
   delete fFinStateLbx;
   delete fANmE;
   delete fAllFinStatesCkb;
   delete fInteractionLFrm;
   delete fInteractionCFrm;
   delete fInteractionRFrm;
   delete fInteractionHFrm;

   delete fMinCutNmE;
   delete fMaxCutNmE;
   delete fMinCutLb;
   delete fMaxCutLb;
   delete fCutVarCbx;
   delete fCutsUFrm;
   delete fCutsLFrm;
   delete fCutsVFrm;

   delete fSumSpacer;
   delete fQelBitMaskCkb;
   delete fResBitMaskCkb;
   delete fDisBitMaskCkb;

   delete fBtnSpacer;
   delete fOkBtn;
   delete fCancelBtn;
   delete fResetBtn;

   delete fPlotTypeFrmLt;
   delete fTotXSecConfFrmLt;
   delete fDifXSecConfFrmLt;
   delete fSFConfFrmLt;
   delete fPlotVarFrmLt;
   delete fInteractionFrmLt;
   delete fCutsSumsFrmLt;
   delete fBtnFrmLt;

   delete fPlotTypeGrpf;
   delete fTotXSecConfGrpf;
   delete fDifXSecConfGrpf;
   delete fSFConfGrpf;
   delete fPlotVarGrpf;
   delete fInteractionGrpf;
   delete fCutsGrpf;
   delete fSumsGrpf;
   delete fBtnHFrm;

   delete fMain;
}
//______________________________________________________________________________
TGGroupFrame * NeuGenInputDialog::BuildXSecTypeFrame(void)
{
  TGGroupFrame * grpf = new TGGroupFrame(fMain, "plot type", kVerticalFrame);

  fPlotTypeHFrm = new TGHorizontalFrame(grpf, 10, 10);

  //-- plot type
  fPlotTypeLb = new TGLabel(fPlotTypeHFrm, new TGString( "quantity: "));

  fPlotTypeCbx = new TGComboBox(fPlotTypeHFrm, 97);
  utils::gui::FillComboBox( fPlotTypeCbx,  k_neugen_plot_type );
  fPlotTypeCbx->Resize(100, 20);

  //-- # of points
  fNpLb = new TGLabel(fPlotTypeHFrm, new TGString( "  npoints: "));
  fNPointsNmE = new TGNumberEntry(fPlotTypeHFrm, 0, 4, 1, TGNumberFormat::kNESInteger);

  //-- flux to scale with
  fFluxLb = new TGLabel(fPlotTypeHFrm, new TGString( "  scaling flux: "));

  fFluxCbx = new TGComboBox(fPlotTypeHFrm, 98);
  utils::gui::FillComboBox( fFluxCbx,  k_scaling_flux );
  fFluxCbx->Resize(60, 20);

  //-- add widgets to horizontal frame
  fPlotTypeHFrm -> AddFrame ( fPlotTypeLb  );
  fPlotTypeHFrm -> AddFrame ( fPlotTypeCbx );
  fPlotTypeHFrm -> AddFrame ( fNpLb        );
  fPlotTypeHFrm -> AddFrame ( fNPointsNmE  );
  fPlotTypeHFrm -> AddFrame ( fFluxLb      );
  fPlotTypeHFrm -> AddFrame ( fFluxCbx     );

  //-- add horizontal frame to group frame
  grpf -> AddFrame ( fPlotTypeHFrm );

  return grpf;
}
//______________________________________________________________________________
TGGroupFrame * NeuGenInputDialog::BuildConfigTotalXSecFrame(void)
{
  TGGroupFrame * grpf = new TGGroupFrame(fMain, "total xsec config", kVerticalFrame);

  fTotXSecConfHFrm = new TGHorizontalFrame(grpf, 10, 10);

  fMinEnergyNmE = new TGNumberEntry(fTotXSecConfHFrm, 0., 6, 1, TGNumberFormat::kNESReal);
  fMaxEnergyNmE = new TGNumberEntry(fTotXSecConfHFrm, 0., 6, 1, TGNumberFormat::kNESReal);

  fMinEnergyLb = new TGLabel(fTotXSecConfHFrm, new TGString( "   Emin:  "));
  fMaxEnergyLb = new TGLabel(fTotXSecConfHFrm, new TGString( "   Emax:  "));

  fTotXSecSpacer = new TGLabel(fTotXSecConfHFrm, new TGString( "   "));

  fTotXSecConfHFrm -> AddFrame ( fMinEnergyLb   );
  fTotXSecConfHFrm -> AddFrame ( fMinEnergyNmE  );
  fTotXSecConfHFrm -> AddFrame ( fMaxEnergyLb   );
  fTotXSecConfHFrm -> AddFrame ( fMaxEnergyNmE  );
  fTotXSecConfHFrm -> AddFrame ( fTotXSecSpacer );

  grpf->AddFrame(fTotXSecConfHFrm);

  return grpf;
}
//______________________________________________________________________________
TGGroupFrame * NeuGenInputDialog::BuildConfigDiffXSecFrame(void)
{
  UInt_t kh = kHorizontalFrame;
  TGNumberFormat::EStyle rstyle = TGNumberFormat::kNESReal;

  TGGroupFrame * grpf = new TGGroupFrame(fMain, "differential xsec config", kh);

  fEnergyLb     = new TGLabel       (grpf, new TGString( "     E:  "));
  fEnergyNmE    = new TGNumberEntry (grpf, 0., 7, 3, rstyle);
  fPlotRangeLb  = new TGLabel       (grpf, new TGString( "  range: "));
  fPlotRangeCbx = new TGComboBox    (grpf, 100);

  utils::gui::FillComboBox( fPlotRangeCbx,  k_neugen_plot_range_option );
  fPlotRangeCbx->Resize(80, 20);

  grpf -> AddFrame ( fEnergyLb     );
  grpf -> AddFrame ( fEnergyNmE    );
  grpf -> AddFrame ( fPlotRangeLb  );
  grpf -> AddFrame ( fPlotRangeCbx );

  return grpf;
}
//______________________________________________________________________________
TGGroupFrame * NeuGenInputDialog::BuildConfigSFFrame(void)
{
  UInt_t kh = kHorizontalFrame;

  TGNumberFormat::EStyle rstyle = TGNumberFormat::kNESReal;

  TGGroupFrame * grpf = new TGGroupFrame(fMain,"structure function config", kh);

  fSFRawDisCkb = new TGCheckButton (grpf, "Raw DIS - ", 133);
  fSFFixVarLb  = new TGLabel       (grpf, new TGString( "  fixed var:    "));
  fSFFixVarNmE = new TGNumberEntry (grpf, 0., 9, 3, rstyle);
  fSFLb        = new TGLabel       (grpf, new TGString( "  SF: "));
  fSFTypeCbx   = new TGComboBox    (grpf, 197);

  utils::gui::FillComboBox( fSFTypeCbx,  k_neugen_sf );
  fSFTypeCbx->Resize(50, 20);

  grpf -> AddFrame ( fSFRawDisCkb );
  grpf -> AddFrame ( fSFFixVarLb  );
  grpf -> AddFrame ( fSFFixVarNmE );
  grpf -> AddFrame ( fSFLb        );
  grpf -> AddFrame ( fSFTypeCbx   );

  return grpf;
}
//______________________________________________________________________________
TGGroupFrame * NeuGenInputDialog::BuildPlotVarFrame(void)
{
  UInt_t kh = kHorizontalFrame;
  TGNumberFormat::EStyle rstyle = TGNumberFormat::kNESReal;

  TGGroupFrame * grpf = new TGGroupFrame(fMain,
                           "kinematic var: x axis for Diff-Xsec/SF plots", kh);

  fPlotVarLb  = new TGLabel       (grpf, new TGString(" plot variable:  "));
  fPlotVarCbx = new TGComboBox    (grpf, 99);
  fMinVarLb   = new TGLabel       (grpf, new TGString( "  min: "));
  fMaxVarLb   = new TGLabel       (grpf, new TGString( "  max: "));
  fMinVarNmE  = new TGNumberEntry (grpf, 0., 6, 1, rstyle);
  fMaxVarNmE  = new TGNumberEntry (grpf, 0., 6, 1, rstyle);

  utils::gui::FillComboBox( fPlotVarCbx,  k_neugen_plot_variable );
  fPlotVarCbx->Resize(100, 20);

  grpf -> AddFrame ( fPlotVarLb    );
  grpf -> AddFrame ( fPlotVarCbx   );
  grpf -> AddFrame ( fMinVarLb     );
  grpf -> AddFrame ( fMinVarNmE    );
  grpf -> AddFrame ( fMaxVarLb     );
  grpf -> AddFrame ( fMaxVarNmE    );

  return grpf;
}
//______________________________________________________________________________
TGGroupFrame * NeuGenInputDialog::BuildInteractionFrame(void)
{
  UInt_t kv = kVerticalFrame;
  TGNumberFormat::EStyle istyle = TGNumberFormat::kNESInteger;

  TGGroupFrame * grpf = new TGGroupFrame(fMain, "interaction", kv);

  //-- GroupFrame contains HorizontalFrame which contains 3 CompositeFrames

  fInteractionHFrm = new TGHorizontalFrame(grpf, 10, 10);

  fInteractionLFrm = new TGCompositeFrame(fInteractionHFrm, 1, 1, kv);
  fInteractionCFrm = new TGCompositeFrame(fInteractionHFrm, 1, 1, kv);
  fInteractionRFrm = new TGCompositeFrame(fInteractionHFrm, 1, 1, kv);

  //-- left composite frame

  fProbeLb     = new TGLabel(fInteractionLFrm, new TGString(" Neutrino:")    );
  fWkCurrLb    = new TGLabel(fInteractionLFrm, new TGString(" Weak current:"));
  fInitStateLb = new TGLabel(fInteractionLFrm, new TGString(" Init. state:") );

  fProbeCbx     = new TGComboBox (fInteractionLFrm, 11);
  fWkCurrCbx    = new TGComboBox (fInteractionLFrm, 12);
  fInitStateCbx = new TGComboBox (fInteractionLFrm, 13);

  utils::gui::FillComboBox ( fProbeCbx,     k_neugen_nu        );
  utils::gui::FillComboBox ( fWkCurrCbx,    k_neugen_wcurr     );
  utils::gui::FillComboBox ( fInitStateCbx, k_neugen_initstate );

  fProbeCbx     -> Resize(150, 22);
  fWkCurrCbx    -> Resize(150, 22);
  fInitStateCbx -> Resize(150, 22);

  fALb  = new TGLabel       (fInteractionLFrm, new TGString( "  Nucl.target mass number "));
  fANmE = new TGNumberEntry (fInteractionLFrm, 0, 10, 1, istyle);

  fInteractionLFrm -> AddFrame ( fProbeLb      );
  fInteractionLFrm -> AddFrame ( fProbeCbx     );
  fInteractionLFrm -> AddFrame ( fWkCurrLb     );
  fInteractionLFrm -> AddFrame ( fWkCurrCbx    );
  fInteractionLFrm -> AddFrame ( fInitStateLb  );
  fInteractionLFrm -> AddFrame ( fInitStateCbx );
  fInteractionLFrm -> AddFrame ( fALb          );
  fInteractionLFrm -> AddFrame ( fANmE         );

  //-- central frame (spacer)

  fInteractionSpacer = new TGLabel(fInteractionCFrm, new TGString("    "));

  fInteractionCFrm->AddFrame(fInteractionSpacer);

  //-- right composite frame

  fFinStateLb  = new TGLabel(fInteractionRFrm, new TGString( " Final state:"));
  fFinStateLbx = new TGListBox(fInteractionRFrm, 21);

  utils::gui::FillListBox(fFinStateLbx,  k_neugen_finstate);

  fFinStateLbx->Resize(165, 90);

  fAllFinStatesCkb = new TGCheckButton(fInteractionRFrm, "INCLUSIVE",  40);

  fInteractionRFrm -> AddFrame ( fAllFinStatesCkb    );
  fInteractionRFrm -> AddFrame ( fInteractionSpacer  );
  fInteractionRFrm -> AddFrame ( fFinStateLb         );
  fInteractionRFrm -> AddFrame ( fFinStateLbx        );

  //-- add left/right composite frame to main frame

  fInteractionHFrm -> AddFrame ( fInteractionLFrm );
  fInteractionHFrm -> AddFrame ( fInteractionCFrm );
  fInteractionHFrm -> AddFrame ( fInteractionRFrm );

  grpf->AddFrame(fInteractionHFrm);

  return grpf;
}
//______________________________________________________________________________
TGGroupFrame * NeuGenInputDialog::BuildCutsFrame(void)
{
  TGGroupFrame * grpf = new TGGroupFrame(fCutsSumsHFrm, "cut variable", kVerticalFrame);

  fCutsVFrm = new TGVerticalFrame(grpf, 10, 10);

  fCutsUFrm = new TGCompositeFrame(fCutsVFrm, 1, 1, kVerticalFrame);
  fCutsLFrm = new TGCompositeFrame(fCutsVFrm, 1, 1, kVerticalFrame);

  //-- upper composite frame

  fCutVarCbx = new TGComboBox(fCutsUFrm, 21);

  utils::gui::FillComboBox( fCutVarCbx,  k_neugen_cut_variable );

  fCutVarCbx -> Resize(175, 20);

  fCutsUFrm -> AddFrame ( fCutVarCbx  );

  //-- lower composite frame

  TGMatrixLayout * cuts_nentries_matrix_layout = new TGMatrixLayout(fCutsLFrm, 0, 2, 1);
  fCutsLFrm->SetLayoutManager( cuts_nentries_matrix_layout );

  fMinCutNmE = new TGNumberEntry(fCutsLFrm, 0.000, 9, 3, TGNumberFormat::kNESReal);
  fMaxCutNmE = new TGNumberEntry(fCutsLFrm, 0.000, 9, 3, TGNumberFormat::kNESReal);

  fMinCutLb = new TGLabel(fCutsLFrm, new TGString( "min:"));
  fMaxCutLb = new TGLabel(fCutsLFrm, new TGString( "max:"));

  fCutsLFrm -> AddFrame ( fMinCutLb  );
  fCutsLFrm -> AddFrame ( fMinCutNmE );
  fCutsLFrm -> AddFrame ( fMaxCutLb  );
  fCutsLFrm -> AddFrame ( fMaxCutNmE );

  //-- add upper/lower composite frames to main group frame

  fCutsVFrm -> AddFrame ( fCutsUFrm );
  fCutsVFrm -> AddFrame ( fCutsLFrm );

  grpf->AddFrame(fCutsVFrm);

  return grpf;
}
//______________________________________________________________________________
TGGroupFrame * NeuGenInputDialog::BuildSumsFrame(void)
{
  TGGroupFrame * grpf = new TGGroupFrame(fCutsSumsHFrm, "sums", kVerticalFrame);

  fSumSpacer = new TGLabel(grpf, new TGString("     "));

  fQelBitMaskCkb = new TGCheckButton(grpf, "Include QEL",  31);
  fResBitMaskCkb = new TGCheckButton(grpf, "Include RES",  32);
  fDisBitMaskCkb = new TGCheckButton(grpf, "Include DIS",  33);

  grpf -> AddFrame ( fSumSpacer  );
  grpf -> AddFrame ( fQelBitMaskCkb );
  grpf -> AddFrame ( fResBitMaskCkb );
  grpf -> AddFrame ( fDisBitMaskCkb );

  return grpf;
}
//______________________________________________________________________________
TGHorizontalFrame * NeuGenInputDialog::BuildButtonsFrame(void)
{
  TGHorizontalFrame * hframe = new TGHorizontalFrame(fMain, 10, 10);

  fOkBtn     = new TGTextButton(hframe, "   &Ok   ",  1);
  fResetBtn  = new TGTextButton(hframe, "  &Reset ",  2);
  fCancelBtn = new TGTextButton(hframe, " &Cancel ",  3);

  fBtnSpacer = new TGLabel(hframe, new TGString("     "));

  fOkBtn     -> Connect("Clicked()",
                              "genie::nuvld::NeuGenInputDialog", this, "OK()");
  fResetBtn  -> Connect("Clicked()",
                           "genie::nuvld::NeuGenInputDialog", this, "Reset()");
  fCancelBtn -> Connect("Clicked()",
                          "genie::nuvld::NeuGenInputDialog", this, "Cancel()");

  hframe -> AddFrame( fOkBtn     );
  hframe -> AddFrame( fBtnSpacer );
  hframe -> AddFrame( fResetBtn  );
  hframe -> AddFrame( fBtnSpacer );
  hframe -> AddFrame( fCancelBtn );

  return hframe;
}
//______________________________________________________________________________
void NeuGenInputDialog::OK(void)
{
  NeuGenCards * cards = NeuGenCards::Instance();

  cards->CurrInputs()->SetNBins        ( this->ReadNPoints()      );
  cards->CurrInputs()->SetPlotType     ( this->ReadPlotType()     );
  cards->CurrInputs()->SetEmin         ( this->ReadEnergyMin()    );
  cards->CurrInputs()->SetEmax         ( this->ReadEnergyMax()    );
  cards->CurrInputs()->SetE            ( this->ReadEnergy()       );
  cards->CurrInputs()->SetPlotVar      ( this->ReadPlotVar()      );
  cards->CurrInputs()->SetFlux         ( this->ReadScalingFlux()  );
  cards->CurrInputs()->SetRange        ( this->ReadPlotRange()    );
  cards->CurrInputs()->SetPlotVarMin   ( this->ReadPlotVarMin()   );
  cards->CurrInputs()->SetPlotVarMax   ( this->ReadPlotVarMax()   );
  cards->CurrInputs()->SetNeutrino     ( this->ReadNeutrino()     );
  cards->CurrInputs()->SetWkCurrent    ( this->ReadWkCurrent()    );
  cards->CurrInputs()->SetFinalState   ( this->ReadFinalState()   );
  cards->CurrInputs()->SetInitialState ( this->ReadInitialState() );
  cards->CurrInputs()->SetCutVar       ( this->ReadCutVar()       );
  cards->CurrInputs()->SetCutVarMin    ( this->ReadCutVarMin()    );
  cards->CurrInputs()->SetCutVarMax    ( this->ReadCutVarMax()    );
  cards->CurrInputs()->SetIncludeQel   ( this->ReadQelBitInMask() );
  cards->CurrInputs()->SetIncludeRes   ( this->ReadResBitInMask() );
  cards->CurrInputs()->SetIncludeDis   ( this->ReadDisBitInMask() );
  cards->CurrInputs()->SetInclusive    ( this->ReadInclusive()    );
  cards->CurrInputs()->SetSFRawDis     ( this->ReadSFRawDis()     );
  cards->CurrInputs()->SetSFFixedVar   ( this->ReadSFFixedVar()   );
  cards->CurrInputs()->SetSF           ( this->ReadSF()           );

  this->Report(); // write out user selections to the GUI

  //cout << *(cards->CurrInputs());

  fMain->SendCloseMessage();
}
//______________________________________________________________________________
void NeuGenInputDialog::Reset(void)
{
// Reset dialog entries

  this->Defaults();
}
//______________________________________________________________________________
void NeuGenInputDialog::Defaults(void)
{
// Set some "default" entries

  fNPointsNmE   -> SetIntNumber (100); //! number of points
  fPlotTypeCbx  -> Select(0);          //! xsec type = 'total'
  fPlotRangeCbx -> Select(0);          //! plot range = 'automatic'

  fEnergyNmE    -> SetNumber (   0.0 );
  fMinEnergyNmE -> SetNumber (   0.1 );
  fMaxEnergyNmE -> SetNumber ( 100.0 );
  fMinVarNmE    -> SetNumber (   0.0 );
  fMaxVarNmE    -> SetNumber (   0.0 );
  fMinCutNmE    -> SetNumber (   0.0 );
  fMaxCutNmE    -> SetNumber (   0.0 );

  fProbeCbx     -> Select(0);
  fWkCurrCbx    -> Select(0);
  fInitStateCbx -> Select(0);
  fFinStateLbx  -> Select(0);
  fANmE         -> SetNumber (1);

  fFluxCbx      -> Select (0);  // scaling flux  = 'none'
  fPlotVarCbx   -> Select (0);  // plot variable = 'none'
  fCutVarCbx    -> Select (0);  // cut variable  = 'none'

  // sums - all check buttons = ON
  fQelBitMaskCkb -> SetOn (kTRUE );
  fResBitMaskCkb -> SetOn (kTRUE );
  fDisBitMaskCkb -> SetOn (kTRUE );

  // SF defaults
  fSFTypeCbx    -> Select (0);
  fSFFixVarNmE  -> SetNumber ( 0.0 );
}
//______________________________________________________________________________
void NeuGenInputDialog::LoadLastEntries(void)
{
// Load the dialog widgets with the last known entries...
// This was proven to be the best dialog initialization tactic, since users
// usually make successive, small changes/corrections after they enter their
// initial input values.

  NeuGenCards * cards = NeuGenCards::Instance();

  NeuGenInputs * inp = cards->CurrInputs();

  fNPointsNmE -> SetIntNumber ( inp->NBins() );

  fPlotTypeCbx->Select(
          utils::gui::ComboBoxSelectionId(k_neugen_plot_type,
                                           inp->PlotTypeString().c_str()) );
  fPlotRangeCbx->Select(
          utils::gui::ComboBoxSelectionId(k_neugen_plot_range_option,
                                          inp->PlotRangeString().c_str()) );
  fProbeCbx->Select(
          utils::gui::ComboBoxSelectionId(k_neugen_nu,
                                              inp->ProbeString().c_str()) );
  fWkCurrCbx->Select(
          utils::gui::ComboBoxSelectionId(k_neugen_wcurr,
                                          inp->WkCurrentString().c_str()) );
  fInitStateCbx->Select(
          utils::gui::ComboBoxSelectionId(k_neugen_initstate,
                                          inp->InitStateString().c_str()) );
  fFinStateLbx->Select(
          utils::gui::ListBoxSelectionId(k_neugen_finstate,
                                         inp->FinalStateString().c_str()) );
  fFluxCbx->Select(
          utils::gui::ComboBoxSelectionId(k_scaling_flux,
                                             inp->FluxIdString().c_str()) );
  fPlotVarCbx->Select(
          utils::gui::ComboBoxSelectionId(k_neugen_plot_variable,
                                            inp->PlotVarString().c_str()) );
  fCutVarCbx->Select(
          utils::gui::ComboBoxSelectionId(k_neugen_cut_variable,
                                             inp->CutVarString().c_str()) );
  fSFTypeCbx->Select(
          utils::gui::ComboBoxSelectionId(k_neugen_sf,
                                                 inp->SFString().c_str()) );

  fEnergyNmE    -> SetNumber ( inp->Energy()       );
  fMinEnergyNmE -> SetNumber ( inp->EnergyMin()    );
  fMaxEnergyNmE -> SetNumber ( inp->EnergyMax()    );
  fMinVarNmE    -> SetNumber ( inp->PlotVarMin()   );
  fMaxVarNmE    -> SetNumber ( inp->PlotVarMax()   );
  fMinCutNmE    -> SetNumber ( inp->CutVarMin()    );
  fMaxCutNmE    -> SetNumber ( inp->CutVarMax()    );
  fSFFixVarNmE  -> SetNumber ( inp->SFFixedVar()   );
  fANmE         -> SetNumber ( inp->A()            );

  fQelBitMaskCkb   -> SetOn ( inp->IncludeQel() );
  fResBitMaskCkb   -> SetOn ( inp->IncludeRes() );
  fDisBitMaskCkb   -> SetOn ( inp->IncludeDis() );
  fAllFinStatesCkb -> SetOn ( inp->Inclusive()  );
  fSFRawDisCkb     -> SetOn ( inp->SFRawDis()   );
}
//______________________________________________________________________________
void NeuGenInputDialog::PositionRelativeToParent(const TGWindow * main)
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
int NeuGenInputDialog::ReadNPoints(void)
{
  return fNPointsNmE->GetIntNumber();
}
//______________________________________________________________________________
bool NeuGenInputDialog::ReadQelBitInMask(void)
{
  return (fQelBitMaskCkb->GetState() == kButtonDown);
}
//______________________________________________________________________________
bool NeuGenInputDialog::ReadResBitInMask(void)
{
  return (fResBitMaskCkb->GetState() == kButtonDown);
}
//______________________________________________________________________________
bool NeuGenInputDialog::ReadDisBitInMask(void)
{
  return (fDisBitMaskCkb->GetState() == kButtonDown);
}
//______________________________________________________________________________
bool NeuGenInputDialog::ReadInclusive(void)
{
  return (fAllFinStatesCkb->GetState() == kButtonDown);
}
//______________________________________________________________________________
bool NeuGenInputDialog::ReadSFRawDis(void)
{
  return (fSFRawDisCkb->GetState() == kButtonDown);
}
//______________________________________________________________________________
float NeuGenInputDialog::ReadEnergy(void)
{
  return (float) fEnergyNmE->GetNumber();
}
//______________________________________________________________________________
float NeuGenInputDialog::ReadEnergyMin(void)
{
  return (float) fMinEnergyNmE->GetNumber();
}
//______________________________________________________________________________
float NeuGenInputDialog::ReadEnergyMax(void)
{
  return (float) fMaxEnergyNmE->GetNumber();
}
//______________________________________________________________________________
float NeuGenInputDialog::ReadPlotVarMin(void)
{
  return (float) fMinVarNmE->GetNumber();
}
//______________________________________________________________________________
float NeuGenInputDialog::ReadPlotVarMax(void)
{
  return (float) fMaxVarNmE->GetNumber();
}
//______________________________________________________________________________
float NeuGenInputDialog::ReadCutVarMin(void)
{
  return (float) fMinCutNmE->GetNumber();
}
//______________________________________________________________________________
float NeuGenInputDialog::ReadCutVarMax(void)
{
  return (float) fMaxCutNmE->GetNumber();
}
//______________________________________________________________________________
float NeuGenInputDialog::ReadSFFixedVar(void)
{
  return (float) fSFFixVarNmE->GetNumber();
}
//______________________________________________________________________________
float NeuGenInputDialog::ReadA(void)
{
  return (float) fANmE->GetNumber();
}
//______________________________________________________________________________
string NeuGenInputDialog::ReadFinalState(void)
{
  return utils::gui::ListBoxSelectionAsString(fFinStateLbx, k_neugen_finstate);
}
//______________________________________________________________________________
string NeuGenInputDialog::ReadInitialState(void)
{
  return utils::gui::ComboBoxSelectionAsString(fInitStateCbx, k_neugen_initstate);
}
//______________________________________________________________________________
string NeuGenInputDialog::ReadPlotType(void)
{
  return utils::gui::ComboBoxSelectionAsString(fPlotTypeCbx, k_neugen_plot_type);
}
//______________________________________________________________________________
string NeuGenInputDialog::ReadScalingFlux(void)
{
  return utils::gui::ComboBoxSelectionAsString(fFluxCbx, k_scaling_flux);
}
//______________________________________________________________________________
string NeuGenInputDialog::ReadPlotVar(void)
{
  return utils::gui::ComboBoxSelectionAsString(fPlotVarCbx, k_neugen_plot_variable);
}
//______________________________________________________________________________
string NeuGenInputDialog::ReadCutVar(void)
{
  return utils::gui::ComboBoxSelectionAsString(fCutVarCbx, k_neugen_cut_variable);
}
//______________________________________________________________________________
string NeuGenInputDialog::ReadPlotRange(void)
{
  return utils::gui::ComboBoxSelectionAsString(fPlotVarCbx, k_neugen_plot_range_option);
}
//______________________________________________________________________________
string NeuGenInputDialog::ReadNeutrino(void)
{
  return utils::gui::ComboBoxSelectionAsString(fProbeCbx, k_neugen_nu);
}
//______________________________________________________________________________
string NeuGenInputDialog::ReadWkCurrent(void)
{
  return utils::gui::ComboBoxSelectionAsString(fWkCurrCbx, k_neugen_wcurr);
}
//______________________________________________________________________________
string NeuGenInputDialog::ReadSF(void)
{
  return utils::gui::ComboBoxSelectionAsString(fSFTypeCbx, k_neugen_sf);
}
//______________________________________________________________________________
void NeuGenInputDialog::Report(void)
{
  SysLogSingleton * syslog = SysLogSingleton::Instance();

  syslog->Log()->AddLine("NeuGEN inputs:");

  syslog->Log()->AddLine( Concat(
                  "npoints............", this->ReadNPoints())              );
  syslog->Log()->AddLine( Concat(
                  "plot type..........", this->ReadPlotType().c_str())     );
  syslog->Log()->AddLine( Concat(
                  "Emin...............", this->ReadEnergyMin())            );
  syslog->Log()->AddLine( Concat(
                  "Emax...............", this->ReadEnergyMax())            );
  syslog->Log()->AddLine( Concat(
                  "E..................", this->ReadEnergy())               );
  syslog->Log()->AddLine( Concat(
                  "Plot variable......", this->ReadPlotVar().c_str())      );
  syslog->Log()->AddLine( Concat(
                  "Scaling flux.......", this->ReadScalingFlux().c_str())  );
  syslog->Log()->AddLine( Concat(
                  "Plot range.........", this->ReadPlotRange().c_str())    );
  syslog->Log()->AddLine( Concat(
                  "Plot var./ min.....", this->ReadPlotVarMin())           );
  syslog->Log()->AddLine( Concat(
                  "Plot var. - max....", this->ReadPlotVarMax())           );
  syslog->Log()->AddLine( Concat(
                  "Neutrino type......", this->ReadNeutrino().c_str())     );
  syslog->Log()->AddLine( Concat(
                  "Weak current.......", this->ReadWkCurrent().c_str())    );
  syslog->Log()->AddLine( Concat(
                  "Final state........", this->ReadFinalState().c_str())   );
  syslog->Log()->AddLine( Concat(
                  "Initial state......", this->ReadInitialState().c_str()) );
  syslog->Log()->AddLine( Concat(
                  "A..................", this->ReadA())                    );
  syslog->Log()->AddLine( Concat(
                  "Cut variable.......", this->ReadCutVar().c_str())       );
  syslog->Log()->AddLine( Concat(
                  "Cut var./ min......", this->ReadCutVarMin())            );
  syslog->Log()->AddLine( Concat(
                  "Cut var./ max......", this->ReadCutVarMax())            );
  syslog->Log()->AddLine( Concat(
                  "Inclusive XSec.....", this->ReadInclusive())            );
  syslog->Log()->AddLine( Concat(
                  "QEL bit/mask.......", this->ReadQelBitInMask())         );
  syslog->Log()->AddLine( Concat(
                  "RES bit/mask.......", this->ReadResBitInMask())         );
  syslog->Log()->AddLine( Concat(
                  "DIS bit/mask.......", this->ReadDisBitInMask())         );
  syslog->Log()->AddLine( Concat(
                  "SF type............", this->ReadSF().c_str())           );
  syslog->Log()->AddLine( Concat(
                  "SF fixed var.......", this->ReadSFFixedVar())           );
  syslog->Log()->AddLine( Concat(
                  "SF Raw DIS.........", this->ReadSFRawDis())             );
}
//______________________________________________________________________________

