//_____________________________________________________________________________
/*!

\class    genie::nuvld::NeuGenInputDialog

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#include <iostream>

#include <TGraph.h>
#include <TCanvas.h>
#include <TRootEmbeddedCanvas.h>

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
using namespace genie::string_utils;
using namespace genie::nuvld::facades;

ClassImp(NeuGenInputDialog)

//______________________________________________________________________________
bool NeuGenInputDialog::_have_no_history = true;
//______________________________________________________________________________
NeuGenInputDialog::NeuGenInputDialog(const TGWindow * p,
                      const TGWindow * main, UInt_t w, UInt_t h, UInt_t options)
{
  _main = new TGTransientFrame(p, main, w, h, options);
  _main->Connect("CloseWindow()",
                      "genie::nuvld::NeuGenInputDialog", this, "CloseWindow()");

  // Build frames
  
  _xsec_type_grpf     = this->BuildXSecTypeFrame();
  _tot_xsec_conf_grpf = this->BuildConfigTotalXSecFrame();
  _dis_xsec_conf_grpf = this->BuildConfigDiffXSecFrame();
  _interaction_grpf   = this->BuildInteractionFrame();
  
  _cuts_sums_hframe   = new TGHorizontalFrame(_main, 10, 10);
  
  _cuts_grpf = this->BuildCutsFrame();
  _sums_grpf = this->BuildSumsFrame();

  _cuts_sums_hframe -> AddFrame( _cuts_grpf );
  _cuts_sums_hframe -> AddFrame( _sums_grpf );

  _button_frame = this->BuildButtonsFrame();

  // Build frame layouts

  _xsec_type_frame_layout     = new TGLayoutHints(kLHintsTop | kLHintsCenterX,  2, 2, 2, 2);
  _tot_xsec_conf_frame_layout = new TGLayoutHints(kLHintsTop | kLHintsCenterX,  2, 2, 2, 2);
  _dis_xsec_conf_frame_layout = new TGLayoutHints(kLHintsTop | kLHintsCenterX,  2, 2, 2, 2);
  _interaction_frame_layout   = new TGLayoutHints(kLHintsTop | kLHintsCenterX,  2, 2, 2, 2);
  _cuts_sums_frame_layout     = new TGLayoutHints(kLHintsTop | kLHintsCenterX,  2, 2, 2, 2);
  _button_frame_layout        = new TGLayoutHints(kLHintsTop | kLHintsCenterX,  2, 2, 2, 2);  
  _button_frame_layout        = new TGLayoutHints(kLHintsTop | kLHintsCenterX,  2, 2, 2, 2);

  // Add Frames to Main Frame

  _main -> AddFrame( _xsec_type_grpf,     _xsec_type_frame_layout     );
  _main -> AddFrame( _tot_xsec_conf_grpf, _tot_xsec_conf_frame_layout );
  _main -> AddFrame( _dis_xsec_conf_grpf, _dis_xsec_conf_frame_layout );
  _main -> AddFrame( _interaction_grpf,   _interaction_frame_layout   );
  _main -> AddFrame( _cuts_sums_hframe,   _cuts_sums_frame_layout     );
  _main -> AddFrame( _button_frame,       _button_frame_layout        );

  _main->MapSubwindows();
  _main->Resize();

  this->PositionRelativeToParent(main); // position relative to the parent's window

  if (_have_no_history) this->Defaults();        // load some default values, or
  else                  this->LoadLastEntries(); // last known set of user inputs
  
  _have_no_history = false;
  
  _main->SetWindowName("NeuGEN Inputs Dialog");

  _main->MapWindow();
}
//______________________________________________________________________________
NeuGenInputDialog::~NeuGenInputDialog()
{
//-- note: the order is significant

   delete _xsec_type_label;
   delete _xsec_type_combox;
   delete _npoints_label;
   delete _npoints;
   delete _flux_label;
   delete _flux_combox;
   delete _xsec_type_hframe;
   
   delete _E_min;
   delete _E_max;
   delete _E_min_label;
   delete _E_max_label;
   delete _tot_xsec_spacer;
   delete _tot_xsec_conf_hframe;

   delete _E;
   delete _E_label;
   delete _plot_var_label;
   delete _plot_var_combox;
   delete _var_min;
   delete _var_max;
   delete _plot_range_label;
   delete _plot_range_combox;
   delete _var_min_label;
   delete _var_max_label;
   delete _diff_xsec_conf_spacer;
   delete _dif_xsec_conf_uframe;
   delete _dif_xsec_conf_cframe;
   delete _dif_xsec_conf_lframe;
   delete _dif_xsec_conf_vframe;

   delete _neutrino_combo;
   delete _wcurrent_combo;
   delete _initstate_combo;
   delete _neutrino_label;
   delete _wcurrent_label;
   delete _initstate_label;
   delete _int_spacer_label;   
   delete _finstate_label;
   delete _finstate_listbox;
   delete _all_finstates;
   delete _int_lframe;
   delete _int_cframe;
   delete _int_rframe;
   delete _int_hframe;
   
   delete _cut_min;
   delete _cut_max;
   delete _cut_min_label;
   delete _cut_max_label;
   delete _cut_variable_combox;
   delete _cuts_uframe;
   delete _cuts_lframe;
   delete _cuts_vframe;

   delete _sum_spacer;
   delete _qel_bit_mask;
   delete _res_bit_mask;
   delete _dis_bit_mask;

   delete _button_spacer;
   delete _ok_button;
   delete _cancel_button;
   delete _reset_button;

   delete _xsec_type_frame_layout;
   delete _tot_xsec_conf_frame_layout;
   delete _dis_xsec_conf_frame_layout;
   delete _interaction_frame_layout;
   delete _cuts_sums_frame_layout;
   delete _button_frame_layout;

   delete _xsec_type_grpf;
   delete _tot_xsec_conf_grpf;
   delete _dis_xsec_conf_grpf;
   delete _interaction_grpf;
   delete _cuts_grpf;
   delete _sums_grpf;
   delete _button_frame;

   delete _main;
}
//______________________________________________________________________________
TGGroupFrame * NeuGenInputDialog::BuildXSecTypeFrame(void)
{
  TGGroupFrame * grpf = new TGGroupFrame(_main, "cross section", kVerticalFrame);

  _xsec_type_hframe = new TGHorizontalFrame(grpf, 10, 10);

  //-- xsec type
  _xsec_type_label = new TGLabel(_xsec_type_hframe, new TGString( "type: "));

  _xsec_type_combox = new TGComboBox(_xsec_type_hframe, 97);
  gui_utils::FillComboBox( _xsec_type_combox,  k_neugen_xsec_type );
  _xsec_type_combox->Resize(80, 20);

  //-- # of points
  _npoints_label = new TGLabel(_xsec_type_hframe, new TGString( "  npoints: "));
  _npoints = new TGNumberEntry(_xsec_type_hframe, 0, 4, 1, TGNumberFormat::kNESInteger);

  //-- flux to scale with
  _flux_label = new TGLabel(_xsec_type_hframe, new TGString( "  scaling flux: "));

  _flux_combox = new TGComboBox(_xsec_type_hframe, 98);
  gui_utils::FillComboBox( _flux_combox,  k_scaling_flux );
  _flux_combox->Resize(60, 20);

  //-- add widgets to horizontal frame
  _xsec_type_hframe -> AddFrame ( _xsec_type_label  );
  _xsec_type_hframe -> AddFrame ( _xsec_type_combox );
  _xsec_type_hframe -> AddFrame ( _npoints_label    );
  _xsec_type_hframe -> AddFrame ( _npoints          );
  _xsec_type_hframe -> AddFrame ( _flux_label       );
  _xsec_type_hframe -> AddFrame ( _flux_combox      );

  //-- add horizontal frame to group frame
  grpf -> AddFrame ( _xsec_type_hframe );

  return grpf;
}
//______________________________________________________________________________
TGGroupFrame * NeuGenInputDialog::BuildConfigTotalXSecFrame(void)
{
  TGGroupFrame * grpf = new TGGroupFrame(_main, "total xsec config", kVerticalFrame);

  _tot_xsec_conf_hframe = new TGHorizontalFrame(grpf, 10, 10);

  _E_min = new TGNumberEntry(_tot_xsec_conf_hframe, 0., 6, 1, TGNumberFormat::kNESReal);
  _E_max = new TGNumberEntry(_tot_xsec_conf_hframe, 0., 6, 1, TGNumberFormat::kNESReal);

  _E_min_label = new TGLabel(_tot_xsec_conf_hframe, new TGString( "   Emin:  "));
  _E_max_label = new TGLabel(_tot_xsec_conf_hframe, new TGString( "   Emax:  "));

  _tot_xsec_spacer = new TGLabel(_tot_xsec_conf_hframe, new TGString( "   "));

  _tot_xsec_conf_hframe -> AddFrame ( _E_min_label     );
  _tot_xsec_conf_hframe -> AddFrame ( _E_min           );
  _tot_xsec_conf_hframe -> AddFrame ( _E_max_label     );
  _tot_xsec_conf_hframe -> AddFrame ( _E_max           );
  _tot_xsec_conf_hframe -> AddFrame ( _tot_xsec_spacer );

  grpf->AddFrame(_tot_xsec_conf_hframe);

  return grpf;
}
//______________________________________________________________________________
TGGroupFrame * NeuGenInputDialog::BuildConfigDiffXSecFrame(void)
{
  TGGroupFrame * grpf = new TGGroupFrame(_main, "differential xsec config", kVerticalFrame);

  //-- main vertical frame + upper/central/lower composite frames
  _dif_xsec_conf_vframe = new TGVerticalFrame(grpf, 10, 10);

  _dif_xsec_conf_uframe = new TGCompositeFrame(_dif_xsec_conf_vframe, 1, 1, kHorizontalFrame);
  _dif_xsec_conf_cframe = new TGCompositeFrame(_dif_xsec_conf_vframe, 1, 1, kHorizontalFrame);
  _dif_xsec_conf_lframe = new TGCompositeFrame(_dif_xsec_conf_vframe, 1, 1, kHorizontalFrame);

  //-- upper frame

  _E_label = new TGLabel(_dif_xsec_conf_uframe, new TGString( "     E:  "));
  _E = new TGNumberEntry(_dif_xsec_conf_uframe, 0., 7, 3, TGNumberFormat::kNESReal);

  _plot_var_label = new TGLabel(_dif_xsec_conf_uframe, new TGString( "     plot variable:  "));
  _plot_var_combox = new TGComboBox(_dif_xsec_conf_uframe, 99);

  gui_utils::FillComboBox( _plot_var_combox,  k_neugen_plot_variable );
  _plot_var_combox->Resize(100, 20);

  _dif_xsec_conf_uframe -> AddFrame ( _E_label         );
  _dif_xsec_conf_uframe -> AddFrame ( _E               );
  _dif_xsec_conf_uframe -> AddFrame ( _plot_var_label  );
  _dif_xsec_conf_uframe -> AddFrame ( _plot_var_combox );

  //-- central frame (spacer)

  _diff_xsec_conf_spacer = new TGLabel(_dif_xsec_conf_cframe, new TGString(" "));

  _dif_xsec_conf_cframe->AddFrame( _diff_xsec_conf_spacer );

  //-- lower frame

  _var_min = new TGNumberEntry(_dif_xsec_conf_lframe, 0., 6, 1, TGNumberFormat::kNESReal);
  _var_max = new TGNumberEntry(_dif_xsec_conf_lframe, 0., 6, 1, TGNumberFormat::kNESReal);

  _plot_range_label = new TGLabel(_dif_xsec_conf_lframe, new TGString( "  range: "));
  _var_min_label    = new TGLabel(_dif_xsec_conf_lframe, new TGString( "  min: "));
  _var_max_label    = new TGLabel(_dif_xsec_conf_lframe, new TGString( "  max: "));

  _plot_range_combox = new TGComboBox(_dif_xsec_conf_lframe, 100);

  gui_utils::FillComboBox( _plot_range_combox,  k_neugen_plot_range_option );
  _plot_range_combox->Resize(80, 20);

  _dif_xsec_conf_lframe -> AddFrame ( _plot_range_label  );
  _dif_xsec_conf_lframe -> AddFrame ( _plot_range_combox );
  _dif_xsec_conf_lframe -> AddFrame ( _var_min_label     );
  _dif_xsec_conf_lframe -> AddFrame ( _var_min           );
  _dif_xsec_conf_lframe -> AddFrame ( _var_max_label     );
  _dif_xsec_conf_lframe -> AddFrame ( _var_max           );

  //-- add upper/lower composite frames to vertical frame

  _dif_xsec_conf_vframe -> AddFrame ( _dif_xsec_conf_uframe );
  _dif_xsec_conf_vframe -> AddFrame ( _dif_xsec_conf_cframe );
  _dif_xsec_conf_vframe -> AddFrame ( _dif_xsec_conf_lframe );

  grpf->AddFrame(_dif_xsec_conf_vframe);

  return grpf;
}
//______________________________________________________________________________
TGGroupFrame * NeuGenInputDialog::BuildInteractionFrame(void)
{
  TGGroupFrame * grpf = new TGGroupFrame(_main, "interaction", kVerticalFrame);

  //-- GroupFrame contains HorizontalFrame which contains 3 CompositeFrames

  _int_hframe = new TGHorizontalFrame(grpf, 10, 10);

  _int_lframe = new TGCompositeFrame(_int_hframe, 1, 1, kVerticalFrame);
  _int_cframe = new TGCompositeFrame(_int_hframe, 1, 1, kVerticalFrame);
  _int_rframe = new TGCompositeFrame(_int_hframe, 1, 1, kVerticalFrame);

  //-- left composite frame

  _neutrino_label  = new TGLabel(_int_lframe, new TGString(" Neutrino:")    );
  _wcurrent_label  = new TGLabel(_int_lframe, new TGString(" Weak current:"));
  _initstate_label = new TGLabel(_int_lframe, new TGString(" Init. state:") );

  _neutrino_combo  = new TGComboBox (_int_lframe, 11);
  _wcurrent_combo  = new TGComboBox (_int_lframe, 12);
  _initstate_combo = new TGComboBox (_int_lframe, 13);

  gui_utils::FillComboBox ( _neutrino_combo,  k_neugen_nu        );
  gui_utils::FillComboBox ( _wcurrent_combo,  k_neugen_wcurr     );
  gui_utils::FillComboBox ( _initstate_combo, k_neugen_initstate );

  _neutrino_combo  -> Resize(150, 22);
  _wcurrent_combo  -> Resize(150, 22);
  _initstate_combo -> Resize(150, 22);

  _int_lframe -> AddFrame ( _neutrino_label  );
  _int_lframe -> AddFrame ( _neutrino_combo  );
  _int_lframe -> AddFrame ( _wcurrent_label  );
  _int_lframe -> AddFrame ( _wcurrent_combo  );
  _int_lframe -> AddFrame ( _initstate_label );
  _int_lframe -> AddFrame ( _initstate_combo );

  //-- central frame (spacer)

  _int_spacer_label = new TGLabel(_int_cframe, new TGString("    "));

  _int_cframe->AddFrame(_int_spacer_label);

  //-- right composite frame

  _finstate_label = new TGLabel(_int_rframe, new TGString( " Final state:"));
  _finstate_listbox  = new TGListBox(_int_rframe, 21);

  gui_utils::FillListBox(_finstate_listbox,  k_neugen_finstate);

  _finstate_listbox->Resize(165, 78);

  _all_finstates = new TGCheckButton(_int_rframe, "INCLUSIVE Cross Section",  40);

  _int_rframe -> AddFrame ( _all_finstates    );
  _int_rframe -> AddFrame ( _int_spacer_label );
  _int_rframe -> AddFrame ( _finstate_label   );
  _int_rframe -> AddFrame ( _finstate_listbox );

  //-- add left/right composite frame to main frame

  _int_hframe -> AddFrame ( _int_lframe );
  _int_hframe -> AddFrame ( _int_cframe );
  _int_hframe -> AddFrame ( _int_rframe );

  grpf->AddFrame(_int_hframe);

  return grpf;
}
//______________________________________________________________________________
TGGroupFrame * NeuGenInputDialog::BuildCutsFrame(void)
{
  TGGroupFrame * grpf = new TGGroupFrame(_cuts_sums_hframe, "cut variable", kVerticalFrame);

  _cuts_vframe = new TGVerticalFrame(grpf, 10, 10);

  _cuts_uframe = new TGCompositeFrame(_cuts_vframe, 1, 1, kVerticalFrame);
  _cuts_lframe = new TGCompositeFrame(_cuts_vframe, 1, 1, kVerticalFrame);

  //-- upper composite frame

  _cut_variable_combox = new TGComboBox(_cuts_uframe, 21);

  gui_utils::FillComboBox( _cut_variable_combox,  k_neugen_cut_variable );

  _cut_variable_combox -> Resize(175, 20);

  _cuts_uframe -> AddFrame ( _cut_variable_combox  );

  //-- lower composite frame

  TGMatrixLayout * cuts_nentries_matrix_layout = new TGMatrixLayout(_cuts_lframe, 0, 2, 1);
  _cuts_lframe->SetLayoutManager( cuts_nentries_matrix_layout );

  _cut_min = new TGNumberEntry(_cuts_lframe, 0.000, 9, 3, TGNumberFormat::kNESReal);
  _cut_max = new TGNumberEntry(_cuts_lframe, 0.000, 9, 3, TGNumberFormat::kNESReal);

  _cut_min_label = new TGLabel(_cuts_lframe, new TGString( "min:"));
  _cut_max_label = new TGLabel(_cuts_lframe, new TGString( "max:"));

  _cuts_lframe -> AddFrame ( _cut_min_label );
  _cuts_lframe -> AddFrame ( _cut_min       );
  _cuts_lframe -> AddFrame ( _cut_max_label );
  _cuts_lframe -> AddFrame ( _cut_max       );

  //-- add upper/lower composite frames to main group frame

  _cuts_vframe -> AddFrame ( _cuts_uframe );
  _cuts_vframe -> AddFrame ( _cuts_lframe );

  grpf->AddFrame(_cuts_vframe);

  return grpf;
}
//______________________________________________________________________________
TGGroupFrame * NeuGenInputDialog::BuildSumsFrame(void)
{
  TGGroupFrame * grpf = new TGGroupFrame(_cuts_sums_hframe, "sums", kVerticalFrame);

  _sum_spacer = new TGLabel(grpf, new TGString("     "));

  _qel_bit_mask = new TGCheckButton(grpf, "Add QEL channel",        31);
  _res_bit_mask = new TGCheckButton(grpf, "Add all RES channels",   32);
  _dis_bit_mask = new TGCheckButton(grpf, "Add all DIS channels",   33);

  grpf -> AddFrame ( _sum_spacer  );
  grpf -> AddFrame ( _qel_bit_mask );
  grpf -> AddFrame ( _res_bit_mask );
  grpf -> AddFrame ( _dis_bit_mask );

  return grpf;
}
//______________________________________________________________________________
TGHorizontalFrame * NeuGenInputDialog::BuildButtonsFrame(void)
{
  TGHorizontalFrame * hframe = new TGHorizontalFrame(_main, 10, 10);

  _ok_button     = new TGTextButton(hframe, "   &Ok   ",  1);
  _reset_button  = new TGTextButton(hframe, "  &Reset ",  2);
  _cancel_button = new TGTextButton(hframe, " &Cancel ",  3);

  _button_spacer = new TGLabel(hframe, new TGString("     "));

  _ok_button     -> Connect("Clicked()",
                              "genie::nuvld::NeuGenInputDialog", this, "OK()");
  _reset_button  -> Connect("Clicked()",
                           "genie::nuvld::NeuGenInputDialog", this, "Reset()");
  _cancel_button -> Connect("Clicked()",
                          "genie::nuvld::NeuGenInputDialog", this, "Cancel()");

  hframe -> AddFrame( _ok_button     );
  hframe -> AddFrame( _button_spacer );
  hframe -> AddFrame( _reset_button  );
  hframe -> AddFrame( _button_spacer );
  hframe -> AddFrame( _cancel_button );

  return hframe;
}
//______________________________________________________________________________
void NeuGenInputDialog::OK(void)
{
  NeuGenCards * cards = NeuGenCards::Instance();

  cards->CurrInputs()->SetNBins        ( this->ReadNPoints()      );
  cards->CurrInputs()->SetXSecType     ( this->ReadXSecType()     );
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
  cards->CurrInputs()->SetQelBitInMask ( this->ReadQelBitInMask() );
  cards->CurrInputs()->SetResBitInMask ( this->ReadResBitInMask() );
  cards->CurrInputs()->SetDisBitInMask ( this->ReadDisBitInMask() );
  cards->CurrInputs()->SetInclusive    ( this->ReadInclusive()    );

  this->Report(); // write out user selections to the GUI

  //cout << *(cards->CurrInputs());
  
  _main->SendCloseMessage();
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

  _npoints -> SetIntNumber (100);      // number of points

  _xsec_type_combox -> Select(0);     // xsec type = 'total'

  _plot_range_combox -> Select(0);    // plot range = 'automatic'

  _E       -> SetNumber (   0.0 );
  _E_min   -> SetNumber (   0.1 );
  _E_max   -> SetNumber ( 100.0 );
  _var_min -> SetNumber (   0.0 );
  _var_max -> SetNumber (   0.0 );
  _cut_min -> SetNumber (   0.0 );
  _cut_max -> SetNumber (   0.0 );

  _neutrino_combo    -> Select(0);
  _wcurrent_combo    -> Select(0);
  _initstate_combo   -> Select(0);
  _finstate_listbox  -> Select(0);

  _flux_combox         -> Select (0);  // scaling flux  = 'none'
  _plot_var_combox     -> Select (0);  // plot variable = 'none'
  _cut_variable_combox -> Select (0);  // cut variable  = 'none'

  // sums - all check buttons = ON
  _qel_bit_mask -> SetOn (kTRUE );
  _res_bit_mask -> SetOn (kTRUE );
  _dis_bit_mask -> SetOn (kTRUE );
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

  _npoints -> SetIntNumber ( inp->NBins() );

  _xsec_type_combox->Select(
          gui_utils::ComboBoxSelectionId(k_neugen_xsec_type,
                                          inp->XSecTypeString().c_str()) );

  _plot_range_combox->Select(
          gui_utils::ComboBoxSelectionId(k_neugen_plot_range_option,
                                          inp->PlotRangeString().c_str()) );

  _neutrino_combo->Select(
          gui_utils::ComboBoxSelectionId(k_neugen_nu,
                                             inp->NuTypeString().c_str()) );

  _wcurrent_combo->Select(
          gui_utils::ComboBoxSelectionId(k_neugen_wcurr,
                                          inp->WkCurrentString().c_str()) );

  _initstate_combo->Select(
          gui_utils::ComboBoxSelectionId(k_neugen_initstate,
                                          inp->InitStateString().c_str()) );

  _finstate_listbox->Select(
          gui_utils::ListBoxSelectionId(k_neugen_finstate,
                                         inp->FinalStateString().c_str()) );

  _flux_combox->Select(
          gui_utils::ComboBoxSelectionId(k_scaling_flux,
                                             inp->FluxIdString().c_str()) );

  _plot_var_combox->Select(
          gui_utils::ComboBoxSelectionId(k_neugen_plot_variable,
                                            inp->PlotVarString().c_str()) );

  _cut_variable_combox->Select(
          gui_utils::ComboBoxSelectionId(k_neugen_cut_variable,
                                             inp->CutVarString().c_str()) );

  _E       -> SetNumber ( inp->Energy()     );
  _E_min   -> SetNumber ( inp->EnergyMin()  );
  _E_max   -> SetNumber ( inp->EnergyMax()  );
  _var_min -> SetNumber ( inp->PlotVarMin() );
  _var_max -> SetNumber ( inp->PlotVarMax() );
  _cut_min -> SetNumber ( inp->CutVarMin()  );
  _cut_max -> SetNumber ( inp->CutVarMax()  );

  _qel_bit_mask -> SetOn ( inp->QelBitInMask() );
  _res_bit_mask -> SetOn ( inp->ResBitInMask() );
  _dis_bit_mask -> SetOn ( inp->DisBitInMask() );

  _all_finstates-> SetOn ( inp->Inclusive()    );
}
//______________________________________________________________________________
void NeuGenInputDialog::PositionRelativeToParent(const TGWindow * main)
{
// position relative to the parent's window

  int ax, ay;
  Window_t wdum;

  gVirtualX->TranslateCoordinates(main->GetId(), _main->GetParent()->GetId(),
             (Int_t)(((TGFrame *) main)->GetWidth() - _main->GetWidth()) >> 1,
             (Int_t)(((TGFrame *) main)->GetHeight() - _main->GetHeight()) >> 1,
             ax, ay, wdum);

  _main->Move(ax, ay);
}
//______________________________________________________________________________
int NeuGenInputDialog::ReadNPoints(void)
{
  return _npoints->GetIntNumber();
}
//______________________________________________________________________________
bool NeuGenInputDialog::ReadQelBitInMask(void)
{
  return (_qel_bit_mask->GetState() == kButtonDown);
}
//______________________________________________________________________________
bool NeuGenInputDialog::ReadResBitInMask(void)
{
  return (_res_bit_mask->GetState() == kButtonDown);
}
//______________________________________________________________________________
bool NeuGenInputDialog::ReadDisBitInMask(void)
{
  return (_dis_bit_mask->GetState() == kButtonDown);
}
//______________________________________________________________________________
bool NeuGenInputDialog::ReadInclusive(void)
{
  return (_all_finstates->GetState() == kButtonDown);
}
//______________________________________________________________________________
float NeuGenInputDialog::ReadEnergy(void)
{
  return (float) _E->GetNumber();
}
//______________________________________________________________________________
float NeuGenInputDialog::ReadEnergyMin(void)
{
  return (float) _E_min->GetNumber();
}
//______________________________________________________________________________
float NeuGenInputDialog::ReadEnergyMax(void)
{
  return (float) _E_max->GetNumber();
}
//______________________________________________________________________________
float NeuGenInputDialog::ReadPlotVarMin(void)
{
  return (float) _var_min->GetNumber();
}
//______________________________________________________________________________
float NeuGenInputDialog::ReadPlotVarMax(void)
{
  return (float) _var_max->GetNumber();
}
//______________________________________________________________________________
float NeuGenInputDialog::ReadCutVarMin(void)
{
  return (float) _cut_min->GetNumber();
}
//______________________________________________________________________________
float NeuGenInputDialog::ReadCutVarMax(void)
{
  return (float) _cut_max->GetNumber();
}
//______________________________________________________________________________
string NeuGenInputDialog::ReadFinalState(void)
{
  return gui_utils::ListBoxSelectionAsString(
                                          _finstate_listbox, k_neugen_finstate);
}    
//______________________________________________________________________________
string NeuGenInputDialog::ReadInitialState(void)
{
  return gui_utils::ComboBoxSelectionAsString(
                                          _initstate_combo, k_neugen_initstate);
}
//______________________________________________________________________________
string NeuGenInputDialog::ReadXSecType(void)
{
  return gui_utils::ComboBoxSelectionAsString(
                                         _xsec_type_combox, k_neugen_xsec_type);
}
//______________________________________________________________________________
string NeuGenInputDialog::ReadScalingFlux(void)
{
  return gui_utils::ComboBoxSelectionAsString(_flux_combox, k_scaling_flux);
}
//______________________________________________________________________________
string NeuGenInputDialog::ReadPlotVar(void)
{
  return gui_utils::ComboBoxSelectionAsString(
                                      _plot_var_combox, k_neugen_plot_variable);
}
//______________________________________________________________________________
string NeuGenInputDialog::ReadCutVar(void)
{
  return gui_utils::ComboBoxSelectionAsString(
                                   _cut_variable_combox, k_neugen_cut_variable);
}
//______________________________________________________________________________
string NeuGenInputDialog::ReadPlotRange(void)
{
  return gui_utils::ComboBoxSelectionAsString(
                                  _plot_var_combox, k_neugen_plot_range_option);
}
//______________________________________________________________________________
string NeuGenInputDialog::ReadNeutrino(void)
{
  return gui_utils::ComboBoxSelectionAsString(_neutrino_combo, k_neugen_nu);
}
//______________________________________________________________________________
string NeuGenInputDialog::ReadWkCurrent(void)
{
  return gui_utils::ComboBoxSelectionAsString(_wcurrent_combo, k_neugen_wcurr);
}
//______________________________________________________________________________
void NeuGenInputDialog::Report(void)
{
  SysLogSingleton * syslog = SysLogSingleton::Instance();

  syslog->Log()->AddLine("NeuGEN inputs:");

  syslog->Log()->AddLine( Concat(
                       "npoints..........", this->ReadNPoints())             );
  syslog->Log()->AddLine( Concat(
                       "xsec type........", this->ReadXSecType().c_str())    );
  syslog->Log()->AddLine( Concat(
                       "Emin.............", this->ReadEnergyMin())           );
  syslog->Log()->AddLine( Concat(
                       "Emax.............", this->ReadEnergyMax())           );
  syslog->Log()->AddLine( Concat(
                       "E................", this->ReadEnergy())              );
  syslog->Log()->AddLine( Concat(
                       "Plot variable....", this->ReadPlotVar().c_str())     );
  syslog->Log()->AddLine( Concat(
                       "Scaling flux.....", this->ReadScalingFlux().c_str()) );
  syslog->Log()->AddLine( Concat(
                       "Plot range.......", this->ReadPlotRange().c_str())   );
  syslog->Log()->AddLine( Concat(
                       "Plot var. - min..", this->ReadPlotVarMin())          );
  syslog->Log()->AddLine( Concat(
                       "Plot var. - max..", this->ReadPlotVarMax())          );
  syslog->Log()->AddLine( Concat(
                       "Neutrino type....", this->ReadNeutrino().c_str())    );
  syslog->Log()->AddLine( Concat(
                       "Weak current.....", this->ReadWkCurrent().c_str())   );
  syslog->Log()->AddLine( Concat(
                       "Final state......", this->ReadFinalState().c_str())  );
  syslog->Log()->AddLine( Concat(
                       "Initial state....", this->ReadInitialState().c_str()));
  syslog->Log()->AddLine( Concat(
                       "Cut variable.....", this->ReadCutVar().c_str())      );
  syslog->Log()->AddLine( Concat(
                       "Cut var. - min...", this->ReadCutVarMin())           );
  syslog->Log()->AddLine( Concat(
                       "Cut var. - max...", this->ReadCutVarMax())           );
  syslog->Log()->AddLine( Concat(
                        "Inclusive XSec..", this->ReadInclusive())           );
  syslog->Log()->AddLine( Concat(
                        "QEL bit/mask....", this->ReadQelBitInMask())        );
  syslog->Log()->AddLine( Concat(
                        "RES bit/mask....", this->ReadResBitInMask())        );
  syslog->Log()->AddLine( Concat(
                        "DIS bit/mask....", this->ReadDisBitInMask())        );  
}
//______________________________________________________________________________

