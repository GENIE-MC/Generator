//_____________________________________________________________________________
/*!

\class    genie::nuvld::NeuGenInputDialog

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#ifndef _NEUGEN_INPUT_DIALOG_H_
#define _NEUGEN_INPUT_DIALOG_H_

#include <string>
#include <TApplication.h>
#include <TVirtualX.h>
#include <TGClient.h>
#include <TGFrame.h>
#include <TGLabel.h>
#include <TGButton.h>
#include <TGComboBox.h>
#include <TGNumberEntry.h>
#include <TGTextEdit.h>
#include <TGStatusBar.h>
#include <RQ_OBJECT.h>

using std::string;

namespace genie {
namespace nuvld {

class NeuGenInputDialog {

RQ_OBJECT("NeuGenInputDialog")

public:
   NeuGenInputDialog(const TGWindow *p, const TGWindow *main,
                           UInt_t w, UInt_t h, UInt_t options = kVerticalFrame);
   virtual ~NeuGenInputDialog();

   void CloseWindow  (void) { delete this;               }
   void Cancel       (void) { _main->SendCloseMessage(); }
   void Reset        (void);
   void OK           (void);

private:

   TGGroupFrame *      BuildXSecTypeFrame         (void);
   TGGroupFrame *      BuildConfigTotalXSecFrame  (void);
   TGGroupFrame *      BuildConfigDiffXSecFrame   (void);
   TGGroupFrame *      BuildInteractionFrame      (void);
   TGGroupFrame *      BuildCutsFrame             (void);
   TGGroupFrame *      BuildSumsFrame             (void);
   TGGroupFrame *      BuildDrawOptionsFrame      (void);
   TGHorizontalFrame * BuildButtonsFrame          (void);

   void   PositionRelativeToParent (const TGWindow * main);
   
   void   Defaults         (void);
   void   LoadLastEntries  (void);
   void   Report           (void);
   
   int    ReadNPoints      (void);
   bool   ReadQelBitInMask (void);
   bool   ReadResBitInMask (void);
   bool   ReadDisBitInMask (void);
   bool   ReadInclusive    (void);
   float  ReadEnergy       (void);
   float  ReadEnergyMin    (void);
   float  ReadEnergyMax    (void);
   float  ReadPlotVarMin   (void);
   float  ReadPlotVarMax   (void);
   float  ReadCutVarMin    (void);
   float  ReadCutVarMax    (void);
   string ReadInitialState (void);
   string ReadFinalState   (void);
   string ReadXSecType     (void);
   string ReadScalingFlux  (void);
   string ReadPlotVar      (void);
   string ReadCutVar       (void);
   string ReadPlotRange    (void);
   string ReadNeutrino     (void);
   string ReadWkCurrent    (void);
   
   TGTransientFrame *    _main;
   TGGroupFrame *        _xsec_type_grpf;
   TGGroupFrame *        _tot_xsec_conf_grpf;
   TGGroupFrame *        _dis_xsec_conf_grpf;
   TGGroupFrame *        _interaction_grpf;
   TGGroupFrame *        _cuts_grpf;
   TGGroupFrame *        _sums_grpf;
   TGHorizontalFrame *   _button_frame;
   TGHorizontalFrame *   _xsec_type_hframe;
   TGHorizontalFrame *   _tot_xsec_conf_hframe;
   TGHorizontalFrame *   _dif_xsec_conf_hframe;
   TGHorizontalFrame *   _int_hframe;
   TGHorizontalFrame *   _cuts_sums_hframe;
   TGVerticalFrame *     _cuts_vframe;
   TGVerticalFrame *     _dif_xsec_conf_vframe;
   TGCompositeFrame *    _dif_xsec_conf_uframe;
   TGCompositeFrame *    _dif_xsec_conf_cframe;
   TGCompositeFrame *    _dif_xsec_conf_lframe;
   TGCompositeFrame *    _int_lframe;
   TGCompositeFrame *    _int_cframe;
   TGCompositeFrame *    _int_rframe;
   TGCompositeFrame *    _cuts_uframe;
   TGCompositeFrame *    _cuts_lframe;
   TGLayoutHints *       _xsec_type_frame_layout;
   TGLayoutHints *       _tot_xsec_conf_frame_layout;
   TGLayoutHints *       _dis_xsec_conf_frame_layout;
   TGLayoutHints *       _interaction_frame_layout;
   TGLayoutHints *       _cuts_sums_frame_layout;
   TGLayoutHints *       _draw_options_frame_layout;
   TGLayoutHints *       _button_frame_layout;
   TGLayoutHints *       _xsec_type_rbutton_layout;
   TGCheckButton *       _all_finstates;
   TGCheckButton *       _qel_bit_mask;
   TGCheckButton *       _res_bit_mask;
   TGCheckButton *       _dis_bit_mask;
   TGComboBox *          _xsec_type_combox;
   TGComboBox *          _flux_combox;
   TGComboBox *          _plot_var_combox;
   TGComboBox *          _plot_range_combox;
   TGComboBox *          _neutrino_combo;
   TGComboBox *          _wcurrent_combo;
   TGComboBox *          _initstate_combo;
   TGComboBox *          _cut_variable_combox;
   TGListBox *           _finstate_listbox;
   TGButton *            _ok_button;
   TGButton *            _reset_button;
   TGButton *            _cancel_button;
   TGNumberEntry *       _npoints;
   TGNumberEntry *       _E;
   TGNumberEntry *       _E_min;
   TGNumberEntry *       _E_max;
   TGNumberEntry *       _var_min;
   TGNumberEntry *       _var_max;      
   TGNumberEntry *       _cut_min;
   TGNumberEntry *       _cut_max;
   TGLabel *             _xsec_type_label;
   TGLabel *             _flux_label;
   TGLabel *             _plot_var_label;
   TGLabel *             _npoints_label;
   TGLabel *             _E_label;
   TGLabel *             _E_min_label;
   TGLabel *             _E_max_label;
   TGLabel *             _plot_range_label;
   TGLabel *             _diff_xsec_conf_spacer;
   TGLabel *             _var_min_label;
   TGLabel *             _var_max_label;
   TGLabel *             _cut_min_label;
   TGLabel *             _cut_max_label;
   TGLabel *             _int_spacer_label;   
   TGLabel *             _neutrino_label;
   TGLabel *             _wcurrent_label;
   TGLabel *             _initstate_label;
   TGLabel *             _finstate_label;
   TGLabel *             _tot_xsec_spacer;
   TGLabel *             _button_spacer;
   TGLabel *             _sum_spacer;
   static bool           _have_no_history;
   
   ClassDef(NeuGenInputDialog, 0)
};

} // nuvld namespace
} // genie namespace

#endif

