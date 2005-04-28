//_____________________________________________________________________________
/*!

\class    genie::nuvld::NeuGenInputDialog

\brief    A GUI dialog for specifying what information to extract from a
          neutrino generator.

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#ifndef _NEUGEN_INPUT_DIALOG_H_
#define _NEUGEN_INPUT_DIALOG_H_

#include <string>

#include <RQ_OBJECT.h>

class TGGroupFrame;
class TGHorizontalFrame;
class TGTransientFrame;
class TGLayoutHints;
class TGComboBox;
class TGListBox;
class TGButton;
class TGNumberEntry;
class TGLabel;

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
   void Cancel       (void) { fMain->SendCloseMessage(); }
   void Reset        (void);
   void OK           (void);

private:

   TGGroupFrame *      BuildXSecTypeFrame         (void);
   TGGroupFrame *      BuildConfigTotalXSecFrame  (void);
   TGGroupFrame *      BuildConfigDiffXSecFrame   (void);
   TGGroupFrame *      BuildConfigSFFrame         (void);
   TGGroupFrame *      BuildPlotVarFrame          (void);  
   TGGroupFrame *      BuildInteractionFrame      (void);
   TGGroupFrame *      BuildCutsFrame             (void);
   TGGroupFrame *      BuildSumsFrame             (void);
   TGGroupFrame *      BuildDrawOptionsFrame      (void);
   TGHorizontalFrame * BuildButtonsFrame          (void);

   void    PositionRelativeToParent (const TGWindow * main);
   
   void    Defaults         (void);
   void    LoadLastEntries  (void);
   void    Report           (void);
   
   int     ReadNPoints      (void);
   bool    ReadQelBitInMask (void);
   bool    ReadResBitInMask (void);
   bool    ReadDisBitInMask (void);
   bool    ReadInclusive    (void);
   bool    ReadSFRawDis     (void);
   float   ReadEnergy       (void);
   float   ReadEnergyMin    (void);
   float   ReadEnergyMax    (void);
   float   ReadPlotVarMin   (void);
   float   ReadPlotVarMax   (void);
   float   ReadCutVarMin    (void);
   float   ReadCutVarMax    (void);
   float   ReadSFFixedVar   (void);
   float   ReadA            (void);
   string  ReadInitialState (void);
   string  ReadFinalState   (void);
   string  ReadPlotType     (void);
   string  ReadScalingFlux  (void);
   string  ReadPlotVar      (void);
   string  ReadCutVar       (void);
   string  ReadPlotRange    (void);
   string  ReadNeutrino     (void);
   string  ReadWkCurrent    (void);
   string  ReadSF           (void);
   
   TGTransientFrame *    fMain;
   TGGroupFrame *        fPlotTypeGrpf;
   TGGroupFrame *        fTotXSecConfGrpf;
   TGGroupFrame *        fDifXSecConfGrpf;
   TGGroupFrame *        fSFConfGrpf;
   TGGroupFrame *        fPlotVarGrpf;
   TGGroupFrame *        fInteractionGrpf;
   TGGroupFrame *        fCutsGrpf;
   TGGroupFrame *        fSumsGrpf;
   TGHorizontalFrame *   fBtnHFrm;
   TGHorizontalFrame *   fPlotTypeHFrm;
   TGHorizontalFrame *   fTotXSecConfHFrm;
   TGHorizontalFrame *   fInteractionHFrm;
   TGHorizontalFrame *   fCutsSumsHFrm;
   TGVerticalFrame *     fCutsVFrm;
   TGCompositeFrame *    fInteractionLFrm;
   TGCompositeFrame *    fInteractionCFrm;
   TGCompositeFrame *    fInteractionRFrm;
   TGCompositeFrame *    fCutsUFrm;
   TGCompositeFrame *    fCutsLFrm;
   TGLayoutHints *       fPlotTypeFrmLt;
   TGLayoutHints *       fTotXSecConfFrmLt;
   TGLayoutHints *       fDifXSecConfFrmLt;
   TGLayoutHints *       fSFConfFrmLt;
   TGLayoutHints *       fPlotVarFrmLt;
   TGLayoutHints *       fInteractionFrmLt;
   TGLayoutHints *       fCutsSumsFrmLt;
   TGLayoutHints *       fBtnFrmLt;
   TGCheckButton *       fAllFinStatesCkb;
   TGCheckButton *       fQelBitMaskCkb;
   TGCheckButton *       fResBitMaskCkb;
   TGCheckButton *       fDisBitMaskCkb;
   TGCheckButton *       fSFRawDisCkb;
   TGComboBox *          fPlotTypeCbx;
   TGComboBox *          fFluxCbx;
   TGComboBox *          fPlotVarCbx;
   TGComboBox *          fPlotRangeCbx;
   TGComboBox *          fProbeCbx;
   TGComboBox *          fWkCurrCbx;
   TGComboBox *          fInitStateCbx;
   TGComboBox *          fCutVarCbx;
   TGComboBox *          fSFTypeCbx;
   TGListBox *           fFinStateLbx;
   TGButton *            fOkBtn;
   TGButton *            fResetBtn;
   TGButton *            fCancelBtn;
   TGNumberEntry *       fNPointsNmE;
   TGNumberEntry *       fEnergyNmE;
   TGNumberEntry *       fMinEnergyNmE;
   TGNumberEntry *       fMaxEnergyNmE;
   TGNumberEntry *       fMinVarNmE;
   TGNumberEntry *       fMaxVarNmE;      
   TGNumberEntry *       fMinCutNmE;
   TGNumberEntry *       fMaxCutNmE;
   TGNumberEntry *       fSFFixVarNmE;
   TGNumberEntry *       fANmE;
   TGLabel *             fPlotTypeLb;
   TGLabel *             fFluxLb;
   TGLabel *             fPlotVarLb;
   TGLabel *             fNpLb;
   TGLabel *             fEnergyLb;
   TGLabel *             fMinEnergyLb;
   TGLabel *             fMaxEnergyLb;
   TGLabel *             fPlotRangeLb;
   TGLabel *             fMinVarLb;
   TGLabel *             fMaxVarLb;
   TGLabel *             fSFFixVarLb;
   TGLabel *             fSFLb;
   TGLabel *             fMinCutLb;
   TGLabel *             fMaxCutLb;
   TGLabel *             fInteractionSpacer;   
   TGLabel *             fProbeLb;
   TGLabel *             fWkCurrLb;
   TGLabel *             fInitStateLb;
   TGLabel *             fFinStateLb;
   TGLabel *             fTotXSecSpacer;
   TGLabel *             fBtnSpacer;
   TGLabel *             fSumSpacer;
   TGLabel *             fALb;
   static bool           fHaveNoHistory;
   
   ClassDef(NeuGenInputDialog, 0)
};

} // nuvld namespace
} // genie namespace

#endif

