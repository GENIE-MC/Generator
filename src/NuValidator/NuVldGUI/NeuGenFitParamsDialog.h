//_____________________________________________________________________________
/*!

\class    genie::nuvld::NeuGenFitParamsDialog

\brief    A GUI dialog for selecting NeuGEN physics parameters to be fitted

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#ifndef _NEUGEN_FIT_PARAMS_DIALOG_H_
#define _NEUGEN_FIT_PARAMS_DIALOG_H_

#include <TROOT.h>
#include <TApplication.h>
#include <TVirtualX.h>
#include <TSystem.h>
#include <TGClient.h>
#include <TGFrame.h>
#include <TGLabel.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <RQ_OBJECT.h>

#include "NuVldGUI/NeuGenFitParams.h"

namespace genie {
namespace nuvld {

class NeuGenFitParamsDialog {

RQ_OBJECT("NeuGenFitParamsDialog")

public:
   NeuGenFitParamsDialog(const TGWindow *p, const TGWindow *main,
                  UInt_t w, UInt_t h, UInt_t options, NeuGenFitParams * ngfp);
   virtual ~NeuGenFitParamsDialog();

   void CloseWindow (void) { delete this;               }
   void Cancel      (void) { fMain->SendCloseMessage(); }
   void Ok          (void);
   void Reset       (void);
   
private:

   void Report(void) const;
   void SetFitParameters (void);
   void LoadLastEntries  (void);
   void PositionRelativeToParent(const TGWindow * main);

   NeuGenFitParams *  fNGFP;
   
   TGTransientFrame * fMain;
   TGCompositeFrame * fBtnCmpFrm;
   TGGroupFrame  *    fFitParamsGrpFrm;
   TGButton *         fOkBtn;
   TGButton *         fResetBtn;
   TGButton *         fCancelBtn;
   TGLayoutHints *    fBtnLt;
   
   TGCheckButton *    fFitParamChkB[22];
   TGNumberEntry *    fFitParamMinNmEV[22];
   TGNumberEntry *    fFitParamMaxNmEV[22];
   TGNumberEntry *    fFitParamStepNmEV[22];
   TGLabel *          fFitParamLb;
   TGLabel *          fFitSpacerLb;
   TGLabel *          fMinFitParamLb;
   TGLabel *          fMaxFitParamLb;
   TGLabel *          fStepFitParamLb;

   static bool        fHaveNoHistory;
   
   ClassDef(NeuGenFitParamsDialog, 0)
};

} // nuvld namespace
} // genie namespace

#endif // _TEXT_ENTRY_DIALOG_H_
