//_____________________________________________________________________________
/*!

\class    genie::nuvld::NeuGenConfigDialog

\brief    A GUI dialog for reading NeuGEN physics configuration parameters

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#ifndef _NEUGEN_CONFIG_DIALOG_H_
#define _NEUGEN_CONFIG_DIALOG_H_

#include <TApplication.h>
#include <TVirtualX.h>
#include <TGClient.h>
#include <TGFrame.h>
#include <TGLabel.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGTextEdit.h>
#include <TGStatusBar.h>
#include <RQ_OBJECT.h>

namespace genie {
namespace nuvld {

class NeuGenConfigDialog {

RQ_OBJECT("NeuGenConfigDialog")

public:
   NeuGenConfigDialog(const TGWindow *p, const TGWindow *main,
                           UInt_t w, UInt_t h, UInt_t options = kVerticalFrame);
   virtual ~NeuGenConfigDialog();

   void CloseWindow    (void) { delete this;               }
   void Cancel         (void) { fMain->SendCloseMessage(); }
   
   void UseNewValues   (void);
   void UseDefaults    (void);
   void RestoreDefaults(void);

private:

   //-- methods used for building the dialog
   
   void BuildAxialMassFrame      (void);
   void BuildQelParamsFrame      (void);
   void BuildResParamsFrame      (void);
   void BuildCohParamsFrame      (void);
   void BuildPdfSetFrame         (void);
   void BuildKnoParamsFrame      (void);
   void BuildDisResTuningFrame   (void);   
   void BuildButtonsFrame        (void);   
   void PositionRelativeToParent (const TGWindow * main);

   //-- methods used for reading out the dialog widgets

   void  Report       (void) const;
   int   ReadPdfGroup (void) const;
   int   ReadPdfSet   (void) const;
   float ReadMaQel    (void) const;
   float ReadMaRes    (void) const;
   float ReadMaCoh    (void) const;
   float ReadQelFa0   (void) const;
   float ReadQelEta   (void) const;
   float ReadResOmega (void) const;
   float ReadResZ     (void) const;
   float ReadCohR0    (void) const;
   float ReadCohREI   (void) const;
   float ReadKnoAvp   (void) const;
   float ReadKnoAvn   (void) const;
   float ReadKnoAvbp  (void) const;
   float ReadKnoAvbn  (void) const;
   float ReadKnoBvp   (void) const;
   float ReadKnoBvn   (void) const;
   float ReadKnoBvbp  (void) const;
   float ReadKnoBvbn  (void) const;
   float ReadDR2vp    (void) const;
   float ReadDR2vn    (void) const;
   float ReadDR2vbp   (void) const;
   float ReadDR2vbn   (void) const;
   float ReadDR3vp    (void) const;
   float ReadDR3vn    (void) const;
   float ReadDR3vbp   (void) const;
   float ReadDR3vbn   (void) const;

   //-- private data members 
   
   TGTransientFrame *  fMain;
   TGGroupFrame *      fMaGrpFrm;
   TGGroupFrame *      fQelPrmGrpFrm;
   TGGroupFrame *      fResPrmGrpFrm;
   TGGroupFrame *      fCohPrmGrpFrm;
   TGGroupFrame *      fPdfGrpFrm;
   TGGroupFrame *      fDisResGrpFrm;
   TGGroupFrame *      fKnoGrpFrm;
   TGHorizontalFrame * fBtnHFrm;
   TGLayoutHints *     fMaFrmLt;
   TGLayoutHints *     fQelPrmFrmLt;
   TGLayoutHints *     fResPrmFrmLt;
   TGLayoutHints *     fCohPrmFrmLt;
   TGLayoutHints *     fPdfFrmLt;      
   TGLayoutHints *     fDisResFrmLt;
   TGLayoutHints *     fKnoFrmLt;
   TGLayoutHints *     fBtnFrmLt;
   TGNumberEntry *     fMaQelNum;
   TGNumberEntry *     fMaResNum;
   TGNumberEntry *     fMaCohNum;
   TGNumberEntry *     fFa0QelNum;
   TGNumberEntry *     fEtaQelNum;
   TGNumberEntry *     fOmegaResNum;
   TGNumberEntry *     fZResNum;
   TGNumberEntry *     fR0CohNum;
   TGNumberEntry *     fREICohNum;
   TGNumberEntry *     fPdfGroupNum;
   TGNumberEntry *     fPdfSetNum;
   TGNumberEntry *     fKnoBvpNum;
   TGNumberEntry *     fKnoBvnNum;
   TGNumberEntry *     fKnoBvbpNum;
   TGNumberEntry *     fKnoBvbnNum;
   TGNumberEntry *     fKnoAvpNum;
   TGNumberEntry *     fKnoAvnNum;
   TGNumberEntry *     fKnoAvbpNum;
   TGNumberEntry *     fKnoAvbnNum;
   TGNumberEntry *     fDisResM2vpNum;
   TGNumberEntry *     fDisResM2vnNum;
   TGNumberEntry *     fDisResM2vbpNum;
   TGNumberEntry *     fDisResM2vbnNum;
   TGNumberEntry *     fDisResM3vpNum;
   TGNumberEntry *     fDisResM3vnNum;
   TGNumberEntry *     fDisResM3vbpNum;
   TGNumberEntry *     fDisResM3vbnNum;
   TGLabel *           fMaQelLbl;
   TGLabel *           fMaResLbl;
   TGLabel *           fMaCohLbl;
   TGLabel *           fFa0QelLbl;
   TGLabel *           fEtaQelLbl;
   TGLabel *           fOmegaResLbl;
   TGLabel *           fZResLbl;
   TGLabel *           fR0CohLbl;
   TGLabel *           fREICohLbl;
   TGLabel *           fPdfGroupLbl;
   TGLabel *           fPdfSetLbl;
   TGLabel *           fMultLbl;
   TGLabel *           fMultEq2Lbl;
   TGLabel *           fMultEq3Lbl;
   TGLabel *           fKnoISvpLbl;
   TGLabel *           fKnoISvnLbl;
   TGLabel *           fKnoISvbpLbl;
   TGLabel *           fKnoISvbnLbl;
   TGLabel *           fKnoALbl;   
   TGLabel *           fKnoBLbl;   
   TGLabel *           fKnoNullLbl;   
   TGLabel *           fISvpLbl;
   TGLabel *           fISvnLbl;
   TGLabel *           fISvbpLbl;
   TGLabel *           fISvbnLbl;
   TGButton *          fUseInputsBtn;
   TGButton *          fCancelBtn;
   TGButton *          fRestoreDefaultBtn;
   TGButton *          fUseDefaultsBtn;

   ClassDef(NeuGenConfigDialog, 0)
};

} // nuvld namespace
} // genie namespace

#endif

