//_____________________________________________________________________________
/*!

\class    genie::nuvld::NuVldMainFrame

\brief    NuValidator GUI prototype - main frame

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#ifndef _NUVLD_MAIN_FRAME_H_
#define _NUVLD_MAIN_FRAME_H_

#include <vector>

#include <TApplication.h>
#include <TVirtualX.h>
#include <TSystem.h>
#include <TGListBox.h>
#include <TGComboBox.h>
#include <TGClient.h>
#include <TGFrame.h>
#include <TGIcon.h>
#include <TGLabel.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGTextEntry.h>
#include <TGMsgBox.h>
#include <TGMenu.h>
#include <TGCanvas.h>
#include <TGTab.h>
#include <TGFileDialog.h>
#include <TGTextEdit.h>
#include <TGStatusBar.h>
#include <TGProgressBar.h>
#include <TGColorSelect.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TRootEmbeddedCanvas.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TF1.h>
#include <RQ_OBJECT.h>

#include "DBUtils/DBTable.h"
#include "DBUtils/vXSecTableRow.h"
#include "DBUtils/eDiffXSecTableRow.h"
#include "Facades/NGCardPairList.h"
#include "NuVldGUI/NeuGenFitParams.h"
#include "Facades/XSecVsEnergy.h"
#include "NuVldGUI/GuiHelpHandler.h"
#include "NuVldGUI/GuiStackHandler.h"
#include "NuVldGUI/GuiDBHandler.h"
#include "NuVldGUI/GuiXmlFileHandler.h"
#include "NuVldGUI/GuiFitKernel.h"
#include "NuVldGUI/vDataSelectionDialog.h"
#include "NuVldGUI/MenuCommandId.h"

class TLatex;

using std::vector;

namespace genie {
namespace nuvld {

class NuVldMainFrame : public TGMainFrame {

public:

   friend class GuiHelpHandler;
   
   NuVldMainFrame(const TGWindow * p, UInt_t w, UInt_t h);
   virtual ~NuVldMainFrame();

   //-- methods for handling GUI events

   void CloseWindow  (void) { gApplication->Terminate(0); }
   void Exit         (void) { this->CloseWindow();        }
   
   void HandleMenu  (Int_t id);

   void ConfigNeugenPhysics   (void);  
   void ConfigNeugenProcess   (void);
   void SelectNeuGenFitParams (void);
   void RunNulook             (void);  
   void LoadExtXSecPrediction (void);

   void HandleSaveCanvas      (void);
   
   void DrawDBTable           (void);
   void PrintDBTable          (void);
   void DrawCurrentDBTable    (void);
   void PrintCurrentDBTable   (void);

   void SetCurrDBTable        (void);
   void RetrieveStackedDBTable(void);

   void DrawNeugenXSecVsEnergy (XSecVsEnergy * xs, TRootEmbeddedCanvas * ecanvas, bool show_titles = true);
   
   bool CheckNeugenCards(void);

   //-- GUI fitter methods

   void  ResetFitterTab       (void);
   void  RunFitter            (void);
   void  RunMcScanner         (void);
   void  Run1dScanner         (void);
   void  Run2dScanner         (void);
   void  RunPostFitProcessor  (void);
   void  PrintFitParameters   (void);
   void  DrawResiduals        (void);
   void  PlotXSecBoundaries   (TCanvas * c, bool clear);

   //-- methods for reseting SQL GUI widgets & viewers

   void ResetSqlSelections    (void);
   void ResetElSqlSelections  (void);
   void ResetNuSqlSelections  (void);
   void ResetCommonSelections (void);
   void SelectAllNuExp        (void);
   void SelectAllNuXSec       (void);
   void SelectAllNuProbes     (void);
   void SelectAllNuTargets    (void);
   void SelectAllElExp        (void);
   void SelectAllElTargets    (void);
   void ClearViewer           (void);

   //-- methods for switching tabs

   void OpenPlotterTab    (void);
   void OpenDataViewerTab (void);
   void OpenFitterTab     (void);
   void OpenSessionLogTab (void);
   
   //-- methods for poping up data selection dialogs

   void PopupNuDataSelectionDialog   (void);
   void PopupNuMeasurementListDialog (void);

private:

   //-- initialization & configuration methods

   void Init                   (void);
   void InitializeHandlers     (void);
   void InitializeSyslog       (void);
   void InitializeBrowser      (void);
   void ConfigHandlers         (void);

   //-- methods for building main frame gui widgets

   void                DefineLayoutHints        (void);
   TGMenuBar *         BuildMenuBar             (void);
   TGTab *             BuildSqlTab              (void);
   TGTab *             BuildDataTab             (void);
   TGGroupFrame *      BuildUpperButtonFrame    (void);
   TGHorizontalFrame * BuildSelectionStackFrame (void);
   TGHorizontalFrame * BuildLowerButtonFrame    (void);
   TGStatusBar *       BuildStatusBar           (void);   
   void                FillNuSqlFrame           (void);
   void                FillElSqlFrame           (void);
   void                AddCommonCheckButtons    (void);
   void                FillFitterFrame          (void);
   void                CreateUpperFrameButtons  (TGGroupFrame * gf);
   void                SetUpperFrameButtonText  (void);
   void                ConnectUpperFrameButtons (void);
   const TGPicture *   Pic  (const char * name, int x, int y);   
   const char *        Icon (const char * name);

   //-- methods for handling data selections
   
   string  ReadXSecSelectionListbox  (void);

   string  NuDataSelections (void);
   string  ElDataSelections (void);

   bool    ScaleWithEnergy (void);
   string  PlotVariable    (void);
      
   //-- methods mimicing the vDataSelectionDialog interface for
   //   the data selection tabs embedded in the main GUI window

   string NuTabBundleSelectionsInString (void);
   string NuTabBundleKeyListInString    (void);
   string NuTabBundleCutsInString       (void);
   string NuTabBundleDrawOptInString    (void);
      
   string ElTabBundleSelectionsInString (void);
   string ElTabBundleKeyListInString    (void);
   string ElTabBundleCutsInString       (void);
   string ElTabBundleDrawOptInString    (void);

   //-- methods for extracting cross section data
   
   DBTable<vXSecTableRow> *     FillNuXSecTable     (void);
   DBTable<eDiffXSecTableRow> * FillElDiffXSecTable (void);

   //-- GUI widgets
   
   TGMainFrame *             fMain;
   TGMenuBar *               fMenu;
   TGPopupMenu *             fMenuFile;
   TGPopupMenu *             fMenuDBase;
   TGPopupMenu *             fMenuNeuGen;
   TGPopupMenu *             fMenuFit;
   TGPopupMenu *             fMenuHelp;
   TGTab *                   fTabSql;
   TGTab *                   fTabData;
   TGCompositeFrame *        fTabPlotter;
   TGCompositeFrame *        fTabDataViewer;
   TGCompositeFrame *        fTabFitter;
   TGCompositeFrame *        fTabLog;
   TGCompositeFrame *        fTabNuSql;
   TGCompositeFrame *        fTabElSql;
   TGCompositeFrame *        fMainFrame;
   TGCompositeFrame *        fMainTopFrame;
   TGCompositeFrame *        fMainBottomFrame;
   TGCompositeFrame *        fMainLeftFrame;
   TGCompositeFrame *        fMainRightFrame;
   TGCompositeFrame *        fFitterLeftFrame;
   TGCompositeFrame *        fFitterRightFrame;
   TGStatusBar *             fStatusBar;
   TGHProgressBar *          fProgressBar;
   TGTextEdit *              fDataViewer;   
   TGTextEdit *              fLog;
   TGTextEdit *              fFitTxtResults;
   TRootEmbeddedCanvas *     fPlotTabEmbCnv;
   TRootEmbeddedCanvas *     fFitTabFuncEmbCnv;
   TRootEmbeddedCanvas *     fFitTabChisqEmbCnv;   
   TGLayoutHints *           fMenuBarLt;
   TGLayoutHints *           fMenuBarItemLt;
   TGLayoutHints *           fMenuBarHelpLt;
   TGLayoutHints *           fPlotterTabLt;
   TGLayoutHints *           fDataViewTabLt;
   TGLayoutHints *           fFitterTabLt;
   TGLayoutHints *           fNuSqlTabLt; 
   TGLayoutHints *           fElSqlTabLt;
   TGLayoutHints *           fLogTabLt;
   TGLayoutHints *           fDataTabLt;
   TGLayoutHints *           fSqlTabLt;
   TGLayoutHints *           fProgressBarLt;
   TGLayoutHints *           fSelStackLt;
   TGLayoutHints *           fExitBtnLt;
   TGLayoutHints *           fLeftBtnLt;
   TGLayoutHints *           fStatusBarLt;
   TGLayoutHints *           fMLeftFrameLt;
   TGLayoutHints *           fMRightFrameLt;
   TGLayoutHints *           fFitLeftFrameLt;
   TGLayoutHints *           fFitRightFrameLt;
   TGMatrixLayout *          fBtnMatrixLt;
   TGMatrixLayout *          fEnergyMatrixLt;
   vector<TGMatrixLayout * > fElVarRangeLt;   
   TGGroupFrame *            fNuXSecErrGrpFrm;
   TGGroupFrame *            fNuExpGrpFrm;
   TGGroupFrame *            fNuXSecGrpFrm;
   TGGroupFrame *            fNuInitStateGrpFrm;
   vector<TGGroupFrame * >   fElVarRangeGrpFrm;   
   TGGroupFrame *            fImgBtnGrpFrm;
   TGGroupFrame *            fEnergyGrpFrm;
   TGGroupFrame *            fFitterGrpFrm;
   TGGroupFrame *            fFitFreeParamGrpFrm;
   TGGroupFrame *            fFitBtnGrpFrm;
   TGGroupFrame *            fElExpGrpFrame;
   TGGroupFrame *            fElTgGrpFrm;
   TGGroupFrame *            fElDrawXGrpFrm;
   TGListBox *               fNuXSecErrLBx;
   TGListBox *               fNuExpLBx;
   TGListBox *               fNuProcLBx;
   TGListBox *               fNuTypeLBx;
   TGListBox *               fNuTgtLBx;
   TGListBox *               fElExpLBx;
   TGListBox *               fElTgtLBx;
   TGComboBox *              fTableStackCBx;
   TGComboBox *              fConfigStackCBx;
   TGComboBox *              fFitterCBx;
   TGComboBox *              fElDrawXCBx;   
   TGPictureButton *         fExitBtn;
   TGPictureButton *         fOpenXmlBtn;
   TGPictureButton *         fParseXmlBtn;
   TGPictureButton *         fDBConnectBtn;
   TGPictureButton *         fDBBootstrapBtn;
   TGPictureButton *         fSqlQInpBtn;
   TGPictureButton *         fSqlQFileBtn;
   TGPictureButton *         fDBUploadBtn;
   TGPictureButton *         fNeugenConfigBtn;
   TGPictureButton *         fNeugenProcBtn;
   TGPictureButton *         fNeugenRunBtn;
   TGPictureButton *         fDrawDataBtn;
   TGPictureButton *         fViewClearBtn;
   TGPictureButton *         fSaveBtn;
   TGPictureButton *         fHelpBtn;
   TGPictureButton *         fDurhamBtn;
   TGPictureButton *         fAboutBtn;  
   TGPictureButton *         fSelResetBtn;
   TGPictureButton *         fPrintDataBtn;   
   TGPictureButton *         fDBCloseBtn;
   TGPictureButton *         fDBCheckBtn;
   TGPictureButton *         fDBInfoBtn;
   TGPictureButton *         fStackTableBtn;
   TGPictureButton *         fStackConfigBtn;
   TGPictureButton *         fLinkStackedBtn;
   TGPictureButton *         fDelStackedBtn;
   TGPictureButton *         fDoFitBtn;
   TGPictureButton *         fPrmScanBtn;
   TGPictureButton *         fPrmScan1dBtn;
   TGPictureButton *         fPrmScan2dBtn;
   TGPictureButton *         fResetFitBtn;   
   TGTextButton *            fShowFullNuDialogTBtn;
   TGTextButton *            fShowExpertNuDialogTBtn;
   TGTextButton *            fSelectNeuGenFitParams;
   TGHorizontalFrame *       fProgressBarHFrm;
   TGHorizontalFrame *       fStackHFrm;   
   TGCheckButton *           fAllNuExpChkB;
   TGCheckButton *           fAllNuProcChkB;
   TGCheckButton *           fAllNuTypesChkB;
   TGCheckButton *           fAllNuTgtChkB;
   TGCheckButton *           fAllElExpChkB;
   TGCheckButton *           fAllElTgtChkB;
   TGCheckButton *           fScaleWithEvChkB;
   TGCheckButton *           fShowColorCodeChkB;
   TGCheckButton *           fShowExtLegendChkB;   
   TGCheckButton *           fUseStackedChkB;   
   TGNumberEntry *           fEMinNmE;
   TGNumberEntry *           fEMaxNmE;
   vector<TGNumberEntry * >  fElVarMinNmEV;
   vector<TGNumberEntry * >  fElVarMaxNmEV;
   TGNumberEntry *           fXMinNmE;
   TGNumberEntry *           fXMaxNmE;
   TGTextEntry *             fStackTableNameTxE;
   TGTextEntry *             fStackConfigNameTxE;
   TGLabel *                 fMinELb;
   TGLabel *                 fMaxELb;
   TGLabel *                 fNuTabBtnSpacerLb;
   TGLabel *                 fXMinLb;
   TGLabel *                 fXMaxLb;
   TGLabel *                 fStackDBTableLb;
   TGLabel *                 fStackConfigLb;
   TGLabel *                 fLinkSelLb;
   TGLabel *                 fLFitSpacerLb;
   TGLabel *                 fRFitSpacerLb;
   TLatex *                  fLtxAuth;
   TLatex *                  fLtxLink;

   //-- other private date members

   DBConnection *            fDBC;
   NeuGenFitParams *         fNGFP;
   
   //-- 'action' objects that handle some classes of GUI events
   
   GuiHelpHandler *          fHelpHandler;
   GuiDBHandler *            fDBaseHandler;
   GuiXmlFileHandler *       fXmlFileHandler;
   GuiStackHandler *         fStackHandler;
   GuiFitKernel *            fFitKernel;

   bool                      fPlotterShowIsOn;
   
   bool                      _v_slc_dialog_requires_attn;
   DataSelectionDialog *     _active_v_slc_dialog;
   XSecVsEnergy *            _xsec_vs_energy;

ClassDef(NuVldMainFrame, 0)
};

} // nuvld namespace
} // genie namespace

#endif

