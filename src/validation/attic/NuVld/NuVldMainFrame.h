//_____________________________________________________________________________
/*!

\class    genie::nuvld::NuVldMainFrame

\brief    NuValidator GUI prototype - main frame

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#ifndef _NUVLD_MAIN_FRAME_H_
#define _NUVLD_MAIN_FRAME_H_

#include <TApplication.h>
#include <TVirtualX.h>
#include <TGFrame.h>
#include <RQ_OBJECT.h>

#include "ValidationTools/NuVld/DBTable.h"
#include "ValidationTools/NuVld/DBNuXSecTableRow.h"
#include "ValidationTools/NuVld/DBElDiffXSecTableRow.h"
#include "ValidationTools/NuVld/DBSFTableRow.h"
#include "ValidationTools/NuVld/GuiMenuCommandId.h"
#include "ValidationTools/NuVld/NuVldConfig.h"

class TSystem;
class TGListBox;
class TGComboBox;
class TGClient;
class TGIcon;
class TGLabel;
class TGButton;
class TGNumberEntry;
class TGPopupMenu;
class TGTextEntry;
class TGGuiMsgBox;
class TGDockableFrame;
class TGMenu;
class TGMenuBar;
class TGCanvas;
class TGTab;
class TGFileDialog;
class TGTextEdit;
class TGStatusBar;
class TGHProgressBar;
class TGColorSelect;
class TGPictureButton;
class TGTextButton;
class TGCheckButton;
class TCanvas;
class TGraphAsymmErrors;
class TRootEmbeddedCanvas;
class TLatex;

using namespace genie;

namespace genie {

class Spline;

namespace nuvld {

class DBConnection;
class GuiHelpHandler;
class GuiStackHandler;
class GuiDBHandler;
class GuiXmlFileHandler;
class NeuGenFitParams;
class GuiDataSelectionDialog;
class GuiNuDataSelectionTab;
class GuiElDataSelectionTab;
class GuiSFDataSelectionTab;

class NuVldMainFrame : public TGMainFrame {

public:

   friend class GuiHelpHandler;

   NuVldMainFrame(const TGWindow * p, UInt_t w, UInt_t h, const NuVldConfig & my_config);
   virtual ~NuVldMainFrame();

   //-- methods for handling GUI events

   void CloseWindow  (void) { gApplication->Terminate(0); }
   void Exit         (void) { this->CloseWindow();        }

   void HandleMenu  (Int_t id);

   void ConfigNeugenPhysics   (void);
   void ConfigNeugenProcess   (void);
   void SelectNeuGenFitParams (void);
   void RunNeuGen             (void);
   void RetrieveNeuGenCards   (void);
   void LoadExtXSecPrediction (void);
   void HandleSaveCanvas      (void);
   void HandleSaveTableAsText (void);
   void HandleSaveSpline      (void);
   void DrawDBTable           (void);
   void PrintDBTable          (void);
   void DrawCurrentDBTable    (void);
   void PrintCurrentDBTable   (void);
   void SetCurrDBTable        (void);
   void RetrieveStackedDBTable(void);
   void UpdateFrameSize       (void);

   void DrawSpline (Spline * xs, TRootEmbeddedCanvas * ecanvas, bool show_titles = true);

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
   void ResetCommonSelections (void);
   void ClearViewer           (void);

   //-- methods for switching tabs

   void OpenPlotterTab        (void);
   void OpenDataViewerTab     (void);
   void OpenFitterTab         (void);
   void OpenSessionLogTab     (void);

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
   void                AddCommonCheckButtons    (void);
   void                FillFitterFrame          (void);
   void                CreateUpperFrameButtons  (TGGroupFrame * gf);
   void                SetUpperFrameButtonText  (void);
   void                ConnectUpperFrameButtons (void);
   void                HandleDisabledNeuGEN     (void);
   void                InitNotify               (void);
   const TGPicture *   Pic  (const char * name, int x, int y);
   const char *        Icon (const char * name);

   //-- methods for handling data selections

   bool    ScaleWithEnergy (void);
   string  PlotVariable    (void);
   string  ErrorOption     (void);

   //-- methods for extracting cross section data

   DBTable<DBNuXSecTableRow> *     FillNuXSecTable     (void);
   DBTable<DBElDiffXSecTableRow> * FillElDiffXSecTable (void);
   DBTable<DBSFTableRow> *         FillSFTable         (void);

   //-- GUI widgets

   TGMainFrame *             fMain;
   TGMenuBar *               fMenu;
   TGDockableFrame *         fMenuDock;
   TGDockableFrame *         fDBSelDock;
   TGDockableFrame *         fBottomFrmDock;
   TGPopupMenu *             fMenuFile;
   TGPopupMenu *             fMenuView;
   TGPopupMenu *             fMenuDBase;
   TGPopupMenu *             fMenuExpImp;
   TGPopupMenu *             fMenuNeuGen;
   TGPopupMenu *             fMenuGENIE;
   TGPopupMenu *             fMenuFit;
   TGPopupMenu *             fMenuHelp;
   TGTab *                   fTabSql;
   TGTab *                   fTabData;
   TGCompositeFrame *        fTabPlotter;
   TGCompositeFrame *        fTabDataViewer;
   TGCompositeFrame *        fTabFitter;
   TGCompositeFrame *        fTabLog;
   TGCompositeFrame *        fTabNuSql; // owned by tab object
   TGCompositeFrame *        fTabElSql; // owned by tab object
   TGCompositeFrame *        fTabSFSql; // owned by tab object
   TGCompositeFrame *        fMainFrame;
   TGCompositeFrame *        fMainTopFrame;
   TGCompositeFrame *        fMainMiddleFrame;
   TGCompositeFrame *        fMainBottomFrame;
   TGCompositeFrame *        fMainDBottomFrame;
   TGCompositeFrame *        fMainLeftFrame;
   TGCompositeFrame *        fMainRightFrame;
   TGCompositeFrame *        fBottomLeftFrame;
   TGCompositeFrame *        fBottomRightFrame;
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
   TGLayoutHints *           fSFSqlTabLt;
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
   TGGroupFrame *            fImgBtnGrpFrm;
   TGGroupFrame *            fEnergyGrpFrm;
   TGGroupFrame *            fFitterGrpFrm;
   TGGroupFrame *            fFitFreeParamGrpFrm;
   TGGroupFrame *            fFitBtnGrpFrm;
   TGComboBox *              fTableStackCBx;
   TGComboBox *              fConfigStackCBx;
   TGComboBox *              fFitterCBx;
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
   TGCheckButton *           fShowColorCodeChkB;
   TGCheckButton *           fShowExtLegendChkB;
   TGCheckButton *           fUseStackedChkB;
   TGNumberEntry *           fEMinNmE;
   TGNumberEntry *           fEMaxNmE;
   TGNumberEntry *           fXMinNmE;
   TGNumberEntry *           fXMaxNmE;
   TGTextEntry *             fStackTableNameTxE;
   TGTextEntry *             fStackConfigNameTxE;
   TGLabel *                 fXMinLb;
   TGLabel *                 fXMaxLb;
   TGLabel *                 fStackDBTableLb;
   TGLabel *                 fStackConfigLb;
   TGLabel *                 fLinkSelLb;
   TGLabel *                 fLFitSpacerLb;
   TGLabel *                 fRFitSpacerLb;
   TLatex *                  fLtxAuth;
   TLatex *                  fLtxLink;

   //-- data selection tabs

   GuiNuDataSelectionTab *   fNuXSecTab;
   GuiElDataSelectionTab *   fElXSecTab;
   GuiSFDataSelectionTab *   fSFTab;

   //-- 'action' objects that handle some classes of GUI events

   GuiHelpHandler *          fHelpHandler;
   GuiDBHandler *            fDBaseHandler;
   GuiXmlFileHandler *       fXmlFileHandler;
   GuiStackHandler *         fStackHandler;
// GuiFitKernel *            fFitKernel;

   //-- other private date members

   DBConnection *            fDBC;
   NeuGenFitParams *         fNGFP;
   bool                      fPlotterShowIsOn;
   Spline *                  fSpline;
   NuVldConfig *             fMyConfig;        // nuvld session config
   UInt_t                    fW0;              // original window width
   UInt_t                    fH0;              // original window height

ClassDef(NuVldMainFrame, 0)
};

} // nuvld namespace
} // genie namespace

#endif

