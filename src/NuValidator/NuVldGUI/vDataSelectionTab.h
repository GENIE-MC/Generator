//_____________________________________________________________________________
/*!

\class    genie::nuvld::vDataSelectionTab

\brief    Neutrino Cross Section Data Selection Graphical Tab

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  November 01, 2004
*/
//_____________________________________________________________________________

#ifndef _NEUTRINO_DATA_SELECTION_TAB_H_
#define _NEUTRINO_DATA_SELECTION_TAB_H_

#include <RQ_OBJECT.h>

#include "NuVldGUI/DataSelectionDialog.h"

class TGFrame;
class TGListBox;
class TGButton;
class TGNumberEntry;
class TGLabel;

namespace genie {
namespace nuvld {

class DBConnection;

class vDataSelectionTab : public DataSelectionDialog {

RQ_OBJECT("vDataSelectionTab")

public:
   vDataSelectionTab(TGMainFrame * main, DBConnection * db = 0);
   virtual ~vDataSelectionTab();

   TGCompositeFrame * Create(TGCompositeFrame * tf, int width, int height);

   string BundleSelectionsInString     (void);
   string BundleKeyListInString        (void);
   string BundleCutsInString           (void);
   string BundleDrawOptInString        (void);
   void   ResetSelections              (void);
   void   SelectAllExp                 (void);
   void   SelectAllXSec                (void);
   void   SelectAllProbes              (void);
   void   SelectAllTargets             (void);
   string ReadXSecErrorListbox         (void);
   void   PopupNuDataSelectionDialog   (void);
   void   PopupNuMeasurementListDialog (void);

private:

   //-- GUI widgets

   TGMainFrame *         fMain;
   TGCompositeFrame *    fTabNuSql;
   TGLayoutHints *       fNuSqlTabLt;
   TGMatrixLayout *      fEnergyMatrixLt;
   TGGroupFrame *        fEnergyGrpFrm;
   TGGroupFrame *        fNuXSecErrGrpFrm;
   TGGroupFrame *        fNuExpGrpFrm;
   TGGroupFrame *        fNuXSecGrpFrm;
   TGGroupFrame *        fNuInitStateGrpFrm;
   TGListBox *           fNuXSecErrLBx;
   TGListBox *           fNuExpLBx;
   TGListBox *           fNuProcLBx;
   TGListBox *           fNuTypeLBx;
   TGListBox *           fNuTgtLBx;
   TGTextButton *        fShowFullNuDialogTBtn;
   TGTextButton *        fShowExpertNuDialogTBtn;
   TGCheckButton *       fAllNuExpChkB;
   TGCheckButton *       fAllNuProcChkB;
   TGCheckButton *       fAllNuTypesChkB;
   TGCheckButton *       fAllNuTgtChkB;
   TGNumberEntry *       fEMinNmE;
   TGNumberEntry *       fEMaxNmE;
   TGLabel *             fMinELb;
   TGLabel *             fMaxELb;
   TGLabel *             fNuTabBtnSpacerLb;
   TGCheckButton *       fScaleWithEvChkB;

   DBConnection *        fDBC;            ///< NuVld database connection

   bool                  fPopupDialogLAM; ///< LAM (Look-At-Me) flag for popup dialog
   DataSelectionDialog * fPopupDialog;    ///< popup delection dialog interface

   ClassDef(vDataSelectionTab, 0)
};

} // nuvld namespace
} // genie namespace

#endif // _NEUTRINO_DATA_SELECTION_TAB_H_

