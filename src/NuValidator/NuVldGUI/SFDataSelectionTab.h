//_____________________________________________________________________________
/*!

\class    genie::nuvld::SFDataSelectionTab

\brief    Structure Function Data Selection Graphical Tab

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  November 01, 2004
*/
//_____________________________________________________________________________

#ifndef _STRUCTURE_FUNCTION_DATA_SELECTION_TAB_H_
#define _STRUCTURE_FUNCTION_DATA_SELECTION_TAB_H_

#include <RQ_OBJECT.h>

#include "NuVldGUI/DataSelectionDialog.h"

class TGFrame;
class TGListBox;
class TGComboBox;
class TGButton;
class TGNumberEntry;

using std::string;

namespace genie {
namespace nuvld {

class DBConnection;

class SFDataSelectionTab : public DataSelectionDialog {

RQ_OBJECT("SFDataSelectionTab")

public:
   SFDataSelectionTab(DBConnection * db = 0);
   virtual ~SFDataSelectionTab();
       
   TGCompositeFrame * Create(TGCompositeFrame * tf, int width, int height);

   string BundleSelectionsInString (void);   
   string BundleKeyListInString    (void);
   string BundleCutsInString       (void);
   string BundleDrawOptInString    (void);
   void   ResetSelections          (void);
   void   SelectAllExp             (void);
   void   SelectAllProbes          (void);
   void   SelectAllTargets         (void);
   //void   SFLoadx                  (void);

private:

   //-- GUI widgets

   TGCompositeFrame *  fTabSFSql;
   TGGroupFrame *      fErrGrpFrm;
   TGGroupFrame *      fExpGrpFrm;
   TGGroupFrame *      fSFGrpFrm;
   TGGroupFrame *      fKineGrpFrm;
   TGGroupFrame *      fInitStateGrpFrm;
   TGGroupFrame *      fPlotVarGrpFrm;
   TGListBox *         fErrLBx;
   TGListBox *         fExpLBx;
   TGListBox *         fSFLBx;
   TGListBox *         fProbeLBx;
   TGListBox *         fTgtLBx;
   TGListBox *         fRLBx;
   TGComboBox *        fPlotVarCBx;
   //TGListBox *         fSFxLBx;
   //TGTextButton *      fSFLoadxTBtn;
   TGCheckButton *     fAllExpChkB;
   TGCheckButton *     fAllProbesChkB;
   TGCheckButton *     fAllTgtChkB;
   TGNumberEntry *     fMinQ2NmE;
   TGNumberEntry *     fMaxQ2NmE;
   TGNumberEntry *     fMinXNmE;
   TGNumberEntry *     fMaxXNmE;
   DBConnection *      fDBC;
      
   ClassDef(SFDataSelectionTab, 0)
};

} // nuvld namespace
} // genie namespace

#endif // _STRUCTURE_FUNCTION_DATA_SELECTION_TAB_H_

