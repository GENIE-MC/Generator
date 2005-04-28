//_____________________________________________________________________________
/*!

\class    genie::nuvld::eDataSelectionTab

\brief    Electron Cross Section Data Selection Graphical Tab

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  November 01, 2004
*/
//_____________________________________________________________________________

#ifndef _ELECTRON_DATA_SELECTION_TAB_H_
#define _ELECTRON_DATA_SELECTION_TAB_H_

#include <vector>

#include <RQ_OBJECT.h>

#include "NuVldGUI/DataSelectionDialog.h"

class TGListBox;
class TGComboBox;
class TGButton;
class TGNumberEntry;
class TGLabel;

using std::vector;

namespace genie {
namespace nuvld {

class DBConnection;

class eDataSelectionTab : public DataSelectionDialog {

RQ_OBJECT("eDataSelectionTab")

public:
   eDataSelectionTab(DBConnection * db = 0);
   virtual ~eDataSelectionTab();
      
   TGCompositeFrame * Create(TGCompositeFrame * tf, int width, int height);

   string BundleSelectionsInString (void);
   string BundleKeyListInString    (void);
   string BundleCutsInString       (void);
   string BundleDrawOptInString    (void);
   void   ResetSelections          (void);
   void   SelectAllExp             (void);
   void   SelectAllTargets         (void);
   
private:

   //-- GUI widgets

   TGCompositeFrame *        fTabElSql;   
   TGGroupFrame *            fElExpGrpFrame;
   TGGroupFrame *            fElTgGrpFrm;
   TGGroupFrame *            fElDrawXGrpFrm;
   TGListBox *               fElExpLBx;
   TGListBox *               fElTgtLBx;
   TGComboBox *              fElDrawXCBx;
   TGCheckButton *           fAllElExpChkB;
   TGCheckButton *           fAllElTgtChkB;
   vector<TGMatrixLayout * > fElVarRangeLt;
   vector<TGGroupFrame * >   fElVarRangeGrpFrm;   
   vector<TGNumberEntry * >  fElVarMinNmEV;
   vector<TGNumberEntry * >  fElVarMaxNmEV;      
   DBConnection *            fDBC;
      
   ClassDef(eDataSelectionTab, 0)
};

} // nuvld namespace
} // genie namespace

#endif // _ELECTRON_DATA_SELECTION_TAB_H_

