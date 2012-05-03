//_____________________________________________________________________________
/*!

\class    genie::nuvld::GuiTablePrinter

\brief    Responds to GUI events for printing DBTable<T>'s

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  August 26, 2004
*/
//_____________________________________________________________________________

#ifndef _GUI_TABLE_PRINTER_H_
#define _GUI_TABLE_PRINTER_H_

#include <string>

#include <TObject.h>

#include  "ValidationTools/NuVld/DBTable.h"
#include  "ValidationTools/NuVld/DBNuXSecTableRow.h"
#include  "ValidationTools/NuVld/DBElDiffXSecTableRow.h"
#include  "ValidationTools/NuVld/DBSFTableRow.h"

using std::string;

namespace genie {
namespace nuvld {

class GuiTablePrinter : public TObject {

  friend class NuVldMainFrame;
  
private:
  GuiTablePrinter();
  virtual ~GuiTablePrinter();

  //-- table printing methods
    
  void PrintTable (genie::nuvld::DBTable<DBNuXSecTableRow> *     table) const;
  void PrintTable (genie::nuvld::DBTable<DBElDiffXSecTableRow> * table) const;
  void PrintTable (genie::nuvld::DBTable<DBSFTableRow> *         table) const;

  //-- options

  void ScaleXSecWithEnergy(bool tf) { fScaleXSecWithEnergy = tf; }

  //-- table formatting methods  
        
  string PrintNuXSecTableHeader          (void) const;
  string PrintElDiffXSecTableHeader      (void) const;
  string PrintSFTableHeader              (void) const;
  string PrintElDiffXSecTableHeaderUnits (void) const;
  string PrintTableRowAsString           (const DBNuXSecTableRow * row) const;
  string PrintTableRowAsString           (const DBElDiffXSecTableRow * row) const;
  string PrintTableRowAsString           (const DBSFTableRow * row) const;
  string PrintTableSeparator             (int n) const;

  bool fScaleXSecWithEnergy;

ClassDef(GuiTablePrinter, 0)
};

} // nuvld namespace
} // genie namespace

#endif // _GUI_TABLE_PRINTER_H_

