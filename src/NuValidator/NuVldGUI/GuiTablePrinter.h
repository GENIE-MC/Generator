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

#include "DBUtils/DBTable.h"
#include "DBUtils/vXSecTableRow.h"
#include "DBUtils/eDiffXSecTableRow.h"
#include "DBUtils/SFTableRow.h"

using std::string;

namespace genie {
namespace nuvld {


class GuiTablePrinter : public TObject {

  friend class NuVldMainFrame;
  
private:

  GuiTablePrinter();
  virtual ~GuiTablePrinter();

  //-- table printing methods
    
  void PrintTable (genie::nuvld::DBTable<vXSecTableRow> *     table) const;
  void PrintTable (genie::nuvld::DBTable<eDiffXSecTableRow> * table) const;
  void PrintTable (genie::nuvld::DBTable<SFTableRow> *        table) const;

  //-- options

  void ScaleXSecWithEnergy(bool tf) { fScaleXSecWithEnergy = tf; }

  //-- table formatting methods  
        
  string PrintNuXSecTableHeader          (void) const;
  string PrintElDiffXSecTableHeader      (void) const;
  string PrintSFTableHeader              (void) const;
  string PrintElDiffXSecTableHeaderUnits (void) const;
  string PrintTableRowAsString           (const vXSecTableRow * row) const;
  string PrintTableRowAsString           (const eDiffXSecTableRow * row) const;
  string PrintTableRowAsString           (const SFTableRow * row) const;
  string PrintTableSeparator             (int n) const;

  bool fScaleXSecWithEnergy;

ClassDef(GuiTablePrinter, 0)
};

} // nuvld namespace
} // genie namespace

#endif // _GUI_TABLE_PRINTER_H_

