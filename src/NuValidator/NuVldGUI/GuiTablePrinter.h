//_____________________________________________________________________________
/*!

\class    genie::nuvld::GuiTablePrinter

\brief    Responds to GUI events for printing DBTable<T>'s

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

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

using std::string;

namespace genie {
namespace nuvld {


class GuiTablePrinter : public TObject {

  friend class NuVldMainFrame;
  
private:

  GuiTablePrinter();
  virtual ~GuiTablePrinter();

  //-- table printing methods
    
  void PrintXSecTable (genie::nuvld::DBTable<vXSecTableRow> *     table) const;
  void PrintXSecTable (genie::nuvld::DBTable<eDiffXSecTableRow> * table) const;

  //-- drawing options

  void ScaleXSecWithEnergy(bool tf) { fScaleXSecWithEnergy = tf; }

  //-- table formatting methods  
        
  string PrintNuXSecTableHeader          (void) const;
  string PrintElDiffXSecTableHeader      (void) const;
  string PrintElDiffXSecTableHeaderUnits (void) const;
  string PrintNuXSecTableSeparator       (void) const;
  string PrintElDiffXSecTableSeparator   (void) const;
  string PrintXSecTableRowAsString       (const vXSecTableRow * row) const;
  string PrintXSecTableRowAsString       (const eDiffXSecTableRow * row) const;

  bool fScaleXSecWithEnergy;

ClassDef(GuiTablePrinter, 0)
};

} // nuvld namespace
} // genie namespace

#endif // _GUI_TABLE_PRINTER_H_

