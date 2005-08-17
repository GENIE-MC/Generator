//____________________________________________________________________________
/*!

\class    genie::GHepPrinter

\brief    Prints the GHep record

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  October 07, 2004

*/
//____________________________________________________________________________

#ifndef _GHEP_PRINTER_H_
#define _GHEP_PRINTER_H_

#include <TGTextEdit.h>

#include "EVGCore/EventRecord.h"

namespace genie {

class GHepPrinter {

public:

   GHepPrinter();
   virtual ~GHepPrinter();

   void SetTextEdit (TGTextEdit *  txedit) { fGHep = txedit; }
   void Print       (EventRecord * ev_rec);

private:

  TGTextEdit * fGHep;
};

}       // genie namespace

#endif  // _GHEP_PRINTER_H_

