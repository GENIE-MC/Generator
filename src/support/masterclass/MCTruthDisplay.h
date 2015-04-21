//____________________________________________________________________________
/*!

\class    genie::masterclass::MCTruthDisplay

\brief    Display MC truth info

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  November 30, 2007

*/
//____________________________________________________________________________

#ifndef _MC_TRUTH_DISPLAY_H_
#define _MC_TRUTH_DISPLAY_H_

#include <TRootEmbeddedCanvas.h>
#include <TGTextEdit.h>

#include "EVGCore/EventRecord.h"

namespace genie {
 namespace masterclass {

   class MCTruthDisplay {
   public:
     MCTruthDisplay(TRootEmbeddedCanvas * ec=0, TGTextEdit * gtx=0);
    ~MCTruthDisplay();
     void DrawDiagram      (EventRecord * event);
     void PrintEventRecord (EventRecord * event);
   private:
     TRootEmbeddedCanvas * fEmbeddedCanvas;
     TGTextEdit *          fGTxt;
   };

 } // masterclass namespace
} // genie namespace

#endif  // _MC_TRUTH_DISPLAY_H_

