//____________________________________________________________________________
/*!

\class    genie::masterclass::MCTruthDisplay

\brief    Display MC truth info

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  November 30, 2007

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _MC_TRUTH_DISPLAY_H_
#define _MC_TRUTH_DISPLAY_H_

#include <TRootEmbeddedCanvas.h>
#include <TGTextEdit.h>

#include "Framework/EventGen/EventRecord.h"

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
