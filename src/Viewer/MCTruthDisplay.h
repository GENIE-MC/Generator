//____________________________________________________________________________
/*!

\class    genie::gview::MCTruthDisplay

\brief    Display MC truth info

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  November 30, 2007

*/
//____________________________________________________________________________

#ifndef _MC_TRUTH_DISPLAY_H_
#define _MC_TRUTH_DISPLAY_H_

class TRootEmbeddedCanvas;
class TGTextEdit;

namespace genie {

 class EventRecord;

 namespace gview {

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

 } // gview namespace
} // genie namespace

#endif  // _MC_TRUTH_DISPLAY_H_

