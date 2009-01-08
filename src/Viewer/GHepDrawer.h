//____________________________________________________________________________
/*!

\class    genie::GHepDrawer

\brief    Draws a generated event

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  November 30, 2007

*/
//____________________________________________________________________________

#ifndef _GHEP_DRAWER_H_
#define _GHEP_DRAWER_H_

class TRootEmbeddedCanvas;

namespace genie {

class EventRecord;

class GHepDrawer {

public:
   GHepDrawer();
  ~GHepDrawer();

   void SetEmbeddedCanvas (TRootEmbeddedCanvas * ec);
   void Draw              (EventRecord * ev_rec);

private:
   TRootEmbeddedCanvas * fEmbeddedCanvas;
};

}       // genie namespace
#endif  // _RENDERER_QEL_H_

