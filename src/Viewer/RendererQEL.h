//____________________________________________________________________________
/*!

\class    genie::RendererQEL

\brief    Draws the Feynman diagram for a QEL event

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  October 07, 2004

*/
//____________________________________________________________________________

#ifndef _RENDERER_QEL_H_
#define _RENDERER_QEL_H_

#include "Viewer/Renderer.h"

namespace genie {

class EventRecord;

class RendererQEL : public Renderer {

public:

   RendererQEL();
   virtual ~RendererQEL();

   void DrawDiagram (EventRecord * ev_rec);
};

}       // genie namespace

#endif  // _RENDERER_QEL_H_

