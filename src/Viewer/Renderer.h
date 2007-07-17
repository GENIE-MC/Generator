//____________________________________________________________________________
/*!

\class    genie::Renderer

\brief    Feynman diagram renderer ABC

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 07, 2004

*/
//____________________________________________________________________________

#ifndef _RENDERER_H_
#define _RENDERER_H_

class TRootEmbeddedCanvas;

namespace genie {

class EventRecord;

class Renderer {

public:

   virtual void DrawDiagram (EventRecord * ev_rec) = 0;

   void         SetEmbeddedCanvas (TRootEmbeddedCanvas * ec);
   const char * P4AsString(double E, double px, double py, double pz);

protected:

   Renderer();
   virtual ~Renderer();

   TRootEmbeddedCanvas * fEmbeddedCanvas;
};

}       // genie namespace

#endif  // _RENDERER_H_

