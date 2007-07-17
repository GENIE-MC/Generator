//____________________________________________________________________________
/*!

\class    genie::Renderer

\brief    Feynman diagram renderer ABC

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 07, 2004

*/
//____________________________________________________________________________

#include <iomanip>
#include <sstream>

#include <TRootEmbeddedCanvas.h>

#include "EVGCore/EventRecord.h"
#include "Viewer/Renderer.h"

using std::ostringstream;
using std::setprecision;
using std::string;

using namespace genie;

//______________________________________________________________________________
Renderer::Renderer() 
{

}
//______________________________________________________________________________
Renderer::~Renderer()
{

}
//______________________________________________________________________________
void Renderer::SetEmbeddedCanvas(TRootEmbeddedCanvas * ec)
{
  fEmbeddedCanvas = ec;
}
//______________________________________________________________________________
const char * Renderer::P4AsString(double E, double px, double py, double pz) 
{
  ostringstream p4;

  p4 << setprecision(3);
  p4 << "(" << E << ", " << px << ", " << py << ", " << pz << ")";

  return p4.str().c_str();
}
//______________________________________________________________________________

