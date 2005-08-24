//____________________________________________________________________________
/*!

\program gViewer

\brief   

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created July 14, 2004
*/
//____________________________________________________________________________

#include <TROOT.h>
#include <TApplication.h>
#include <TGClient.h>

#include "Viewer/GenieViewer.h"

using namespace genie;

int main(int argc, char ** argv)
{
  TApplication genie_gui("GENIE", &argc, argv);

  GenieViewer main_window(gClient->GetRoot(), 700, 350);
  genie_gui.Run();

  return 0;
}
