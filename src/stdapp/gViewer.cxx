//____________________________________________________________________________
/*!

\program gViewer

\brief   

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created July 14, 2004

\cpright Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
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
