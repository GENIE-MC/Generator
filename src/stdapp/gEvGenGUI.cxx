//____________________________________________________________________________
/*!

\program gevgen_gui

\brief   

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created July 14, 2004

\cpright Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <TROOT.h>
#include <TApplication.h>
#include <TGClient.h>

#include "Viewer/GViewerMainFrame.h"

using namespace genie;
using namespace genie::gview;

int main(int argc, char ** argv)
{
  TApplication ggui("GENIE", &argc, argv);

  GViewerMainFrame gvmf(gClient->GetRoot(), 700, 350);
  ggui.Run();

  return 0;
}

