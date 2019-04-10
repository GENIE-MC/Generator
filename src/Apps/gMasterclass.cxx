//____________________________________________________________________________
/*!

\program gmstcl

\brief   GENIE Masterclass GUI

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created July 14, 2004

\cpright Copyright (c) 2003-2019, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <TROOT.h>
#include <TApplication.h>
#include <TGClient.h>

#include "Tools/Masterclass/GNuMcMainFrame.h"

using namespace genie;
using namespace genie::masterclass;

int main(int argc, char ** argv)
{
  TApplication ggui("GENIE - Neutrino Masterclass App", &argc, argv);
  GNuMcMainFrame gvmf(gClient->GetRoot(), 700, 350);
  ggui.Run();
  return 0;
}


