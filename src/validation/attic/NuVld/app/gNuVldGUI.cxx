//____________________________________________________________________________
/*!

\program  gNuVldGUI

\brief    Runs the NuValidator Graphical User Interface 

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created October 05, 2004
*/
//____________________________________________________________________________

#include <TROOT.h>
#include <TApplication.h>
#include <TGClient.h>

#include "Messenger/Messenger.h"
#include "ValidationTools/NuVld/NuVldMainFrame.h"
#include "ValidationTools/NuVld/NuVldConfig.h"

using genie::Messenger;
using namespace genie::nuvld;

int main(int argc, char ** argv)
{
  Messenger * msg = Messenger::Instance();
  msg->SetPriorityLevel("NuVld", pINFO);

  TApplication vld_gui("VLD", &argc, argv);

  NuVldConfig config;
  config.AutoDetect();

  NuVldMainFrame main_window(gClient->GetRoot(), 800, 1100, config);
  vld_gui.Run();

  return 0;
}
