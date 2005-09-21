//____________________________________________________________________________
/*!

\program  nuvld

\brief    Runs the NuValidator Graphical User Interface (GUI)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created October 05, 2004
*/
//____________________________________________________________________________

#include <TROOT.h>
#include <TApplication.h>
#include <TGClient.h>

#include "Messenger/Messenger.h"
#include "NuVldGUI/NuVldMainFrame.h"
#include "NuVldGUI/NuVldConfig.h"

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
