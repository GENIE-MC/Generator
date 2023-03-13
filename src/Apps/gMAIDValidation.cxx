//____________________________________________________________________________
/*!

\program gMAIDValidation

\brief   Generates MAID Form Factors for all fitted resonances

\syntax  gMAIDValidation [-o output_name.root]

         -o :
          Specifies a name to be used in the output files.
          Default: maid_validation.root

\author  Julia Tena Vidal <jtenavidal \at tauex.tau.ac.il>
Tel Aviv University

\created 13 Mar 2023

\cpright Copyright (c) 2023-2033, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         
*/
//____________________________________________________________________________

#include <vector>
#include <string>

#include <TFile.h>
#include <TMath.h>
#include <TPostScript.h>
#include <TPavesText.h>
#include <TText.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TPaletteAxis.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/PartonDistributions/PDF.h"
#include "Framework/Utils/CmdLnArgParser.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/Style.h"

using namespace std;
using namespace genie;
using namespace genie::utils;

// globals
string        gOptOutFile = "pdf_comp"; // -o argument

// function prototypes
void GetCommandLineArgs (int argc, char ** argv);
void MakePlots     (void);

//___________________________________________________________________
int main(int argc, char ** argv)
{
  utils::style::SetDefaultStyle();

  GetCommandLineArgs (argc,argv);   // Get command line arguments
  MakePlots();   // Produce all output plots and fill output n-tuple

  return 0;
}
//_________________________________________________________________________________
void MakePlots (void)
{

}
//_________________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  // necessary for setting from whence it gets ModelConfiguration.xml
  RunOpt::Instance()->ReadFromCommandLine(argc,argv);

  CmdLnArgParser parser(argc,argv);

  if(parser.OptionExists('o')){
    gOptOutFile = parser.Arg('o');
  }

}
//_________________________________________________________________________________

