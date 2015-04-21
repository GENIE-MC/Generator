//____________________________________________________________________________
/*!

\program gvld_hadronz_test

\brief   Hadronization validation program
         (comparing GENIE models with neutrino bubble chamber data).

\syntax  gvld_hadronz_test -g genie_inputs.xml [-f fmt]

         Options:

          [] Denotes an optional argument.
          
          -f Output plot format (either `ps', `eps', `gif' or `root').
             If `ps' is selected, all plots are saved in a single document.
             If `eps', `gif' or `root' is selected, then each canvas is 
             stored separately. Default: `ps'.
          -g An input XML file for specifying a GHEP event list for each
             model to be considered in the hadronization benchmark test.
             For info on the XML file format see the GSimFiles class documentation.

\author  Tingjun Yang, Hugh Gallagher, Pauli Kehayias, Costas Andreopoulos 

\created March 1, 2009

\cpright Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <cstdlib>
#include <cassert>
#include <map>
#include <vector>
#include <string>

#include "Messenger/Messenger.h"
#include "Utils/StringUtils.h"
#include "Utils/Style.h"
#include "Utils/GSimFiles.h"
#include "Utils/CmdLnArgParser.h"
#include "validation/Hadronization/HadPlots.h"
#include "validation/Hadronization/HadPlotter.h"

using namespace std;
using namespace genie;
using namespace genie::utils;
using namespace genie::mc_vs_data;

// prototypes
void LoadFilesAndBookPlots (void);
void Analyze               (void);
void Plot                  (void);
void End                   (void);
void GetCommandLineArgs    (int argc, char** argv);
void PrintSyntax           (void);

// globals & user inputs
vector<HadPlots *> gHP;
GSimFiles          gOptGenieInputs(false,10);
string             gFmt = "ps";

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs(argc, argv);
  style::SetDefaultStyle();
  LoadFilesAndBookPlots();
  Analyze();
  Plot();
  End();

  LOG("gvldtest", pNOTICE) << "Done!";

  return 0;
}
//____________________________________________________________________________
void LoadFilesAndBookPlots(void)
{
  vector<HadPlots *>::iterator hpvit;
  int nfiles = 0;

  // loop over models considered at the current validation test
  for(int imodel=0; imodel < gOptGenieInputs.NModels(); imodel++) {

    string model_name = gOptGenieInputs.ModelTag(imodel);
    vector<string> & event_filenames = gOptGenieInputs.EvtFileNames(imodel);

      LOG("gvldtest", pNOTICE)
            << "Booking plots for model: " << model_name;

      HadPlots * hadplot = new HadPlots(model_name);
      
      // include all data files for current model
      vector<string>::const_iterator file_iter = event_filenames.begin();
      for( ; file_iter != event_filenames.end(); ++file_iter) {
        
        if(nfiles==kMaxFiles)  {
          LOG("gvldtest",pFATAL) 
              << "Number of Input Files greater than Maximum: " << 
              kMaxFiles;
          gAbortingInErr=true;
          exit(1) ; 
        }

         string filename = *file_iter;
         LOG("gvldtest", pNOTICE)
            << " Loading data from file:.....: " << filename;
         hadplot->LoadData(filename);
         nfiles++;
      }// file_iter

      // store
      gHP.push_back(hadplot);
  }
}
//____________________________________________________________________________
void Analyze(void)
{
  vector<HadPlots *>::iterator hpvit = gHP.begin();
  for( ; hpvit != gHP.end(); ++hpvit) {
    HadPlots * curr_model = *hpvit;
    curr_model->Analyze();
  }
}
//____________________________________________________________________________
void Plot(void)
{
  // plot all bubble chamber data and start superimposing model predictions
  HadPlotter plotter(gFmt);
  vector<HadPlots *>::iterator hpvit = gHP.begin();
  for( ; hpvit != gHP.end(); ++hpvit) {
    HadPlots * curr_model = *hpvit;
    plotter.AddPlots(*curr_model);
  }
  plotter.ShowPlots();
}
//____________________________________________________________________________
void End(void)
{
  vector<HadPlots *>::iterator hpvit = gHP.begin();
  for( ; hpvit != gHP.end(); ++hpvit) {
    HadPlots * curr_model = *hpvit;
    if(curr_model) {
       delete curr_model;
       curr_model = 0;
    }
  }
  gHP.clear();
}
//____________________________________________________________________________
void GetCommandLineArgs(int argc, char** argv)
{
  CmdLnArgParser parser(argc,argv);

  // help?
  bool help = parser.OptionExists('h');
  if(help) {
      PrintSyntax();
      exit(0);
  }

  // get GENIE inputs
  if(parser.OptionExists('g')) {
     string inputs = parser.ArgAsString('g');
     bool ok = gOptGenieInputs.LoadFromFile(inputs);
     if(!ok) { 
        LOG("gvldtest", pFATAL) 
          << "Could not read validation program inputs from: " << inputs;
        PrintSyntax();
        gAbortingInErr=true;
        exit(1);
     }
  } 

  // output plot format
  if(parser.OptionExists('f')) {
    gFmt = parser.ArgAsString('f');
  } else {
    gFmt = "ps";
  }

  if(gOptGenieInputs.NModels()==0) {
    LOG("gvldtest", pFATAL) << "** No input model data to analyze";
    gAbortingInErr=true;
    exit(1);
  }

  LOG("gvldtest", pFATAL) << "Input data: ";  
  LOG("gvldtest", pFATAL) << gOptGenieInputs;
}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("vldtest", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << " gvld_hadronz_test -g genie_inputs.xml [-f format] \n";
}
//____________________________________________________________________________
