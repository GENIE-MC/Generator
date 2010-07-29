//____________________________________________________________________________
/*!

\program gvld_hadronz_test

\brief   Hadronization validation program
         (comparing GENIE models with neutrino bubble chamber data).

\syntax  gvld_hadronz_test -g genie_inputs.xml [-f fmt]

         Options:

          [] Denotes an optional argument.
          
          -f Output plot format (0: eps, 1: gif) [default: gif]          
          -g An input XML file for specifying a GHEP event list for each
             model to be considered in the hadronization benchmark test.

             The input file should look like:

             <?xml version="1.0" encoding="ISO-8859-1"?>
             <vld_inputs>
               <model name="a_model_name">
                 <evt_file format="ghep"> /model_1/evtfile0.root </evt_file>
                 <evt_file format="ghep"> /model_1/evtfile1.root </evt_file>
                 <evt_file format="ghep"> /model_1/evtfile2.root </evt_file>
                 ...
               </model>

               <model name="another_model_name">
                 <evt_file format="ghep"> /model_2/evtfile0.root </evt_file>
                 <evt_file format="ghep"> /model_2/evtfile1.root </evt_file>
                 <evt_file format="ghep"> /model_2/evtfile2.root </evt_file>
                 ...
               </model>
               ...
             </vld_inputs>


\author  Tingjun Yang 

\created March 1, 2009

\cpright Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
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
#include "Utils/VldTestInputs.h"
#include "Utils/CmdLnArgParser.h"
#include "ValidationTools/Hadronization/HadPlots.h"
#include "ValidationTools/Hadronization/HadPlotter.h"

using namespace std;
using namespace genie;
using namespace genie::utils;
using namespace genie::utils::vld;
using namespace genie::vld_hadronization;

// prototypes
void LoadFilesAndBookPlots (void);
void Analyze               (void);
void Plot                  (void);
void End                   (void);
void GetCommandLineArgs    (int argc, char** argv);
void PrintSyntax           (void);

// globals & user inputs
vector<HadPlots *> gHP;
VldTestInputs      gOptGenieInputs(false,10);;
int                gFmt = 1;

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs(argc, argv);
  style::SetDefaultStyle();
  LoadFilesAndBookPlots();
  Analyze();
  Plot();
  End();

  LOG("vldtest", pNOTICE) << "Done!";

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

      LOG("vldtest", pNOTICE)
            << "Booking plots for model: " << model_name;

      HadPlots * hadplot = new HadPlots(model_name);
      
      // include all data files for current model
      vector<string>::const_iterator file_iter = event_filenames.begin();
      for( ; file_iter != event_filenames.end(); ++file_iter) {
        
        if(nfiles==kMaxFiles)  {
          LOG("vldtest",pFATAL) 
              << "Number of Input Files greater than Maximum: " << 
              kMaxFiles;
          gAbortingInErr=true;
          exit(1) ; 
        }

         string filename = *file_iter;
         LOG("vldtest", pNOTICE)
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
  bool in_eps = (gFmt==0);
  HadPlotter plotter(in_eps);
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
     int format = parser.ArgAsInt('f');
     if(format==0 || format==1) {
       gFmt=format;
     }
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
    << " gvld_hadronz_test -g genie_inputs.xml\n [-f format]";
}
//____________________________________________________________________________
