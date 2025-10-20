//____________________________________________________________________________
/*!
\program gcalchedisdiffxsec
\brief   GENIE utility program calculating differential cross sections from 
         HEDIS package.
         Syntax :
           gcalchedisdiffxsec -p nu -t tgt -o root_file
                              -x table_type
                              --tune genie_tune
                              --event-generator-list list_name
         Note :
           [] marks optional arguments.
           <> marks a list of arguments out of which only one can be
              selected at any given time.
         Options :
           -p
              the neutrino pdg code
           -t
              the target pdg code (format: 10LZZZAAAI)
           -x
              path to ascii file with x values
           -y
              path to ascii file with y values
           -e
              path to ascii file with energy values
           -o
              output ROOT file name
           --tune
              Specifies a GENIE comprehensive neutrino interaction model tune.
           --event-generator-list
              List of event generators to load in event generation drivers.
              [default: "Default"].
           --message-thresholds
              Allows users to customize the message stream thresholds.
              The thresholds are specified using an XML file.
              See $GENIE/config/Messenger.xml for the XML schema.
        ***  See the User Manual for more details and examples. ***
\author  Alfonso Garcia <alfonsog \at nikhef.nl>
         NIKHEF (Amsterdam)
\created May 26, 2021
\cpright Copyright (c) 2003-2025, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include "TMath.h"
#include "TSystem.h"
#include "TTree.h"
#include "TFile.h"

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/EventGen/InteractionList.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/GEVGDriver.h"
#include "Framework/Utils/AppInit.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/CmdLnArgParser.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Interaction/InitialState.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Conventions/Units.h"

#include <fstream>
#include <iostream>

using namespace genie;

void   PrintSyntax          (void);
void   DecodeCommandLine    (int argc, char * argv[]);

vector<double> ReadListFromPath (string path);

const double epsilon = 1e-5;

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//INPUT ARGUMENTS
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

int fNu             = -1;
int fTgt            = -1;
string fPathXlist   = "";
string fPathYlist   = "";
string fPathElist   = "";
string fOutFileName = "";
bool fSaveAll       = false;

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//MAIN PROGRAM
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

int main(int argc, char ** argv) {

  DecodeCommandLine(argc,argv); 

  RunOpt::Instance()->BuildTune();
  utils::app_init::MesgThresholds(RunOpt::Instance()->MesgThresholdFiles());

  string fChannel = RunOpt::Instance()->EventGeneratorList();

  GEVGDriver evg_driver;
  InitialState init_state(fTgt, fNu);
  evg_driver.SetEventGeneratorList(fChannel);
  evg_driver.Configure(init_state);

  vector<double> ve = ReadListFromPath(fPathElist);
  vector<double> vx = ReadListFromPath(fPathXlist);
  vector<double> vy = ReadListFromPath(fPathYlist);

  double widthx = (TMath::Log10(vx[1])-TMath::Log10(vx[0]));

  LOG("gcalchedisdiffxsec", pINFO) << widthx;

  int quark;
  double ei;
  double bx;
  double by;
  double dxsec;
  
  TString treename = Form("diffxsec_nu%d_tgt%d_%s",fNu,fTgt,fChannel.c_str());

  LOG("gcalchedisdiffxsec", pINFO) << treename;

  TTree * tree = new TTree(treename,treename);
  tree->Branch( "Quark",    &quark, "Quark/I"    );
  tree->Branch( "Ei",       &ei,    "Ei/D"       );
  tree->Branch( "Bx",       &bx,    "Bx/D"       );
  tree->Branch( "By",       &by,    "By/D"       );
  tree->Branch( "DiffXsec", &dxsec, "DiffXsec/D" );

  const InteractionList * intlst   = evg_driver.Interactions();

  InteractionList::const_iterator intliter;
  for(intliter = intlst->begin(); intliter != intlst->end(); ++intliter) {

    const Interaction * interaction = *intliter;

    quark = 10000 * interaction->InitState().Tgt().HitQrkPdg();
    if (!interaction->InitState().Tgt().HitSeaQrk()) {
      if ( quark>0 ) quark += 100;
      else           quark -= 100;
    }
    quark += 1 * interaction->ExclTag().FinalQuarkPdg();

    LOG("gcalchedisdiffxsec", pINFO) << "Current interaction: " << interaction->AsString();

    const XSecAlgorithmI * xsec_alg = evg_driver.FindGenerator(interaction)->CrossSectionAlg();

    for ( unsigned i=0; i<ve.size(); i++ ) {

      ei = ve[i];

      LOG("gcalchedisdiffxsec", pINFO) << "Energy: " << ei << " [GeV]";

      TLorentzVector p4(0,0,ei,ei);
      interaction->InitStatePtr()->SetProbeP4(p4);

      for ( unsigned j=0; j<vy.size(); j++ ) {

        double z = vy[j]; // z = (eo-em)/(ei-em)

        by = (1-z)*(ei-10.)/ei; // y = 1 - eo/ei;

        LOG("gcalchedisdiffxsec", pDEBUG) << "  z: " << z << "  y: " << by;

        if      ( by==1. ) by -= 1e-4; 
        else if ( by==0. ) {
          by = 5./ei;
          double by_prev = (1-vy[j-1])*(ei-10.)/ei;
          LOG("gcalchedisdiffxsec", pDEBUG) << "  by: " << by << "  by_prev: " << by_prev << "  " << vy[j-1];
          if (by>by_prev) by = 1e-4;
        }
        
        if ( by<0 || by>1 ) continue;
        
        interaction->KinePtr()->Sety(by);
        
        if (fSaveAll) {
          for ( unsigned k=0; k<vx.size(); k++ ) {
            bx = vx[k];
            interaction->KinePtr()->Setx(bx);
            utils::kinematics::UpdateWQ2FromXY(interaction);
            dxsec = xsec_alg->XSec(interaction, kPSxyfE) / units::cm2;
            LOG("gcalchedisdiffxsec", pDEBUG) << "x: " << bx << "  y: " << by << " -> d2sigmadxdy[E=" << ei << "GeV] = " << dxsec << " cm2";
            tree->Fill();
          }
        }
        else {
          dxsec = 0.;
          for ( unsigned k=0; k<vx.size(); k++ ) {
            bx = vx[k];
            interaction->KinePtr()->Setx(bx);
            utils::kinematics::UpdateWQ2FromXY(interaction);
            dxsec += xsec_alg->XSec(interaction, kPSxyfE) * bx;
          }
          dxsec *=  widthx*TMath::Log(10) / units::cm2;
          LOG("gcalchedisdiffxsec", pDEBUG) << "  y: " << by << " -> d2sigmady[E=" << ei << "GeV] = " << dxsec << " cm2";
          tree->Fill();
        }
        
      }
    
    }

  }

  TFile * outfile = new TFile(fOutFileName.c_str(),"RECREATE");
  tree->Write(tree->GetName());
  outfile->Close();

}


//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//READ LIST
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
vector<double> ReadListFromPath(string path) {

  vector<double> list;

  std::ifstream infile(path.c_str());

  //check to see that the file was opened correctly:
  if (!infile.is_open()) {
    LOG("gcalchedisdiffxsec", pFATAL) << "There was a problem opening the input file!";
    exit(1);//exit or do additional error checking
  }

  double val = 0.;
  while (infile >> val) list.push_back(val);

  return list;

}


//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//INPUT PARSER
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

void DecodeCommandLine(int argc, char * argv[]) {

  // Common run options. Set defaults and read.
  RunOpt::Instance()->ReadFromCommandLine(argc,argv);

  // Parse run options for this app
  CmdLnArgParser parser(argc,argv);

  if( parser.OptionExists('p') ){ 
    fNu = parser.ArgAsInt('p');
    LOG("gcalchedisdiffxsec", pINFO) << "Probe = " << fNu;
  }          
  else {
    LOG("gcalchedisdiffxsec", pFATAL) << "Unspecified input neutrino type!";
    PrintSyntax();
    exit(1);
  }

  if( parser.OptionExists('t') ){ 
    fTgt = parser.ArgAsInt('t');
    LOG("gcalchedisdiffxsec", pINFO) << "Target = " << fTgt;
  }          
  else {
    LOG("gcalchedisdiffxsec", pFATAL) << "Unspecified input target type!";
    PrintSyntax();
    exit(1);
  }
  
  if( parser.OptionExists('x') ){ 
    fPathXlist = parser.ArgAsString('x');
    LOG("gcalchedisdiffxsec", pINFO) << "PathXlist = " << fPathXlist;
  }          
  else {
    LOG("gcalchedisdiffxsec", pFATAL) << "Unspecified input path to Xlist!";
    PrintSyntax();
    exit(1);
  }

  if( parser.OptionExists('y') ){ 
    fPathYlist = parser.ArgAsString('y');
    LOG("gcalchedisdiffxsec", pINFO) << "PathYlist = " << fPathYlist;
  }          
  else {
    LOG("gcalchedisdiffxsec", pFATAL) << "Unspecified input path to Ylist!";
    PrintSyntax();
    exit(1);
  }

  if( parser.OptionExists('e') ){ 
    fPathElist = parser.ArgAsString('e');
    LOG("gcalchedisdiffxsec", pINFO) << "PathElist = " << fPathElist;
  }          
  else {
    LOG("gcalchedisdiffxsec", pFATAL) << "Unspecified input path to Elist!";
    PrintSyntax();
    exit(1);
  }

  if( parser.OptionExists('s') ) {
    fSaveAll = true;
    LOG("gcalchedisdiffxsec", pINFO) << "SaveAll = " << fSaveAll;
  }          

  if( parser.OptionExists('o') ){ 
    fOutFileName = parser.ArgAsString('o');
    LOG("gcalchedisdiffxsec", pINFO) << "OutFileName = " << fOutFileName;
  }          
  else {
    LOG("gcalchedisdiffxsec", pFATAL) << "Unspecified output file name!";
    PrintSyntax();
    exit(1);
  }

}

//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gcalchedisdiffxsec", pNOTICE)
      << "\n\n" << "Syntax:" << "\n"
      << "   gcalchedisdiffxsec -p nu -t tgt -o root_file\n"
      << "            -x pathxlist\n"
      << "            -y pathylist\n"
      << "            -e pathelist\n"
      << "            --tune genie_tune\n"
      << "            --event-generator-list list_name\n"
      << "            [--message-thresholds xml_file]\n";
}