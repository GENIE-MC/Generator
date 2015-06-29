//____________________________________________________________________________
/*!

\program gtestDISSF

\brief   Program used for testing / debugging the GENIE DIS SF models

         Syntax :
          gtestDISSF -a model -c config [-m mode] [-x x] [-q Q2]

         Options :
           -a  DIS SF model (algorithm name, eg genie::BYStructureFuncModel)
           -c  DIS SF model configuration
           -m  mode (1: make std SF ntuple, 2: vertical slice)
               [default:1]
           -x  Specify Bjorken x to be used at the vertical slice
           -q  Specify mom. transfer Q2(>0) to be used at the vertical slice
         

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created June 20, 2004

\cpright Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <string>

#include <TFile.h>
#include <TTree.h>

#include "Algorithm/Algorithm.h"
#include "Algorithm/AlgFactory.h"
#include "Base/DISStructureFunc.h"
#include "Base/DISStructureFuncModelI.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "Utils/CmdLnArgParser.h"

using namespace genie;
using std::string;

void GetCommandLineArgs(int argc, char ** argv);
void PrintSyntax(void);

void BuildStdNtuple (void);
void VerticalSlice  (void);

int    gMode        = 1;
double gX           = 0;
double gQ2          = 0;
string gDISSFAlg    = "";
string gDISSFConfig = "";

//__________________________________________________________________________
int main(int argc, char ** argv)
{
  // -- parse the command line arguments (model choice,...)
  GetCommandLineArgs(argc,argv);

  if(gMode==1) BuildStdNtuple();
  if(gMode==2) VerticalSlice ();

  return 0;
}
//__________________________________________________________________________
void BuildStdNtuple(void)
{
  // -- define initial states & x,Q2 values to compute the DIS SFs 
  const int kNNu   = 6;
  const int kNNuc  = 2;
  const int kNTgt  = 3;
  const int kNX    = 12;
  const int kNQ2   = 60;

  int neutrino    [kNNu]   = { kPdgNuE, kPdgAntiNuE, kPdgNuMu, kPdgAntiNuMu, kPdgNuTau, kPdgAntiNuTau };
  int hit_nucleon [kNNuc]  = { kPdgProton, kPdgNeutron };
  int target      [kNTgt]  = { kPdgTgtFreeP, kPdgTgtFreeN, kPdgTgtFe56 };

  double xbj[kNX] = {
        0.0150, 0.0450, 0.0800, 0.1250, 0.1750, 0.2250,
        0.2750, 0.3500, 0.4500, 0.5500, 0.6500, 0.7500
  };
  double Q2[kNQ2] = {
        0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1.0,  2.0, 
        3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0, 10.0, 11.0, 12.0, 
       13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0,
       23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 42.0,
       33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 52.0,
       43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 51.0, 62.0
  };

  // -- request the specified DIS SF model
  AlgFactory * algf = AlgFactory::Instance();
  const DISStructureFuncModelI * dissf_model =
      dynamic_cast<const DISStructureFuncModelI *> (
                              algf->GetAlgorithm(gDISSFAlg, gDISSFConfig));
  assert(dissf_model);

  // -- create a structure functions objects and attach the specified model
  DISStructureFunc dissf;
  dissf.SetModel(dissf_model);

  // -- define the outout TTree
  TTree * dissf_nt = new TTree("dissf_nt","DIS SF");
  int    brNu;    // event number
  int    brNuc;   // hit nucleon PDG code
  int    brTgt;   // neutrino PDG code
  double brXbj;   //
  double brQ2;    //
  double brF1;    //
  double brF2;    //
  double brF3;    //
  double brF4;    //
  double brF5;    //
  dissf_nt->Branch("nu",  &brNu,  "nu/I");
  dissf_nt->Branch("nuc", &brNuc, "nuc/I");
  dissf_nt->Branch("tgt", &brTgt, "tgt/I");
  dissf_nt->Branch("x",   &brXbj, "x/D");
  dissf_nt->Branch("Q2",  &brQ2,  "Q2/D");
  dissf_nt->Branch("F1",  &brF1,  "F1/D");
  dissf_nt->Branch("F2",  &brF2,  "F2/D");
  dissf_nt->Branch("F3",  &brF3,  "F3/D");
  dissf_nt->Branch("F4",  &brF4,  "F4/D");
  dissf_nt->Branch("F5",  &brF5,  "F5/D");

  // -- loop over the initial states
  for(int inu=0; inu<kNNu; inu++) {
     for(int inuc=0; inuc<kNNuc; inuc++) {
        for(int itgt=0; itgt<kNTgt; itgt++) {

           // -- get a DIS interaction object & access its kinematics
           Interaction * interaction = Interaction::DISCC(
                                 target[itgt],hit_nucleon[inuc],neutrino[inu]);
           Kinematics * kine = interaction->KinePtr();

           // -- loop over x,Q2
           for(int ix=0; ix<kNX; ix++) {
             for(int iq=0; iq<kNQ2; iq++) {

                // update the interaction kinematics
                kine->Setx(xbj[ix]);
                kine->SetQ2(Q2[iq]);

                // calculate the structure functions for the input interaction
                dissf.Calculate(interaction);

                // update all tree branches & save
                brNu  = neutrino[inu];
                brNuc = hit_nucleon[inuc];
                brTgt = target[itgt];
                brXbj = xbj[ix];
                brQ2  = Q2[iq];
                brF1  = dissf.F1();
                brF2  = dissf.F2();
                brF3  = dissf.F3();
                brF4  = dissf.F4();
                brF5  = dissf.F5();
                dissf_nt->Fill();
             }//Q2
           }//x

           delete interaction;

        }//tgt
     }//nuc
  }//nu

  TFile f("./dissf.root","recreate");
  dissf_nt->Write();
  f.Close();
}
//__________________________________________________________________________
void VerticalSlice(void)
{
  // manual override of Registry mesg level for more verbose output
  Messenger * msg = Messenger::Instance();
  msg->SetPriorityLevel("DISSF",     pDEBUG);
  msg->SetPriorityLevel("BodekYang", pDEBUG);
  msg->SetPriorityLevel("PDFLIB",    pDEBUG);

  // -- request the specified DIS SF model
  AlgFactory * algf = AlgFactory::Instance();
  const DISStructureFuncModelI * dissf_model =
      dynamic_cast<const DISStructureFuncModelI *> (
                              algf->GetAlgorithm(gDISSFAlg, gDISSFConfig));
  assert(dissf_model);

  // -- create a structure functions objects and attach the specified model
  DISStructureFunc dissf;
  dissf.SetModel(dissf_model);

  const int kNNu   = 2;
  const int kNNuc  = 2;
  const int kNTgt  = 1;

  int neutrino    [kNNu]   = { kPdgNuMu,   kPdgAntiNuMu };
  int hit_nucleon [kNNuc]  = { kPdgProton, kPdgNeutron  };
  int target      [kNTgt]  = { kPdgTgtFe56 };

  for(int inu=0; inu<kNNu; inu++) {
     for(int inuc=0; inuc<kNNuc; inuc++) {
        for(int itgt=0; itgt<kNTgt; itgt++) {

           // -- get a DIS interaction object & access its kinematics
           Interaction * interaction = Interaction::DISCC(
                              target[itgt],hit_nucleon[inuc],neutrino[inu]);
           Kinematics * kine = interaction->KinePtr();
           kine->Setx(gX);
           kine->SetQ2(gQ2);

           LOG("test", pNOTICE) << "*** Vertical SF slice for: ";
           LOG("test", pNOTICE) << *interaction;

           // calculate the structure functions for the input interaction
           dissf.Calculate(interaction);

           LOG("test", pNOTICE) << dissf;

           delete interaction;
        }
     }
  }

  // print algorithms & heir configurations
  LOG("test", pNOTICE) << *algf;
}
//__________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
// Parse the command line arguments

  CmdLnArgParser parser(argc,argv);

  // DIS SF alg:
  if( parser.OptionExists('a') ) {
    LOG("test", pINFO) << "Reading DIS SF algorithm name";
    gDISSFAlg = parser.ArgAsString('a');
  } else {
    LOG("test", pINFO) << "No DIS SF algorithm was specified";
    PrintSyntax();
    exit(1);
  }

  // DIS SF config:
  if( parser.OptionExists('c') ) {
    LOG("test", pINFO) << "Reading DIS SF algorithm config name";
    gDISSFConfig = parser.ArgAsString('c');
  } else {
    LOG("test", pINFO) << "No DIS SF algorithm config was specified";
    PrintSyntax();
    exit(1);
  }

  // testDISSF mode:
  if( parser.OptionExists('m') ) {
    LOG("test", pINFO) << "Reading testDISSF mode";
    gMode = parser.ArgAsInt('m');
  } else {
    LOG("test", pINFO) << "No testDISSF was specified. Using default";
  }

  // x,Q2 for vertical slice mode
  if(gMode==2) {
    // x:
    if( parser.OptionExists('x') ) {
      LOG("test", pINFO) << "Reading x";
      gX = parser.ArgAsDouble('x');
    } else {
      LOG("test", pINFO) 
        << "No Bjorken x was specified for vertical slice";
      PrintSyntax();
      exit(1);
    }
    if( parser.OptionExists('q') ) {
      LOG("test", pINFO) << "Reading Q2";
      gQ2 = parser.ArgAsDouble('q');
    } else {
      LOG("test", pINFO) 
        << "No momentum transfer Q2 was specified for vertical slice";
      PrintSyntax();
      exit(1);
    }
  }//mode=2
}
//__________________________________________________________________________
void PrintSyntax(void)
{
  LOG("test", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
          << "  testDISSF -a model -c config [-m mode] [-x x] [-q Q2]\n";
}
//____________________________________________________________________________


