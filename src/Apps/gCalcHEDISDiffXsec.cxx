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
              differential cross section as function of:
              1 = Energy of outgoing lepton (Eo)
              2 = Inelasticity (y = 1 - Eo/Ei)
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
\cpright Copyright (c) 2003-2021, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include "TMath.h"
#include "TSystem.h"
#include "TH3D.h"
#include "TH2D.h"
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

#include <iostream>

using namespace genie;

void   PrintSyntax          (void);
void   DecodeCommandLine    (int argc, char * argv[]);
void   WriteDiffXSecDEDX    (GEVGDriver evg_driver, TH3D * hist);
void   WriteDiffXSecDYDX    (GEVGDriver evg_driver, TH3D * hist);

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//GLOBAL VARIABLES
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

const int    nlog10e   = 500;
const double log10emin = 2.;
const double log10emax = 10.;

const int    nlog10y   = 500;
const double log10ymin = -10.;
const double log10ymax = 0.;

const int    nlog10x   = 500;
const double log10xmin = -10.;
const double log10xmax = 0.;

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//INPUT ARGUMENTS
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

int fNu             = -1;
int fTgt            = -1;
int fTableType      = -1;
string fOutFileName = "";

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

  TString sufix = Form("%d_%d_%s",fNu,fTgt,fChannel.c_str());

  LOG("gcalchedisdiffxsec", pDEBUG) << sufix;

  TH3D * h3d;
  TH2D * h2d;
  if (fTableType==1) {
    
    h3d = new TH3D("dxsec_deodx_"+sufix,";log_{10}(E_{#nu}[GeV]);log_{10}(E_{out}[GeV]);log_{10}(x)",nlog10e,log10emin,log10emax,nlog10e,log10emin,log10emax,nlog10x,log10xmin,log10xmax);
    WriteDiffXSecDEDX(evg_driver,h3d);
    
    h2d = new TH2D("dxsec_deo_"+sufix,";log_{10}(E_{#nu}[GeV]);log_{10}(E_{out}[GeV]);d#sigma/dE_{out}",nlog10e,log10emin,log10emax,nlog10e,log10emin,log10emax);
    for ( int ix=1; ix<=h3d->GetNbinsX(); ix++ ) {
      double ei = TMath::Power(10.,h3d->GetXaxis()->GetBinCenter(ix));
      for ( int iy=1; iy<=h3d->GetNbinsY(); iy++ ) {
        double dxsec = 0.;
        for ( int iz=1; iz<=h3d->GetNbinsZ(); iz++ ) {
          double x = TMath::Power(10.,h3d->GetZaxis()->GetBinCenter(iz));
          dxsec += h3d->GetBinContent(ix,iy,iz) * x;
        }          
        dxsec *= h3d->GetZaxis()->GetBinWidth(1) * TMath::Log(10.) / ei;
        h2d->SetBinContent(ix,iy,dxsec);
      }
    }
    
  
  }
  else if (fTableType==2) {
    
    h3d = new TH3D("dxsec_dydx_"+sufix,";log_{10}(E_{#nu}[GeV]);log_{10}(y);log_{10}(x)",nlog10e,log10emin,log10emax,nlog10y,log10ymin,log10ymax,nlog10x,log10xmin,log10xmax);
    WriteDiffXSecDYDX(evg_driver,h3d);
    
    h2d = new TH2D("dxsec_dy_"+sufix,";log_{10}(E_{#nu}[GeV]);log_{10}(y);d#sigma/dy",nlog10e,log10emin,log10emax,nlog10y,log10ymin,log10ymax);
    for ( int ix=1; ix<=h3d->GetNbinsX(); ix++ ) {
      for ( int iy=1; iy<=h3d->GetNbinsY(); iy++ ) {
        double dxsec = 0.;
        for ( int iz=1; iz<=h3d->GetNbinsZ(); iz++ ) {
          double x = TMath::Power(10.,h3d->GetZaxis()->GetBinCenter(iz));
          dxsec += h3d->GetBinContent(ix,iy,iz) * x;
        }          
        dxsec *= h3d->GetZaxis()->GetBinWidth(1) * TMath::Log(10.);
        h2d->SetBinContent(ix,iy,dxsec);
      }
    }
  }

  TFile * outfile = new TFile(fOutFileName.c_str(),"RECREATE");
  h3d->Write(h3d->GetName());
  h2d->Write(h2d->GetName());
  outfile->Close();

}


//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//DIFFERENTIAL CROSS SECTION AS FUNCTION OF THE ENERGY OF OUTGOING LEPTON
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

void WriteDiffXSecDEDX(GEVGDriver evg_driver, TH3D * hist) {

  const InteractionList * intlst   = evg_driver.Interactions();

  LOG("gcalchedisdiffxsec", pDEBUG) << *intlst;

  InteractionList::const_iterator intliter;
  for(intliter = intlst->begin(); intliter != intlst->end(); ++intliter) {

    const Interaction * interaction = *intliter;
    const XSecAlgorithmI * xsec_alg = evg_driver.FindGenerator(interaction)->CrossSectionAlg();

    LOG("gcalchedisdiffxsec", pINFO) << "Current interaction: " << interaction->AsString();

    for ( int ix=1; ix<=hist->GetNbinsX(); ix++ ) {

      double ei = TMath::Power(10.,hist->GetXaxis()->GetBinCenter(ix));

      LOG("gcalchedisdiffxsec", pINFO) << "Energy: " << ei << " [GeV]";

      TLorentzVector p4(0,0,ei,ei);
      interaction->InitStatePtr()->SetProbeP4(p4);

      for ( int iy=1; iy<=hist->GetNbinsY(); iy++ ) {

        double eo = TMath::Power(10.,hist->GetYaxis()->GetBinCenter(iy));
     
        double y = 1. - eo / ei;
        if (y<0 || y>1) continue;
        interaction->KinePtr()->Sety(y);
        
        for ( int iz=1; iz<=hist->GetNbinsZ(); iz++ ) {

          double x = TMath::Power(10.,hist->GetZaxis()->GetBinCenter(iz));
          interaction->KinePtr()->Setx(x);
          utils::kinematics::UpdateWQ2FromXY(interaction);
          double dxsec = xsec_alg->XSec(interaction, kPSxyfE);

          LOG("gcalchedisdiffxsec", pDEBUG) << "x: " << x << "  eo: " << eo << " -> dsdxdeo[E=" << ei << "GeV] = " << dxsec;

          hist->SetBinContent(ix,iy,iz,hist->GetBinContent(ix,iy,iz)+dxsec);

        }

      }
    
    }

  }

}


//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//DIFFERENTIAL CROSS SECTION AS FUNCTION OF THE INELASTICITY
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

void WriteDiffXSecDYDX(GEVGDriver evg_driver, TH3D * hist) {

  const InteractionList * intlst   = evg_driver.Interactions();

  InteractionList::const_iterator intliter;
  for(intliter = intlst->begin(); intliter != intlst->end(); ++intliter) {

    const Interaction * interaction = *intliter;
    const XSecAlgorithmI * xsec_alg = evg_driver.FindGenerator(interaction)->CrossSectionAlg();

    LOG("gcalchedisdiffxsec", pINFO) << "Current interaction: " << interaction->AsString();

    for ( int ix=1; ix<=hist->GetNbinsX(); ix++ ) {

      double ei = TMath::Power(10.,hist->GetXaxis()->GetBinCenter(ix));

      LOG("gcalchedisdiffxsec", pINFO) << "Energy: " << ei << " [GeV]";

      TLorentzVector p4(0,0,ei,ei);
      interaction->InitStatePtr()->SetProbeP4(p4);

      for ( int iy=1; iy<=hist->GetNbinsY(); iy++ ) {

        double y = TMath::Power(10.,hist->GetYaxis()->GetBinCenter(iy));
     
        interaction->KinePtr()->Sety(y);

        for ( int iz=1; iz<=hist->GetNbinsZ(); iz++ ) {

          double x = TMath::Power(10.,hist->GetZaxis()->GetBinCenter(iz));
          interaction->KinePtr()->Setx(x);
          utils::kinematics::UpdateWQ2FromXY(interaction);
          double dxsec = xsec_alg->XSec(interaction, kPSxyfE);

          LOG("gcalchedisdiffxsec", pDEBUG) << "x: " << x << "  y: " << y << " -> dsdxdy[E=" << ei << "GeV] = " << dxsec;

          hist->SetBinContent(ix,iy,iz,hist->GetBinContent(ix,iy,iz)+dxsec);

        }
        
      }
    
    }

  }

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
    fTableType = parser.ArgAsInt('x');
    LOG("gcalchedisdiffxsec", pINFO) << "TableType = " << fTableType;
  }          
  else {
    LOG("gcalchedisdiffxsec", pFATAL) << "Unspecified input table type!";
    PrintSyntax();
    exit(1);
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
      << "   gcalchedisdiffxsec -p nu -t tgt -o root_file -x table_type\n"
      << "            --tune genie_tune\n"
      << "            --event-generator-list list_name\n"
      << "            [--message-thresholds xml_file]\n";
}