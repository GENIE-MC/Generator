//____________________________________________________________________________
/*!

\program gvld_merenyi_test

\brief   Merenyi test (validating GENIE intranuclear rescattering models with
         neutrino bubble chamber data)

\syntax  ./gvld_merenyi_test

\help     Specify D2 or Ne analysis and the event file.
          The script will look for the event file in the same directory.
          The included examples are "d2_v5.dat" and "neon_v15.dat".
          The event files use the `ghad' format.

          The program analyzes the data and displays the fractional populations 
          for the various reactions on the screen. Neutrino energy, muon angle 
          cosine, and charged pion momentum plots are created and saved in the 
          same directory (pl1.gif, pl2.gif, pl3.gif).
           
          The program also writes the fractional populations and the histograms 
          drawn in the above plots to text files according to this table:
          (i)   Fractional populations: "populations.txt"    
                First line includes the number of simulated events and the 
                number of events kept.
          (ii)  Neutrino energies: "nu_energies.txt"    
                First line includes the number of simulated events.
          (iii) Muon angle cosines: "muon_angles.txt"    
                First line includes the number of simulated events.
          (iv)  Charged pion momenta: "pion_momenta.txt"   
                First line includes the number of charged pions.
          Fractional populations are listed as in Table 1 of the Merenyi 
          et al. paper.

\author  Pauli Kehayias
         Hugh Gallagher

\created January 13, 2009

\cpright Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>

#include <TParticle.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <Riostream.h>

#include "validation/Merenyi/MerenyiNuEvent.h"
#include "validation/Merenyi/MerenyiAnalysisD2.h"
#include "validation/Merenyi/MerenyiAnalysisNe.h"

using namespace std;
using namespace genie;
using namespace genie::vld_merenyi;

// D2 constants
//
const int numDatEvents1_d2 = 1826;
const int numDatEvents2_d2 = 785;     //estimated d2 pi+ and pi- counts

// Ne constants
//
const int numDatEvents1_ne = 490;
const int numDatEvents2_ne = 186;     //estimated ne pi+ and pi- counts

// globals
//
int numDatEvents1, numDatEvents2;
string nuE_file, nuEFrac_file, muAng_file, piMom_file;

// function prototypes
//
string getType   (void);
string getFile   (void);
void   drawPlot1 (MerenyiAnalysis* ptr, string type);	
void   drawPlot2 (MerenyiAnalysis* ptr, string type);
void   drawPlot3 (MerenyiAnalysis* ptr, string type); 

//____________________________________________________________________________
int main()
{
  string type, file;
  MerenyiAnalysis *anaPtr;
  MerenyiAnalysisD2 d2Ana; MerenyiAnalysisNe neAna;

  // pick target
  type = getType();		

  if (type == "d2")
  {
    anaPtr = &d2Ana;
    numDatEvents1 = numDatEvents1_d2;
    numDatEvents2 = numDatEvents2_d2;
    nuE_file = "original_energies_norm_d2.txt";
    nuEFrac_file = "merenyi_pl1_frac_d2.txt";
    muAng_file = "original_angles_d2.txt";
    piMom_file = "original_piMomenta_d2.txt";

  }

  else if (type == "ne")
  {
    anaPtr = &neAna;
    numDatEvents1 = numDatEvents1_ne;
    numDatEvents2 = numDatEvents2_ne;
    nuE_file = "original_energies_norm_ne.txt";
    nuEFrac_file = "merenyi_pl1_frac_ne.txt";
    muAng_file = "original_angles_ne.txt";
    piMom_file = "original_piMomenta_ne.txt";
  }

  else if (type == "q")
    exit(0);

  else {
    cout << "Error in target selection" << endl;
    exit(0);
  }

  file = getFile();		//set simulated data file
  cout << nuE_file << endl;

  //run the analysis
  anaPtr->setFile(file);
  anaPtr->readData();
  anaPtr->dispPops();

  //plot the results
  drawPlot1(anaPtr, type);
  drawPlot2(anaPtr, type);
  drawPlot3(anaPtr, type);

  //print out lists of collected fractional populations, neutrino
  //energies, muon angle cosines, and charged pion momenta
  anaPtr->printPops();
  anaPtr->printNuEHist();
  anaPtr->printMuAngHist();
  anaPtr->printPiMomHist();

  cout << "\nSimulated Merenyi analysis complete\n";
  return 0;
}

//____________________________________________________________________________
string getType(void)
{
  string anaType;

  do {
    cout << "Enter target (d2 = deuterium, ne = neon, q = quit): ";
    cin >> anaType;
  } while (anaType != "d2" && anaType != "ne" && anaType != "q");

  return anaType;
}
//____________________________________________________________________________
string getFile(void)
{
  string inputFile;

  cout << "Enter file name: ";
  cin >> inputFile;

  return inputFile;
}
//____________________________________________________________________________
//
// plots neutrino energies
//
void drawPlot1(MerenyiAnalysis* ptr, string type)	
{
  ifstream inDat, inCorr;

  double x[28], yMc[28], yDat[28], exMc[28], eyMc[28], exDat[28], eyDat[28];
  double numEvents = ptr->getNumCCEvents();
  double p, q;	//scratch variables

  TCanvas *canvas = new TCanvas("c1", "Neutrino Energies");
  canvas->SetFillColor(10);

  //read in the data
  inDat.open(nuE_file.c_str());		//original data (normalized to 1)

  inDat >> p;	//skip two lines
  inDat >> p;

  for (int i = 0; i < 28; i++)
  {
    p = ptr->getNuEi(i+2);
    x[i] = 0.2*i + 0.5;
    yMc[i] = p;
    exMc[i] = 0.1;
    eyMc[i] = pow(p*(1-p) / numEvents, 0.5);

    inDat >> q;
    yDat[i] = q;
    exDat[i] = 0;

    inCorr >> q;
    if (q == 0)
      eyDat[i] = 0;
    else
      eyDat[i] = yDat[i] * pow((1-q) / (q*numDatEvents1), 0.5);
  }

  inDat.close();  inCorr.close();

  //make the plot
  TGraphErrors *dataGr = new TGraphErrors(28, x, yDat, exDat, eyDat);
  dataGr->SetMarkerStyle(21);

  TGraphErrors *mcGr = new TGraphErrors(28, x, yMc, exMc, eyMc);
  mcGr->SetFillColor(2);

  if (type == "d2")
    mcGr->SetTitle("#nu_{#mu} Energies (D target)");
  else
    mcGr->SetTitle("#nu_{#mu} Energies (Ne target)");

  mcGr->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
  mcGr->GetXaxis()->CenterTitle();
  mcGr->GetXaxis()->SetLimits(0, 6);
  mcGr->GetYaxis()->SetTitle("Fraction of Events / (0.2 GeV)");
  mcGr->GetYaxis()->CenterTitle();
  mcGr->GetYaxis()->SetTitleOffset(1.2);

  mcGr->Draw("a2");
  dataGr->Draw("p");

  canvas->SaveAs("pl1.gif");
  cout << "plot 1 done" << endl;
}
//____________________________________________________________________________
//
// plots muon angle cosines
//
void drawPlot2(MerenyiAnalysis* ptr, string type)      
{
  ifstream inDat;

  double x[20], yMc[20], yDat[20], exMc[20], eyMc[20], exDat[20], eyDat[20];
  double numEvents = ptr->getNumEvents();
  double p;	//scratch variable

  TCanvas *canvas = new TCanvas("c2", "Muon Angles");
  canvas->SetFillColor(10);

  //read in the data
  inDat.open(muAng_file.c_str());         //original data

  for (int i = 0; i < 20; i++)
  {
    p = ptr->getMuAngi(i);
    x[i] = 0.1 * i - 0.95;
    yMc[i] = p;
    exMc[i] = 0.05;
    eyMc[i] = pow(p*(1-p) / numEvents, 0.5);

    inDat >> p;
    yDat[i] = p;
    exDat[i] = 0;
    eyDat[i] = pow(p*(1-p) / numDatEvents1, 0.5);
  }

  inDat.close();

  TGraphErrors *dataGr = new TGraphErrors(20, x, yDat, exDat, eyDat);
  dataGr->SetMarkerStyle(21);

  TGraphErrors *mcGr = new TGraphErrors(20, x, yMc, exMc, eyMc);
  mcGr->SetFillColor(2);

  if (type == "d2")
    mcGr->SetTitle("#mu^{-} Production Angle Cosines (D target)");
  else
    mcGr->SetTitle("#mu^{-} Production Angle Cosines (Ne target)");

  mcGr->GetXaxis()->SetTitle("Cosine(#nu_{#mu} - #mu^{-})");
  mcGr->GetXaxis()->CenterTitle();
  mcGr->GetXaxis()->SetLimits(-1, 1);
  mcGr->GetYaxis()->SetTitle("Fraction of Events / (0.1)");
  mcGr->GetYaxis()->CenterTitle();
  mcGr->GetYaxis()->SetTitleOffset(1.2);

  mcGr->Draw("a2");
  dataGr->Draw("p");

  canvas->SaveAs("pl2.gif");
  cout << "plot 2 done" << endl;
}
//____________________________________________________________________________
//
// plots charged pion momenta
//
void drawPlot3(MerenyiAnalysis* ptr, string type)      
{
  ifstream inDat;

  double x[20], yMc[20], yDat[20], exMc[20], eyMc[20], exDat[20], eyDat[20];
  double numEvents = ptr->getNumEvents();
  double p;     //scratch variable

  TCanvas *canvas = new TCanvas("c3", "Charged Pion Momenta");
  canvas->SetFillColor(10);

  //read in the data
  inDat.open(piMom_file.c_str());         //original data

  for (int i = 0; i < 20; i++)
  {
    p = ptr->getPiMomi(i);
    x[i] = 0.1 * i + 0.05;
    yMc[i] = p;
    exMc[i] = 0.05;
    eyMc[i] = pow(p*(1-p) / numEvents, 0.5);

    inDat >> p;
    yDat[i] = p;
    exDat[i] = 0;
    eyDat[i] =  pow(p*(1-p) / numDatEvents2_d2, 0.5);
  }

  inDat.close();

  TGraphErrors *dataGr = new TGraphErrors(20, x, yDat, exDat, eyDat);
  dataGr->SetMarkerStyle(21);

  TGraphErrors *mcGr = new TGraphErrors(20, x, yMc, exMc, eyMc);

  mcGr->GetYaxis()->SetTitleOffset(1.2);
  mcGr->SetFillColor(2);

  if (type == "d2")
    mcGr->SetTitle("Final State Pion Momenta (D target)");
  else
    mcGr->SetTitle("Final State Pion Momenta (Ne target)");

  mcGr->GetXaxis()->SetTitle("Pion Momentum (GeV/c)");
  mcGr->GetXaxis()->CenterTitle();
  mcGr->GetXaxis()->SetLimits(0, 2);
  mcGr->GetYaxis()->SetTitle("Fraction of #pi^{#pm} / (0.1 GeV)");
  mcGr->GetYaxis()->CenterTitle();
  mcGr->GetYaxis()->SetTitleOffset(1.2);

  mcGr->Draw("a2");
  dataGr->Draw("p");

  canvas->SaveAs("pl3.gif");
  cout << "plot 3 done" << endl;
}
//____________________________________________________________________________
