//____________________________________________________________________________
/*!

\program testGiBUUData

\brief   Program used for testing and debugging the input GiBUU data tables 
         and interpolation.

\syntax  gtestGiBUUData

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created May 31, 2009

\cpright Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <TFile.h>
#include <TTree.h>

#include "BaryonResonance/BaryonResUtils.h"
#include "GiBUU/GiBUUData.h"
#include "Messenger/Messenger.h"
#include "Numerical/Spline.h"
#include "PDG/PDGCodes.h"

using namespace genie;
using namespace genie::utils;

void SaveToRootFile(void);

int main(int /*argc*/, char ** /*argv*/)
{
  SaveToRootFile();
  return 0;
}

void SaveToRootFile(void)
{
  //
  // get GiBUU data
  //

  GiBUUData * gd = GiBUUData::Instance();

  const GiBUUData::FormFactors & ff = gd->FF();

  //
  // define inputs
  //

  const int kNRes    =  13;
  const int kNNuc    =   2;
  const int kNIntTyp =   3;
  const int kNQ2     = 200;

  Resonance_t resonance[kNRes] = { 
      kP33_1232, kS11_1535, kD13_1520, kS11_1650,
      kD15_1675, kS31_1620, kD33_1700, kP11_1440,
      kP13_1720, kF15_1680, kP31_1910, kF35_1905, kF37_1950  };

  int nucleon[kNNuc] = { 
      kPdgProton, kPdgNeutron };

  InteractionType_t interaction[kNIntTyp] = {
      kIntWeakCC, kIntWeakNC, kIntEM };

  double Q2min = 0;
  double Q2max = 4;
  double dQ2   = (Q2max-Q2min)/(kNQ2-1);

  //
  // define the output tree
  //

  int    bResId;        // resonance id
  int    bNucleon;      // nucleon pdg
  int    bIntType;      // interaction type (cc, nc, em)
  int    bIsoX2;        // resonance isospin x 2
  double bQ2;           // Q2 
  double bC3V;          // C3V f/f
  double bC4V;          // C4V f/f
  double bC5V;          // C5V f/f
  double bC6V;          // C6V f/f
  double bC3A;          // C3A f/f
  double bC4A;          // C4A f/f
  double bC5A;          // C5A f/f
  double bC6A;          // C6A f/f
  double bF1V;          // F1V f/f
  double bF2V;          // F2V f/f
  double bFA;           // FA  f/f
  double bFP;           // FP  f/f

  TTree * gibuu_ffres = new TTree("gibuu_ffres","GiBUU resonance form factors");

  gibuu_ffres -> Branch ("res",    &bResId,    "res/I"   );
  gibuu_ffres -> Branch ("nuc",    &bNucleon,  "nuc/I"   );
  gibuu_ffres -> Branch ("intyp",  &bIntType,  "intyp/I" );
  gibuu_ffres -> Branch ("isoX2",  &bIsoX2,    "isoX2/I" );
  gibuu_ffres -> Branch ("Q2",     &bQ2,       "Q2/D"    );
  gibuu_ffres -> Branch ("C3V",    &bC3V,      "C3V/D"   );
  gibuu_ffres -> Branch ("C4V",    &bC4V,      "C4V/D"   );
  gibuu_ffres -> Branch ("C5V",    &bC5V,      "C5V/D"   );
  gibuu_ffres -> Branch ("C6V",    &bC6V,      "C6V/D"   );
  gibuu_ffres -> Branch ("C3A",    &bC3A,      "C3A/D"   );
  gibuu_ffres -> Branch ("C4A",    &bC4A,      "C4A/D"   );
  gibuu_ffres -> Branch ("C5A",    &bC5A,      "C5A/D"   );
  gibuu_ffres -> Branch ("C6A",    &bC6A,      "C6A/D"   );
  gibuu_ffres -> Branch ("F1V",    &bF1V,      "F1V/D"   );
  gibuu_ffres -> Branch ("F2V",    &bF2V,      "F2V/D"   );
  gibuu_ffres -> Branch ("FA",     &bFA,       "FA/D"    );
  gibuu_ffres -> Branch ("FP",     &bFP,       "FP/D"    );

  //
  // calculate(interpolate) resonance f/f and fill-in the tree
  //

  for(int r=0; r<kNRes; r++) {
    for(int i=0; i<kNNuc; i++) {
      for(int j=0; j<kNIntTyp; j++) {

         bResId   = (int) resonance[r];
         bNucleon = nucleon[i];
         bIntType = interaction[j];
         bIsoX2   = res::IsDelta(resonance[r]) ? 3 : 1;

         for(int iq2=0; iq2<kNQ2; iq2++) {
             double Q2 = Q2min + iq2 * dQ2;

             bQ2  = Q2;
             bC3V = ff.C3V(Q2, resonance[r], nucleon[i], interaction[j]);
             bC4V = ff.C4V(Q2, resonance[r], nucleon[i], interaction[j]);
             bC5V = ff.C5V(Q2, resonance[r], nucleon[i], interaction[j]);
             bC6V = ff.C6V(Q2, resonance[r], nucleon[i], interaction[j]);
             bC3A = ff.C3A(Q2, resonance[r], nucleon[i], interaction[j]);
             bC4A = ff.C4A(Q2, resonance[r], nucleon[i], interaction[j]);
             bC5A = ff.C5A(Q2, resonance[r], nucleon[i], interaction[j]);
             bC6A = ff.C6A(Q2, resonance[r], nucleon[i], interaction[j]);
             bF1V = ff.F1V(Q2, resonance[r], nucleon[i], interaction[j]);
             bF2V = ff.F2V(Q2, resonance[r], nucleon[i], interaction[j]);
             bFA  = ff.FA (Q2, resonance[r], nucleon[i], interaction[j]);
             bFP  = ff.FP (Q2, resonance[r], nucleon[i], interaction[j]);

             gibuu_ffres->Fill();
         }//Q2
      }//interaction
    }//nucleon
  }//res

  //
  // save
  //
  TFile f("./gibuu_data.root", "recreate");
  gibuu_ffres->Write();
  f.Close();
}
