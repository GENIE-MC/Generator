//____________________________________________________________________________
/*!

\program gtestCOHNCGammaXSec.C

\brief   Macro to use for testing/debugging the COH NCGamma model and compare to external
         cross section results.

         Run with: (make sure to enclose the macro in "" )
         `genie -l "COHNCGammaXSec.C(target, probe_pdg, max_nu_energy)"`

          optional additional arguements are:
         `genie -l "COHNCGammaXSec.C(target, probe_pdg, max_nu_energy, theta_incoming_nu, theta_gamma, phi_gamma)"`

          For example: an incoming neutrino (14) with energy of 1GeV interacting with a 40Ar (1000180400) nucleus
          `genie -l "COHNCGammaXSec.C(14,1000180400,1)"`
      
          See below in "Main parameters for COH NC Single Gamma Cross Section for defintion of the parameters.

\author  Jon Sensenig inspired by gtestDISSF test file

\created July 23, 2020

*/
//____________________________________________________________________________

#include <string>
#include <vector>
#include <array>
#include <cstdio>

#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>


#include "Framework/Algorithm/Algorithm.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Conventions/Units.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Utils/CmdLnArgParser.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/Coherent/XSection/COHFormFactorI.h"
#include "Physics/Coherent/XSection/COHFormFactorMap.h"
#include "Physics/Coherent/XSection/DeVriesFormFactor.h"
#include "Physics/Coherent/XSection/FourierBesselFFCalculator.h"
#include "Physics/Coherent/XSection/AlvarezRusoSalaCOHGammaPXSec.h"

#include "COHNCGamma_ext/NCgamma_Diff_Cross_Section.h"

// Load the shared library for the external NC COH Gamma xsection
R__LOAD_LIBRARY(COHNCGamma_ext/libGenie_XSection.so)

using namespace genie;
using namespace std;

// Constructor for external X-Sec code
using namespace NC_gamma;
//Diff_Cross_Section CS("nu", "40Ar");

// Define so cleaner code later
typedef AlvarezRusoSalaCOHGammaPXSec ARSXSec;

const unsigned int steps  = 100; // Number of plotted xsec points
const unsigned int angles = 7;   // Number of plotted gamma theta angles

void gtestCOHNCGammaXSec(int tgt, int prb, double E, double theta_nu=1., double phi=0.);
void XSectionTest(const TFile & f);
void SigmaEg(const Interaction *i, const ARSXSec * xsec_alg, Kinematics * kine); 
std::pair<std::string,std::string> ExtXSec(int nu_probe, int tgt);
void PlotXSecs(double Eg_arr[steps], double gxsec[][steps], double exsec[][steps]);

////////////////////////////
//
// Main parameters for COH NC Single Gamma Cross Section 
//
////////////////////////////

// Require set by main function
int target       = -1000;// Target nucleus [pdg]
int probe        = -1000;// Incoming particle pdg [pdg]
double probe_E   = -1.;  // Incoming nu energy [GeV]
// Defaults
int prod         = 22; // Produced FS particle = gamma [pdg]
double theta_lep = 1.; // Outgoing nu wrt incoming nu [deg]
double phi_lep   = 0.; // Outgoing nu in plane transverse to incoming nu [deg]
double theta_g   = 1.; // Gamma angle wrt to incoming nu [deg]
double phi_g     = 0.; // Gamma angle in plane transverse to incoming nu [deg]

array<double, angles> theta_gamma = {3,5,10,30,-3,-5,-10}; // FS gamma wrt z-axis in deg
array<short, angles>  root_color  = {600,1,632,416,616,800,432};

/////////////////////////////

double Eg_arr [steps];
double xsec_arr[angles][steps];
double ext_xsec_arr[angles][steps];

//__________________________________________________________________________
void gtestCOHNCGammaXSec (int tgt, int prb, double E, double theta_nu, double phi)
{

  target = tgt;
  probe = prb;
  probe_E = E;
  theta_lep = theta_nu;
  phi_g = phi;

  TFile file("./coh_ncgamma.root","recreate");
  
  XSectionTest(file);  // Top level cross section

  file.Close();
}
//__________________________________________________________________________
void XSectionTest(const TFile & file)
{

  cout << "\n ============= Beginning COH NC 1gamma Test ============== \n" << endl;

  // Instantiate external cross section
  std::pair<std::string,std::string> params = ExtXSec(probe, target);
  Diff_Cross_Section CS(params.first, params.second); // 1st = probe 2nd = target

  double range   = probe_E;
  double E_lep   = probe_E;

  // -- request the COH NC Gamma model
  AlgFactory * algf = AlgFactory::Instance();
  AlgId id("AlvarezRusoSalaCOHGammaPXSec","Default");
  const Algorithm * algXsec = algf->GetAlgorithm(id);
  const XSecAlgorithmI * xsec_alg = dynamic_cast<const XSecAlgorithmI *>(algXsec);

  // Set up the interaction and kinematics
  Interaction * i = Interaction::COHNC(target, probe, prod, probe_E);
  Kinematics * kine = i->KinePtr();

  // Set up the interaction kinematics
  TVector3       p3_g;
  TVector3       p3_lep;
  TLorentzVector p4_g;
  TLorentzVector p4_lep;
  TLorentzVector p4_nu(0., 0., probe_E, probe_E); // probe along z

  KinePhaseSpace_t t = kPSEgOlOgfE;

  for( unsigned int ang = 0; ang < angles; ang++ ) {
    theta_g = theta_gamma[ang]*TMath::DegToRad();
    E_lep = probe_E;
    for( unsigned int eg = 0; eg < steps; eg++ ) {
                                                                                                         
      double E_g = probe_E - E_lep;
      // Set FS gamma 4-mom.
      p3_g.SetMagThetaPhi( E_g, theta_g*TMath::DegToRad(), 
			   phi_g*TMath::DegToRad() );
      p4_g.SetVect( p3_g );
      p4_g.SetE( E_g );
      
      // Set FS lepton (nu) 4-mom.
      p3_lep.SetMagThetaPhi(  E_lep, theta_lep*TMath::DegToRad(), 0. );
      p4_lep.SetVect( p3_lep );
      p4_lep.SetE( E_lep );
                                                                                                        
      kine->SetHadSystP4( p4_g );
      kine->SetFSLeptonP4( p4_lep );
                 
      TLorentzVector q = p4_nu-p4_lep ;
      double Q2 = -q.Mag2();
      double x = Q2/(2 * E_g * constants::kNucleonMass );
      double y = E_g/probe_E;
      kine->Setx( x );
      kine->Sety( y );
      utils::kinematics::UpdateWQ2FromXY(i);
      double coh_t = (q - p4_g).Mag() ;
      kine->Sett( coh_t ) ;

      // Check the validity
      cout << " Gamma theta angle " << theta_gamma[ang]
           << " Gamma Energy " << E_g << " Q2 " << Q2 << " x " << x << " y " << y
           << " Valid process: " << xsec_alg->ValidProcess(i) 
           << " Valid Kinematics: " <<  xsec_alg->ValidKinematics(i) << endl;
      
      // Calculate Genie Xsection
      xsec_arr[ang][eg] = (xsec_alg->XSec(i, t) * 1.e41)/units::cm2;
      
      // Calculate external Xsection
      // If invalid process or kinematics set cross section to 0
      if ( !xsec_alg->ValidProcess(i) || !xsec_alg->ValidKinematics(i) ) {
         ext_xsec_arr[ang][eg] = 0.;
      } else {
        ext_xsec_arr[ang][eg] = ( CS.getDiffCrossSection( probe_E, E_lep, theta_lep*TMath::DegToRad(),
                                                          theta_g*TMath::DegToRad(), phi_g*TMath::DegToRad()
                                                          )*1e41 ) / units::cm2;
      } 

      cout << "---> GENIE Cross Section: " << xsec_arr[ang][eg] 
           << "  External Cross Section: " << ext_xsec_arr[ang][eg] << endl;

      Eg_arr[eg] = E_g;
      E_lep -= range / double( steps-1 );

    } // gamma E
  } // gamma angle

  PlotXSecs(Eg_arr, xsec_arr, ext_xsec_arr);
}
//__________________________________________________________________________
std::pair<std::string,std::string> ExtXSec(int nu_probe, int tgt)
{
  std::pair<std::string,std::string> ext_pdgs("","");

  if (nu_probe == kPdgNuMu)          ext_pdgs.first = "nu";
  else if (nu_probe == kPdgAntiNuMu) ext_pdgs.first = "anti-nu";
  else std::cout << "Unknown probe pdg!" << std::endl;

  if (tgt == 1000060120) ext_pdgs.second = "12C";
  else if (tgt == 1000180400) ext_pdgs.second = "40Ar";
  else std::cout << "Unknown nucleus pdg!" << std::endl;

  return ext_pdgs; 
}
//__________________________________________________________________________
void PlotXSecs(double Eg_arr[steps], double gxsec[][steps], double exsec[][steps]) {


  TLegend *legend = new TLegend(0.5,0.7,0.9,0.9);
  TCanvas *c1 = new TCanvas("c1","COH NC Gamma",10,10,1000,450);
  c1->Divide(2,1);

  auto mg = new TMultiGraph();
  auto ext_mg = new TMultiGraph();

  char text[50];
  for( unsigned int k = 0; k < angles; k++ ) {
    auto gr = new TGraph( steps, Eg_arr, gxsec[k] ); gr->SetLineColor( root_color[k] ); gr->SetLineWidth(2);
    auto egr = new TGraph( steps, Eg_arr, exsec[k] ); egr->SetLineColor( root_color[k] ); egr->SetLineWidth(2);

    snprintf( text, 50, "#theta_{#gamma} = %d",(int)theta_gamma[k] );
    legend->AddEntry( gr, (char*)text, "l" );

    mg->Add(gr);
    ext_mg->Add(egr);
  }

  stringstream description ; 

  TParticlePDG * tprobe = PDGLibrary::Instance() -> Find( probe ) ;
  TParticlePDG * ttgt = PDGLibrary::Instance() -> Find( target ) ;

  description << tprobe -> Title() << " on " << ttgt -> Title() ;
  description << " #theta_{l}=" << theta_l << "^{#circ}" ;
  description << " #phi_{#gamma}=" << phi_g << "^{#circ}" ;
  description << ";E_{#gamma} [GeV];#frac{d^{5}#sigma}{dE_{#gamma}d#Omega_{l}d#Omega_{#gamma}} [10^-41 #frac{cm^2}{GeV}]" ;

  std::string genie_title = "GENIE " + description.str() ;
  std::string ext_title = "Eduardo " + description.str() ;

  mg->SetTitle( genie_title.c_str() );
  ext_mg->SetTitle( ext_title.c_str() );

  c1->cd(1);
  mg->Draw("ACP");
  legend->Draw("same");
  c1 -> Update() ;
  c1->cd(2);
  ext_mg->Draw("ACP");

  mg->Write();
  ext_mg->Write();
                                                                                                                                 
}
