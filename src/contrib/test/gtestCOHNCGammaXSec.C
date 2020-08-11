//____________________________________________________________________________
/*!

\program gtestCOHNCGammaXSec

\brief   Crude v0 program used for testing/debugging the COH NCGamma model

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

#include "Framework/Algorithm/Algorithm.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Utils/CmdLnArgParser.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/Coherent/XSection/COHFormFactorI.h"
#include "Physics/Coherent/XSection/COHFormFactorMap.h"
#include "Physics/Coherent/XSection/DeVriesFormFactor.h"
#include "Physics/Coherent/XSection/FourierBesselFFCalculator.h"

#include "COHNCGamma_ext/NCgamma_Diff_Cross_Section.h"

// Load the shared library for the external NC COH Gamma xsection
R__LOAD_LIBRARY(COHNCGamma_ext/libGenie_XSection.so)

using namespace genie;
using namespace std;

using namespace NC_gamma;
Diff_Cross_Section CS("nu", "12C");

struct SigVars {
  double E_nu_probe;
  double E_lep;
  double theta_lep;
  double phi_lep;
  double E_g;
  double theta_g;
  double phi_g;
  double range;
} sigma_vars;
  
typedef AlvarezRusoSalaCOHGammaPXSec ARSXSec;
const unsigned int steps = 50;

array<double,7> theta_gamma {3,5,10,30,-3,-5,-10}; // FS gamma wrt z-axis in deg
array<short,7> root_color = {600,1,632,416,616,800,432};

void gtestCOHNCGammaXSec();
void XSectionTest(const TFile & f);
void GetCLIArgs();
void SigmaEg(const Interaction *i, const ARSXSec * xsec_alg, Kinematics * kine, SigVars sig_param); 
//void SigmaEnu(const Interaction *i, const ARSXSec * xsec_alg, Kinematics * kine, SigVars sig_param); 
double ExtXSec(int nu_probe, int tgt, double Enu, double Enu_final, double theta_l, double theta_g, double phi_g);
void SigmaCosTh(const Interaction *i, const ARSXSec * xsec_alg, Kinematics * kine, SigVars sig_param); 


//__________________________________________________________________________
void gtestCOHNCGammaXSec ()
{
  TFile file("./coh_ncgamma.root","recreate");

  XSectionTest(file);  // Top level cross section

  file.Close();
}

//__________________________________________________________________________
void GetCLIArgs()
{
// TODO implement
}
//__________________________________________________________________________
void XSectionTest(const TFile & file)
{

  cout << "\n ============= Beginning Delta Current Trace Test ============== \n" << endl;

  SigVars sigma_vars;
  int target = 1000180400;   // 40Ar
  //int target = 1000060120; // 12C
  int probe  = kPdgNuMu; // nu mu
  int prod   = 22;           // gamma

  // Fill the struct
  sigma_vars.E_nu_probe = 1.;
  sigma_vars.theta_lep  = 1.;
  sigma_vars.phi_lep    = 0.;
  sigma_vars.theta_g    = 10.; 
  sigma_vars.phi_g      = 10.;
  sigma_vars.range      = sigma_vars.E_nu_probe;
  sigma_vars.E_lep      = sigma_vars.E_nu_probe;

  // -- request the COH NC Gamma model
  AlgFactory * algf = AlgFactory::Instance();
  AlgId id("AlvarezRusoSalaCOHGammaPXSec","Default");
  const Algorithm * algXsec = algf->GetAlgorithm(id);
  const ARSXSec * xsec_alg = dynamic_cast<const ARSXSec *>(algXsec);

  // Set up the interaction and kinematics
  Interaction * i = Interaction::COHNC(target, probe, prod, sigma_vars.E_nu_probe);
  Kinematics * kine = i->KinePtr();
  sigma_vars.E_nu_probe = i->InitState().ProbeE(kRfLab);

  // Variable xsection is calculated wrt
  string plot = "E_gamma";
  if ( plot == "E_gamma" ) {
    SigmaEg(i, xsec_alg, kine, sigma_vars);
  } else if ( plot == "E_nu" ) {
    //SigmaEnu(i, xsec_alg, kine, sigma_vars);
  } else if ( plot == "Cos_th" ) {
    SigmaCosTh(i, xsec_alg, kine, sigma_vars);
  } else {
    cout << "Unknown cross section variable! Choose E_gamma, E_nu, or Cos_th" << endl;
  }

  cout << "End of Delta Current Trace test!" << endl;
}
//__________________________________________________________________________
void SigmaEg( const Interaction *i, const ARSXSec * xsec_alg, Kinematics * kine, SigVars sig_param)
{

  TLegend *legend = new TLegend(0.5,0.7,0.9,0.9);
  auto mg = new TMultiGraph();
  auto ext_mg = new TMultiGraph();
  TCanvas *c1 = new TCanvas("c1","COH NC Gamma",10,10,1000,450);
  c1->Divide(2,1);

  double Eg_arr [steps];
  double xsec_arr[theta_gamma.size()][steps];
  double ext_xsec_arr[theta_gamma.size()][steps];

  // Set up the interaction kinematics
  TVector3       p3_g;
  TVector3       p3_lep;
  TLorentzVector p4_g;
  TLorentzVector p4_lep;
  TLorentzVector p4_nu(0., 0., sig_param.E_nu_probe, sig_param.E_nu_probe); // probe along z

  KinePhaseSpace_t t = kPSEgTlTgPgfE;

  for( unsigned int ang = 0; ang < theta_gamma.size(); ang++ ) {
    sig_param.theta_g = theta_gamma[ang]*TMath::DegToRad();
    sig_param.E_lep = sig_param.E_nu_probe;
    for( unsigned int eg = 0; eg<steps; eg++ ) {
                                                                                                         
      double E_g = sig_param.E_nu_probe - sig_param.E_lep;
      // Set FS gamma 4-mom.
      p3_g.SetMagThetaPhi( E_g, sig_param.theta_g*TMath::DegToRad(), 
                                sig_param.phi_g*TMath::DegToRad() );
      p4_g.SetVect( p3_g );
      p4_g.SetE( E_g );
      
      // Set FS lepton (nu) 4-mom.
      p3_lep.SetMagThetaPhi(  sig_param.E_lep, sig_param.theta_lep*TMath::DegToRad(), 
                                               sig_param.phi_lep*TMath::DegToRad() );
      p4_lep.SetVect( p3_lep );
      p4_lep.SetE( sig_param.E_lep );
                                                                                                        
      kine->SetHadSystP4( p4_g );
      kine->SetFSLeptonP4( p4_lep );
                 
      TLorentzVector q = p4_nu-p4_lep ;
      double Q2 = -q.Mag2();
      double x = Q2/(2 * E_g * constants::kNucleonMass );
      double y = E_g/sig_param.E_nu_probe;
      kine->Setx( x );
      kine->Sety( y );
      utils::kinematics::UpdateWQ2FromXY(i);
      double t = (q - p4_g).Mag() ;
      kine->Sett( t ) ;
      cout << "E_g " << E_g << " sig_param.E_nu_probe " << sig_param.E_nu_probe << endl;                                                                                            
      // Check the validity
      cout << " Gamma theta angle " << theta_gamma[ang]
           << " Gamma Energy " << E_g << " Q2 " << Q2 << " x " << x << " y " << y
           << " Valid process: " << xsec_alg->ValidProcess(i) 
           << " Valid Kinematics: " <<  xsec_alg->ValidKinematics(i) << endl;
      
      // Calculate Genie Xsection
      xsec_arr[ang][eg] = (xsec_alg->XSec(i, t) * 1.e41)/units::cm2;
      cout << "---> Cross section " << xsec_arr[ang][eg] << endl;
      LOG("gtestCOHNCGammaXSec", pINFO ) << "---> Cross section " << xsec_arr[ang][eg];
      // Calculate external Xsection
      // If invalid process or kinematics cross section is 0
      if ( !xsec_alg->ValidProcess(i) || !xsec_alg->ValidKinematics(i) ) {
         ext_xsec_arr[ang][eg] = 0.;
      } else {
        ext_xsec_arr[ang][eg] = ( ExtXSec( 14, 1000180400, sig_param.E_nu_probe, sig_param.E_lep,  
                                  sig_param.theta_lep*TMath::DegToRad(), 
                                  sig_param.theta_g*TMath::DegToRad(),
                                  sig_param.phi_g*TMath::DegToRad() )*1.e41 ) / units::cm2 ;
      } 

      Eg_arr[eg] = E_g;
      sig_param.E_lep -= sig_param.range / double( steps-1 );
    } // gamma E
  } // gamma angle

  char text[50];
  for( unsigned int k = 0; k < theta_gamma.size(); k++ ) {
    auto gr = new TGraph( steps, Eg_arr, xsec_arr[k] ); gr->SetLineColor( root_color[k] ); gr->SetLineWidth(2);
    snprintf( text, 50, "#theta_{#gamma} = %d",(int)theta_gamma[k] );
    legend->AddEntry( gr, (char*)text, "l" );
    mg->Add(gr);
  }

  mg->SetTitle("Genie ^{40}Ar #bar{#nu}_{#mu} #theta_{l}=1 #phi_{l}=0 #phi_{#gamma}=10; E_{#gamma} (GeV); XSec #sigma (x10^{-41})");
  c1->cd(1);
  mg->Draw("ACP");
  legend->Draw();
  mg->Write();

  for( unsigned int m = 0; m < theta_gamma.size(); m++ ) {
    auto egr = new TGraph( steps, Eg_arr, ext_xsec_arr[m] ); egr->SetLineColor( root_color[m] ); egr->SetLineWidth(2);
    ext_mg->Add(egr);
  }
  ext_mg->SetTitle("Ext ^{40}Ar #bar{#nu}_{#mu} #theta_{l}=1 #phi_{l}=0 #phi_{#gamma}=10; E_{#gamma} (GeV); XSec #sigma (x10^{-41})");
  c1->cd(2);
  ext_mg->Draw("ACP");
}
//__________________________________________________________________________
double ExtXSec(int nu_probe, int tgt, double Enu, double Enu_final, double theta_l, double theta_g, double phi_g)
{
  std::string probe;
  std::string nucleus;

  if (nu_probe == kPdgNuMu)          probe = "nu";
  else if (nu_probe == kPdgAntiNuMu) probe = "anti-nu";
  else return 0.;

  if (tgt == 1000060120) nucleus = "12C";
  if (tgt == 1000180400) nucleus = "40Ar";
  else return 0.;

  CS.setMode(probe); // Set nu or anti-nu
  CS.setNucleus(nucleus); // Set target nucleus
  double xsec = CS.getDiffCrossSection(Enu, Enu_final, theta_l, theta_g, phi_g);
  return xsec; // Return the xsection from external C++ code
}
//__________________________________________________________________________
void SigmaCosTh(const Interaction *i, const ARSXSec * xsec_alg, Kinematics * kine, SigVars sig_param)
{
}
