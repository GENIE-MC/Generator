#include "TMath.h"
#include "TSystem.h"
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
using namespace std;

#include "TStyle.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLorentzVector.h"

#include "FluxDrivers/GNuMIFlux.h"
#include "PDG/PDGCodeList.h"
#include "PDG/PDGUtils.h"

using namespace genie::flux;

#include "FluxDrivers/GNuMIFlux.h"
#include "PDG/PDGCodeList.h"
#include "PDG/PDGUtils.h"

// typedefs and forward declarations
typedef genie::flux::GNuMIFluxPassThroughInfo fluxentry_t;

double fabserrX(double a, double b) 
{ return TMath::Abs(a-b)/TMath::Max(TMath::Abs(b),1.0e-30); }

//////////////////////////////////////////////////////////////////////////////

int calc_enuwgt(const fluxentry_t& fluxentry,
                double x, double y, double z,
                double& enu, double& wgt_xy, xypartials& partials);

//////////////////////////////////////////////////////////////////////////////

void test_gnumi(int mxpull = 0)
{
  // mxpull:  max number of entries to process 0=all, -N=until N mismatches

  // currently TFile not TChain, so better resolve to a single file
  string fluxpatt = "";
  fluxpatt = "$FLUXPATH/v19/fluka05_le010z185i/root/fluka05_le010z185i_1.root";
  //fluxpatt = "$FLUXPATH/g4numi/g4fluka_le010z185i_leak_v10_0001.root";
  //fluxpatt = "$GENIE/src/FluxDrivers/gnumi_test/*g3*.root";
  //fluxpatt = "*g3*.root";
  string fluxfilename(gSystem->ExpandPathName(fluxpatt.c_str()));

  genie::flux::GNuMIFlux* gnumi = new genie::flux::GNuMIFlux();
  gnumi->LoadBeamSimData(fluxfilename);
  gnumi->SetGenWeighted(true);  // generated entries w/ fWeight != 1.0

  const genie::PDGCodeList& pdgList = gnumi->FluxParticles();
  int nfluxtypes = pdgList.size();
  cout << "GNuMIFlux supplies " << nfluxtypes << " PDG particles of type: ";
  // don't try iterators .. apparently they aren't usable in CINT
  for (int indx=0; indx < nfluxtypes ; ++indx) 
    cout << pdgList[indx] << " ";
  cout << endl;

  cout << "MaxEnergy " << gnumi->MaxEnergy() << endl;

  cout << "scan for max weight at NearDet" << endl;
  //gnumi->UseFluxAtNearDetCenter();
  gnumi->SetLengthUnits(genie::flux::GNuMIFlux::kcm);
  gnumi->SetBeamFluxWindow(genie::flux::GNuMIFlux::kMinosNearCenter);
  gnumi->ScanForMaxWeight();

  int npull, nxy, file_debug;
  double xbeam, ybeam, zbeam = 103935.; // distance in cm to near "window"

  ifstream xypoints("XYPOINTS.TXT");
  
  xypoints >> npull >> nxy >> zbeam >> file_debug;

  TH1D* hist_ferr_enu = new TH1D("ferr_enu","ferr_enu",1000,0.0,0.000005);
  TH1D* hist_ferr_wgt = new TH1D("ferr_wgt","ferr_wgt",1000,0.0,0.01);
  hist_ferr_enu->Sumw2();
  hist_ferr_wgt->Sumw2();
  TH1D* hist_ferr_enu0 = new TH1D("ferr_enu0","ferr_enu0",1000,0.0,0.000005);
  TH1D* hist_ferr_wgt0 = new TH1D("ferr_wgt0","ferr_wgt0",1000,0.0,0.01);
  hist_ferr_enu0->Sumw2();
  hist_ferr_wgt0->Sumw2();
  gStyle->SetOptStat(111111);

  if ( mxpull > 0 && npull > mxpull ) npull = mxpull;

  for ( int ipull = 1; ipull <= npull ; ++ipull ) {

    bool genok = gnumi->GenerateNext();
    const fluxentry_t& fluxentry = gnumi->PassThroughInfo();

    double wgtprd_gen = gnumi->Weight();
    double enu_gen = gnumi->Momentum().E();
    if ( wgtprd_gen != 1.0 ) {
      // weighted flux ... we can validate "center" values
      // assume used kMinosNearCenter in SetBeamFluxWindow
      double wgtprd_ntuple = fluxentry.nimpwt * fluxentry.nwtnear;
      double enu_ntuple    = fluxentry.nenergyn;
      double ferr_enu0     = fabserrX(enu_gen,enu_ntuple);
      double ferr_wgtprd0  = fabserrX(wgtprd_gen,wgtprd_ntuple);
      const double eps0h = 1.0e-6; // 2.5e-5;
      if ( ferr_enu0 > eps0h || ferr_wgtprd0 > eps0h ) {
        hist_ferr_enu0->Fill(ferr_enu0,1.0);
        hist_ferr_wgt0->Fill(ferr_wgtprd0,1.0);
      }
      const double eps0 = 2.5e-5; // 2.5e-5;
      if ( ferr_enu0 > eps0 || ferr_wgtprd0 > eps0 ) {
        char xenu0 = ( ferr_enu0 > eps0 ) ? '#' : ' ';
        char xwgt0 = ( ferr_wgtprd0 > eps0 ) ? '#' : ' ';
        cout << "Checking \"center\" energy/weight " 
             << " enu " << setw(8) << enu_gen << " " << setw(8) << enu_ntuple
             << " err " << setw(11) << ferr_enu0 << xenu0
             << " wgtprd " << setw(8) << wgtprd_gen << " " << setw(8) << wgtprd_ntuple
             << " err " << setw(11) << ferr_wgtprd0 << xwgt0
             << endl;
      }

    }


    //bool atend = gnumi->End();
    //cout << "gnumi->GenerateNext() returned " << (genok?"true":"false") 
    //     << " End() returns " << (atend?"true":"false")
    //     << endl;

    //cout << "Current entry: "
    //      << " PdgCode=" << gnumi->PdgCode()
    //     << " Weight=" << gnumi->Weight()
    //     << endl;
    // gnumi->PrintCurrent();

    for ( int ixy = 0; ixy <= nxy ; ++ixy ) {
      xypartials part_file, part_calc;
      if ( file_debug ) part_file.ReadStream(xypoints);

      double enu, wgt_xy;
      double enu_in, wgt_xy_in;
      int ipull_in, ixy_in;
      xypoints >> ipull_in >> ixy_in >> xbeam >> ybeam >> enu_in >> wgt_xy_in;
      if ( ! xypoints.good() || ( ipull_in != ipull ) || ( ixy_in != ixy) ) {
        cout << "HEY!!! XYPOINTS.TXT "
             << " file status = " << (int)xypoints.good() 
             << " ipull " << ipull_in << " " << ipull 
             << " ixy " << ixy_in << " " << ixy 
             << endl;
        return;
      }
      fluxentry.CalcEnuWgt(xbeam,ybeam,zbeam,enu,wgt_xy,part_calc);

      double ferr_enu    = fabserrX(enu,enu_in);
      double ferr_wgt_xy = fabserrX(wgt_xy,wgt_xy_in);
      const double eps_enu = 2.5e-5; // 0.003; // 6.0e-4; // 2.5e-4; // 1.5e-5; // 1.0e-6;
      const double eps_wgt = 2.5e-5; // 0.003; // 6.0e-4; // 2.5e-4; // 2.5e-5; //
      if ( ferr_enu > eps_enu || ferr_wgt_xy > eps_wgt ) {
        cout << endl;
        char xenu = ( ferr_enu > eps_enu ) ? '#':' ';
        char xwgt = ( ferr_wgt_xy > eps_wgt ) ? '#':' ';
        cout << "pull " << setw(6) << ipull << " ixy " << setw(2) << ixy 
             << " enu " << setw(8) << enu << " " << setw(8) << enu_in 
             << " err " << setw(11) << ferr_enu << xenu
             << " wgt_xy " << setw(11) << wgt_xy << " " << setw(11) << wgt_xy_in 
             << " err " << setw(11) << ferr_wgt_xy << xwgt
             << "  ptype " << setw(4) << part_file.ptype 
             << " ntype " << setw(3) << part_file.ntype << endl;
        hist_ferr_enu->Fill(ferr_enu,1.0);
        hist_ferr_wgt->Fill(ferr_wgt_xy,1.0);
      } 
      //else {
      //  cout << "pull " << ipull << " ixy " << ixy 
      //       << " okay ptype " << part_file.ptype << " ntype " << part_file.ntype << endl;
      //}

      if (file_debug) {
        int np = part_file.Compare(part_calc);
        if (np>0) {
          std::cout << fluxentry << endl;
          part_file.Print();
          part_calc.Print();

          if ( mxpull < 0 ) {
            cout << "COUNTDOWN mismatch error count @ " << mxpull << endl;
            mxpull += 1;
            if ( mxpull == 0 ) {
              cout << "REACHED LIMIT ON MISMATCHES ... exit" << endl;
              return; // reached limit on errors, asked to stop on error
            }
          }
        }
      }


      //if ( part_calc.ptype == 13 || part_calc.ptype == -13 ) {
      //  cout << "STOPPING HERE FOR INSPECTION" << endl;
      //  return;
      //}

    } //loop over xy points for given entryç
  } // loop over entries (pulls)

  TCanvas* c1 = new TCanvas("c1","all positions");
  c1->Divide(1,2);
  c1->cd(1);
  gPad->SetLogy(true);
  hist_ferr_enu->Draw();
  c1->cd(2);
  gPad->SetLogy(true);
  hist_ferr_wgt->Draw();

  TCanvas* c2 = new TCanvas("c2","center weight");
  c2->Divide(1,2);
  c2->cd(1);
  gPad->SetLogy(true);
  hist_ferr_enu0->Draw();
  c2->cd(2);
  gPad->SetLogy(true);
  hist_ferr_wgt0->Draw();

}

//////////////////////////////////////////////////////////////////////////////
int calc_enuwgt(const fluxentry_t& fluxentry, 
                double xpos, double ypos, double zpos,
                double& enu, double& wgt_xy, xypartials& partials)
{

  // Neutrino Energy and Weigth at arbitrary point
  // based on:
  //   NuMI-NOTE-BEAM-0109 (MINOS DocDB # 109)
  //   Title:   Neutrino Beam Simulation using PAW with Weighted Monte Carlos
  //   Author:  Rick Milburn
  //   Date:    1995-10-01

  // history:
  // jzh  3/21/96 grab R.H.Milburn's weighing routine
  // jzh  5/ 9/96 substantially modify the weighting function use dot product 
  //              instead of rotation vecs to get theta get all info except 
  //              det from ADAMO banks neutrino parent is in Particle.inc
  //              Add weighting factor for polarized muon decay
  // jzh  4/17/97 convert more code to double precision because of problems 
  //              with Enu>30 GeV
  // rwh 10/ 9/08 transliterate function from f77 to C++

  // original function description:
  //   Real function for use with PAW Ntuple To transform from destination
  //   detector geometry to the unit sphere moving with decaying hadron with
  //   velocity v, BETA=v/c, etc..  For (pseudo)scalar hadrons the decays will
  //   be isotropic in this  sphere so the fractional area (out of 4-pi) is the
  //   fraction of decays that hit the target.  For a given target point and 
  //   area, and given x-y components of decay transverse location and slope,
  //   and given decay distance from target ans given decay GAMMA and 
  //   rest-frame neutrino energy, the lab energy at the target and the 
  //   fractional solid angle in the rest-frame are determined.
  //   For muon decays, correction for non-isotropic nature of decay is done.

  // Arguments:
  //    genie::flux::GNuMIFluxPassThroughInfo fluxentry :: ntuple entry info
  //    double x, y, z :: position to evaluate for (enu,wgt_xy) 
  //        in *beam* frame coordinates  (?units)
  //    double enu, wgt_xy :: resulting energy and weight
  // Return:
  //    int :: error code
  // Assumptions:
  //    Energies given in GeV
  //    Particle codes have been translated from GEANT into PDG codes

  // for now ... these _should_ come from DB
  // but use these hard-coded values to "exactly" reproduces old code
  const double kPIMASS = 0.13957;
  const double kKMASS  = 0.49368;
  const double kK0MASS = 0.49767;
  const double kMUMASS = 0.105658389;

  const int kpdg_nue       =   12;  // extended Geant 53
  const int kpdg_nuebar    =  -12;  // extended Geant 52
  const int kpdg_numu      =   14;  // extended Geant 56
  const int kpdg_numubar   =  -14;  // extended Geant 55

  const int kpdg_muplus    =  -13;  // Geant  5
  const int kpdg_muminus   =   13;  // Geant  6
  const int kpdg_pionplus  =  211;  // Geant  8
  const int kpdg_pionminus = -211;  // Geant  9
  const int kpdg_k0long    =  130;  // Geant 10  ( K0=311, K0S=310 )
  const int kpdg_k0short   =  310;  // Geant 16
  const int kpdg_k0mix     =  311;  
  const int kpdg_kaonplus  =  321;  // Geant 11
  const int kpdg_kaonminus = -321;  // Geant 12

  const double kRDET = 100.0;   // set to flux per 100 cm radius

  enu    = 0.0;  // don't know what the final value is
  wgt_xy = 0.0;  // but set these in case we return early due to error


  // in principle we should get these from the particle DB
  // but for consistency testing use the hardcoded values
  double parent_mass = kPIMASS;
  switch ( fluxentry.ptype ) {
  case kpdg_pionplus:
  case kpdg_pionminus:
    parent_mass = kPIMASS;
    break;
  case kpdg_kaonplus:
  case kpdg_kaonminus:
    parent_mass = kKMASS;
    break;
  case kpdg_k0long:
  case kpdg_k0short:
  case kpdg_k0mix:
    parent_mass = kK0MASS;
    break;
  case kpdg_muplus:
  case kpdg_muminus:
    parent_mass = kMUMASS;
    break;
  default:
    std::cerr << "NU_REWGT unknown particle type " << fluxentry.ptype
              << std::endl;
    return 1;
  }

  double parentp2 = ( fluxentry.pdpx*fluxentry.pdpx +
                      fluxentry.pdpy*fluxentry.pdpy +
                      fluxentry.pdpz*fluxentry.pdpz );
  double parent_energy = TMath::Sqrt( parentp2 +
                                     parent_mass*parent_mass);
  double parentp = TMath::Sqrt( parentp2 );

  double gamma     = parent_energy / parent_mass;
  double gamma_sqr = gamma * gamma;
  double beta_mag  = TMath::Sqrt( ( gamma_sqr - 1.0 )/gamma_sqr );

  // Get the neutrino energy in the parent decay CM
  double enuzr = fluxentry.necm;
  // Get angle from parent line of flight to chosen point in beam frame
  double rad = TMath::Sqrt( (xpos-fluxentry.vx)*(xpos-fluxentry.vx) +
                            (ypos-fluxentry.vy)*(ypos-fluxentry.vy) +
                            (zpos-fluxentry.vz)*(zpos-fluxentry.vz) );

  double costh_pardet = ( fluxentry.pdpx*(xpos-fluxentry.vx) +
                          fluxentry.pdpy*(ypos-fluxentry.vy) +
                          fluxentry.pdpz*(zpos-fluxentry.vz) ) 
                        / ( parentp * rad);
  if ( costh_pardet >  1.0 ) costh_pardet =  1.0;
  if ( costh_pardet < -1.0 ) costh_pardet = -1.0;
  double theta_pardet = TMath::ACos(costh_pardet);

  // Weighted neutrino energy in beam, approx, good for small theta
  //RWH//double emrat = 1.0 / ( gamma * ( 1.0 - beta_mag * TMath::Cos(theta_pardet)));
  double emrat = 1.0 / ( gamma * ( 1.0 - beta_mag * costh_pardet ));
  enu = emrat * enuzr;  // ! the energy ... normally

  // RWH-debug
  bool debug = false;
  if (debug) {
    cout << setprecision(15);
    cout << "ptype " << fluxentry.ptype << " m " << parent_mass 
         << " p " << parentp << " e " << parent_energy << " gamma " << gamma
         << " beta " << beta_mag << endl;

    cout << " enuzr " << enuzr << " rad " << rad << " costh " << costh_pardet
         << " theta " << theta_pardet << " emrat " << emrat 
         << " enu " << enu << endl;
  }
  partials.parent_mass   = parent_mass;
  partials.parentp       = parentp;
  partials.parent_energy = parent_energy;
  partials.gamma         = gamma;
  partials.beta_mag      = beta_mag;
  partials.enuzr         = enuzr;
  partials.rad           = rad;
  partials.costh_pardet  = costh_pardet;
  partials.theta_pardet  = theta_pardet;
  partials.emrat         = emrat;
  partials.eneu          = enu;

  // Get solid angle/4pi for detector element
  double sangdet = ( kRDET*kRDET / ( (zpos-fluxentry.vz)*(zpos-fluxentry.vz)))/4.0;

  // Weight for solid angle and lorentz boost
  wgt_xy = sangdet * ( emrat * emrat );  // ! the weight ... normally

  partials.sangdet       = sangdet;
  partials.wgt           = wgt_xy;
  partials.ptype         = fluxentry.ptype; // assume already PDG

  // Done for all except polarized muon decay
  // in which case need to modify weight 
  // (must be done in double precision)
  if ( fluxentry.ptype == kpdg_muplus || fluxentry.ptype == kpdg_muminus) {
    double beta[3], p_dcm_nu[4], p_nu[3], p_pcm_mp[3], partial;

    // Boost neu neutrino to mu decay CM
    beta[0] = fluxentry.pdpx / parent_energy;
    beta[1] = fluxentry.pdpy / parent_energy;
    beta[2] = fluxentry.pdpz / parent_energy;
    p_nu[0] = (xpos-fluxentry.vx)*enu/rad;
    p_nu[1] = (ypos-fluxentry.vy)*enu/rad;
    p_nu[2] = (zpos-fluxentry.vz)*enu/rad;
    partial = gamma * 
      (beta[0]*p_nu[0] + beta[1]*p_nu[1] + beta[2]*p_nu[2] );
    partial = enu - partial/(gamma+1.0);
    p_dcm_nu[0] = p_nu[0] - beta[0]*gamma*partial;
    p_dcm_nu[1] = p_nu[1] - beta[1]*gamma*partial;
    p_dcm_nu[2] = p_nu[2] - beta[2]*gamma*partial;
    p_dcm_nu[3] = TMath::Sqrt( p_dcm_nu[0]*p_dcm_nu[0] +
                               p_dcm_nu[1]*p_dcm_nu[1] +
                               p_dcm_nu[2]*p_dcm_nu[2] );

    partials.betanu[0]     = beta[0];
    partials.betanu[1]     = beta[1];
    partials.betanu[2]     = beta[2];
    partials.p_nu[0]       = p_nu[0];
    partials.p_nu[1]       = p_nu[1];
    partials.p_nu[2]       = p_nu[2];
    partials.partial1      = partial;
    partials.p_dcm_nu[0]   = p_dcm_nu[0];
    partials.p_dcm_nu[1]   = p_dcm_nu[1];
    partials.p_dcm_nu[2]   = p_dcm_nu[2];
    partials.p_dcm_nu[3]   = p_dcm_nu[3];

    // Boost parent of mu to mu production CM
    double particle_energy = fluxentry.ppenergy;
    gamma = particle_energy/parent_mass;
    beta[0] = fluxentry.ppdxdz * fluxentry.pppz / particle_energy;
    beta[1] = fluxentry.ppdydz * fluxentry.pppz / particle_energy;
    beta[2] =                    fluxentry.pppz / particle_energy;
    partial = gamma * ( beta[0]*fluxentry.muparpx + 
                        beta[1]*fluxentry.muparpy + 
                        beta[2]*fluxentry.muparpz );
    partial = fluxentry.mupare - partial/(gamma+1.0);
    p_pcm_mp[0] = fluxentry.muparpx - beta[0]*gamma*partial;
    p_pcm_mp[1] = fluxentry.muparpy - beta[1]*gamma*partial;
    p_pcm_mp[2] = fluxentry.muparpz - beta[2]*gamma*partial;
    double p_pcm = TMath::Sqrt ( p_pcm_mp[0]*p_pcm_mp[0] +
                                 p_pcm_mp[1]*p_pcm_mp[1] +
                                 p_pcm_mp[2]*p_pcm_mp[2] );

    //cout << " muparpxyz " << fluxentry.muparpx << " "
    //     << fluxentry.muparpy << " " << fluxentry.muparpz << endl;
    //cout << " beta " << beta[0] << " " << beta[1] << " " << beta[2] << endl;
    //cout << " gamma " << gamma << " partial " << partial << endl;
    //cout << " p_pcm_mp " << p_pcm_mp[0] << " " << p_pcm_mp[1] << " " << p_pcm_mp[2] << " " << p_pcm << endl;
    partials.muparent_px   = fluxentry.muparpx;
    partials.muparent_py   = fluxentry.muparpy;
    partials.muparent_pz   = fluxentry.muparpz;
    partials.gammamp       = gamma;
    partials.betamp[0]     = beta[0];
    partials.betamp[1]     = beta[1];
    partials.betamp[2]     = beta[2];
    partials.partial2      = partial;
    partials.p_pcm_mp[0]   = p_pcm_mp[0];
    partials.p_pcm_mp[1]   = p_pcm_mp[1];
    partials.p_pcm_mp[2]   = p_pcm_mp[2];
    partials.p_pcm         = p_pcm;

    const double eps = 1.0e-30;  // ? what value to use
    if ( p_pcm < eps || p_dcm_nu[3] < eps ) {
      return 3; // mu missing parent info?
    }
    // Calc new decay angle w.r.t. (anti)spin direction
    double costh = ( p_dcm_nu[0]*p_pcm_mp[0] +
                     p_dcm_nu[1]*p_pcm_mp[1] +
                     p_dcm_nu[2]*p_pcm_mp[2] ) /
                   ( p_dcm_nu[3]*p_pcm );
    if ( costh >  1.0 ) costh =  1.0;
    if ( costh < -1.0 ) costh = -1.0;
    // Calc relative weight due to angle difference
    double wgt_ratio = 0.0;
    switch ( fluxentry.ntype ) {
    case kpdg_nue:
    case kpdg_nuebar:
      wgt_ratio = 1.0 - costh;
      break;
    case kpdg_numu:
    case kpdg_numubar:
      double xnu = 2.0 * enuzr / kMUMASS;
      wgt_ratio = ( (3.0-2.0*xnu )  - (1.0-2.0*xnu)*costh ) / (3.0-2.0*xnu);
      break;
    default:
      return 2; // bad neutrino type
    }
    wgt_xy = wgt_xy * wgt_ratio;

    partials.ntype     = fluxentry.ntype; // assume converted to PDG
    partials.costhmu   = costh;
    partials.wgt_ratio = wgt_ratio;

  }

  return 0;
}
