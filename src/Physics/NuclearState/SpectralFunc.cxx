//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

         Steven Gardiner <gardiner \at fnal.gov>
         Fermi National Accelerator Laboratory

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ May 01, 2012 - CA
   Pick spectral function data from $GENIE/data/evgen/nucl/spectral_functions

 @ April 29, 2019 - SG
   Full rewrite to use finely-binned histograms (tightly-peaked distributions led
   to very slow rejection sampling), XML-based configuration of spectral function
   data files with lazy initialization
*/
//____________________________________________________________________________

#include <fstream>
#include <sstream>
#include <string>

#include "TF2.h"
#include "TGraph2D.h"
#include "TH2D.h"
#include "TNtupleD.h"
#include "TRandom.h"
#include "TSystem.h"

#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/NuclearState/SpectralFunc.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Numerical/RandomGen.h"

namespace {

  // Replace this with std::to_string when we switch to C++11
  std::string replace_with_std_to_string(int an_integer) {
    std::ostringstream oss;
    oss << an_integer;
    return oss.str();
  }

}

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;

//____________________________________________________________________________
SpectralFunc::SpectralFunc() : NuclearModelI("genie::SpectralFunc")
{

}
//____________________________________________________________________________
SpectralFunc::SpectralFunc(string config) :
  NuclearModelI("genie::SpectralFunc", config)
{

}
//____________________________________________________________________________
SpectralFunc::~SpectralFunc()
{
  // Delete the TH2D objects from the spectral function map
  std::map<int, TH2D*>::iterator begin = fSpectralFunctionMap.begin();
  std::map<int, TH2D*>::iterator end = fSpectralFunctionMap.end();
  for (std::map<int, TH2D*>::iterator iter = begin; iter != end; ++iter) {
    TH2D* hist = iter->second;
    if ( hist ) delete hist;
  }
}
//____________________________________________________________________________
bool SpectralFunc::GenerateNucleon(const Target& target) const
{
  TH2D* sf = this->SelectSpectralFunction( target );

  if ( !sf ) {
    fCurrRemovalEnergy = 0.;
    fCurrMomentum.SetXYZ(0., 0., 0.);
    return false;
  }

  double kmin    = sf->GetXaxis()->GetXmin(); // momentum range
  double kmax    = sf->GetXaxis()->GetXmax();
  double wmin    = sf->GetYaxis()->GetXmin(); // removal energy range
  double wmax    = sf->GetYaxis()->GetXmax();

  LOG("SpectralFunc", pDEBUG) << "Momentum range = ["   << kmin << ", " << kmax << "]";
  LOG("SpectralFunc", pDEBUG) << "Rmv energy range = [" << wmin << ", " << wmax << "]";

  // Temporarily replace ROOT's gRandom RNG with GENIE's so that our call to
  // TH2::GetRandom2() will use GENIE's random numbers
  TRandom* old_gRandom = gRandom;
  RandomGen* rnd = RandomGen::Instance();
  gRandom = &rnd->RndGen();

  // Pick a random nucleon momentum and removal energy using the spectral
  // function histogram
  double kc = 0.;
  double wc = 0.;
  sf->GetRandom2( kc, wc );

  // Restore ROOT's old gRandom now that we've had our fun
  gRandom = old_gRandom;

  // Log our choices for the nucleon momentum and removal energy
  LOG("SpectralFunc", pINFO) << "|p,nucleon| = " << kc;
  LOG("SpectralFunc", pINFO) << "|w,nucleon| = " << wc;

  // Generate momentum components
  double costheta = -1. + 2. * rnd->RndGen().Rndm();
  double sintheta = TMath::Sqrt(1.-costheta*costheta);
  double fi       = 2 * kPi * rnd->RndGen().Rndm();
  double cosfi    = TMath::Cos(fi);
  double sinfi    = TMath::Sin(fi);

  double kx = kc * sintheta * cosfi;
  double ky = kc * sintheta * sinfi;
  double kz = kc * costheta;

  // set generated values
  fCurrRemovalEnergy = wc;
  fCurrMomentum.SetXYZ(kx,ky,kz);

  return true;
}
//____________________________________________________________________________
double SpectralFunc::Prob(double p, double w, const Target & target) const
{
  TH2D* sf = this->SelectSpectralFunction( target );
  if ( !sf ) return 0.;

  // TODO: check whether we should return the probability mass or density here.
  // The previous implementation used the density. For now, we return the mass.
  int bin_num = sf->FindBin( p, w );
  double prob_mass = sf->GetBinContent( bin_num );

  return prob_mass;

  //int bx, by, bz;
  //sf->GetBinXYZ( bin_num, bx, by, bz );

  //double prob_density = prob_mass / sf->GetXaxis()->GetBinWidth( bx )
  //  / sf->GetYaxis()->GetBinWidth( by )
  //  / std::pow(sf->GetXaxis()->GetBinLowEdge( bx ), 2)
  //  / (4. * kPi);

  //return prob_density;
}
//____________________________________________________________________________
void SpectralFunc::Configure(const Registry & config)
{
  Algorithm::Configure( config );
  this->LoadConfig();
}
//____________________________________________________________________________
void SpectralFunc::Configure(string param_set)
{
  Algorithm::Configure( param_set );
  this->LoadConfig();
}
//____________________________________________________________________________
void SpectralFunc::LoadConfig(void)
{
  // Determine the folder that should be searched to find the spectral function
  // data files
  std::string data_path_type;
  this->GetParamDef( "DataPathType", data_path_type, std::string("relative") );

  std::string data_path;
  this->GetParam( "DataPath", data_path );
  // Ensure that the path ends with '/', just in case
  data_path += '/';

  if ( data_path_type == "relative" ) {
    std::string genie_path = gSystem->Getenv("GENIE");
    genie_path += '/';
    data_path = genie_path + data_path;
  }

  fDataPath = data_path;
}
//____________________________________________________________________________
TH2D* SpectralFunc::SelectSpectralFunction(const Target& t) const
{
  // First check whether we've already built the requested spectral function
  int target_pdg = t.Pdg();
  if ( fSpectralFunctionMap.count(target_pdg) ) {
    return fSpectralFunctionMap.at( target_pdg );
  }

  // If not, attempt to build it
  std::string target_pdg_string = replace_with_std_to_string( target_pdg );
  RgKey sf_key( "SpectFuncTable@Pdg=" + target_pdg_string );
  std::string data_filename;
  this->GetParamDef( sf_key, data_filename, std::string() );

  if ( data_filename.empty() ) {
    LOG("SpectralFunc", pERROR) << "** The spectral function for target "
      << target_pdg << " isn't available";
    std::exit( 1 );
  }

  // Prepend the data path to the file name
  data_filename = fDataPath + data_filename;

  TH2D* sf_hist = this->LoadSFDataFile(data_filename);

  LOG("SpectralFunc", pNOTICE) << "Loaded spectral function data"
    " for target with PDG code " << target_pdg << " from the file "
    << data_filename;

  sf_hist->SetName( ("sf_" + target_pdg_string).c_str() );

  // Set the directory to NULL so that this histogram is never auto-deleted
  // (the genie::SpectralFunc object will take ownership)
  sf_hist->SetDirectory( NULL );

  // Store the new spectral function histogram for easy retrieval later
  fSpectralFunctionMap[ target_pdg ] = sf_hist;

  return sf_hist;
}
//____________________________________________________________________________
TH2D* SpectralFunc::LoadSFDataFile(const std::string& full_file_name) const
{
  int num_E_bins, num_p_bins;
  std::ifstream in_file( full_file_name );

  // TODO: add I/O error handling

  // First pass: get the number of bins on each axis
  // and the set of bin edges
  // TODO: consider changing the format so that this can
  // be conveniently done in a single read-through of the
  // input file
  in_file >> num_E_bins >> num_p_bins;
  double p, E, prob_density;
  std::vector<double> ps, Es;
  for ( int ip = 0; ip < num_p_bins; ++ip ) {

    in_file >> p;
    ps.push_back( p * genie::units::MeV );

    for ( int ie = 0; ie < num_E_bins; ++ie ) {
      in_file >> E >> prob_density;
      if ( ip == 0 ) Es.push_back( E * genie::units::MeV );
    }
  }

  // Assume that the last bin along each axis has the same width
  // as the first bin (Noemi's SF tables do not give the last bin's
  // high edge explicitly). Defining a TH2 with variable-size bins
  // requires these high edge values.
  // TODO: revisit this as needed
  double p_width = ps.at(1) - ps.front();
  ps.push_back( ps.back() + p_width );

  double E_width = Es.at(1) - Es.front();
  Es.push_back( Es.back() + E_width );

  in_file.close();

  // We've built the grid, now make a histogram using it
  TH2D* sf_hist = new TH2D("temp_sf_hist",
    "Spectral Function; nucleon momentum (GeV); removal energy (GeV)",
    num_p_bins, &(ps.front()), num_E_bins, &(Es.front()));

  // Do a second pass and set the contents of each bin
  in_file.open( full_file_name );
  in_file >> num_E_bins >> num_p_bins;
  for ( int ip = 0; ip < num_p_bins; ++ip ) {

    in_file >> p;

    // Convert from MeV to GeV
    p *= genie::units::MeV;

    for ( int ie = 0; ie < num_E_bins; ++ie ) {

      in_file >> E >> prob_density;

      // Convert from MeV to GeV
      E *= genie::units::MeV;

      // Convert bin contents from probability density to probability
      // mass for easier sampling
      TAxis* p_axis = sf_hist->GetXaxis();
      int p_idx = p_axis->FindBin( p );
      double p_bin_width = p_axis->GetBinWidth( p_idx );

      TAxis* E_axis = sf_hist->GetYaxis();
      int E_idx = E_axis->FindBin( E );
      double E_bin_width = E_axis->GetBinWidth( E_idx );

      double prob_mass = prob_density * std::pow(p, 2) * 4. * kPi
        * p_bin_width * E_bin_width;

      sf_hist->SetBinContent( sf_hist->FindBin(p, E), prob_mass );
      LOG("SpectralFunc", pDEBUG) << "p = " << p << ", E = " << E << ", dens = " << prob_density
        << ", mass = " << prob_mass << ", p_bin_width = " << p_bin_width << ", E_bin_width = "
        << E_bin_width;
    }
  }

  return sf_hist;
}
