//____________________________________________________________________________
/*
  Copyright (c) 2003-2025, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org
  
  Author: John Plows <komninos-john.plows \at physics.ox.ac.uk>
          University of Oxford

	  Costas Andreopoulos <c.andreopoulos \at cern.ch>
	  University of Liverpool
*/
//____________________________________________________________________________

#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/PrintUtils.h"
#include "Physics/BeamHNL/HNLDecayUtils.h"

using namespace genie;
using namespace genie::hnl;

//____________________________________________________________________________
string genie::utils::hnl::ProdAsString(HNLProd_t hnlprod)
{
  switch(hnlprod) {

  case (kHNLProdPion2Muon):
    return "pi+- --> N + mu+-";
    break;
    
  case (kHNLProdPion2Electron):
    return "pi+- --> N + e+-";
    break;

  case (kHNLProdKaon2Muon):
    return "K+- --> N + mu+-";
    break;

  case (kHNLProdKaon2Electron):
    return "K+- --> N + e+-";
    break;

  case (kHNLProdKaon3Muon):
    return "K+- --> N + mu+- + pi0";
    break;

  case (kHNLProdKaon3Electron):
    return "K+- --> N + e+- + pi0";
    break;

  case (kHNLProdNeuk3Muon):
    return "K0L --> N + mu+- + pi-+";
    break;

  case (kHNLProdNeuk3Electron):
    return "K0L --> N + e+- + pi-+";
    break;

  case (kHNLProdMuon3Numu):
    return "mu-+ --> N + numu(bar) + e-+";
    break;

  case (kHNLProdMuon3Nue):
    return "mu-+ --> N + nue(bar) + e-+";
    break;

  case (kHNLProdMuon3Nutau):
    return "mu-+ --> N + nutau(bar) + e-+";
    break;

  default:
    return "Invalid HNL production mode!";
  }
}
//____________________________________________________________________________
string genie::utils::hnl::AsString(HNLDecayMode_t hnldm)
{
  switch(hnldm) {

  case (kHNLDcyNull):
    return "Unspecified";
    break;

  case (kHNLDcyPiMu):
    return "N -> pi+- mu-+";
    break;

  case (kHNLDcyPiE):
    return "N -> pi+- e-+";
    break;

  case (kHNLDcyPi0Nu):
    return "N -> pi0 v";
    break;

  case (kHNLDcyNuNuNu):
    return "N -> v v v";
    break;

  case (kHNLDcyNuMuMu):
    return "N -> v mu+ mu-";
    break;

  case (kHNLDcyNuEE):
    return "N -> v e+ e-";
    break;

  case (kHNLDcyNuMuE):
    return "N -> v mu+- e-+";
    break;

  case (kHNLDcyPiPi0E):
    return "N -> pi+- pi0 e-+";
    break;

  case (kHNLDcyPiPi0Mu):
    return "N -> pi+- pi0 mu-+";
    break;

  case (kHNLDcyPi0Pi0Nu):
    return "N -> v pi0 pi0";
    break;

  case (kHNLDcyTEST):
    return "N -> v v";
    break;

  default:
    return "Invalid HNL decay mode!";
  }
  //return "Invalid HNL decay mode!";
}
//____________________________________________________________________________
bool genie::utils::hnl::IsProdKinematicallyAllowed(HNLProd_t hnlprod)
{
// Check the input mass of the HNL (M) against the sum of the masses of the
// particles to be produced in the input production mode

  PDGCodeList decay_products = ProductionProductList(hnlprod);

  PDGLibrary * pdglib = PDGLibrary::Instance();

  double Msum = 0.; double Mpar = 0.;
  PDGCodeList::const_iterator it = decay_products.begin();
  for ( ; it != decay_products.end(); ++it)
  {
    int pdg_code = *it;
    TParticlePDG * p = pdglib->Find(pdg_code);
    if(p) {
      if( it == decay_products.begin() ) Mpar = p->Mass();
      else Msum += p->Mass();
    } else {
       LOG("HNL", pERROR)
        << "Decay list includes particle with unrecognised PDG code: "
        << pdg_code;
    }
  }

  return (Mpar > Msum);
}
//____________________________________________________________________________
bool genie::utils::hnl::IsKinematicallyAllowed(HNLDecayMode_t hnldm, double M)
{
// Check the input mass of the HNL (M) against the sum of the masses of the
// particles to be produced in the input decay

  PDGCodeList decay_products = DecayProductList(hnldm);

  PDGLibrary * pdglib = PDGLibrary::Instance();

  double Msum = 0.;
  PDGCodeList::const_iterator it = decay_products.begin();
  for ( ; it != decay_products.end(); ++it)
  {
    int pdg_code = *it;
    TParticlePDG * p = pdglib->Find(pdg_code);
    if(p) {
       Msum += p->Mass();
    } else {
       LOG("HNL", pERROR)
        << "Decay list includes particle with unrecognised PDG code: "
        << pdg_code;
    }
  }

  return (M > Msum);
}
//____________________________________________________________________________
PDGCodeList genie::utils::hnl::ProductionProductList(HNLProd_t hnldm)
{
  // 0th element is parent
  // 1st is HNL
  // the rest are the HNL's siblings
  bool allow_duplicate = true;
  PDGCodeList decay_products(allow_duplicate);

  switch(hnldm) {
    
  case(kHNLProdPion2Muon):
    decay_products.push_back(kPdgPiP);
    decay_products.push_back(kPdgHNL);
    decay_products.push_back(kPdgMuon);
    break;
    
  case(kHNLProdPion2Electron):
    decay_products.push_back(kPdgPiP);
    decay_products.push_back(kPdgHNL);
    decay_products.push_back(kPdgElectron);
    break;
    
  case(kHNLProdKaon2Muon):
    decay_products.push_back(kPdgKP);
    decay_products.push_back(kPdgHNL);
    decay_products.push_back(kPdgMuon);
    break;

  case(kHNLProdKaon2Electron):
    decay_products.push_back(kPdgKP);
    decay_products.push_back(kPdgHNL);
    decay_products.push_back(kPdgElectron);
    break;

  case(kHNLProdKaon3Muon):
    decay_products.push_back(kPdgKP);
    decay_products.push_back(kPdgHNL);
    decay_products.push_back(kPdgMuon);
    decay_products.push_back(kPdgPi0);
    break;
    
  case(kHNLProdKaon3Electron):
    decay_products.push_back(kPdgKP);
    decay_products.push_back(kPdgHNL);
    decay_products.push_back(kPdgElectron);
    decay_products.push_back(kPdgPi0);
    break;

  case(kHNLProdNeuk3Muon):
    decay_products.push_back(kPdgK0L);
    decay_products.push_back(kPdgHNL);
    decay_products.push_back(kPdgMuon);
    decay_products.push_back(kPdgPiP);
    break;
    
  case(kHNLProdNeuk3Electron):
    decay_products.push_back(kPdgK0L);
    decay_products.push_back(kPdgHNL);
    decay_products.push_back(kPdgElectron);
    decay_products.push_back(kPdgPiP);
    break;

  case(kHNLProdMuon3Numu):
    decay_products.push_back(kPdgMuon);
    decay_products.push_back(kPdgHNL);
    decay_products.push_back(kPdgElectron);
    decay_products.push_back(kPdgAntiNuMu);
    break;

  case(kHNLProdMuon3Nue):
    decay_products.push_back(kPdgMuon);
    decay_products.push_back(kPdgHNL);
    decay_products.push_back(kPdgElectron);
    decay_products.push_back(kPdgAntiNuE);
    break;

  case(kHNLProdMuon3Nutau):
    decay_products.push_back(kPdgMuon);
    decay_products.push_back(kPdgHNL);
    decay_products.push_back(kPdgElectron);
    decay_products.push_back(kPdgAntiNuTau);
    break;

  default :
    break;
  }

  return decay_products;
}
//____________________________________________________________________________
PDGCodeList genie::utils::hnl::DecayProductList(HNLDecayMode_t hnldm)
{
  bool allow_duplicate = true;
  PDGCodeList decay_products(allow_duplicate);

  switch(hnldm) {

  case(kHNLDcyPiMu):
    decay_products.push_back(kPdgPiP);
    decay_products.push_back(kPdgMuon);
    break;

  case(kHNLDcyPiE):
    decay_products.push_back(kPdgPiP);
    decay_products.push_back(kPdgElectron);
    break;

  case(kHNLDcyPi0Nu):
    decay_products.push_back(kPdgPi0);
    decay_products.push_back(kPdgNuMu); // could be nue or nutau too!
    break;

  case(kHNLDcyNuNuNu):
    decay_products.push_back(kPdgNuE);
    decay_products.push_back(kPdgNuMu);
    decay_products.push_back(kPdgNuTau); // again, any permutation of {e,mu,tau}^3 works
    break;

  case(kHNLDcyNuMuMu):
    decay_products.push_back(kPdgNuMu);
    decay_products.push_back(kPdgMuon);
    decay_products.push_back(kPdgAntiMuon);
    break;

  case(kHNLDcyNuEE):
    decay_products.push_back(kPdgNuE);
    decay_products.push_back(kPdgElectron);
    decay_products.push_back(kPdgPositron);
    break;

  case(kHNLDcyNuMuE):
    decay_products.push_back(kPdgNuMu);
    decay_products.push_back(kPdgMuon);
    decay_products.push_back(kPdgPositron);
    break;

  case(kHNLDcyPiPi0E):
    decay_products.push_back(kPdgPiP);
    decay_products.push_back(kPdgPi0);
    decay_products.push_back(kPdgElectron);
    break;

  case(kHNLDcyPiPi0Mu):
    decay_products.push_back(kPdgPiP);
    decay_products.push_back(kPdgPi0);
    decay_products.push_back(kPdgMuon);
    break;

  case(kHNLDcyPi0Pi0Nu):
    decay_products.push_back(kPdgPi0);
    decay_products.push_back(kPdgPi0);
    decay_products.push_back(kPdgNuMu);
    break;

  case(kHNLDcyTEST):
    decay_products.push_back(kPdgNuMu);
    decay_products.push_back(kPdgAntiNuMu);
    break;

  default :
    break;
  }

  return decay_products;
}
//____________________________________________________________________________
// see e.g. Physics/HadronTensors/TabulatedLabFrameHadronTensor.cxx for this code
int genie::utils::hnl::GetCfgInt( string file_id, string set_name, string par_name )
{
  const genie::Registry * tmpReg = genie::AlgConfigPool::Instance()
    ->CommonList( file_id, set_name );
  int param = tmpReg->GetInt( par_name );
  return param;
}
//____________________________________________________________________________
// see e.g. Tools/Flux/GNuMIFlux.cxx for this code
std::vector<int> genie::utils::hnl::GetCfgIntVec( string file_id, string set_name, string par_name )
{
  // get vector as string first
  const genie::Registry * tmpReg = genie::AlgConfigPool::Instance()
    ->CommonList( file_id, set_name );
  string str = tmpReg->GetString( par_name );

  // turn string into vector<int>
  // be liberal about separators, users might punctuate for clarity
  std::vector<string> strtokens = genie::utils::str::Split( str, " ,;:()[]=" );
  std::vector<int> vect;
  size_t ntok = strtokens.size();

  for( size_t i=0; i < ntok; ++i ){
    std::string trimmed = utils::str::TrimSpaces( strtokens[i] );
    if( " " == trimmed || "" == trimmed ) continue; // skip empty strings
    int val = strtod( trimmed.c_str(), (char **)NULL );
    vect.push_back( val );
  }

  return vect;
}
//____________________________________________________________________________
double genie::utils::hnl::GetCfgDouble( string file_id, string set_name, string par_name )
{
  const genie::Registry * tmpReg = genie::AlgConfigPool::Instance()
    ->CommonList( file_id, set_name );
  double param = tmpReg->GetDouble( par_name );
  return param;
}
//____________________________________________________________________________
std::vector<double> genie::utils::hnl::GetCfgDoubleVec( string file_id, string set_name, string par_name )
{
  // get vector as string first
  const genie::Registry * tmpReg = genie::AlgConfigPool::Instance()
    ->CommonList( file_id, set_name );
  string str = tmpReg->GetString( par_name );

  // turn string into vector<double>
  // be liberal about separators, users might punctuate for clarity
  std::vector<string> strtokens = genie::utils::str::Split( str, " ,;:()[]=" );
  std::vector<double> vect;
  size_t ntok = strtokens.size();

  for( size_t i=0; i < ntok; ++i ){
    std::string trimmed = utils::str::TrimSpaces( strtokens[i] );
    if( " " == trimmed || "" == trimmed ) continue; // skip empty strings
    double val = strtod( trimmed.c_str(), (char **)NULL );
    vect.push_back( val );
  }

  return vect;
}
//____________________________________________________________________________
bool genie::utils::hnl::GetCfgBool( string file_id, string set_name, string par_name )
{
  const genie::Registry * tmpReg = genie::AlgConfigPool::Instance()
    ->CommonList( file_id, set_name );
  bool param = tmpReg->GetBool( par_name );
  return param;
}
//____________________________________________________________________________
std::vector<bool> genie::utils::hnl::GetCfgBoolVec( string file_id, string set_name, string par_name )
{
  // get vector as string first
  const genie::Registry * tmpReg = genie::AlgConfigPool::Instance()
    ->CommonList( file_id, set_name );
  string str = tmpReg->GetString( par_name );

  // turn string into vector<bool>
  // be liberal about separators, users might punctuate for clarity
  std::vector<string> strtokens = genie::utils::str::Split( str, " ,;:()[]=" );
  std::vector<bool> vect;
  size_t ntok = strtokens.size();

  for( size_t i=0; i < ntok; ++i ){
    std::string trimmed = utils::str::TrimSpaces( strtokens[i] );
    if( " " == trimmed || "" == trimmed ) continue; // skip empty strings
    bool val = strtod( trimmed.c_str(), (char **)NULL );
    vect.push_back( val );
  }

  return vect;
}
//____________________________________________________________________________
string genie::utils::hnl::GetCfgString( string file_id, string set_name, string par_name )
{
  const genie::Registry * tmpReg = genie::AlgConfigPool::Instance()
    ->CommonList( file_id, set_name );
  string param = tmpReg->GetString( par_name );
  return param;
}
//____________________________________________________________________________
