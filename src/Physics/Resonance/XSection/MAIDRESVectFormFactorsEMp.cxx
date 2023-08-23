//____________________________________________________________________________
/*
  Copyright (c) 2023-2025, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org

  Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
  University of Liverpool & STFC Rutherford Appleton Laboratory
*/
//____________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/Resonance/XSection/MAIDRESVectFormFactorsEMp.h"
#include "Framework/Utils/StringUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
MAIDRESVectFormFactorsEMp::MAIDRESVectFormFactorsEMp() :
  RESVectFormFactorsI("genie::MAIDRESVectFormFactorsEMp")
{

}
//____________________________________________________________________________
MAIDRESVectFormFactorsEMp::MAIDRESVectFormFactorsEMp(string config) :
  RESVectFormFactorsI("genie::MAIDRESVectFormFactorsEMp", config)
{

}
//____________________________________________________________________________
MAIDRESVectFormFactorsEMp::~MAIDRESVectFormFactorsEMp()
{

}
//____________________________________________________________________________
RESVectFFAmplitude MAIDRESVectFormFactorsEMp::Compute( const Interaction interaction ) {
  RESVectFFAmplitude ampl ; 

  const InitialState & init_state = interaction.InitState();
  const Target & target = init_state.Tgt();
  const Resonance_t res = interaction.ExclTag().Resonance();
  
  //Get kinematical parameters
  const Kinematics & kinematics = interaction.Kine();
  double q2 = kinematics.q2();
  double q4 = TMath::Power( q2, 2 ) ; 
  double q6 = TMath::Power( q2, 3 ) ; 
  double q8 = TMath::Power( q2, 4 ) ; 

  double W2     = TMath::Power( kinematics.W(), 2);
  double Mnuc2  = TMath::Power(target.HitNucMass(), 2);
  double MR  = utils::res::Mass(res);
  double kR = (MR*MR - Mnuc2)/(2.*MR); 

  if( res == kP33_1232 ) { 
    // Dipole form facor
    double GD = 1./TMath::Power(1-q2/fDipoleMass,2) ;
    // Transition form factors
    double GM = fAM0_P33_1232 * (1. - fBetaM_P33_1232 * q2) * TMath::Exp( fGammaM_P33_1232 * q2 ) * GD;
    double GE = fAE0_P33_1232 * (1. - fBetaE_P33_1232 * q2) * TMath::Exp( fGammaE_P33_1232 * q2 ) * GD;
    double GC = fAC0_P33_1232 * (1. - fBetaC_P33_1232 * q2) / (1. - fDC_P33_1232 * q2/(4.*Mnuc2))*2*MR/kR*TMath::Exp( fGammaC_P33_1232 * q2 ) * GD; 
    ampl.SetAmplA12( (3.*GE+GM)/2. );
    ampl.SetAmplA32( TMath::Sqrt(3.)/2.*(GE-GM) );
    ampl.SetAmplS12( TMath::Sqrt(2.)*GC );
  } else { 
    double A120 = fA120P[res] ;
    double A12Alpha1 = fA12Alpha1P[res] ;
    double A12Alpha2 = fA12Alpha2P[res] ;
    double A12Alpha4 = fA12Alpha4P[res] ;
    double A12Beta = fA12BetaP[res] ;
    double A320 = fA320P[res] ;
    double A32Alpha1 = fA32Alpha1P[res] ;
    double A32Alpha2 = fA32Alpha2P[res] ;
    double A32Alpha4 = fA32Alpha4P[res] ;
    double A32Beta = fA32BetaP[res] ;
    double S120 = fS120P[res] ;
    double S12Alpha1 = fS12Alpha1P[res] ;
    double S12Alpha2 = fS12Alpha2P[res] ;
    double S12Alpha4 = fS12Alpha4P[res] ;
    double S12Beta = fS12BetaP[res] ;

    ampl.SetAmplA12( A120 * ( 1 - A12Alpha1 * q2 + A12Alpha2 * q4 + A12Alpha4 * q8 ) * TMath::Exp( A12Beta * q2 ) ) ; 
    ampl.SetAmplA32( A320 * ( 1 - A32Alpha1 * q2 + A32Alpha2 * q4 + A32Alpha4 * q8 ) * TMath::Exp( A32Beta * q2 ) ) ; 
    ampl.SetAmplS12( S120 * ( 1 - S12Alpha1 * q2 + S12Alpha2 * q4 + S12Alpha4 * q8 ) * TMath::Exp( S12Beta * q2 ) ) ; 
  }

  return ampl;
}
//____________________________________________________________________________
void MAIDRESVectFormFactorsEMp::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void MAIDRESVectFormFactorsEMp::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
void MAIDRESVectFormFactorsEMp::LoadConfig(void)
{
  bool good_config = true ; 

  this->GetParam( "AM@P33(1232)", fAM0_P33_1232 ) ; 
  this->GetParam( "AE@P33(1232)", fAE0_P33_1232 ) ;
  this->GetParam( "AC@P33(1232)", fAC0_P33_1232 ) ; 
  
  this->GetParam( "BetaM@P33(1232)", fBetaM_P33_1232 ) ;
  this->GetParam( "BetaE@P33(1232)", fBetaE_P33_1232 ) ;
  this->GetParam( "BetaC@P33(1232)", fBetaC_P33_1232 ) ;

  this->GetParam( "GammaM@P33(1232)", fGammaM_P33_1232 ) ;
  this->GetParam( "GammaE@P33(1232)", fGammaE_P33_1232 ) ;
  this->GetParam( "GammaC@P33(1232)", fGammaC_P33_1232 ) ;

  this->GetParam( "DeltaC@P33(1232)", fDC_P33_1232 ) ;
  this->GetParam( "DipoleMass", fDipoleMass ) ;
 
  auto kres_list_A12_0_p = GetConfig().FindKeys("A120P@") ;
  if( kres_list_A12_0_p.size() == 0 ) good_config = false ; 
  for( auto kiter = kres_list_A12_0_p.begin(); kiter != kres_list_A12_0_p.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"@");
    assert(kv.size()==2);
    Resonance_t res_id = utils::res::FromString( (kv[1]).c_str() );
    this->GetParam( key, fA120P[res_id] ) ; 
  }

  auto kres_list_A12_alpha1_p = GetConfig().FindKeys("A12Alpha1P@") ;
  if( kres_list_A12_alpha1_p.size() == 0 ) good_config = false ; 
  for( auto kiter = kres_list_A12_alpha1_p.begin(); kiter != kres_list_A12_alpha1_p.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"@");
    assert(kv.size()==2);
    Resonance_t res_id = utils::res::FromString( (kv[1]).c_str() );
    this->GetParam( key, fA12Alpha1P[res_id] ) ; 
  }

  auto kres_list_A12_alpha2_p = GetConfig().FindKeys("A12Alpha2P@") ;
  if( kres_list_A12_alpha2_p.size() == 0 ) good_config = false ; 
  for( auto kiter = kres_list_A12_alpha2_p.begin(); kiter != kres_list_A12_alpha2_p.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"@");
    assert(kv.size()==2);
    Resonance_t res_id = utils::res::FromString( (kv[1]).c_str() );
    this->GetParam( key, fA12Alpha2P[res_id] ) ; 
  }

  auto kres_list_A12_alpha4_p = GetConfig().FindKeys("A12Alpha4P@") ;
  if( kres_list_A12_alpha4_p.size() == 0 ) good_config = false ; 
  for( auto kiter = kres_list_A12_alpha4_p.begin(); kiter != kres_list_A12_alpha4_p.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"@");
    assert(kv.size()==2);
    Resonance_t res_id = utils::res::FromString( (kv[1]).c_str() );
    this->GetParam( key, fA12Alpha4P[res_id] ) ; 
  }

  auto kres_list_A12_beta_p = GetConfig().FindKeys("A12BetaP@") ;
  if( kres_list_A12_beta_p.size() == 0 ) good_config = false ; 
  for( auto kiter = kres_list_A12_beta_p.begin(); kiter != kres_list_A12_beta_p.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"@");
    assert(kv.size()==2);
    Resonance_t res_id = utils::res::FromString( (kv[1]).c_str() );
    this->GetParam( key, fA12BetaP[res_id] ) ; 
  }

  auto kres_list_A32_0_p = GetConfig().FindKeys("A320P@") ;
  if( kres_list_A32_0_p.size() == 0 ) good_config = false ; 
  for( auto kiter = kres_list_A32_0_p.begin(); kiter != kres_list_A32_0_p.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"@");
    assert(kv.size()==2);
    Resonance_t res_id = utils::res::FromString( (kv[1]).c_str() );
    this->GetParam( key, fA320P[res_id] ) ; 
  }

  auto kres_list_A32_alpha1_p = GetConfig().FindKeys("A32Alpha1P@") ;
  if( kres_list_A32_alpha1_p.size() == 0 ) good_config = false ; 
  for( auto kiter = kres_list_A32_alpha1_p.begin(); kiter != kres_list_A32_alpha1_p.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"@");
    assert(kv.size()==2);
    Resonance_t res_id = utils::res::FromString( (kv[1]).c_str() );
    this->GetParam( key, fA32Alpha1P[res_id] ) ; 
  }

  auto kres_list_A32_alpha2_p = GetConfig().FindKeys("A32Alpha2P@") ;
  if( kres_list_A32_alpha2_p.size() == 0 ) good_config = false ; 
  for( auto kiter = kres_list_A32_alpha2_p.begin(); kiter != kres_list_A32_alpha2_p.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"@");
    assert(kv.size()==2);
    Resonance_t res_id = utils::res::FromString( (kv[1]).c_str() );
    this->GetParam( key, fA32Alpha2P[res_id] ) ; 
  }

  auto kres_list_A32_alpha4_p = GetConfig().FindKeys("A32Alpha4P@") ;
  if( kres_list_A32_alpha4_p.size() == 0 ) good_config = false ; 
  for( auto kiter = kres_list_A32_alpha4_p.begin(); kiter != kres_list_A32_alpha4_p.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"@");
    assert(kv.size()==2);
    Resonance_t res_id = utils::res::FromString( (kv[1]).c_str() );
    this->GetParam( key, fA32Alpha4P[res_id] ) ; 
  }

  auto kres_list_A32_beta_p = GetConfig().FindKeys("A32BetaP@") ;
  if( kres_list_A32_beta_p.size() == 0 ) good_config = false ; 
  for( auto kiter = kres_list_A32_beta_p.begin(); kiter != kres_list_A32_beta_p.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"@");
    assert(kv.size()==2);
    Resonance_t res_id = utils::res::FromString( (kv[1]).c_str() );
    this->GetParam( key, fA32BetaP[res_id] ) ; 
  }


  auto kres_list_S12_0_p = GetConfig().FindKeys("S120P@") ;
  if( kres_list_S12_0_p.size() == 0 ) good_config = false ; 
  for( auto kiter = kres_list_S12_0_p.begin(); kiter != kres_list_S12_0_p.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"@");
    assert(kv.size()==2);
    Resonance_t res_id = utils::res::FromString( (kv[1]).c_str() );
    this->GetParam( key, fS120P[res_id] ) ; 
  }

  auto kres_list_S12_alpha1_p = GetConfig().FindKeys("S12Alpha1P@") ;
  if( kres_list_S12_alpha1_p.size() == 0 ) good_config = false ; 
  for( auto kiter = kres_list_S12_alpha1_p.begin(); kiter != kres_list_S12_alpha1_p.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"@");
    assert(kv.size()==2);
    Resonance_t res_id = utils::res::FromString( (kv[1]).c_str() );
    this->GetParam( key, fS12Alpha1P[res_id] ) ; 
  }

  auto kres_list_S12_alpha2_p = GetConfig().FindKeys("S12Alpha2P@") ;
  if( kres_list_S12_alpha2_p.size() == 0 ) good_config = false ; 
  for( auto kiter = kres_list_S12_alpha2_p.begin(); kiter != kres_list_S12_alpha2_p.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"@");
    assert(kv.size()==2);
    Resonance_t res_id = utils::res::FromString( (kv[1]).c_str() );
    this->GetParam( key, fS12Alpha2P[res_id] ) ; 
  }

  auto kres_list_S12_alpha4_p = GetConfig().FindKeys("S12Alpha4P@") ;
  if( kres_list_S12_alpha4_p.size() == 0 ) good_config = false ; 
  for( auto kiter = kres_list_S12_alpha4_p.begin(); kiter != kres_list_S12_alpha4_p.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"@");
    assert(kv.size()==2);
    Resonance_t res_id = utils::res::FromString( (kv[1]).c_str() );
    this->GetParam( key, fS12Alpha4P[res_id] ) ; 
  }

  auto kres_list_S12_beta_p = GetConfig().FindKeys("S12BetaP@") ;
  if( kres_list_S12_beta_p.size() == 0 ) good_config = false ; 
  for( auto kiter = kres_list_S12_beta_p.begin(); kiter != kres_list_S12_beta_p.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"@");
    assert(kv.size()==2);
    Resonance_t res_id = utils::res::FromString( (kv[1]).c_str() );
    this->GetParam( key, fS12BetaP[res_id] ) ; 
  }

  if( ! good_config ) { 
    LOG("MAIDRESVectFormFactorsEMp", pERROR ) << " Configuration failed.";
    exit(78) ;
  }

}
