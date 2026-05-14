//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool

         Changes required to implement the GENIE Boosted Dark Matter module
         were installed by Josh Berger (Univ. of Wisconsin)

         Changes required to implement the GENIE Dark Neutrino module
         were installed by Iker de Icaza (Univ. of Sussex)
	 
	 Changes required to implement the GENIE BeamHNL module
	 were installed by John Plows (Univ. of Oxford)
*/
//____________________________________________________________________________

#include <iostream>
#include <string>

#include <TSystem.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Conventions/Constants.h"


using std::string;

using namespace genie;

//____________________________________________________________________________
PDGLibrary * PDGLibrary::fInstance = 0;
//____________________________________________________________________________
PDGLibrary::PDGLibrary()
{
  if( ! LoadDBase() ) LOG("PDG", pERROR) << "Could not load PDG data";

#ifdef __GENIE_DARK_NEUTRINO_ENABLED__
  LOG("PDG", pINFO) << "Loading Dark sector Info";
  if ( ! AddDarkSector() ) { 
    LOG("PDG", pFATAL) << "Could not load Dark Neutrino data";
    exit(78);
  }
#endif // __GENIE_DARK_NEUTRINO_ENABLED__

#ifdef __GENIE_HEAVY_NEUTRAL_LEPTON_ENABLED__
  LOG("PDG", pINFO) << "Loading Heavy Neutral Lepton data";
  if( ! AddHNL() ){
    LOG("PDG", pFATAL) << "Could not load Heavy Neutral Lepton data";
    exit(78);
  }
#endif // #ifdef __GENIE_HEAVY_NEUTRAL_LEPTON_ENABLED__
  
  fInstance =  0;
}
//____________________________________________________________________________
PDGLibrary::~PDGLibrary()
{
  fInstance = 0;
}
//____________________________________________________________________________
PDGLibrary * PDGLibrary::Instance()
{
  if(fInstance == 0) {
    LOG("PDG", pINFO) << "PDGLibrary late initialization";

    static PDGLibrary::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();

    fInstance = new PDGLibrary;
  }
  return fInstance;
}
//____________________________________________________________________________
TDatabasePDG * PDGLibrary::DBase(void)
{
  return fDatabasePDG;
}
//____________________________________________________________________________
TParticlePDG * PDGLibrary::Find(int pdgc, bool must_exist )
{

  TParticlePDG * ret = fDatabasePDG->GetParticle(pdgc);
  if(ret) return ret;

  if ( must_exist ) {
    LOG("PDG", pERROR) << "Requested missing particle with PDG: " << pdgc ;
  }

  return ret ;
}

//____________________________________________________________________________
bool PDGLibrary::LoadDBase(void)
{
  fDatabasePDG = TDatabasePDG::Instance();

  // loading PDG data from $GENIE/config/
  const char* altpdgtable = gSystem->Getenv("GENIE_PDG_TABLE");
  if ( altpdgtable ) {
    if ( ! (gSystem->AccessPathName(altpdgtable) ) ) {
        LOG("PDG", pINFO) << "Load PDG data from $GENIE_PDG_TABLE: "
                          << altpdgtable;
        fDatabasePDG->ReadPDGTable( altpdgtable );
        return true;
    }
  }

  if ( gSystem->Getenv("GENIE") ) {
    string base_dir = string( gSystem->Getenv("GENIE") );
    base_dir += string("/data/evgen/catalogues/pdg/") ; 

    string file_name = "genie_pdg_table.txt" ; 
    const Registry * reg = AlgConfigPool::Instance()->CommonList("Param", "PDG");
    if( reg ) {
      file_name = reg -> GetString("PDG-TableName") ;
      LOG("PDG", pINFO) << "Found file name specification: " << file_name ;

    }
    
    string path = base_dir + file_name ;

    if ( ! (gSystem->AccessPathName(path.c_str()) ) ) {
        LOG("PDG", pINFO) << "Load PDG data from: " << path;
        fDatabasePDG->ReadPDGTable( path.c_str() );
        return true;
    }
  }

  // no PDG data in $GENIE/config/ - Try $ROOTSYS/etc/

  if(gSystem->Getenv("ROOTSYS")) {
    string base_dir  = string( gSystem->Getenv("ROOTSYS") );
    string path = base_dir  + string("/etc/pdg_table.txt");

    if ( !(gSystem->AccessPathName(path.c_str())) ) {
        LOG("PDG", pINFO) << "Load PDG data from: " << path;
        fDatabasePDG->ReadPDGTable( path.c_str() );
        return true;
     }
  }

  LOG("PDG", pERROR) << " *** The PDG extensions will not be loaded!! ***";
  return false;
};
//____________________________________________________________________________
void PDGLibrary::AddDarkMatter(double mass, double med_ratio)
{
// Add dark matter particle to PDG database

  double med_mass = mass*med_ratio;
  TParticlePDG * dm_particle = fDatabasePDG->GetParticle(kPdgDarkMatter);
  TParticlePDG * med_particle = fDatabasePDG->GetParticle(kPdgMediator);
  if (!dm_particle) {
    // Name Title Mass Stable Width Charge Class PDG
    fDatabasePDG->AddParticle("chi_dm","chi_dm",mass,true,0.,0,"DarkMatter",kPdgDarkMatter);
  }
  else {
    assert(dm_particle->Mass() == mass);
  }
  if (!med_particle) {
    // Name Title Mass Stable Width Charge Class PDG
    fDatabasePDG->AddParticle("Z_prime","Z_prime",med_mass,true,0.,0,"DarkMatter",kPdgMediator);
  }
  else {
    assert(med_particle->Mass() == med_mass);
  }
}
//____________________________________________________________________________
bool PDGLibrary::AddHNL()
{
  // Add HNL to PDG database
  const Registry * reg = AlgConfigPool::Instance()->CommonList("HNL", "ParameterSpace");
  if (!reg) {
    LOG("PDG", pERROR) << "Cannot find HNL ParameterSpace param_set";
    return false;
  }
  TParticlePDG * hnl = fDatabasePDG->GetParticle(kPdgHNL);
  if (!hnl) {
    // Name Title Mass Stable Width Charge Class PDG
    fDatabasePDG->AddParticle("HNL","HNL",reg->GetDouble("HNL-Mass"),true,0.,0,"HNL",kPdgHNL);
    fDatabasePDG->AddParticle("HNLBar","HNLBar",reg->GetDouble("HNL-Mass"),true,0.,0,"HNL",-1*kPdgHNL);
  }
  return true;
}
//____________________________________________________________________________
bool PDGLibrary::AddDarkSector()
{
  // Add dark neutrino particles to PDG database

  const Registry * reg = AlgConfigPool::Instance()->CommonList("Dark", "Masses");
  if(!reg) {
    LOG("PDG", pERROR) << "The Dark Sector masses not available.";
    return false;
  }
  TParticlePDG * dnu_particle = fDatabasePDG->GetParticle(kPdgDarkNeutrino);
  TParticlePDG * anti_dnu_particle = fDatabasePDG->GetParticle(kPdgAntiDarkNeutrino);
  TParticlePDG * med_particle = fDatabasePDG->GetParticle(kPdgDNuMediator);
  if (!dnu_particle) {
    // Name Title Mass Stable Width Charge Class PDG
    fDatabasePDG->AddParticle("nu_D","#nu_{D}",reg->GetDouble("Dark-NeutrinoMass"),
                              true,0.,0,"DarkNeutrino",kPdgDarkNeutrino);
  }
  if (!anti_dnu_particle) {
    // Name Title Mass Stable Width Charge Class PDG
    fDatabasePDG->AddParticle("nu_D_bar","#bar{#nu}_{D}",reg->GetDouble("Dark-NeutrinoMass"),
                              true,0.,0,"DarkNeutrino",kPdgAntiDarkNeutrino);
  }
  if (!med_particle) {
    // Name Title Mass Stable Width Charge Class PDG
    fDatabasePDG->AddParticle("Z_D","Z_{D}",reg->GetDouble("Dark-MediatorMass"),
                              true,0.,0,"DarkNeutrino",kPdgDNuMediator);
  }
  return true;
}
//____________________________________________________________________________
#ifdef __GENIE_INCL_ENABLED__
void PDGLibrary::AddHypernucleus(int pdg_hypernucleus){
  if (fDatabasePDG->GetParticle(pdg_hypernucleus)) {
    LOG("PDG", pINFO) << "Hyper-nucleus (" << pdg_hypernucleus
                      << ") already in PDG library";
    return;
  }
  int abs_pdg = std::abs(pdg_hypernucleus);
  int S = (abs_pdg / 10000000) % 10;
  if (S != 1) {
    LOG("PDG", pFATAL)
      << "Hypernucleus with strangeness S=" << S
      << " (PDG=" << pdg_hypernucleus << ") not supported. Only S=1 is handled.";
    exit(1);
  }
  int sign = (pdg_hypernucleus < 0) ? -1 : 1;     // sign for particle and anti-particle
  int pdg_ion = sign * (abs_pdg - S * 10000000);  // strip strangeness digit(s)

  TParticlePDG *ion = fDatabasePDG->GetParticle(pdg_ion);
  if (!ion) {
    LOG("PDG", pFATAL) << "Base ion " << pdg_ion
                       << " not found for hypernucleus " << pdg_hypernucleus;
    exit(1);
  }

  // m(hypernucleus) ≈ m(ion) + S*(m_Λ - m_N); ignores Λ binding energy.
  const double dm_lambda_n = 0.17608; // GeV, m_Λ - m_n
  double hypernucl_mass = ion->Mass() + S * dm_lambda_n;

  std::string name = std::string(ion->GetName()) + "_L" + std::to_string(S);
  LOG("PDG", pINFO) << "Adding hyper-nucleus: " << name << " => " << pdg_hypernucleus;
  fDatabasePDG->AddParticle(name.c_str(), name.c_str(), hypernucl_mass,
                            true, 0, ion->Charge(), "HyperIon", pdg_hypernucleus);
}
void PDGLibrary::AddVirtualCluster(int pdg_virtual){
  if (fDatabasePDG->GetParticle(pdg_virtual)) {
    LOG("PDG", pINFO) << "Virtual cluster (" << pdg_virtual
                      << ") already in PDG library";
    return;
  }
  int A = (pdg_virtual/10)%1000;
  int Z = (pdg_virtual/10000)%1000;
  double virtual_cluster_mass = A*genie::constants::kNucleonMass;
  int charge = Z;
  std::string name = "VCluster_Z" + std::to_string(Z) + "_A" + std::to_string(A);
  LOG("PDG", pINFO) << "Adding virtual cluster: " << name << " => " << pdg_virtual;
  fDatabasePDG->AddParticle(name.c_str(), name.c_str(), virtual_cluster_mass,
                            true, 0, charge, "VirtualCluster", pdg_virtual);
}
#else
// you don't need this without INCL
void PDGLibrary::AddHypernucleus(int pdg) {
  LOG("PDG", pERROR) << "INCL support not compiled in";
}
void PDGLibrary::AddVirtualCluster(int pdg) {
  LOG("PDG", pERROR) << "INCL support not compiled in";
}
#endif
//____________________________________________________________________________
// EDIT: need a way to clear and then reload the PDG database
void PDGLibrary::ReloadDBase(void)
{
  if(fDatabasePDG) {
    delete fDatabasePDG;
  }

  if( ! LoadDBase() ) LOG("PDG", pERROR) << "Could not load PDG data";
}
//____________________________________________________________________________
