//____________________________________________________________________________
/*!

\program gtestFGPauliBlockSuppr

\brief   Plot suppression factor due to Pauli-blocking.

\author  Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool

\created June 20, 2004

\cpright Copyright (c) 2003-2025, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         
*/
//____________________________________________________________________________

#include <TFile.h>
#include <TNtuple.h>

#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Physics/NuclearState/FermiMomentumTablePool.h"
#include "Physics/NuclearState/FermiMomentumTable.h"

using namespace genie;

int main(int /*argc*/, char ** /*argv*/)
{
  TNtuple * nt = new TNtuple("nt","","Z:A:nucl:Q2:Rdef:Rm10p:Rp10p");

  const int ntgt = 4;
  const int nnuc = 4;

  int target  [ntgt] = { kPdgTgtDeuterium, kPdgTgtC12, kPdgTgtO16, kPdgTgtFe56 };
  int nucleon [nnuc] = { kPdgProton, kPdgNeutron };

  const int    N        = 3001;
  const double Q2min    = 0.000001*units::GeV2;
  const double Q2max    = 10*units::GeV2;
  const double logQ2min = TMath::Log(Q2min);
  const double logQ2max = TMath::Log(Q2max);
  const double dlogQ2   = (logQ2max-logQ2min)/(N-1);
  const double kMN      = constants::kNucleonMass;
  const double pmax     = 0.5*units::GeV;

  FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
  const FermiMomentumTable * kft  = kftp->GetTable("Default");

  for(int itgt=0; itgt<ntgt; itgt++) {
   for(int inuc=0; inuc<nnuc; inuc++) {

     double kFi = kft->FindClosestKF(target[itgt], nucleon[inuc]);

     int Z = pdg::IonPdgCodeToZ(target[itgt]);
     int A = pdg::IonPdgCodeToA(target[itgt]);

     LOG("test", pNOTICE) 
       << "nuclear target = " << target[itgt] << ", nucleon = " << nucleon[inuc];

     for(int i=0; i<N; i++) {
       double Q2 = TMath::Exp(logQ2min + i*dlogQ2); 
       double Rdef  = utils::nuclear::RQEFG_generic(-1*Q2,kMN,    kFi,    kFi,pmax);
       double Rm10p = utils::nuclear::RQEFG_generic(-1*Q2,kMN,0.9*kFi,0.9*kFi,pmax); // kF -> +10%
       double Rp10p = utils::nuclear::RQEFG_generic(-1*Q2,kMN,1.1*kFi,1.1*kFi,pmax); // kF -> -10%

       LOG("test", pNOTICE) 
         << "Q2 = " << Q2 << " GeV, Rdef = " << Rdef 
         << ", R(kF->0.9kF) = " << Rm10p << ", R(kF->1.1*kF) = " << Rp10p;

       nt->Fill(Z,A,nucleon[inuc],Q2,Rdef,Rm10p,Rp10p);

    }//i
   }//inuc
  }//itgt

  TFile f("./fg_pauli_suppression_factors.root","recreate");
  nt->Write();
  f.Close();

  return 0;
}
