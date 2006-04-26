//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - November 25, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <sstream>

#include <TGenPhaseSpace.h>
#include <TMCParticle6.h>
#include <TF1.h>

#include "Algorithm/AlgFactory.h"
#include "Fragmentation/KNOHadronization2.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/MathUtils.h"
#include "Utils/PrintUtils.h"
#include "Utils/FragmRecUtils.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using std::ostringstream;

using namespace genie;
using namespace genie::utils::print;

const UInt_t kIsForward = 1<<17; // forwardness bit in TMCParticle bit field.

//____________________________________________________________________________
KNOHadronization2::KNOHadronization2() :
HadronizationModelI("genie::KNOHadronization2")
{

}
//____________________________________________________________________________
KNOHadronization2::KNOHadronization2(string config) :
HadronizationModelI("genie::KNOHadronization2", config)
{

}
//____________________________________________________________________________
KNOHadronization2::~KNOHadronization2()
{

}
//____________________________________________________________________________
void KNOHadronization2::Initialize(void) const
{

}
//____________________________________________________________________________
TClonesArray * KNOHadronization2::Hadronize(
                                        const Interaction * interaction) const
{
  LOG("KNOHad2", pDEBUG) << *fConfig;

  //-- Retrieve & run a version of the original KNO Hadronization model

  const HadronizationModelI * base_kno  = this->OriginalHadronizationModel();

  TClonesArray * particle_list = base_kno->Hadronize(interaction);

  int multiplicity = particle_list->GetEntries();

  LOG("KNOHad2", pINFO)
      << "\n Original KNO Hadronizer returned a list with "
                                             << multiplicity << " particles";

  assert(multiplicity > 1); // if everything was ok at previous step

  if(multiplicity == 2) return particle_list; // take no futhrer action for now

  utils::fragmrec::Print(particle_list);

  //-- Compute the 'forwardness' based on the expected Forard (xF>0) and
  //   Backward (xF<0) tot/+/-/neutral multiplicities

  bool ok = this->ComputeForwardness(particle_list, interaction);

  if(!ok) return particle_list;

  //-- Get the number of forward & backward going particles that we start with

  SLOG("KNOHad2", pINFO) << "NF = " << this->NParticles(particle_list, true )
                       << ", NB = " << this->NParticles(particle_list, false);

  //-- Compute the total fwd & bkw mass and set the total 4-momenta
  //   of the backward & foward system to decay
  //   Generate momentum kicks as a tunable gaussian fraction of the
  //   availabe momentum

  double W  = interaction->GetKinematics().W();

/*
  double PL = this->PKick(particle_list, interaction);

  if(PL<0) return particle_list;

  TLorentzVector p4_fwd(0,0, PL, W/2.);
  TLorentzVector p4_bkw(0,0,-PL, W/2.);

  LOG("KNOHad2", pINFO)
                << "\n 4P(fwd) = " << P4AsString(&p4_fwd)
                                    << "\n 4P(bkw) = " << P4AsString(&p4_bkw);

  //-- get foward & backward hadronic system

  TClonesArray * fwd_system = this->DecayExclusive(p4_fwd, particle_list, true );
  TClonesArray * bkw_system = this->DecayExclusive(p4_fwd, particle_list, false);

  //-- merge them into a single particle list

  TClonesArray * new_particle_list = new TClonesArray("TMCParticle", multiplicity);

  TMCParticle * particle = 0;

  TIter fwd_piter(fwd_system);
  TIter bkw_piter(bkw_system);

  int i = 0;

  while( (particle = (TMCParticle * ) fwd_piter.Next()) )
                    new ( (*new_particle_list)[i++] ) TMCParticle(*particle);

  while( (particle = (TMCParticle * ) bkw_piter.Next()) )
                    new ( (*new_particle_list)[i++] ) TMCParticle(*particle);

  utils::fragmrec::Print(new_particle_list);

  // pT...
  this->Cooling(new_particle_list);
*/

  this->Cooling2(particle_list, W);

//  utils::fragmrec::Print(new_particle_list);
  utils::fragmrec::Print(particle_list);

  // the cooling gets you out of the HCMS - boost back

//  this->BoostToHCMS(new_particle_list);
  this->BoostToHCMS(particle_list);

//  utils::fragmrec::Print(new_particle_list);
  utils::fragmrec::Print(particle_list);

  // delete the old particle list & the temporary bkw/fwd particle lists

//  particle_list -> Delete();
//  fwd_system    -> Delete();
//  bkw_system    -> Delete();

//  delete particle_list;
//  delete fwd_system;
//  delete bkw_system;

  // set the container as its element 'owner' & return

//  new_particle_list->SetOwner(true);

//  return new_particle_list;

  return particle_list;
}
//____________________________________________________________________________
const HadronizationModelI *
                     KNOHadronization2::OriginalHadronizationModel(void) const
{
// Retrieve a pre-configured version of the original KNO Hadronization model
// that this algorithm improves.

  assert(
      fConfig->Exists("base-kno-alg-name")  &&
      fConfig->Exists("base-kno-param-set")
  );

  string alg_name  = fConfig->GetString("base-kno-alg-name");
  string param_set = fConfig->GetString("base-kno-param-set");

  AlgFactory * algf = AlgFactory::Instance();

  const Algorithm * algbase = algf->GetAlgorithm(alg_name, param_set);

  const HadronizationModelI * base_kno =
                        dynamic_cast<const HadronizationModelI *> (algbase);

  assert(base_kno);

  return base_kno;
}
//____________________________________________________________________________
bool KNOHadronization2::ComputeForwardness(
          TClonesArray * particle_list, const Interaction * interaction) const
{
  //-- Get the expected Forward (xF>0) & Backward <xF<0), tot/+/-/neutral
  //   multipliciies according to the external (XML) configuration params
  //   of this hadronization model

  map<Multiplicity_t, double> mlt = this->ExpectedMultiplicities(interaction);

  double nF     = mlt[kMltFwdTotal];
  double nB     = mlt[kMltBkwTotal];
  double nF_neg = mlt[kMltFwdNegative];
  double nB_neg = mlt[kMltBkwNegative];
  double nF_pos = mlt[kMltFwdPositive];
  double nB_pos = mlt[kMltBkwPositive];
  double nF_0   = mlt[kMltFwdNeutral];
  double nB_0   = mlt[kMltBkwNeutral];

  //-- Get tot/+/-/neutral multiplicity in this event

  int multiplicity = particle_list->GetEntries();

  double n    = (double) multiplicity;
  double npos = (double) utils::fragmrec::NPositives(particle_list);
  double nneg = (double) utils::fragmrec::NNegatives(particle_list);
  double n0   = utils::math::NonNegative(n - npos - nneg);

  SLOG("KNOHad2", pINFO)
     << "KNO: N = " << n << ", N+ = " << npos << ", N- = " << nneg << ", N0 = " << n0;

  //-- Normalize bkw+fwd sums of the expected tot/+/-/neutral to the
  //   multiplicities of this event.
  //   Since normalized F/B multiplicities will be used to form
  //   forward-ness/backward-ness probabilities correct them for
  //   the leading correctly signed pion (put by hand in fwd hemisphere)
  //   and the leading nucleon (put by hand in fwd hemisphere);

  //----- specifically for vp

  if(n<2 || npos<2 || nF<1 || nB<1 || nF_pos<1 || nB_pos<1) return false;

  n    -= 2;
  npos -= 2;

  nF--;
  nB--;
  nF_pos--;
  nB_pos--;

  double norm_tot = (nF     + nB     > 0) ? n    / (nF     + nB    ) : 0.;
  double norm_pos = (nF_pos + nB_pos > 0) ? npos / (nF_pos + nB_pos) : 0.;
  double norm_neg = (nF_neg + nB_neg > 0) ? nneg / (nF_neg + nB_neg) : 0.;
  double norm_0   = (nF_0   + nB_0   > 0) ? n0   / (nF_0   + nB_0  ) : 0.;

  SLOG("KNOHad2", pINFO)
         << "NORMS: tot = " << norm_tot << ", + = "
                  << norm_pos << ", - = " << norm_neg << ", 0 = " << norm_0;

  nF     *= norm_tot;
  nB     *= norm_tot;
  nF_neg *= norm_neg;
  nB_neg *= norm_neg;
  nF_pos *= norm_pos;
  nB_pos *= norm_pos;
  nF_0   *= norm_0;
  nB_0   *= norm_0;

  SLOG("KNOHad2", pINFO)
     << "TOT: N = " << n << ", N+ = " << npos << ", N- = " << nneg << ", N0 = " << n0;
  SLOG("KNOHad2", pINFO)
     << "norm/FWD: N = " << nF << ", N+ = " << nF_pos << ", N- = " << nF_neg << ", N0 = " << nF_0;
  SLOG("KNOHad2", pINFO)
     << "norm/BKW: N = " << nB << ", N+ = " << nB_pos << ", N- = " << nB_neg << ", N0 = " << nB_0;


  //-- Form forward-ness "probabilities" for +/-/0's

  double fwd_probability_0   = (n0>0)     ? nF_0   /n0   : -1;
  double fwd_probability_pos = (npos > 0) ? nF_pos /npos : -1;
  double fwd_probability_neg = (nneg > 0) ? nF_neg /nneg : -1;

  SLOG("KNOHad2", pINFO) << "Probabilities:"
                         << " fwd[0] = " << fwd_probability_0
                         << " fwd[+] = " << fwd_probability_pos
                         << " fwd[-] = " << fwd_probability_neg;

  RandomGen * rnd = RandomGen::Instance();

  TIter piter(particle_list);

  TMCParticle * p = 0;

  double W  = interaction->GetKinematics().W();
  double WF = W/2;
  double WB = W/2;

  bool is_first = true;

  while( (p = (TMCParticle *) piter.Next()) ) {

      int    pdgc = p->GetKF();
      double mass = p->GetMass();

      double Q    = PDGLibrary::Instance()->Find(pdgc)->Charge();

      bool is_nucleon    = pdg::IsNeutronOrProton(pdgc);
      bool is_pion       = (pdgc==211);
      bool is_first_pion = false;

      if(is_pion && is_first) {

        is_first      = false;
        is_first_pion = true;
      }

      //---- assign forward-ness
      if(is_nucleon) {

         p->SetBit(kIsForward, false);
         WB -= mass;

      } else if (is_first_pion) {

         p->SetBit(kIsForward, true);
         WF -= mass;

      } else {
          double x = rnd->Random1().Rndm();
          double Pfwd = 0;
          if      ( Q < 0 ) Pfwd = fwd_probability_neg;
          else if ( Q > 0 ) Pfwd = fwd_probability_pos;
          else              Pfwd = fwd_probability_0;

          if (x < Pfwd) {
                  p->SetBit(kIsForward);
                  WF -= mass;
          }
          else {
                if(WB - mass > 0.010) {
                       p->SetBit(kIsForward, false);
                       WB -= mass;
                } else {
                       p->SetBit(kIsForward, true);
                       WF -= mass;
                }
          }
      } // not leading F/B hadrons
  } // p

  if(WF < 0 || WB < 0) return false;

  return true;
}
//____________________________________________________________________________
double KNOHadronization2::PKick(
          TClonesArray * particle_list, const Interaction * interaction) const
{
  double fwd_mass = this->TotalMass(particle_list, true );
  double bkw_mass = this->TotalMass(particle_list, false);

  double W = interaction->GetKinematics().W();

  double PAvail = W/2. - TMath::Max(fwd_mass, bkw_mass);

  RandomGen * rnd = RandomGen::Instance();

  double lambda = rnd->Random1().Gaus(0.7, 0.2);

  if( TMath::Abs(lambda) >= 1. ) lambda = 1;

  if(PAvail < 0.050) lambda = 0;

  double PL = lambda*PAvail;

  double M = 0;

  TMCParticle * p = 0;

  TIter piter(particle_list);

  if( this->NParticles(particle_list, true ) == 1 ) {

     while( (p = (TMCParticle *) piter.Next()) )
                          if( p->TestBit(kIsForward) ) M = p->GetMass();

  } else if ( this->NParticles(particle_list, false ) == 1 ) {

     while( (p = (TMCParticle *) piter.Next()) )
                        if( ! p->TestBit(kIsForward) ) M = p->GetMass();
  }

  if( M>0 ) {
     PL = TMath::Sqrt(W*W/4 - M*M);
  }

  if(PL > PAvail) return -1;

  LOG("KNOHad2", pINFO)
      << "\n W = " << W
      << ", M(fwd) = " << fwd_mass << ", M(bkw) = " << bkw_mass
      << "\n Pavail = " << PAvail << " / lambda = " << lambda;

  return PL;
}
//____________________________________________________________________________
int KNOHadronization2::NParticles(TClonesArray * part_list, bool fwd) const
{
  TIter piter(part_list);

  TMCParticle * p = 0;

  int N = 0;

  while( (p = (TMCParticle *) piter.Next()) ) {

      if( fwd &&  p->TestBit(kIsForward) ) N++;
      if(!fwd && !p->TestBit(kIsForward) ) N++;
  }
  return N;
}
//____________________________________________________________________________
double KNOHadronization2::TotalMass(TClonesArray * part_list, bool fwd) const
{
  TIter piter(part_list);

  TMCParticle * p = 0;

  double total_mass = 0;

  while( (p = (TMCParticle *) piter.Next()) ) {

     double m = p->GetMass();

     if( fwd &&  p->TestBit(kIsForward)) total_mass += m;
     if(!fwd && !p->TestBit(kIsForward)) total_mass += m;

  }
  return total_mass;
}
//____________________________________________________________________________
void KNOHadronization2::BoostToHCMS(TClonesArray * particle_list) const
{
   LOG("KNOHad2", pINFO) << "Boosting particles to HCMS";

   double sum_px = 0, sum_py = 0, sum_pz = 0, sum_E = 0;

   TMCParticle * particle = 0;

   TIter piter(particle_list);

   while( (particle = (TMCParticle * )piter.Next()) ) {

       sum_px += (particle->GetPx());
       sum_py += (particle->GetPy());
       sum_pz += (particle->GetPz());
       sum_E  += (particle->GetEnergy());
   }

   piter.Reset();

   LOG("KNOHad2", pINFO) << "Sum{Energy} = " << sum_E;

   assert(sum_E>0);

   TVector3 beta(-sum_px/sum_E, -sum_py/sum_E, -sum_pz/sum_E);

   while( (particle = (TMCParticle * )piter.Next()) ) {

       double px = (particle->GetPx());
       double py = (particle->GetPy());
       double pz = (particle->GetPz());
       double E  = (particle->GetEnergy());

       TLorentzVector p4( px, py, pz, E );

       p4.Boost(beta);

       particle->SetPx     ( p4.Px()     );
       particle->SetPy     ( p4.Py()     );
       particle->SetPz     ( p4.Pz()     );
       particle->SetEnergy ( p4.Energy() );
   }
}
//____________________________________________________________________________
TClonesArray * KNOHadronization2::DecayExclusive(
           TLorentzVector & sp4, TClonesArray * particle_list, bool fwd) const
{
  int N = this->NParticles(particle_list, fwd);

  double mass[N];
  int    pdgc[N];

  TIter piter(particle_list);

  TMCParticle * p = 0;

  int i = 0;

  while( (p = (TMCParticle *) piter.Next()) ) {

     double m  = p->GetMass();
     int    kf = p->GetKF();

     if( p->TestBit(kIsForward) == fwd) {

         mass[i] = m;
         pdgc[i] = kf;
         i++;
     }
  }

  TClonesArray * new_particle_list = new TClonesArray("TMCParticle", N);

  if(N > 1) {

    LOG("KNOHad2", pINFO) << "Setting the decay to the phase space generator";

    TGenPhaseSpace phgen;

    bool is_permitted = phgen.SetDecay(sp4, N, mass);

    assert(is_permitted);

    LOG("KNOHad2", pINFO) << "Generating phase space";

    phgen.Generate();

    for(int j = 0; j < N; j++) {

       LOG("KNOHad2", pINFO)
           << "Adding final state particle PDGC = " << pdgc[j]
                                    << " with mass = " << mass[j] << " GeV";

       TLorentzVector * p4 = phgen.GetDecay(j);

       int    code = pdgc[j];
       double M    = mass[j];
       double E    = p4->Energy();
       double px   = p4->Px();
       double py   = p4->Py();
       double pz   = p4->Pz();

       new ( (*new_particle_list)[j] )
                         TMCParticle(1,code,0,0,0,px,py,pz,E,M,0,0,0,0,0);

       if(fwd) (*new_particle_list)[j]->SetBit(kIsForward);
    }

  }
  else {

   LOG("KNOHad2", pINFO)
           << "Adding final state particle PDGC = " << pdgc[0]
                                    << " with mass = " << mass[0] << " GeV";
    int    code = pdgc[0];
    double px   = sp4.Px();
    double py   = sp4.Py();
    double pz   = sp4.Pz();
    double E    = sp4.Energy();
    double M    = mass[0];

    new ( (*new_particle_list)[0] )
                             TMCParticle(1,code,0,0,0,px,py,pz,E,M,0,0,0,0,0);

    if(fwd) (*new_particle_list)[0]->SetBit(kIsForward);
  }

  // set the container as its element 'owner' & return

  new_particle_list->SetOwner(true);

  return new_particle_list;
}
//____________________________________________________________________________
void KNOHadronization2::Cooling(TClonesArray * particle_list) const
{
  // pT2 is selected from an exponential distribution

  TF1 pt2func("pt2func","exp(-x/[0])",0.,2.);

  pt2func.SetParameter(0,0.65);

  RandomGen * rnd = RandomGen::Instance();

  TMCParticle * particle = 0;

  TIter piter(particle_list);

  double sum_px  = 0.;
  double sum_py  = 0.;

  double sgn_px  = -1. + 2. * rnd->Random1().Rndm();
  double sgn_py  = -1. + 2. * rnd->Random1().Rndm();

  int j = 0;
  int N = particle_list->GetEntries();

  while( (particle = (TMCParticle * ) piter.Next()) ) {

    LOG("KNOHad2", pINFO)
           << "Cooling particle PDGC = " << particle->GetKF();

    double px   = particle->GetPx();
    double py   = particle->GetPy();
    double pz   = particle->GetPz();

    LOG("KNOHad2", pINFO) << "Getting random pT2";

    double pT2  = pt2func.GetRandom();

    LOG("KNOHad2", pINFO) << "pT^2 = " << pT2;

    double pT2o = px*px + py*py;

    if( pT2 > pT2o ) pT2 = pT2o;

    double dpT2 = utils::math::NonNegative(pT2o - pT2);

    assert( dpT2 >= 0 );

    double dpz  = TMath::Sqrt( pz*pz + dpT2 ) - pz;

    double px2  = pT2 * rnd->Random1().Rndm();
    double py2  = pT2 - px2;

    px   = TMath::Sign( TMath::Sqrt(px2), -sum_px );
    py   = TMath::Sign( TMath::Sqrt(py2), -sum_py );
    pz   = TMath::Sign( pz+dpz, pz );

    if(j == N-1) {
       px = -sum_px;
       py = -sum_py;
    }

    particle->SetPx(px);
    particle->SetPy(py);
    particle->SetPz(pz);

    sum_px += px;
    sum_py += py;

    sgn_px  = -sum_px;
    sgn_py  = -sum_py;

    j++;
  }
}
//____________________________________________________________________________
void KNOHadronization2::Cooling2(TClonesArray * particle_list, double W) const
{
  // pT2 is selected from an exponential distribution

  TF1 pt2func("pt2func","exp(-x/[0])",0.,2.);

  pt2func.SetParameter(0,0.65);

  // xF is selected from a function fitted to experimental data

  TF1 fwd_xf_func("fwd_xf_func","0.10-0.09*x", 0.,1.);
  TF1 bkw_xf_func("bkw_xf_func","0.10+0.9*x", -1.,0.);

  RandomGen * rnd = RandomGen::Instance();

  TMCParticle * particle = 0;

  TIter piter(particle_list);

  double sum_px  = 0.;
  double sum_py  = 0.;
  double sum_pz  = 0.;

  double sgn_px  = -1. + 2. * rnd->Random1().Rndm();
  double sgn_py  = -1. + 2. * rnd->Random1().Rndm();

  int j = 0;
  int N = particle_list->GetEntries();

  while( (particle = (TMCParticle * ) piter.Next()) ) {

    LOG("KNOHad2", pINFO)
           << "Cooling particle PDGC = " << particle->GetKF();

    double px   = particle->GetPx();
    double py   = particle->GetPy();
    double pz   = particle->GetPz();

    double p2   = px*px + py*py + pz*pz;

    LOG("KNOHad2", pINFO) << "p2 = " << p2;

    LOG("KNOHad2", pINFO) << "Getting random xF";

    double xF = 0, pL = 0;

    do {
      if ( particle->TestBit(kIsForward)) xF = fwd_xf_func.GetRandom();
      else                                xF = bkw_xf_func.GetRandom();

      pL = xF*W/2;

      LOG("KNOHad2", pINFO) << "xF = " << xF;
      LOG("KNOHad2", pINFO) << "pL = " << pL;

    } while(pL*pL >= p2);


    double pT2  = p2 - pL*pL;

    LOG("KNOHad2", pINFO) << "pT^2 = " << pT2;

    assert( pT2 >= 0 );

    double px2  = pT2 * rnd->Random1().Rndm();
    double py2  = pT2 - px2;

    px   = TMath::Sign( TMath::Sqrt(px2), -sum_px );
    py   = TMath::Sign( TMath::Sqrt(py2), -sum_py );
    pz   = pL;

    if(j == N-1) {
    //   px = -sum_px;
    //   py = -sum_py;
    }

    particle->SetPx(px);
    particle->SetPy(py);
    particle->SetPz(pz);

    sum_px += px;
    sum_py += py;
    sum_pz += pz;

    sgn_px  = -sum_px;
    sgn_py  = -sum_py;

    j++;
  }
}
//____________________________________________________________________________
map<Multiplicity_t, double>
             KNOHadronization2::ExpectedMultiplicities(
                                        const Interaction * interaction) const
{
// Returns the expected Forward (xF>0) & Backward <xF<0), tot/+/-/neutral
// average multipliciies according to the external (XML) configuration params

  double W = interaction->GetKinematics().W();

  double lnW2 = TMath::Log(W*W);

  int v = interaction->GetInitialState().GetProbePDGCode();
  int N = interaction->GetInitialState().GetTarget().StruckNucleonPDGCode();

  double nF_tot = this->MultiplicityParam('A', kMltFwdTotal,    v,N) +
                  this->MultiplicityParam('B', kMltFwdTotal,    v,N) * lnW2;
  double nB_tot = this->MultiplicityParam('A', kMltBkwTotal,    v,N) +
                  this->MultiplicityParam('B', kMltBkwTotal,    v,N) * lnW2;
  double nF_neg = this->MultiplicityParam('A', kMltFwdNegative, v,N) +
                  this->MultiplicityParam('B', kMltFwdNegative, v,N) * lnW2;
  double nB_neg = this->MultiplicityParam('A', kMltBkwNegative, v,N) +
                  this->MultiplicityParam('B', kMltBkwNegative, v,N) * lnW2;
  double nF_pos = this->MultiplicityParam('A', kMltFwdPositive, v,N) +
                  this->MultiplicityParam('B', kMltFwdPositive, v,N) * lnW2;
  double nB_pos = this->MultiplicityParam('A', kMltBkwPositive, v,N) +
                  this->MultiplicityParam('B', kMltBkwPositive, v,N) * lnW2;

  nF_tot = utils::math::NonNegative( nF_tot );
  nB_tot = utils::math::NonNegative( nB_tot );
  nF_neg = utils::math::NonNegative( nF_neg );
  nB_neg = utils::math::NonNegative( nB_neg );
  nF_pos = utils::math::NonNegative( nF_pos );
  nB_pos = utils::math::NonNegative( nB_pos );

  double nF_0   = utils::math::NonNegative(nF_tot - nF_neg - nF_pos);
  double nB_0   = utils::math::NonNegative(nB_tot - nB_neg - nB_pos);

  nF_0  += TMath::Min(nF_neg, nF_pos);
  nB_0  += TMath::Min(nB_neg, nB_pos);

  LOG("KNOHad2", pINFO)
      << "\n init/FWD: N = " << nF_tot << ", N+ = "
                    << nF_pos << ", N- = " << nF_neg << ", N0 = " << nF_0;
  LOG("KNOHad2", pINFO)
      << "\n init/BKW: N = " << nB_tot << ", N+ = "
                    << nB_pos << ", N- = " << nB_neg << ", N0 = " << nB_0;

  map<Multiplicity_t, double> mlt;

  mlt[kMltFwdTotal]    = nF_tot;
  mlt[kMltBkwTotal]    = nB_tot;
  mlt[kMltFwdNegative] = nF_neg;
  mlt[kMltBkwNegative] = nB_neg;
  mlt[kMltFwdPositive] = nF_pos;
  mlt[kMltBkwPositive] = nB_pos;
  mlt[kMltFwdNeutral]  = nF_0;
  mlt[kMltBkwNeutral]  = nB_0;

  return mlt;
}
//____________________________________________________________________________
double KNOHadronization2::MultiplicityParam(
                          char param, Multiplicity_t mult, int v, int N) const
{
// Returns the parameter 'A' or 'B' (for computing the average multiplicity
// <N> = A + B * ln(W^2) for the input type of multiplicity, neutrino pdg code
// and nucleon pdg-code
// The parameter is retrieved from the algorithm's configuration registry that
// is created from its XML configuration file.
// The registry key for the requested parameter is constructed as:
// [vp/vbp/vn/vbn]-[A/B]-[fwd_tot/bkw-tot/fwd_neg/bkw-neg/fwd_pos/bkw-pos]

  ostringstream key;

  if      ( pdg::IsNeutrino    (v) && pdg::IsProton (N) ) key << "vp-";
  else if ( pdg::IsAntiNeutrino(v) && pdg::IsProton (N) ) key << "vbp-";
  else if ( pdg::IsNeutrino    (v) && pdg::IsNeutron(N) ) key << "vn-";
  else if ( pdg::IsAntiNeutrino(v) && pdg::IsNeutron(N) ) key << "vbn-";
  else
  {
   LOG("KNOHad2", pERROR)
        << "**** Unknown state with v-pdg = " << v << "and nuc-pdg = " << N;

   return 0;
  }

  if      (param == 'A') key << "A-";
  else if (param == 'B') key << "B-";
  else
  {
   LOG("KNOHad2", pERROR) << "**** Unknown parameter = " << param;
   return 0;
  }

  if      ( mult == kMltFwdTotal    ) key << "fwd-tot";
  else if ( mult == kMltBkwTotal    ) key << "bkw-tot";
  else if ( mult == kMltFwdNegative ) key << "fwd-neg";
  else if ( mult == kMltBkwNegative ) key << "bkw-neg";
  else if ( mult == kMltFwdPositive ) key << "fwd-pos";
  else if ( mult == kMltBkwPositive ) key << "bkw-pos";
  else
  {
   LOG("KNOHad2", pERROR) << "**** Unknown multiplicity code = " << mult;
   return 0;
  }

  assert( fConfig->Exists(key.str()) );

  double parameter = fConfig->GetDouble( key.str() );

  return parameter;
}
//____________________________________________________________________________
