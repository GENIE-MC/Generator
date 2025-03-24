//____________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Steven Gardiner <gardiner \at fnal.gov>
 Fermi National Accelerator Laboratory
*/
//____________________________________________________________________________

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include "TMatrixD.h"
 #include "TDecompSVD.h"

#include "Framework/Conventions/Constants.h"
#include "Physics/Common/PrimaryLeptonUtils.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/RunningThreadInfo.h"
#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/ParticleData/PDGLibrary.h"

using namespace genie;
using namespace genie::utils;
using namespace genie::constants;
using namespace std::complex_literals;

//___________________________________________________________________________
void genie::utils::SetPrimaryLeptonPolarization( GHepRecord * ev )
{
// Moved out of the PrimaryLeptonGenerator class to make the same treatment
// accessible for generators that use a more unified approach (e.g.,
// QELEventGenerator and MECGenerator). -- S. Gardiner

  // get the final state primary lepton
  GHepParticle * fsl = ev->FinalStatePrimaryLepton();
  if ( !fsl ) {
    LOG("LeptonicVertex", pERROR)
      << "Final state lepton not set yet! \n" << *ev;
    return;
  }
  //-- Get the interaction
  Interaction * interaction = ev->Summary();
  //-- Access cross section algorithm for running thread
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  const XSecAlgorithmI * xsec_alg = evg->CrossSectionAlg();
  fsl->SetPolarization(xsec_alg->FinalLeptonPolarization(interaction));
    
  LOG("LeptonicVertex", pINFO)
    << "Setting polarization angles for particle: " << fsl->Name();

  if ( fsl->PolzIsSet() ) {
    LOG("LeptonicVertex", pINFO)
      << "Polarization (rad): Polar = "  << fsl->PolzPolarAngle()
      << ", Azimuthal = " << fsl->PolzAzimuthAngle();
  }
  
}
//___________________________________________________________________________
void genie::utils::CalculatePolarizationVectorWithNuclearTensor(
                        TVector3 & polarization,
                        const TLorentzVector & neutrinoMom,
                        const TLorentzVector & leptonMom, 
                        bool isLeftPolarized,
                        const HermitianMatrix & NTensor
)
{
  double k[4], l[4], s[4], eskl[4];
  std::complex<double> jp[4], jm[4];
    
  k[0] = neutrinoMom.E();
  k[1] = -neutrinoMom.Px();
  k[2] = -neutrinoMom.Py();
  k[3] = -neutrinoMom.Pz();
  
  l[0] = leptonMom.E();
  l[1] = -leptonMom.Px();
  l[2] = -leptonMom.Py();
  l[3] = -leptonMom.Pz();
  
  double ml = leptonMom.M();
  
  s[0] = leptonMom.P()/ml;
  s[1] = -leptonMom.Vect().Unit().X()*leptonMom.E()/ml;
  s[2] = -leptonMom.Vect().Unit().Y()*leptonMom.E()/ml;
  s[3] = -leptonMom.Vect().Unit().Z()*leptonMom.E()/ml;
  

  
  // epsilon_\alpha\beta\gamma\delta s^\beta k^\gamma l^\delta
  for (int a = 0; a < 4; a++)
  {
    eskl[a] = 0;
    for (int b = 0; b < 4; b++)
    {
        if (b == a) continue;
        for (int g = 0; g < 4; g++)
        {
            if (g == b || g == a) continue;
            for (int d = 0; d < 4; d++)
            {
                if (d == g || d == b || d == a) continue;
                double sb = s[b]*genie::utils::g(b,b);
                double kg = k[g]*genie::utils::g(g,g);
                double ld = l[d]*genie::utils::g(d,d);
                eskl[a] += e(a,b,g,d)*sb*kg*ld;
            }
        }
    }
  }
    
  double kl = k[0]*l[0] - k[1]*l[1] - k[2]*l[2] - k[3]*l[3];
  double ks = k[0]*s[0] - k[1]*s[1] - k[2]*s[2] - k[3]*s[3];
        
  for (int a = 0; a < 4; a++)
  {
     double aux_plus  = kl + ml*ks;
     double aux_minus = kl - ml*ks;
     if (isLeftPolarized)
     {
        jp[a] = aux_plus  > 0 ? (l[a]*ks - s[a]*kl - 1i*eskl[a] + ml*k[a])/sqrt(aux_plus)   : 0;   //jp_\alpha
        jm[a] = aux_minus > 0 ? (-l[a]*ks + s[a]*kl + 1i*eskl[a] + ml*k[a])/sqrt(aux_minus) : 0;   //jm_\alpha
     }
     else
     {
        jp[a] =  aux_minus > 0 ? (l[a]*ks - s[a]*kl + 1i*eskl[a] - ml*k[a])/sqrt(aux_minus) : 0;   //jp_\alpha
        jm[a] =  aux_plus  > 0 ? (l[a]*ks - s[a]*kl + 1i*eskl[a] + ml*k[a])/sqrt(aux_plus)  : 0;   //jm_\alpha
     }
  }

  std::complex<double> LWpp(0, 0), LWpm(0, 0), LWmp(0, 0), LWmm(0, 0);
  for(int mu = 0; mu < 4; mu++)
  {
     for(int nu = mu;nu < 4; nu++)
     {
        LWpp += jp[mu]*std::conj(jp[nu])*NTensor(mu,nu); // Lpp_\mu\nu*W^\mu\nu
        LWpm += jp[mu]*std::conj(jm[nu])*NTensor(mu,nu); // Lpm_\mu\nu*W^\mu\nu
        LWmp += jm[mu]*std::conj(jp[nu])*NTensor(mu,nu); // Lmp_\mu\nu*W^\mu\nu
        LWmm += jm[mu]*std::conj(jm[nu])*NTensor(mu,nu); // Lmm_\mu\nu*W^\mu\nu
        if (mu != nu)
        {
            LWpp += jp[nu]*std::conj(jp[mu])*NTensor(nu,mu); // Lpp_\mu\nu*W^\mu\nu
            LWpm += jp[nu]*std::conj(jm[mu])*NTensor(nu,mu); // Lpm_\mu\nu*W^\mu\nu
            LWmp += jm[nu]*std::conj(jp[mu])*NTensor(nu,mu); // Lmp_\mu\nu*W^\mu\nu
            LWmm += jm[nu]*std::conj(jm[mu])*NTensor(nu,mu); // Lmm_\mu\nu*W^\mu\nu
        }
    }
  }
  
  std::complex<double> LWppmm = LWpp + LWmm;
  if (LWppmm.real() == 0 && LWppmm.imag() == 0)
  {
     polarization = TVector3(0, 0, 0);
     return;
  } 
  std::complex<double> rhopp = LWpp/LWppmm;
  std::complex<double> rhopm = LWpm/LWppmm;
  std::complex<double> rhomp = LWmp/LWppmm;
  std::complex<double> rhomm = LWmm/LWppmm;
  double PL = std::real(rhopp - rhomm);
  double PP = std::real(rhopm + rhomp);
  double PT = std::imag(rhomp - rhopm);
  
  TVector3 neutrinoMom3 = neutrinoMom.Vect();                                          
  TVector3 leptonMom3 = leptonMom.Vect();
  TVector3 Pz = leptonMom3.Unit();
  TVector3 Px = neutrinoMom3.Cross(leptonMom3).Unit();
  TVector3 Py = Pz.Cross(Px);
  polarization = PT*Px + PP*Py + PL*Pz;
//  std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
//  std::cout << "PL = " << PL << ", PP = " << PP << ", PT = " << PT << "\n";
//  std::cout << "UT@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n" << std::endl;
}
//____________________________________________________________________________
void  genie::utils::CalculatePolarizationVectorWithStructureFunctions(
                                TVector3 & polarization,
                                const TLorentzVector & neutrinoMom,
                                const TLorentzVector & leptonMom, 
                                const TLorentzVector & inNucleonMom,
                                const TLorentzVector & q4,
                                bool isLeftPolarized,
                                double M,
                                double W1,
                                double W2,
                                double W3,
                                double W4,
                                double W5,
                                double W6
)
{
  double M2 = M*M;
  double p[4], q[4], epq[4][4];
  p[0] = inNucleonMom.E();
  p[1] = inNucleonMom.Px();
  p[2] = inNucleonMom.Py();
  p[3] = inNucleonMom.Pz();
  
  q[0] = q4.E();
  q[1] = q4.Px();
  q[2] = q4.Py();
  q[3] = q4.Pz();
  // epsilon^\alpha\beta\gamma\delta p_\gamma q_\delta
  for (int a = 0; a < 4; a++)
  {
    for (int b = 0; b < 4; b++)
    {
        epq[a][b] = 0;
        if (b == a) continue;
        for (int g = 0; g < 4; g++)
        {
            if (g == b || g == a) continue;
            for (int d = 0; d < 4; d++)
            {
                if (d == g || d == b || d == a) continue;
                epq[a][b] += e(a,b,g,d)*genie::utils::g(a,a)*genie::utils::g(b,b)*p[g]*q[d];
            }
        }
    }
  }
  HermitianMatrix NucleonTensor(4);
  for(int mu = 0; mu < 4; mu++)
  {
     for(int nu = mu; nu < 4; nu++)
     {
        double Wreal = -g(mu,nu)*W1 + p[mu]*p[nu]*W2/M2 + q[mu]*q[nu]*W4/M2 + (p[mu]*q[nu] + q[mu]*p[nu])*W5/2/M2;
        double Wimag = epq[mu][nu]*W3/2/M2 + (q[mu]*p[nu] - p[mu]*q[nu])*W6/2/M2;
        NucleonTensor.set(mu, nu, Wreal - 1i*Wimag);  // W^\mu\nu
        if (mu != nu) NucleonTensor.set(nu, mu, Wreal + 1i*Wimag);
    }
  }
  
  CalculatePolarizationVectorWithNuclearTensor(
                                    polarization,
                                    neutrinoMom,
                                    leptonMom, 
                                    isLeftPolarized,
                                    NucleonTensor);
  
}
//____________________________________________________________________________
void  genie::utils::CalculatePolarizationVectorInTargetRestFrame(
                      TVector3 & polarization,
                      const TLorentzVector & neutrinoMomTRF,
                      const TLorentzVector & leptonMomTRF, 
                      bool isLeftPolarized,
                      double M,
                      double W1,
                      double W2,
                      double W3,
                      double W4,
                      double W5,
                      double W6
)
{
  double ml = leptonMomTRF.M();
  double ml2 = ml*ml;
  double M2 = M*M;
  double Ev = neutrinoMomTRF.E();
  double El = leptonMomTRF.E();
  double Pl = leptonMomTRF.P();
  double cost = TMath::Cos( neutrinoMomTRF.Angle(leptonMomTRF.Vect()) );
  double sint = TMath::Sqrt(1 - cost*cost);
  int sign = isLeftPolarized?-1:1;
  double auxm  = (El - Pl*cost)/2/M;
  double auxp  = (El + Pl*cost)/2/M;
  double aux1m = (Pl - El*cost)/2/M;
  double aux1p = (Pl + El*cost)/2/M;
  double aux1  = ml2/2/M2;
  double aux2  = (Ev + El)/M;
  double R     = 2*auxm*(W1 + aux1*W4) + auxp*W2 - sign*(aux2*auxm - aux1)*W3 - aux1*W5;
  if (R == 0)
  {
      polarization = TVector3(0, 0, 0);
      return;
  }
  double PL    = sign*(2*aux1m*(W1 - aux1*W4) + aux1p*W2 - sign*(aux2*aux1m + aux1*cost)*W3 - aux1*cost*W5)/R;
  double PP    = sign*ml*sint*(2*W1 - W2 -sign*Ev*W3/M - ml2*W4/M2 + El*W5/M)/2/M/R;
  double PT    = - ml*Pl*sint*W6/2/M2/R;
  
  
  TVector3 neutrinoMomTRF3 = neutrinoMomTRF.Vect();                                          
  TVector3 leptonMomTRF3 = leptonMomTRF.Vect();
  TVector3 Pz = leptonMomTRF3.Unit();
  TVector3 Px = neutrinoMomTRF3.Cross(leptonMomTRF3).Unit();
  TVector3 Py = Pz.Cross(Px);
  polarization = PT*Px + PP*Py + PL*Pz;
//  std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
//  std::cout << "PL = " << PL  << ", PP = " << PP << ", PT = " << PT << ", R = " << R << "\n";
//  std::cout << "UT@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n" << std::endl;
}
//____________________________________________________________________________


