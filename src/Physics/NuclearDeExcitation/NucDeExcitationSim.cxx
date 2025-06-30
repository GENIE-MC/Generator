//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool
*/
//____________________________________________________________________________

#include <cstdlib>
#include <sstream>

#include <TFile.h>
#include <TH1D.h>
#include <TMath.h>
#include <TParameter.h>
#include <TSystem.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Physics/NuclearDeExcitation/NucDeExcitationSim.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"

using std::ostringstream;

using namespace genie;
using namespace genie::utils;
using namespace genie::constants;
using namespace genie::controls;

//___________________________________________________________________________
NucDeExcitationSim::NucDeExcitationSim() :
EventRecordVisitorI("genie::NucDeExcitationSim")
{

}
//___________________________________________________________________________
NucDeExcitationSim::NucDeExcitationSim(string config) :
EventRecordVisitorI("genie::NucDeExcitationSim", config)
{

}
//___________________________________________________________________________
NucDeExcitationSim::~NucDeExcitationSim()
{

}
//___________________________________________________________________________
void NucDeExcitationSim::ProcessEventRecord(GHepRecord * evrec) const
{
  LOG("NucDeEx", pNOTICE)
     << "Simulating nuclear de-excitation gamma rays";

  GHepParticle * nucltgt = evrec->TargetNucleus();
  if (!nucltgt) {
    LOG("NucDeEx", pINFO)
      << "No nuclear target found - Won't simulate nuclear de-excitation";
    return;
  }

  if(nucltgt->Z()==8) this->OxygenTargetSim(evrec);

  // if oxygen, keep the existing genie behavior

  // if not oxygen, only simulate these for QE events
  // double nucleon knockout produces fewer photons.
  // too much energy and it all leaves as nucleon ejection.
  // either set them to zero, or set them to 1/10 probability.
  // proc_info.IsQuasiElastic();
  // RIK TODO  still need to test that this works.
  // Then move the photon spectrom to p-hole too
  // Then activate for Argon.
  // Then make plots of the spectra.
  // do not do this for coherent ?  MEC, RES, DIS only ?
  // ask if the exotic processes are used and should get it.

  RandomGen * rnd = RandomGen::Instance();

  if(evrec->Summary()->ProcInfo().IsQuasiElastic()){ // || rand < 0.1}
    //suppress_probability_factor = 1.0;
    if(nucltgt->Z()==6) this->CarbonTargetSim(evrec);
    // simulate argon with a dedicated calculation
    if(nucltgt->Z()==18) this->ArgonTargetSim(evrec);
  } else if(evrec->Summary()->ProcInfo().IsResonant() ||
	    evrec->Summary()->ProcInfo().IsMEC() ||
	    evrec->Summary()->ProcInfo().IsDeepInelastic()){

    double suppress_probability_factor = 0.20;
    // deexcitation is less likely after multiple nucleon knockout
    // because there is too much energy.
    // the decay will happen via nucleon or alpha ejection
    // and the KE will be carried away in that fashion, not gamma.
    if(rnd->RndDec().Rndm() < suppress_probability_factor){
      if(nucltgt->Z()==6) this->CarbonTargetSim(evrec);
      // simulate argon with the same probabilities.
      // Argoneut using FLUKA says
      if(nucltgt->Z()==18) this->ArgonTargetSim(evrec);
    } else {
      //std::cout << "Rik made none nonQE " << std::endl;
    }

  }


	   // else if rand < 0.1

  LOG("NucDeEx", pINFO)
     << "Done with this event";
}
//___________________________________________________________________________
void NucDeExcitationSim::OxygenTargetSim(GHepRecord * evrec) const
{
  LOG("NucDeEx", pNOTICE)
     << "Simulating nuclear de-excitation gamma rays for Oxygen target";

  //LOG("NucDeEx", pNOTICE) << *evrec;

  GHepParticle * hitnuc = evrec->HitNucleon();
  if(!hitnuc) return;

  bool p_hole = (hitnuc->Pdg() == kPdgProton);
  double dt   = -1;

  RandomGen * rnd = RandomGen::Instance();

  //
  // ****** P-Hole
  //
  if (p_hole) {
    //
    // * Define all the data required for simulating deexcitations of p-hole states
    //

    // > probabilities for creating a p-hole in the P1/2, P3/2, S1/2 shells
    double Pp12 = 0.25;              // P1/2
    double Pp32 = 0.47;              // P3/2
    double Ps12 = 1. - Pp12 - Pp32;  // S1/2

    // > excited state energy levels & probabilities for P3/2-shell p-holes
    const int np32 = 3;
    double p32Elv[np32] = { 0.00632, 0.00993, 0.01070 };
    double p32Plv[np32] = { 0.872,   0.064,   0.064   };
    // - probabilities for deexcitation modes of P3/2-shell p-hole state '1'
    double p32Plv1_1gamma  = 0.78;  // prob to decay via 1 gamma
    double p32Plv1_cascade = 0.22;  // prob to decay via gamma cascade

    // > excited state energy levels & probabilities for S1/2-shell p-holes
    const int ns12 = 11;
    double s12Elv[ns12] = {
               0.00309, 0.00368, 0.00385, 0.00444, 0.00492,
               0.00511, 0.00609, 0.00673, 0.00701, 0.00703, 0.00734 };
    double s12Plv[ns12] = {
               0.0625,  0.1875,  0.075,   0.1375,  0.1375,
               0.0125,  0.0125,  0.075,   0.0563,  0.0563,  0.1874  };
    // - gamma energies and probabilities for S1/2-shell p-hole excited
    //   states '2','7' and '10' with >1 deexcitation modes
    const int ns12lv2 = 3;
    double s12Elv2[ns12lv2]    = { 0.00309, 0.00369, 0.00385 };
    double s12Plv2[ns12lv2]    = { 0.013,   0.360,   0.625   };
    const int ns12lv7 = 2;
    double s12Elv7[ns12lv7]    = { 0.00609, 0.00673 };
    double s12Plv7[ns12lv7]    = { 0.04,    0.96    };
    const int ns12lv10 = 3;
    double s12Elv10[ns12lv10]  = { 0.00609, 0.00673, 0.00734 };
    double s12Plv10[ns12lv10]  = { 0.050,   0.033,   0.017   };

    // Select one of the P1/2, P3/2 or S1/2
    double rshell = rnd->RndDec().Rndm();
    //
    // >> P1/2 shell
    //
    if(rshell < Pp12) {
        LOG("NucDeEx", pNOTICE)
          << "Hit nucleon left a P1/2 shell p-hole. Remnant is at g.s.";
        return;
    }
    //
    // >> P3/2 shell
    //
    else
    if(rshell < Pp12 + Pp32) {
        LOG("NucDeEx", pNOTICE)
            << "Hit nucleon left a P3/2 shell p-hole";
        // Select one of the excited states
        double rdecmode  = rnd->RndDec().Rndm();
        double prob_sum  = 0;
        int    sel_state = -1;
        for(int istate=0; istate<np32; istate++) {
            prob_sum += p32Plv[istate];
            if(rdecmode < prob_sum) {
              sel_state = istate;
              break;
            }
        }
        LOG("NucDeEx", pNOTICE)
            << "Selected P3/2 excited state = " << sel_state;

        // Decay that excited state
        // >> 6.32 MeV state
        if(sel_state==0) {
            this->AddPhoton(evrec, p32Elv[0], dt);
        }
        // >> 9.93 MeV state
        else
        if(sel_state==1) {
            double r = rnd->RndDec().Rndm();
            // >>> emit a single gamma
            if(r < p32Plv1_1gamma) {
               this->AddPhoton(evrec, p32Elv[1], dt);
            }
            // >>> emit a cascade of gammas
            else
            if(r < p32Plv1_1gamma + p32Plv1_cascade) {
               this->AddPhoton(evrec, p32Elv[1],           dt);
               this->AddPhoton(evrec, p32Elv[1]-p32Elv[0], dt);
            }
        }
        // >> 10.7 MeV state
        else
        if(sel_state==2) {
           // Above the particle production threshold - need to emit
           // a 0.5 MeV kinetic energy proton.
           // Will neglect that given that it is a very low energy
           // kinetic energy nucleon and the intranuke break-up nucleon
           // cross sections are already tuned.
           return;
        }
    } //p3/2
    //
    // >> S1/2 shell
    //
    else if (rshell < Pp12 + Pp32 + Ps12) {
        LOG("NucDeEx", pNOTICE)
            << "Hit nucleon left an S1/2 shell p-hole";
        // Select one of the excited states caused by a S1/2 shell hole
        double rdecmode  = rnd->RndDec().Rndm();
        double prob_sum  = 0;
        int    sel_state = -1;
        for(int istate=0; istate<ns12; istate++) {
            prob_sum += s12Plv[istate];
            if(rdecmode < prob_sum) {
              sel_state = istate;
              break;
            }
        }
        LOG("NucDeEx", pNOTICE)
            << "Selected S1/2 excited state = " << sel_state;

        // Decay that excited state
        bool multiple_decay_modes =
              (sel_state==2 || sel_state==7 || sel_state==10);
        if(!multiple_decay_modes) {
          this->AddPhoton(evrec, s12Elv[sel_state], dt);
        } else {
          int ndec = -1;
          double * pdec = 0, * edec = 0;
          switch(sel_state) {
           case(2) :
              ndec = ns12lv2;  pdec = s12Plv2;  edec = s12Elv2;
              break;
           case(7) :
              ndec = ns12lv7;  pdec = s12Plv7;  edec = s12Elv7;
              break;
           case(10) :
              ndec = ns12lv10; pdec = s12Plv10; edec = s12Elv10;
              break;
           default:
             return;
          }
          double r = rnd->RndDec().Rndm();
          double decmode_prob_sum = 0;
          int sel_decmode = -1;
          for(int idecmode=0; idecmode < ndec; idecmode++) {
             decmode_prob_sum += pdec[idecmode];
             if(r < decmode_prob_sum) {
                 sel_decmode = idecmode;
                 break;
             }
          }
          if(sel_decmode == -1) return;
          this->AddPhoton(evrec, edec[sel_decmode], dt);
        }//mult.dec.ch

    } // s1/2
    else {
    }
  } // p-hole

  //
  // ****** n-hole
  //
  else {
    //
    // * Define all the data required for simulating deexcitations of n-hole states
    //

    // > probabilities for creating a n-hole in the P1/2, P3/2, S1/2 shells
    double Pp12 = 0.25;  // P1/2
    double Pp32 = 0.44;  // P3/2
    double Ps12 = 0.09;  // S1/2
    //>
    double p32Elv = 0.00618;
    //>
    double s12Elv = 0.00703;
    double s12Plv = 0.222;

    // Select one of the P1/2, P3/2 or S1/2
    double rshell = rnd->RndDec().Rndm();
    //
    // >> P1/2 shell
    //
    if(rshell < Pp12) {
        LOG("NucDeEx", pNOTICE)
          << "Hit nucleon left a P1/2 shell n-hole. Remnant is at g.s.";
        return;
    }
    //
    // >> P3/2 shell
    //
    else
    if(rshell < Pp12 + Pp32) {
        LOG("NucDeEx", pNOTICE)
            << "Hit nucleon left a P3/2 shell n-hole";
        this->AddPhoton(evrec, p32Elv, dt);
    }
    //
    // >> S1/2 shell
    //
    else
    if(rshell < Pp12 + Pp32 + Ps12) {
        LOG("NucDeEx", pNOTICE)
            << "Hit nucleon left a S1/2 shell n-hole";
        // only one of the deexcitation modes involve a (7.03 MeV) photon
        double r = rnd->RndDec().Rndm();
        if(r < s12Plv) this->AddPhoton(evrec, s12Elv,dt);
    }
    else {
    }
  } //n-hole
}

//___________________________________________________________________________
//  Carbon
//
//  As of October 2022, we are not prepared to turn on new GENIE functionality
//  that includes full simulation of remnant nucleus deexcitation
//  such as INCL++ or Marley (or FLUKA).
//
//  Based on some MINERvA effort at Duluth, undergrad Brandon Reed and Rik Gran
//  We have coded in the neutron hole decay spectrum predicted by
//  Kamyshkov and Kolbe @article{Kamyshkov:2002wp, arXiv:nucl-th/0206030
//  doi = "10.1103/PhysRevD.67.076007",
//
//  This paper proposes a decay spectrum using a similar concept to the one
//  for neutron-hole and proton-hole used for oxygen that is coded above.
//  They do not present a proton hole calculation.
//  The P3/2 neutron-hole decay leads only to a 2 MeV photon.
//  The S1/2 neutron-hole decay leads to a wide range of outcomes
//      of which this code only gives the photon spectrum, not ejected nucleons.
//  They also present a discussion of two-neutron hole spectra
//      which preferentially lead to nucleon ejection and less often a photon
//
void NucDeExcitationSim::CarbonTargetSim(GHepRecord * evrec) const
{
  // If we've disabled carbon de-excitations, then this function is a no-op
  if ( !fDoCarbon ) return;

  LOG("NucDeEx", pNOTICE)
     << "Simulating nuclear de-excitation gamma rays for Carbon target";

  GHepParticle * hitnuc = evrec->HitNucleon();
  if(!hitnuc) return;

  //  bool p_hole = (hitnuc->Pdg() == kPdgProton);



  double dt   = -1;

  RandomGen * rnd = RandomGen::Instance();


    //
    // * Define all the data required for simulating deexcitations of n-hole states
    //

    //std::cout << "Rik simulating n-hole in carbon " << std::endl;

    // > probabilities for creating a n-hole in the P1/2, P3/2, S1/2 shells
    // Kamyshkov gives it a different way than the oxygen folks did.
    // A probability for the Pstates combined, and branching fractions for S1/2
    double Pp12 = 0.0; //  0.75/6.0; // 0.25;  // P1/2  Rik says set to zero.
    double Pp32 = 4.0/6.0; // 0.44;  // P3/2
    double Ps12 = 2.00/6.0;  //0.09;  // S1/2
    //>
    double p32Elv = 0.0020;
    //>
    const int ns12 = 8;
    double s12Elv[ns12] = {0.0005, 0.0007, 0.0017, 0.0021, 0.0033, 0.0035, 0.0047, 0.0063};
    //double s12Plv[ns12] = {0.21, 0.295, 0.14, 0.26, 0.14, 0.2, 0.03, 0.03};
    // the above multiply by 0.2
    double s12Plv[ns12] = {0.042, 0.059, 0.028, 0.052, 0.028, 0.04, 0.006, 0.006};
    // the above multiply by 0.2 and by 2/6.
    //double s12Plv[ns12] = {0.0140, 0.01967, 0.0933, 0.01733, 0.00933, 0.0133, 0.0020, 0.0020};


    //double s12Elv = 0.00703;
    //double s12Plv = 0.222;

    // Select one of the P1/2, P3/2 or S1/2
    double rshell = rnd->RndDec().Rndm();
    //
    // >> P1/2 shell
    //
    if(rshell < Pp12) {
      // Rik says this probability is already set to zero for Kamyshkov
        LOG("NucDeEx", pNOTICE)
          << "Hit nucleon left a P1/2 shell n-hole. Remnant is at g.s.";
        return;
    }
    //
    // >> P3/2 shell
    //
    else
    if(rshell < Pp12 + Pp32) {
        LOG("NucDeEx", pNOTICE)
            << "Hit nucleon left a P3/2 shell n-hole";

	double myrand  = rnd->RndDec().Rndm();
	if(myrand < 0.2){
	  this->AddPhoton(evrec, p32Elv, dt);
	  //std::cout << "Rik made p32Elv " << p32Elv << std::endl;
	} else {
	  //std::cout << "Rik made none p32" << std::endl;
	}
    }
    //
    // >> S1/2 shell
    //
    else
    if(rshell < Pp12 + Pp32 + Ps12) {
       LOG("NucDeEx", pNOTICE)
            << "Hit nucleon left an S1/2 shell p-hole";
        // Select one of the excited states caused by a S1/2 shell hole
        double rdecmode  = rnd->RndDec().Rndm();
        double prob_sum  = 0.;
        int    sel_state = -1;
        for(int istate=0; istate<ns12; istate++) {
            prob_sum += s12Plv[istate];
            if(rdecmode < prob_sum) {
              sel_state = istate;
              break;
            }
        }
        LOG("NucDeEx", pNOTICE)
            << "Selected S1/2 excited state = " << sel_state;
	if(sel_state >= 0){
	  this->AddPhoton(evrec, s12Elv[sel_state], dt);
	  //std::cout << "Rik made s12Elv " << s12Elv[sel_state] << std::endl;
	} else {
	  //std::cout << "Rik made none s12" << std::endl;
	}
    }
    else {
      //std::cout << "Rik made none at all" << std::endl;
    }


}
//___________________________________________________________________________
void NucDeExcitationSim::AddPhoton(
                         GHepRecord * evrec, double E0, double dt) const
{
// Add a photon at the event record & recoil the remnant nucleus so as to
// conserve energy/momenta
//
  double E = (dt>0) ? this->PhotonEnergySmearing(E0, dt) : E0;

  LOG("NucDeEx", pNOTICE)
    << "Adding a " << E/units::MeV << " MeV photon from nucl. deexcitation";

  GHepParticle * target  = evrec->Particle(1);
  GHepParticle * remnant = 0;
  for(int i = target->FirstDaughter(); i <= target->LastDaughter(); i++) {
    remnant  = evrec->Particle(i);
    if(pdg::IsIon(remnant->Pdg())) break;
  }

  TLorentzVector x4(0,0,0,0);
  TLorentzVector p4 = this->Photon4P(E);
  GHepParticle gamma(kPdgGamma, kIStStableFinalState,1,-1,-1,-1, p4, x4);  // note that this assigns the parent of the photon as the initial-state nucleon/nucleus.  (do we want that??)
  evrec->AddParticle(gamma);


  remnant->SetPx     ( remnant->Px() - p4.Px() );
  remnant->SetPy     ( remnant->Py() - p4.Py() );
  remnant->SetPz     ( remnant->Pz() - p4.Pz() );
  remnant->SetEnergy ( remnant->E()  - p4.E()  );
}
//___________________________________________________________________________
double NucDeExcitationSim::PhotonEnergySmearing(double E0, double dt) const
{
// Returns the smeared energy of the emitted gamma
// E0 : energy of the excited state (GeV)
// dt: excited state lifetime (sec)
//
  double dE = kPlankConstant / (dt*units::s);

  RandomGen * rnd = RandomGen::Instance();
  double E = rnd->RndDec().Gaus(E0 /*mean*/, dE /*sigma*/);

  LOG("NucDeEx", pNOTICE)
     << "<E> = " << E0 << ", dE = " << dE << " -> E = " << E;

  return E;
}
//___________________________________________________________________________
TLorentzVector NucDeExcitationSim::Photon4P(double E) const
{
// Generate a photon 4p

  RandomGen * rnd = RandomGen::Instance();

  double costheta = -1. + 2. * rnd->RndDec().Rndm();
  double sintheta = TMath::Sqrt(TMath::Max(0., 1.-TMath::Power(costheta,2)));
  double phi      = 2*kPi * rnd->RndDec().Rndm();
  double cosphi   = TMath::Cos(phi);
  double sinphi   = TMath::Sin(phi);

  double px = E * sintheta * cosphi;
  double py = E * sintheta * sinphi;
  double pz = E * costheta;

  TLorentzVector p4(px,py,pz,E);
  return p4;
}
//___________________________________________________________________________
// Added on 27 January 2023 by S. Gardiner based on a preliminary MARLEY
// calculation. Only the leading de-excitation photon is generated for now.
void NucDeExcitationSim::ArgonTargetSim( GHepRecord* evrec ) const
{
  // If we've disabled argon de-excitations, then this function is a no-op
  if ( !fDoArgon ) return;

  LOG( "NucDeEx", pNOTICE ) << "Simulating nuclear de-excitation gamma-rays"
    " for an argon target";

  // Load the leading gamma-ray spectra and emission probabilities from prior
  // simulation results
  static bool loaded_hists = false;
  static TH1D* hist_n = nullptr;
  static TH1D* hist_p = nullptr;
  static double prob_gamma_n = 0.;
  static double prob_gamma_p = 0.;

  if ( !loaded_hists ) {
    std::string data_file_name = std::string( gSystem->Getenv("GENIE") )
      + "/data/evgen/nucl/marley_argon_sf_lead_gamma_hists.root";

    TFile data_file( data_file_name.c_str(), "read" );

    TParameter< double >* prob_n = nullptr;
    TParameter< double >* prob_p = nullptr;

    data_file.GetObject( "hist_n", hist_n );
    data_file.GetObject( "hist_p", hist_p );

    hist_n->SetDirectory( nullptr );
    hist_p->SetDirectory( nullptr );

    data_file.GetObject( "prob_gamma_n", prob_n );
    data_file.GetObject( "prob_gamma_p", prob_p );

    assert( hist_n && hist_p && prob_n && prob_p );

    prob_gamma_n = prob_n->GetVal();
    prob_gamma_p = prob_p->GetVal();

    loaded_hists = true;
  }

  GHepParticle* hitnuc = evrec->HitNucleon();
  if ( !hitnuc ) return;

  // Choose the appropriate gamma-ray distribution based on the struck nucleon
  int hit_nuc_pdg = hitnuc->Pdg();
  if ( !genie::pdg::IsNucleon(hit_nuc_pdg) ) return;

  TH1D* lead_gamma_hist = hist_n;
  double gamma_prob = prob_gamma_n;
  if ( hit_nuc_pdg == kPdgProton ) {
    lead_gamma_hist = hist_p;
    gamma_prob = prob_gamma_p;
  }

  // Throw a random number to see if this de-excitation event contains at least
  // one gamma-ray. If it doesn't, just return without doing anything.
  RandomGen* rnd = RandomGen::Instance();
  double r_gamma = rnd->RndDec().Rndm();
  if ( r_gamma > gamma_prob ) {
    LOG( "NucDeEx", pNOTICE ) << "No gamma-ray emitted";
    return;
  }
  // If it does, then sample the leading gamma from the pre-saved distribution
  // and add it to the event. Note that the MARLEY histogram is in MeV, so we
  // also convert to GeV here for consistency with GENIE conventions.
  double E_gamma = lead_gamma_hist->GetRandom() / 1e3;
  this->AddPhoton( evrec, E_gamma, -1. );

  LOG( "NucDeEx", pNOTICE ) << "Added gamma with energy " << E_gamma << " GeV";
}
//_________________________________________________________________________
void NucDeExcitationSim::Configure( const Registry& config )
{
  Algorithm::Configure( config );
  this->LoadConfig();
}
//____________________________________________________________________________
void NucDeExcitationSim::Configure( std::string config )
{
  Algorithm::Configure( config );
  this->LoadConfig();
}
//_________________________________________________________________________
void NucDeExcitationSim::LoadConfig()
{
  // Boolean flag that enables/disables de-excitation handling for carbon
  GetParamDef( "DoCarbon", fDoCarbon, false );

  // Boolean flag that enables/disables de-excitation handling for argon
  GetParamDef( "DoArgon", fDoArgon, false );
}
