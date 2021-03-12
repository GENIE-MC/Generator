#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"

#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Ntuple/NtpMCTreeHeader.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"


using namespace genie ;


void dnu_plots( TString in_file_name  = "gntp.0.ghep.root" ,
                TString out_file_name = "" ) {


  if ( out_file_name == "" ) {
    out_file_name = in_file_name ;
    out_file_name.ReplaceAll( ".ghep.root",
                              ".plots.root" ) ;
  }

  LOG("DNuPlots", pINFO) << "Reading file: " << in_file_name;
  TFile in_file( in_file_name ) ;

  NtpMCTreeHeader * header = dynamic_cast<NtpMCTreeHeader*> (in_file.Get("header"));

  // Get the GENIE GHEP tree and set its branch address
  TTree * in_tree = dynamic_cast<TTree*> (in_file.Get("gtree"));
  NtpMCEventRecord * mcrec = 0;
  in_tree->SetBranchAddress("gmcrec", & mcrec);


  std::map<std::string, TH1*> hists ;

  // plots we need
  TH1* h_E_probe = hists["E_probe"] = new TH1D("h_E_probe", "E_{probe};E_{probe} [GeV]",
                                               80, 0., 6. ) ;
  TH1* h_E_N = hists["E_N"] = new TH1D("h_E_N", "E_{N};E_{N} [GeV]",
                                       80, 0., 6. ) ;
  TH1* h_T_T = hists["T_T"] = new TH1D("h_T_T", "T_{T};T_{T} [GeV]",
                                       100, 0., 0.001 ) ;
  TH1* h_E_vis = hists["E_vis"] = new TH1D("h_E_vis", "E_{vis};E_{vis} [GeV]",
                                           100, 0., 6. ) ;
  TH1* h_vis_dec = hists["vis_dec"] = new TH1D("h_vis_dec", "Visible Decay;visible",
                                               2, -0.5, 1.5 ) ;
  TH1* h_E_prod = hists["E_prod"] = new TH1D("h_E_prod", "E_{prod};E_{prod} [GeV]",
                                             100, 0., 1. ) ;
  TH1* h_theta_N = hists["theta_N"] = new TH1D("h_theta_n", "#theta_{N};#theta_{N} [rad]",
                                               80, 0., TMath::Pi() ) ;
  TH1* h_N_prod = hists["N_prod"] = new TH1D("h_N_prod", "N_{Prod};N_{Prod}",
                                             10, -0.5, 9.5 ) ;
  TH1* h_decay_length = hists["decay_length"] = new TH1D("h_decay_length", "Decay Length;|#Delta x| [fm]",
                                                         500, 0., 1e13 ) ; // maximum is 1 cm
  TH1* h_theta_Med = hists["theta_Med"] = new TH1D("h_theta_Med", "#theta_{Z};#theta_{Z} [rad]",
                                                   80, 0., TMath::Pi() ) ;
  TH1* h_phi_ee = hists["phi_ee"] = new TH1D("h_phi_ee", "#phi_{ee};#phi_{ee} [rad]",
                                             80, 0., TMath::Pi() ) ;
  auto pr_E_ee = new TProfile("pr_E_ee", "Profile of E_{e+} versus E_{e-}",
                              100,0,1.,0,1.);
  auto pr2_E_ee = new TProfile2D("pr2_E_ee", "2D Profile of E_{e+} versus E_{e-}",
                                 100,0,1., 100,0,1.);


  // Event loop
  for(Long64_t i=0; i < in_tree->GetEntries(); i++) {

    in_tree->GetEntry(i);

    EventRecord & event = *(mcrec->event);

    const Interaction & inter = *( event.Summary() ) ;

    const ProcessInfo & proc_info = inter.ProcInfo() ;

    if ( !proc_info.IsCoherentElastic() ) { continue; }
    if ( !proc_info.IsDarkNeutralCurrent() ) { continue; }

    const GHepParticle & N = * event.Particle(2) ;
    if(N.Pdg() != kPdgDarkNeutrino && 
       N.Pdg() != kPdgAntiDarkNeutrino ){
      LOG("DNuPlots", pERROR)
        << "WARNING!\n Dark Neutrino wasn't found.\n"
        << "N.Name(): " << N.Name() << "\n";
      exit(1);
    };

    const GHepParticle * final_neutrino = nullptr ;
    const GHepParticle * mediator = nullptr ;
    std::vector<const GHepParticle*> final_products ;

    const GHepParticle * temp = event.Particle( N.FirstDaughter() ) ;

    // setting the final products might not be easy for ever
    // let's try to get it right once and for all
    if ( N.LastDaughter() - N.FirstDaughter() == 1 ) {
      // the dark neutrino decays only in two bodies
      // so one is the neutrino
      // the other is the dark mediator
      if ( pdg::IsNeutrino( TMath::Abs( temp -> Pdg() ) ) ) {
        final_neutrino = temp ;
        mediator = event.Particle( N.LastDaughter() ) ;
      }
      else {
        mediator = temp ;
        final_neutrino = event.Particle( N.LastDaughter() ) ;
      }

      // the final products are then the daughters of the mediator
      for ( unsigned int i = mediator -> FirstDaughter() ;
            i <= mediator -> LastDaughter() ; ++i ) {
        final_products.push_back( event.Particle( i ) ) ;
      }

    }
    else {
      // otherwise there is a neutrino and there is no mediator
      // so the neutrino is among the decay products fo the dark neutrino,
      // all the rest is products
      // impossible to say which one is the main neutrino, just pick the first in that case

      for ( unsigned int i = N.FirstDaughter() ;
            i <= N.LastDaughter() ; ++i ) {

        temp = event.Particle( i ) ;
        if ( ! final_neutrino ) {
          if ( pdg::IsNeutrino( TMath::Abs( temp -> Pdg() ) ) ) {
            final_neutrino = temp ;
            continue ;
          }
        }
        final_products.push_back( temp ) ;
      }

    }  // the neutrino decays in other than 2 daughters

    // now we have all the particles identified and we can fill the hists
    // note that the mediator pointer might be 0 as it does not necessarily exist

    h_N_prod -> Fill( final_products.size() ) ;

    const TLorentzVector & probe = * event.Probe() -> P4() ;

    const TLorentzVector & p4_N =  * N.P4()  ;
    const TLorentzVector & p4_recoil = * event.Particle(3) ->P4() ;


    h_E_probe -> Fill( probe.E() ) ;
    h_E_N -> Fill( p4_N.E() ) ;

    double t_t = p4_recoil.E() - p4_recoil.Mag() ;
    h_T_T -> Fill( t_t ) ;

    h_theta_N -> Fill( p4_N.Angle( probe.Vect() ) ) ;


    double vis_e = t_t ;
    double prod_e = 0. ;
    bool vis_decay = false ;
    std::vector<const GHepParticle*> charged_particles ;
    for ( const auto & p : final_products ) {
      // photon
      if ( p -> Pdg() == 22 )  vis_e += p -> P4() -> E() ;
      else if ( p -> Charge() != 0. ) {
        vis_decay = true ;
        vis_e += p -> P4() -> E() ;
        charged_particles.push_back(p) ;
      }

      prod_e += p -> P4() -> E() ;

    } // final_products loop

    // if the event has two charged particles it also has a
    // mediator initialised
    if ( charged_particles.size() == 2 ) {
      h_phi_ee -> Fill(
        charged_particles[0]->P4()->Angle(
          charged_particles[1]->P4()->Vect() ) ) ;
      if (charged_particles[0]->Charge()<0){
        pr_E_ee->Fill(charged_particles[0]->P4()->E(), charged_particles[1]->P4()->E(), 1);
        pr2_E_ee->Fill(charged_particles[0]->P4()->E(), charged_particles[1]->P4()->E(), 1);
      }
      else{
        pr_E_ee->Fill(charged_particles[1]->P4()->E(), charged_particles[0]->P4()->E(), 1);
        pr2_E_ee->Fill(charged_particles[1]->P4()->E(), charged_particles[0]->P4()->E(), 1);
      }
    }

    if ( mediator ) {
      // evaluate the lenght of the first and second decay
      TLorentzVector Delta = (* mediator -> X4()) - (* N.X4()) ;
      double total_distance = Delta.Vect().Mag() ;
      double delay = Delta.T() ;

      // std::cout << "Dark neutrino: " << total_distance << " fm in " << delay << " 10^-24 s and beta = " << N.P4() -> Beta() << std::endl ;
      Delta = (*final_products[0] -> X4()) - (* mediator -> X4()) ;
      total_distance = Delta.Vect().Mag() ;
      delay = Delta.T() ;

      // std::cout << "mediator: " << total_distance << " fm in " << delay << " 10^-24 s and beta = " << mediator -> P4() -> Beta() << std::endl ;
      h_theta_Med -> Fill(
        mediator->P4()->Angle( probe.Vect() ) ) ;
    }

    TLorentzVector Delta = (*final_products[0] -> X4()) - (* N.X4() ) ;
    double total_distance = Delta.Vect().Mag() ;
    double delay = Delta.T() ;

    h_decay_length -> Fill( total_distance ) ;
    //std::cout << "Total: " << total_distance << " fm in " << delay << " 10^-24 s" << std::endl ;

    h_vis_dec -> Fill( vis_decay ? 1 : 0 ) ;
    h_E_vis -> Fill( vis_e ) ;
    h_E_prod -> Fill( prod_e ) ;

    mcrec->Clear() ;
  } // event loop


  LOG("DNuPlots", pINFO) << "Writing to file: " << out_file_name;
  TFile out_file ( out_file_name, "RECREATE" ) ;
  out_file.cd() ;

  for ( auto & h : hists ) {
    SLOG("DNuPlots", pINFO) << "Plot: " << h.second->GetTitle();
    h.second -> Write() ;
    delete h.second ;
  }
  SLOG("DNuPlots", pINFO) << "Plot: " << pr_E_ee->GetTitle();
  pr_E_ee -> Write() ;
  delete pr_E_ee ;
  SLOG("DNuPlots", pINFO) << "Plot: " << pr2_E_ee->GetTitle();
  pr2_E_ee -> Write() ;
  delete pr2_E_ee ;

}
