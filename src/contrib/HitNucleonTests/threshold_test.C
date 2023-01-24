#include <iostream>
#include <memory>
#include <cmath>
#include <limits>

#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"

#include "Framework/Interaction/Interaction.h"

#include "Math/Vector4D.h"
#include "Math/Boost.h" 
#include "TFile.h" 
#include "TH2.h"


using namespace genie;
using namespace std;


double find_threshold( const Interaction & i, double prec = 0.00001, double max = 10., double min = 1E-5 ) {

  Interaction in(i);

  // max - min > prec * 0.5 * ( max + min)
  while ( (2-prec)*max > (2+prec)*min ) {

    double centre = 0.5*(max+min);
    in.InitStatePtr()->SetProbeE(centre);
    
    if ( in.PhaseSpace().IsAboveThreshold() ) {
      max = centre;
    }
    else { 
      min = centre; 
    }
    
  }
  
  return 0.5*(min+max);

}


double new_threshold( const Interaction & i ) {

  int finalNucPDG = i.InitState().Tgt().HitNucPdg();
  if ( i.ProcInfo().IsWeakCC() ) finalNucPDG = pdg::SwitchProtonNeutron( finalNucPDG );

  double m_l = i.FSPrimLepton()->Mass();
  double m_p = i.InitState().Probe()->Mass();
  double E_T = i.InitState().Tgt().HitNucP4().E();
  double p_T = i.InitState().Tgt().HitNucP4().P();
  double W = PDGLibrary::Instance()->Find( finalNucPDG )->Mass();
  
  double alpha = m_l*m_l + W*W + 2*W*m_l - m_p*m_p - E_T*E_T + p_T*p_T ;

  double p_T_p = i.InitState().Tgt().HitNucP4().Pz();

  double chi = 2*m_p*E_T - alpha;

  auto a2 = alpha*alpha;
  auto Et2 = E_T*E_T;
  auto ptp2 = p_T_p*p_T_p;
  auto mp2 = m_p*m_p;

  auto delta = a2*Et2 - 4*(Et2-ptp2)*(ptp2*mp2+a2/4);

  double threshold = -1. ;

  if ( E_T > abs(p_T_p) ) {

    if ( p_T_p <= 0. ) {
      //cout << "Simple case, backward target" << endl;
      if ( chi <= 0. ) { 
        threshold = alpha*E_T;
        auto sqrt_delta = sqrt( delta );
        threshold -= sqrt_delta;
        threshold /= 2*(Et2 - ptp2);
      }
 
    } else { // p_t_p > 0  
      //cout << "Simple case, forward target" << endl;
      if ( delta >= 0. ) {
	threshold = alpha*E_T;
	auto sqrt_delta = sqrt( delta );
	threshold += sqrt_delta;
	threshold /= 2*(Et2 - ptp2);
      }
    }
  } else { // E_T <= abs(p_T_p)
    
    if ( p_T_p >= 0. ) {
      threshold = numeric_limits<double>::infinity();
    } else { 
      if ( chi <= 0. ) {
	threshold = alpha*E_T;
        auto sqrt_delta = sqrt( delta );
        threshold -= sqrt_delta;
        threshold /= 2*(Et2 - ptp2);
      }
    }
  }
  
  return max( m_p, threshold );
  
}


void print_initial_conditions( TLorentzVector tgt ) {

  cout << "Target: ( " << tgt.X() << ", " << tgt.Y() << ", " << tgt.Z() << ", " << tgt.E() << " ) with beta: " << tgt.Beta() << endl ;

}


void compare_threshold( double old_th, double new_th ) {

  cout << "Old threshold: " << old_th << " GeV" << endl;
  cout << "New threshold: " << new_th << " GeV" << endl;
  cout << "Relative difference:  " << (new_th - old_th)/old_th*100 << " %" << endl;


}



void threshold_test() {


  unique_ptr<Interaction> i_p {Interaction::QELCC( 1000180400, kPdgNeutron, kPdgNuMu, 2. )} ;

  cout << "--------- Neutrino on Target on shell, forward nucleon ------------" << endl;

  TLorentzVector on_shell_target;
  on_shell_target.SetXYZM( 0.1, -0.2, 0.25, PDGLibrary::Instance()->Find(kPdgNeutron)->Mass() );
  i_p -> InitStatePtr() -> TgtPtr() -> SetHitNucP4( on_shell_target ) ;

  print_initial_conditions( on_shell_target );
  
  auto old_th = find_threshold(*i_p);
  auto new_th = new_threshold(*i_p); 

  compare_threshold( old_th, new_th );
  

  cout << "--------- Neutrino on Target on shell, backward nucleon ------------" << endl;

  on_shell_target.SetXYZM( 0.1, -0.2, -0.25, PDGLibrary::Instance()->Find(kPdgNeutron)->Mass() );
  i_p -> InitStatePtr() -> TgtPtr() -> SetHitNucP4( on_shell_target ) ;

  print_initial_conditions( on_shell_target );
  
  old_th = find_threshold(*i_p);
  new_th = new_threshold(*i_p); 

  compare_threshold( old_th, new_th );

 
  cout << "--------- Neutrino on Target off shell, small beta, forward nucleon ------------" << endl;
  TLorentzVector off_shell_target;
  off_shell_target.SetXYZM( 0.1, -0.2, 0.25, PDGLibrary::Instance()->Find(kPdgNeutron)->Mass() - 0.2 );
  i_p -> InitStatePtr() -> TgtPtr() -> SetHitNucP4( off_shell_target ) ;

  print_initial_conditions( off_shell_target );
  
  old_th = find_threshold(*i_p);
  new_th = new_threshold(*i_p); 

  compare_threshold( old_th, new_th );

  cout << "--------- Neutrino on Target off shell, small beta, backward nucleon ------------" << endl;
  off_shell_target.SetXYZM( 0.1, -0.2, -0.25, PDGLibrary::Instance()->Find(kPdgNeutron)->Mass() - 0.2 );
  i_p -> InitStatePtr() -> TgtPtr() -> SetHitNucP4( off_shell_target ) ;

  print_initial_conditions( off_shell_target );
  
  old_th = find_threshold(*i_p);
  new_th = new_threshold(*i_p); 

  compare_threshold( old_th, new_th );
  
  cout << endl << "--------- Neutrino on Target off shell, beta > 1, forward nucleon --------" << endl;
  off_shell_target.SetXYZT( 0.2, 0.8, 0.1, PDGLibrary::Instance()->Find(kPdgNeutron)->Mass() - 0.25 );
  i_p -> InitStatePtr() -> TgtPtr() -> SetHitNucP4( off_shell_target ) ;

  print_initial_conditions( off_shell_target );

  new_th = new_threshold(*i_p); 

  cout << "New threshold: " << new_th << " GeV" << endl;

  cout << endl << "--------- Neutrino on Target off shell, beta > 1, backward nucleon --------" << endl;
  off_shell_target.SetXYZT( 0.2, 0.8, -0.1, PDGLibrary::Instance()->Find(kPdgNeutron)->Mass() - 0.25 );
  i_p -> InitStatePtr() -> TgtPtr() -> SetHitNucP4( off_shell_target ) ;

  print_initial_conditions( off_shell_target );

  new_th = new_threshold(*i_p); 

  cout << "New threshold: " << new_th << " GeV" << endl;

  cout << endl << "--------- Neutrino on Target off shell, beta > 1, forward nucleon --------" << endl;
  off_shell_target.SetXYZT( 0.2, 0.1, 0.8, PDGLibrary::Instance()->Find(kPdgNeutron)->Mass() - 0.25 );
  i_p -> InitStatePtr() -> TgtPtr() -> SetHitNucP4( off_shell_target ) ;

  print_initial_conditions( off_shell_target );

  new_th = new_threshold(*i_p); 

  cout << "New threshold: " << new_th << " GeV" << endl;

  cout << endl << "--------- Neutrino on Target off shell, beta > 1, backward nucleon --------" << endl;
  off_shell_target.SetXYZT( 0.2, 0.1, -0.8, PDGLibrary::Instance()->Find(kPdgNeutron)->Mass() - 0.25 );
  i_p -> InitStatePtr() -> TgtPtr() -> SetHitNucP4( off_shell_target ) ;

  print_initial_conditions( off_shell_target );

  new_th = new_threshold(*i_p); 

  cout << "New threshold: " << new_th << " GeV" << endl;



}


void threshold_exploration( int nucleon_pdg = kPdgNeutron, 
                            int probe_pdg = kPdgNuMu, 
			    double p_T_L = 0.2,
                            int target_pdg = 1000180400 ) {

  string file_name = PDGLibrary::Instance()->Find(probe_pdg)->GetName();
  file_name += "_on_" ;
  file_name +=  PDGLibrary::Instance()->Find(nucleon_pdg)->GetName() ;
  file_name += ".root" ;

  TFile out_file(file_name.c_str(), "RECREATE");

  constexpr double max_p = 0.9;
  double target_mass =PDGLibrary::Instance()->Find(nucleon_pdg)->Mass() ;
  double max_Et = sqrt( max_p*max_p + target_mass*target_mass );  
  string title = PDGLibrary::Instance()->Find(probe_pdg)->GetTitle();
  title += " on " ;
  title += PDGLibrary::Instance()->Find(nucleon_pdg)->GetTitle() ;
  TH2D new_th_hist( "new_threshold", (title + ";E_{T} (GeV);p_{T #parallel} (GeV/c)").c_str(), 500, .4, max_Et, 100, -max_p, max_p );
  TH2D old_th_hist( "old_threshold", (title + ";E_{T} (GeV);p_{T #parallel} (GeV/c)").c_str(), 500, .4, max_Et, 100, -max_p, max_p );

  unique_ptr<Interaction> i_p {Interaction::QELCC( target_pdg, nucleon_pdg, probe_pdg, 2. )} ;

  for ( int i = 1; i <= new_th_hist.GetNbinsX() ; ++i ) {
    double E_T = new_th_hist.GetXaxis()->GetBinCenter(i);
    for ( int j = 1; j <=  new_th_hist.GetNbinsY() ; ++j ) {
      double p_t_p = new_th_hist.GetYaxis()->GetBinCenter(j);

      TLorentzVector p_T;
      p_T.SetXYZT( 0., p_T_L, p_t_p, E_T );
      i_p -> InitStatePtr() -> TgtPtr() -> SetHitNucP4( p_T ) ;

      auto bin = new_th_hist.GetBin(i,j);

      auto new_th = new_threshold(*i_p);
      if ( new_th < numeric_limits<double>::infinity() )
	new_th_hist.SetBinContent(bin, new_th);
      
      auto old_th = find_threshold(*i_p);
      old_th_hist.SetBinContent(bin, old_th);

    }

  }

  new_th_hist.Write();
  old_th_hist.Write();

}
