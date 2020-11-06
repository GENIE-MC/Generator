
#include "TString.h" 
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/Algorithm.h"
#include "Physics/Coherent/XSection/COHFormFactorI.h"

#include "Framework/Ntuple/NtpMCTreeHeader.h" 
#include "Framework/Ntuple/NtpMCEventRecord.h" 
#include "Framework/EventGen/EventRecord.h" 
#include "Framework/ParticleData/BaryonResUtils.h"

#include "Framework/GHEP/GHepParticle.h"       
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/GHEP/GHepStatus.h"

#include "Framework/ParticleData/PDGCodes.h"

using namespace genie ;


void formfactor_test( std::string algo_name  = "genie::DeVriesFormFactorMap" , 
		      std::string algo_par_set = "Default",
		      TString out_file_name = "" ) {

  if ( out_file_name == "" ) {
    
    std::string name = algo_name.substr( algo_name.find("::") + 2 ) ;
    name += "_ff_graphs.root" ;
        
    out_file_name = name.c_str() ;

  }
    
  std::map<int, TGraph*> proton_graphs, neutron_graphs ; // int = pdg 
  
  AlgId id( algo_name, algo_par_set );
  const Algorithm * algo = AlgFactory::Instance()->GetAlgorithm( id ) ; 
  
  const COHFormFactorI * form_factor = dynamic_cast<const COHFormFactorI*>( algo ) ; 

  for ( unsigned int z = 1 ; z < 100 ; ++z ) {  // a better limit coudl be nice, don't even know what z is
    for ( unsigned int n = z/2 ; n <= 2*z ; ++n ) {
      
      int pdg = pdg::IonPdgCode( n+z, z ) ;

      if ( form_factor -> HasNucleus( pdg ) ) {
	
        std::cout << "FF Has Nucleus " << pdg <<  std::endl;
	Range1D_t q_range = form_factor -> QRange( pdg ) ;

        int nQ = 200;
        double Q = q_range.min;
        double ff_p_arr[nQ];
        double ff_n_arr[nQ];
        double Q_arr[nQ];

	double delta_q = (q_range.max - q_range.min ) / (nQ-1) ; 

        for ( int Qstep = 0; Qstep < nQ; Qstep++ ) {
          ff_p_arr[Qstep] = form_factor -> ProtonFF( Q, pdg ) ;
          ff_n_arr[Qstep] = form_factor -> NeutronFF( Q, pdg ) ;
          Q_arr[Qstep] = Q;
	  Q += delta_q ; 
        }

        std::string nucleus( PDGLibrary::Instance() -> Find( pdg ) -> GetTitle() ) ;
        std::string p_title = "Proton Form Factor for " + nucleus + ";Q [GeV];FF";
        std::string n_title = "Neutron Form Factor for " + nucleus + ";Q [GeV];FF";

        proton_graphs[pdg] = new TGraph( nQ, Q_arr, ff_p_arr );
        proton_graphs[pdg]->SetName( ("proton_FF_" + nucleus).c_str() );
        proton_graphs[pdg]->SetTitle( p_title.c_str() );

        neutron_graphs[pdg] = new TGraph( nQ, Q_arr, ff_n_arr );
        neutron_graphs[pdg]->SetName( ("neutron_FF_" + nucleus).c_str() );
        neutron_graphs[pdg]->SetTitle( n_title.c_str() );

      }
    }
  }


  TFile out_file ( out_file_name, "RECREATE" ) ;
  out_file.cd() ;
  
  for ( auto & g : proton_graphs ) {
    g.second -> Write() ;
    delete g.second ;
  }

    for ( auto & g : neutron_graphs ) {
    g.second -> Write() ;
    delete g.second ;
  }


}



