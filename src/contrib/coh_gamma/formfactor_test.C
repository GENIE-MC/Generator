
#include "TString.h" 
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"

#include "Framework/Algorithm/AlgFactory.h"


#include "Framework/Ntuple/NtpMCTreeHeader.h" 
#include "Framework/Ntuple/NtpMCEventRecord.h" 
#include "Framework/EventGen/EventRecord.h" 
#include "Framework/ParticleData/BaryonResUtils.h"

#include "Framework/GHEP/GHepParticle.h"       
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/GHEP/GHepStatus.h"

#include "Framework/ParticleData/PDGCodes.h"

using namespace genie ;


void formfactor_test( TString algo_name  = "genie::DeVriesFormFactorMap" , 
		      TString algo_par_set = "Default",
		      TString out_file_name = "" ) {

  if ( out_file_name == "" ) {

    /// somethign clelver to get a good output name
  }
    
  std::map<int, TGraph*> proton_graphs, neutrons_graphs ; // int = pdg 
  
  const Algorithm * algo = AlgFactory::Instance()->GetAlgorithm( algo_name, alog_par_set ) ; 
  
  const COHFormFactorI * form_factor = dynamic_cast<const COHFormFactorI*>( algo ) ; 
    
  // loop over nuclei 

  for ( unsigned int z = 2 ; z < 100 ; ++z ) {  // a better limit coudl be nice, don't even know what z is
    for ( unsigned int n = z/2 ; n < 2*z ; ++n ) {
      
      int pdg = pdg::IonPdgCode( n+z, z ) ;

      if ( form_factor -> HasNucleus( pdg ) ) {
	
	Range1D_t q_range = form_factor -> QRange( pdg ) ;

	// build graph or hist for proton and neutron

	// fill them
	
	// save them into the maps

      }
    }
    
  }


  TFile out_file ( out_file_name, "RECREATE" ) ;
  out_file.cd() ;
  
  for ( auto & g : graphs ) {
    g.second -> Write() ;
    delete g.second ;
  }

}



