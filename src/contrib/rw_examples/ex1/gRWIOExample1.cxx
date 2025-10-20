//____________________________________________________________________________
/*!

\program gRWIOExample1

\brief   Simple example to demostrate use of GENIE ReWeight I/O infrastructure.
         It is a 2-step procedure.
	 First, it processes ALL events in the input file, but performs re-weighting
	 only where it is applicable.
	 An "empty" RW IO record will still be generated and written out for those 
	 events where re-weighting does not apply.
	 Second, it opens the GENIE file again, this time in the "READ" mode, and
	 extracts and prints the RW information. 
	 WARNING-1: this example is NOT entirely "fool proof" !!!
	 WARNING-2: this example can be quite verbose, due to all the printouts at step-2 !!!

         Syntax :
           gRWIOExample1 -f input_ghep_file
	   
          -f : Input GENIE event file

\author(s) Julia Yarba (FNAL)

\created May 11, 2016

\cpright Copyright (c) 2003-2025, The GENIE Collaboration
         All rights reserved.
         For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#include <string>
#include <sstream>
#include <cassert>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>

#include "EVGCore/EventRecord.h"
#include "Ntuple/NtpMCEventRecord.h"

#include "Interaction/Interaction.h"

#include "ReWeight/GReWeightI.h"
#include "ReWeight/GSystSet.h"
#include "ReWeight/GSyst.h"
#include "ReWeight/GReWeight.h"
#include "ReWeight/GSystUncertainty.h"

// Modules to calc weights/uncertainties
// by varying specific variables
//
// This one right below is for MaCCQE
//
#include "ReWeight/GReWeightNuXSecCCQE.h"

// I/O for re-weighting info 
//
#include "ReWeight/GReWeightIOBranchDesc.h"
#include "ReWeight/GReWeightIORecord.h"

#include "Utils/CmdLnArgParser.h"
  
int main( int argc, char ** argv )
{
      
   // GENIE event file on input
   //
   std::string isample = "";
   
   genie::CmdLnArgParser* lp = new genie::CmdLnArgParser( argc, argv );
   if ( lp->OptionExists('f') )
   {
      isample = lp->ArgAsString('f');
   }
   
   if ( isample.empty() )
   {
      // Nothing to re-weight, bail out... may also a warning will be useful...
      std::cout << " Missing GENIE event file on input " << std::endl;
      return 1;
   }
      
   // Define what parameter we want to vary/tweak
   //
   // In principle, it should be run-time configurable
   // But for simplicity we will define it here in this particular case
   //
   genie::rew::GSyst_t param_to_tweak = genie::rew::kXSecTwkDial_MaCCQE ;

   // Create a GReWeight object and add to it a set of weight calculators
   //
   genie::rew::GReWeight rw;
   
   // Add weight calculator for MaCCQE  
   // NOTE: will add other weight calculators later
   //
   rw.AdoptWghtCalc( "xsec_ccqe", new genie::rew::GReWeightNuXSecCCQE() );

   // Get GSystSet and include the (single) input systematic parameter
   // NOTE: in this case we use kXSecTwkDial_MaCCQE (for "MaCCQE")
   //
   genie::rew::GSystSet& syst = rw.Systematics();
   // syst.Init( genie::rew::kXSecTwkDial_MaCCQE );
   syst.Init( param_to_tweak );
   
   // By default GReWeightNuXSecCCQE is in `NormAndMaShape' mode 
   // where Ma affects the shape of dsigma/dQ2 and a different param affects the normalization
   // If the input is MaCCQE, switch the weight calculator to `Ma' mode
   //
   genie::rew::GReWeightNuXSecCCQE* rwccqe = dynamic_cast<genie::rew::GReWeightNuXSecCCQE*>( rw.WghtCalc("xsec_ccqe") );  
   rwccqe->SetMode( genie::rew::GReWeightNuXSecCCQE::kModeMa );

   // Now open input GENIE sample (fetched at the beginning of the job)
   //
   TFile* file = new TFile( isample.c_str(), "UPDATE" );
      
   // if invalid input file, bail out
   //
   if ( !file ) 
   {
       std::cout << " Can NOT open input GENIE file " << isample << std::endl;
       return 1;
   }

   //
   // Fetch or create a tree for RW records 
   //
   TTree* rwtree = 0;
   rwtree = dynamic_cast<TTree*>( file->Get("reweighting") );
   if ( !rwtree )
   {
      rwtree = new TTree( "reweighting", "GENIE weights tree" );
      TTree::SetBranchStyle(1);
      rwtree->SetAutoSave( 200000000 );  // autosave when 0.2 Gbyte written 
                                         // - it's the same for "gtree" but I need to double check 
				         // how to *get* autosave from "gtree" 
   }
   
   // Now create a branch to correspond to a specific parameter/variable to vary
   //
   // FIXME !!!
   // In principle, one should also check if a branch, and the corresponding metadata already exist, etc. 
   //
   genie::rew::GReWeightIORecord* rwrec = 0;
   std::string param_name = genie::rew::GSyst::AsString( param_to_tweak ); 
   TBranch* rwbr = rwtree->Branch( param_name.c_str(), 
                                   "genie::rew::GReWeightIORecord", 
				   &rwrec, 32000, 1 ); // FIXME !!! also check more "sophisticated" options  
   assert(rwbr); 
   rwbr->SetAutoDelete(kFALSE);  

   // Add meta-data (UserInfo) to the RW tree
   //
   // MaCCQE=0.99 has been extracted in the course of run under gdb from GReWeightNuXSecCCQE class 
   // (member data fMaDef).
   // In principle, it depends on the physics model - in the case of GReWeightNuXSecCCQE the model 
   // is based on the "LwlynSmithQELCCPXSec" algorithm (see GReWeightNuXSecCCQE::Init() method)
   // A more uniform machinery to access such information would be useful, but details of it need
   // to be discussed additionally.
   //
   // Sigma's (+/-) can be extracted from GSystUncertainty
   //
   genie::rew::GSystUncertainty* syser = genie::rew::GSystUncertainty::Instance();
   double sigpls = syser->OneSigmaErr( param_to_tweak,  1 );
   double sigmin = syser->OneSigmaErr( param_to_tweak, -1 );
   //
   rwtree->GetUserInfo()->Add( new genie::rew::GReWeightIOBranchDesc( param_name, 0.99, sigpls, sigmin ) ); 
      
   // Fetch the Evt tree
   //
   TTree* evtree = dynamic_cast<TTree*>( file->Get("gtree")  );
   
   // Connect Evt record (branch)
   //
   genie::NtpMCEventRecord* mcrec = 0;
   evtree->SetBranchAddress( "gmcrec", &mcrec );
   
   // "Tie" together these trees, Evt & RW !
   //
   evtree->AddFriend( rwtree );
   
   // now loop over events and see what needs to be re-weighted
   //
   int nevt_total = evtree->GetEntries();   
      
   double twk = 0.;
   double wt = 1.;
      
   for ( int iev=0; iev<nevt_total; ++iev )
   {
   
       evtree->GetEntry(iev);
       genie::EventRecord& evt = *(mcrec->event);
       
       //
       // Select events to be re-weighted
       //
       // Specifically, check if it's QEL && WeakCC process
       // because that's what we want to reweight (MaCCQE).
       // Also skip charm events (although if those are quite rare in this case)
       //
       genie::Interaction* interaction = evt.Summary();    
       const genie::ProcessInfo& prinfo = interaction->ProcInfo();   
       const genie::XclsTag&     xclsv  = interaction->ExclTag();
       bool accept = ( prinfo.IsQuasiElastic() && prinfo.IsWeakCC() && !xclsv.IsCharmEvent() );
       if ( !accept ) 
       {
          rwrec = new genie::rew::GReWeightIORecord();
          rwrec->SetOriginalEvtNumber(iev);
          rwtree->Fill();
          delete rwrec;
          rwrec=0;
          continue ;
       }
          
       rwrec = new genie::rew::GReWeightIORecord();
       
       rwrec->SetOriginalEvtNumber( iev );

       twk = -0.5; // in the units of MaCCQE SIGMA !!! (that's how the weigh calculator "understands" it)
       syst.Set( param_to_tweak, twk );
       rw.Reconfigure();
       wt = rw.CalcWeight(evt);
       rwrec->Insert( twk, wt );
       
       twk = 0.5;
       syst.Set( param_to_tweak, twk );
       rw.Reconfigure();
       wt = rw.CalcWeight(evt);
       rwrec->Insert( twk, wt );
       
       rwtree->Fill();
              
       // Clear mc evt record before the next one
       //
       mcrec->Clear();
       if ( rwrec ) delete rwrec;
       rwrec = 0;
         
   }
  
   file->cd();
   rwtree->Write("",TObject::kOverwrite);
   
   delete rwtree;
   rwtree=0;

   file->Close();
   
   delete lp; // destroy input args line parser 
    

   // RE_TEST NOW !!!
   //
   // Open up Genie event file (now also with the RW tree in it), 
   // this time in the READ mode
   //   
   TFile* tfile = new TFile( isample.c_str(), "READ" );

   //
   // Fetch the RW tree
   //
   TTree* rwtree_test = dynamic_cast<TTree*>( tfile->Get("reweighting") );
   
   TList* hdr = rwtree_test->GetUserInfo();
   assert(hdr);
      
   int nentries = rwtree_test->GetEntries();
   int nrw = 0;
   std::cout << " num of entries of re-test: " << nentries << std::endl;   
   genie::rew::GReWeightIORecord* rwrec_test = 0;
   rwtree_test->SetBranchAddress( param_name.c_str(), &rwrec_test );   
   for ( int i=0; i<nentries; ++i )
   {
       rwtree_test->GetEntry(i);
       // genie::EventRecord& evt = *(mcrec->event);
       int nres = rwrec_test->GetNumOfRWResults();
       if ( nres <= 0 ) continue;
       for ( int ir=0; ir<nres; ++ir )
       {
          double twk_test = rwrec_test->GetTweak( ir );
	  double wt_test  = rwrec_test->GetWeight( ir );
	  std::cout << " twk = " << twk_test << " wt = " << wt_test << std::endl;
       }
       nrw++;
   }

   std::cout << " num of re-weighted results of re-test: " << nrw << std::endl;   

   // now print meta-data
   //
   int nhdr = hdr->GetEntries();
   for ( int i=0; i<nhdr; ++i )
   {
      genie::rew::GReWeightIOBranchDesc* brdesc = dynamic_cast<genie::rew::GReWeightIOBranchDesc*>( hdr->At(i) );
      std::cout << " branch name: " << brdesc->GetParameterName() << std::endl;
      std::cout << " parameter: " << brdesc->GetParameterMean() << " " 
                << brdesc->GetParameterSigmaPlus() << " " << brdesc->GetParameterSigmaMinus() << std::endl; 
   }

   return 0;

}

