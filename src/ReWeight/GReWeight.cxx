//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Aug 01, 2009 - CA
   Was adapted from Jim's and Costas' T2K-specific GENIE reweighting code.
   First included in v2.5.1.

*/
//____________________________________________________________________________

#include <vector>

#include <TMath.h>
#include <TString.h>

#include "EVGCore/EventRecord.h"
#include "Messenger/Messenger.h"
#include "ReWeight/GReWeight.h"
#include "ReWeight/GReWeightNuXSec.h"
#include "ReWeight/GReWeightAGKY.h"
#include "ReWeight/GReWeightFZone.h"
#include "ReWeight/GReWeightINuke.h"

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::string;

using namespace genie;
using namespace genie::rew;

//____________________________________________________________________________
GReWeight::GReWeight()
{
  this->Init();
}
//____________________________________________________________________________
GReWeight::~GReWeight()
{
  this->CleanUp();
}
//____________________________________________________________________________
GSystSet & GReWeight::Systematics(void)
{ 
  return fSystSet; 
}
//____________________________________________________________________________
void GReWeight::Reconfigure(void)
{
  LOG("ReW", pNOTICE) << "Reconfiguring ...";

  fSystSet.PrintSummary();

  vector<genie::rew::GSyst_t> svec = fSystSet.AllIncluded();

  vector<genie::rew::GSyst_t>::const_iterator it = svec.begin();
  for( ; it != svec.end(); ++it) {
    GSyst_t syst = *it;
    double val = fSystSet.CurValue(syst);

    fReWeightNuXSec -> SetSystematic(syst, val);
    fReWeightAGKY   -> SetSystematic(syst, val);
    fReWeightFZone  -> SetSystematic(syst, val);
    fReWeightINuke  -> SetSystematic(syst, val);
  }

  fReWeightNuXSec -> Reconfigure();
  fReWeightAGKY   -> Reconfigure();
  fReWeightFZone  -> Reconfigure();
  fReWeightINuke  -> Reconfigure();

  LOG("ReW", pNOTICE) << "Done reconfiguring";
}
//____________________________________________________________________________
double GReWeight::CalcWeight(const genie::EventRecord & event) 
{
// calculate weight for all tweaked physics parameters
//
  double weight_xsec  = fReWeightNuXSec -> CalcWeight(event);  // cross sections
  double weight_agky  = fReWeightAGKY   -> CalcWeight(event);  // hadronization
  double weight_fzone = fReWeightFZone  -> CalcWeight(event);  // form. zone
  double weight_inuke = fReWeightINuke  -> CalcWeight(event);  // intranuke

  double weight = weight_xsec * 
                  weight_agky *
                  weight_fzone *
                  weight_inuke;

  return weight;
}
//____________________________________________________________________________
double GReWeight::CalcChisq(void) 
{
// calculate the sum of penalty terms for all tweaked physics parameters

  double chisq_xsec  = fReWeightNuXSec -> CalcChisq(); 
  double chisq_agky  = fReWeightAGKY   -> CalcChisq();  
  double chisq_fzone = fReWeightFZone  -> CalcChisq();  
  double chisq_inuke = fReWeightINuke  -> CalcChisq(); 

  double chisq = TMath::Max(0., chisq_xsec ) +
                 TMath::Max(0., chisq_agky ) +
                 TMath::Max(0., chisq_fzone) +
                 TMath::Max(0., chisq_inuke);

  return chisq;
}
//____________________________________________________________________________
void GReWeight::Init(void)
{
  fReWeightNuXSec = new GReWeightNuXSec;
  fReWeightAGKY   = new GReWeightAGKY;
  fReWeightFZone  = new GReWeightFZone;
  fReWeightINuke  = new GReWeightINuke;
}
//____________________________________________________________________________
void GReWeight::CleanUp(void)
{
  delete fReWeightNuXSec; 
  delete fReWeightAGKY; 
  delete fReWeightFZone; 
  delete fReWeightINuke; 
}
//____________________________________________________________________________
void GReWeight::Print()
{
  vector<genie::rew::GSyst_t> syst_vec = this->Systematics().AllIncluded();
  int vec_size = syst_vec.size();
	
  cout << "\nCurrent list of tweaking parameter dials:"<< endl;
  for(int i = 0 ; i < vec_size ; i ++){
     LOG("ReW", pNOTICE) 
         << "    "     << GSyst::AsString(syst_vec[i])
         << " set at " << this->Systematics().CurValue(syst_vec[i]);
  }		       	        

  double chi2val = this->CalcChisq();

  LOG("ReW", pNOTICE) << "Chisq_{penalty} = " << chi2val;
}
//____________________________________________________________________________


