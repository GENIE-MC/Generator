//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Aug 01, 2009 - CA
   Was adapted from Jim's and Costas' T2K-specific GENIE reweighting code.
   First included in v2.5.1.
 @ May 18, 2010 - CA
   AdoptWghtCalc(string,GReWeightI*) allows user to decide which weight 
   calculator to include. Weight calculators are owned by GReWeight and are 
   identified by a name. Weight calculators can be retrieved via the 
   WghtCalc(string) method and their reweighting options can be fine-tuned.
*/
//____________________________________________________________________________

#include <vector>

#include <TMath.h>
#include <TString.h>

#include "EVGCore/EventRecord.h"
#include "Messenger/Messenger.h"
#include "Utils/RunOpt.h"
#include "ReWeight/GReWeight.h"

using std::vector;

using namespace genie;
using namespace genie::rew;

//____________________________________________________________________________
GReWeight::GReWeight()
{
  // Disable cacheing that interferes with event reweighting
  RunOpt::Instance()->EnableBareXSecPreCalc(false);
}
//____________________________________________________________________________
GReWeight::~GReWeight()
{
  this->CleanUp();
}
//____________________________________________________________________________
void GReWeight::AdoptWghtCalc(string name, GReWeightI* wcalc)
{
  if(!wcalc) return;

  fWghtCalc.insert(map<string, GReWeightI*>::value_type(name,wcalc));
}
//____________________________________________________________________________
GReWeightI* GReWeight::WghtCalc(string name)
{ 
  map<string, GReWeightI*>::iterator iter = fWghtCalc.find(name);
  if(iter != fWghtCalc.end()) return iter->second;
  
  return 0;
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

  vector<genie::rew::GSyst_t> svec = fSystSet.AllIncluded();

  map<string, GReWeightI *>::iterator it = fWghtCalc.begin();
  for( ; it != fWghtCalc.end(); ++it) {

      GReWeightI * wcalc = it->second;

      vector<genie::rew::GSyst_t>::const_iterator parm_iter = svec.begin();
      for( ; parm_iter != svec.end(); ++parm_iter) {
          GSyst_t syst = *parm_iter;
          double val = fSystSet.Info(syst)->CurValue;
          wcalc->SetSystematic(syst, val);
      }//params

      wcalc->Reconfigure();

  }//weight calculators

  LOG("ReW", pDEBUG) << "Done reconfiguring";
}
//____________________________________________________________________________
double GReWeight::CalcWeight(const genie::EventRecord & event) 
{
// calculate weight for all tweaked physics parameters
//
  double weight = 1.0;
  map<string, GReWeightI *>::iterator it = fWghtCalc.begin();
  for( ; it != fWghtCalc.end(); ++it) {
    GReWeightI * wcalc = it->second;
    double w = wcalc->CalcWeight(event); 
    LOG("ReW", pNOTICE) 
       << "Calculator: " << it->first << " => wght = " << w;	
    weight *= w;
  }
  return weight;
}
//____________________________________________________________________________
double GReWeight::CalcChisq(void) 
{
// calculate the sum of penalty terms for all tweaked physics parameters
//
  double tot_chisq = 0.0;

  map<string, GReWeightI *>::iterator it = fWghtCalc.begin();
  for( ; it != fWghtCalc.end(); ++it) {
    GReWeightI * wcalc = it->second;
    double chisq = wcalc->CalcChisq(); 
    LOG("ReW", pNOTICE) 
       << "Calculator: " << it->first << " => chisq = " << chisq;	
    tot_chisq *= chisq;
  }
  return tot_chisq;
}
//____________________________________________________________________________
void GReWeight::CleanUp(void)
{
  map<string, GReWeightI *>::iterator it = fWghtCalc.begin();
  for( ; it != fWghtCalc.end(); ++it) {
    GReWeightI * rw = it->second;
    if(rw) {
      delete rw;
      rw=0;
    }
  }
  fWghtCalc.clear();
}
//____________________________________________________________________________
void GReWeight::Print()
{
  vector<genie::rew::GSyst_t> syst_vec = this->Systematics().AllIncluded();
  int vec_size = syst_vec.size();

  LOG("ReW", pNOTICE) << "Current set of systematic params:";	
  for(int i = 0 ; i < vec_size ; i ++){
     LOG("ReW", pNOTICE) 
        << " --o "  << GSyst::AsString(syst_vec[i])
        << " is set at " << this->Systematics().Info(syst_vec[i])->CurValue;
  }		       	        

  double chi2val = this->CalcChisq();

  LOG("ReW", pNOTICE) << "Chisq_{penalty} = " << chi2val;
}
//____________________________________________________________________________


