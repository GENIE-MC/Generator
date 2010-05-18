//____________________________________________________________________________
/*
 Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
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
#include "ReWeight/GReWeightNuXSecCCQE.h"
#include "ReWeight/GReWeightAGKY.h"
#include "ReWeight/GReWeightFZone.h"
#include "ReWeight/GReWeightINuke.h"
#include "ReWeight/GReWeightFGM.h"
#include "ReWeight/GReWeightResonanceDecay.h"

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

}
//____________________________________________________________________________
GReWeight::~GReWeight()
{
  this->CleanUp();
}
//____________________________________________________________________________
void GReWeight::AdoptWghtCalc(GReWeightI* wcalc)
{
  if(!wcalc) return;
  fWghtCalc.push_back(wcalc);
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

  vector<GReWeightI *>::iterator wcalc_iter = fWghtCalc.begin();
  for( ; wcalc_iter != fWghtCalc.end(); ++wcalc_iter) {

      GReWeightI * wcalc = *wcalc_iter;

      vector<genie::rew::GSyst_t>::const_iterator parm_iter = svec.begin();
      for( ; parm_iter != svec.end(); ++parm_iter) {
          GSyst_t syst = *parm_iter;
          double val = fSystSet.CurValue(syst);
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
  vector<GReWeightI *>::iterator it = fWghtCalc.begin();
  for( ; it != fWghtCalc.end(); ++it) {
    GReWeightI * wcalc = *it;
    double w = wcalc->CalcWeight(event); 
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

  vector<GReWeightI *>::iterator it = fWghtCalc.begin();
  for( ; it != fWghtCalc.end(); ++it) {
    GReWeightI * wcalc = *it;
    double chisq = wcalc->CalcChisq(); 
    tot_chisq *= chisq;
  }
  return tot_chisq;
}
//____________________________________________________________________________
void GReWeight::CleanUp(void)
{
  vector<GReWeightI *>::iterator it = fWghtCalc.begin();
  for( ; it != fWghtCalc.end(); ++it) {
    GReWeightI * rw = *it;
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


