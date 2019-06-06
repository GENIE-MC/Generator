//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "T2KEvGenMetaData.h"

using std::endl;

ClassImp(genie::utils::T2KEvGenMetaData);

//____________________________________________________________________________
namespace genie
{
 namespace utils
 {
  ostream & operator << (ostream & stream, const T2KEvGenMetaData & md)
  {
     md.Print(stream);
     return stream;
  }
 }
}
//____________________________________________________________________________
void genie::utils::T2KEvGenMetaData::Print(ostream & stream) const
{
  stream << endl;

  if(this->jnubeam_version.size() > 0) {
    stream << "jnubeam version = " << this->jnubeam_version << endl;          
  }
  if(this->jnubeam_file.size() > 0) {
    stream << "flux ntuple filename = " << this->jnubeam_file << endl;
  }
  if(this->detector_location.size() > 0) {
    stream << "detector location = " << this->detector_location << endl;          
  }
  if(this->geom_file.size() > 0) {
    stream << "detector geometry file = " << this->geom_file << endl;          
  }

  map<int, TH1D*> fluxhists = this->flux_hists;
  if(fluxhists.size()>0) {
    stream << "found flux histograms:" << endl;
  }
  map<int, TH1D*>::const_iterator hist_iter = fluxhists.begin();
  while(hist_iter != fluxhists.end()){
    TH1D * curr_hist = (TH1D*) hist_iter->second;
    if(curr_hist){
      stream << " - name = "  << curr_hist->GetName() 
             << " (entries: " << curr_hist->GetEntries() 
             << ", mean: "    << curr_hist->GetMean() 
             << ") --> neutrino pdg = " << hist_iter->first << endl;
      ++hist_iter;
    }//curr_hist
  }//hist_iter

  map<int, double> targetmix = this->target_mix;
  if(targetmix.size()>0) {
    stream << "found target mix:" << endl;
  }
  map<int, double>::const_iterator target_iter = targetmix.begin();
  while(target_iter != targetmix.end()){
     stream << " - target pdg = " << target_iter->first 
            << ", weight fraction = " << target_iter->second << endl;
      ++target_iter;
  }
}
//____________________________________________________________________________

