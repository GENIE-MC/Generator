//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "validation/NuXSec/NuXSecComparison.h"

using namespace genie;
using namespace genie::mc_vs_data;

//____________________________________________________________________________
NuXSecComparison::NuXSecComparison(
    string id, string label, string dataset_keys, NuXSecFunc * xsec_func,
    double Emin,  double Emax, 
    bool in_logx, bool in_logy,  bool scale_with_E
) :
fID          (id),
fLabel       (label),
fDataSetKeys (dataset_keys),
fXSecFunc    (xsec_func),   
fEmin        (Emin),
fEmax        (Emax),
fInLogX      (in_logx),
fInLogY      (in_logy),
fScaleWithE  (scale_with_E)
{

}
//____________________________________________________________________________
NuXSecComparison::~NuXSecComparison()
{

}
//____________________________________________________________________________
bool NuXSecComparison::SameMCPrediction(const NuXSecComparison * comp)
{
  if(!comp) return false;

  if(!comp->XSecFunc()) return false;
  if(!this->XSecFunc()) return false;

  bool same_func = 
     (this->XSecFunc()->Name() == comp->XSecFunc()->Name());
  if(!same_func) return false;

  bool wider_range =
     (comp->Emin() <= this->Emin() && comp->Emax() >= this->Emax());
  if(!wider_range) return false;

  return true;
}
//____________________________________________________________________________

