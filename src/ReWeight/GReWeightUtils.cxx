//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author:  Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 11, 2009 - CA
   Was first added in v2.5.1. Adapted from the T2K-specific version of the
   GENIE reweighting tool.

*/
//____________________________________________________________________________

#include <TMath.h>

#include "ReWeight/GReWeightUtils.h"

//____________________________________________________________________________
double genie::utils::rew::WeightTwkMeanFreePath(
   double prob_def, double prob_twk, bool interacted)
{
// Returns a weight to account for a change in hadron mean free path inside
// insidea nuclear medium.
//
// Inputs:
//   prob_def   : nominal survival probability
//   prob_twk   : survival probability for the tweaked mean free path
//   interacted : flag indicating whether the hadron interacted or escaped
//
// See utils::intranuke::ProbSurvival() for the calculation of probabilities.
//
  double w_mfp = 1.;

  if(interacted) {
     w_mfp = (1-prob_def>0) ?  (1-prob_twk)  /   (1-prob_def)  : 1;
  } else {
     w_mfp = (prob_def>0) ? prob_twk / prob_def : 1;
  }
  w_mfp = TMath::Max(0.,w_mfp);

  return w_mfp;
}
//____________________________________________________________________________
