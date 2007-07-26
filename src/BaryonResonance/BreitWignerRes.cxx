//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - November 22, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "BaryonResonance/BreitWignerRes.h"
#include "BaryonResonance/BaryonResDataSetI.h"
#include "BaryonResonance/BaryonResParams.h"
#include "BaryonResonance/BaryonResUtils.h"
#include "Messenger/Messenger.h"
#include "Utils/BWFunc.h"

using namespace genie;

//______________________________________________________________________
BreitWignerRes::BreitWignerRes() :
BreitWignerI("genie::BreitWignerRes")
{

}
//______________________________________________________________________
BreitWignerRes::BreitWignerRes(string config) :
BreitWignerI("genie::BreitWignerRes", config)
{

}
//______________________________________________________________________
BreitWignerRes::~BreitWignerRes()
{

}
//______________________________________________________________________
double BreitWignerRes::Eval(Resonance_t res, double W) const
{
  //-- instantiate a BaryonResParams object & set the table to lookup
  BaryonResParams res_params;
  res_params.SetDataSet(fBaryonResDataSet);
  res_params.RetrieveData(res);

  //-- get mass, width and norm
  double mass  = res_params.Mass();
  double width = res_params.Width();
  double norm  = res_params.BreitWignerNorm();

  //-- call the actual Breit-Wigner function
  double bw = utils::bwfunc::BreitWigner(W, mass, width, norm);
  return bw;
}
//______________________________________________________________________
void BreitWignerRes::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BreitWignerRes::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BreitWignerRes::LoadConfig(void)
{
// Load the "baryon resonance table" sub-algorithm specified at the algorithm
// configuration

  fBaryonResDataSet =
     dynamic_cast<const BaryonResDataSetI *> (this->SubAlg("BaryonResData"));
  assert(fBaryonResDataSet);
}
//____________________________________________________________________________

