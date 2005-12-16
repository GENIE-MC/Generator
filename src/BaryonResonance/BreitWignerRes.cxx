//____________________________________________________________________________
/*!

\class    genie::BreitWignerRes

\brief    Concrete implementation of the BreitWignerI interface:
          Simple Breit-Wigner distribution with no L-dependent width.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  November 22, 2004

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
  this->LoadSubAlg();
}
//____________________________________________________________________________
void BreitWignerRes::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadSubAlg();
}
//____________________________________________________________________________
void BreitWignerRes::LoadSubAlg(void)
{
// Load the "baryon resonance table" sub-algorithm specified at the algorithm
// configuration

  fBaryonResDataSet =
         dynamic_cast<const BaryonResDataSetI *> (this->SubAlg(
                              "baryon-res-alg-name", "baryon-res-param-set"));
}
//____________________________________________________________________________

