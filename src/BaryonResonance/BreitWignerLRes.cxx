//____________________________________________________________________________
/*!

\class    genie::BreitWignerLRes

\brief    Concrete implementation of the BreitWignerI interface:
          A realistic Breit-Wigner distribution with L-dependent width.

          It is similar with the breit_wigner_L function but rather than
          specifying the Breit-Wigner parameters directly, you specify a
          resonance name and the concrete implementation of BaryonResDataSetI
          to be looked up for extracting those parameters.

          Pre-configured instances can be obtained from the AlgFactory

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  November 22, 2004

*/
//____________________________________________________________________________

#include "Algorithm/AlgFactory.h"
#include "BaryonResonance/BreitWignerLRes.h"
#include "BaryonResonance/BaryonResDataSetI.h"
#include "BaryonResonance/BaryonResParams.h"
#include "BaryonResonance/BaryonResUtils.h"
#include "BaryonResonance/breit_wigner_func.h"
#include "Messenger/Messenger.h"

using namespace genie;

//______________________________________________________________________
BreitWignerLRes::BreitWignerLRes() : BreitWignerI() 
{
  fName = "genie::BreitWignerLRes";
} 
//______________________________________________________________________
BreitWignerLRes::BreitWignerLRes(const char * param_set) :
BreitWignerI(param_set)
{ 
  fName = "genie::BreitWignerLRes";
}
//______________________________________________________________________
BreitWignerLRes::~BreitWignerLRes()
{

} 
//______________________________________________________________________
double BreitWignerLRes::Eval(double W) const
{
  //-- get the specified baryon resonance table

 assert(fConfig->Exists("baryon-res-alg-name") &&
                               fConfig->Exists("baryon-res-param-set"));

  string alg_name  = fConfig->GetString("baryon-res-alg-name");
  string param_set = fConfig->GetString("baryon-res-param-set");

  AlgFactory * algf = AlgFactory::Instance();
 
  const Algorithm * algbase = algf->GetAlgorithm(alg_name, param_set);

  const BaryonResDataSetI * dataset =
                      dynamic_cast<const BaryonResDataSetI *> (algbase);


  //-- instantiate a BaryonResParams object & set the table to lookup

  BaryonResParams res_params;

  res_params.SetDataSet(dataset);

    
  //-- get the specified resonance from the configuration
   
  Resonance_t resonance;

  if ( fConfig->Exists("resonance") ) {

      resonance = res_utils::FromString(
                           fConfig->GetString("resonance").c_str() );
                              
  } else if ( fConfig->Exists("resonance-id") ) {

      resonance = (Resonance_t) fConfig->GetInt("resonance-id");

  } else {

      LOG("BreitWigner", pFATAL) << "Unspecified resonance";

      assert(false);
  }  

   
  //-- retrieve data for the input resonance

  res_params.RetrieveData(resonance);

              
  //-- get mass, width, norm and orbital angular momentum 
  
  Registry config;
    
  config.Set("Res-Orb-Angular-Mom", res_params.OrbitalAngularMom() );
  config.Set("Res-Mass",            res_params.Mass()              );
  config.Set("Res-Width",           res_params.Width()             );
  config.Set("Breit-Wigner-Norm",   res_params.BreitWignerNorm()   );


  //-- run the actual Breit-Wigner resonance code
  
  double bw = breit_wigner_L(W, config);
  
  return bw;
}
//______________________________________________________________________

