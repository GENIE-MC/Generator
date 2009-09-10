//____________________________________________________________________________
/*!

\class   genie::rew::GReWeightINukeHelper

\brief   Helper class for cross section model reweighting

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created October 22, 2005

*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_INTRANUKE_HELPER_H_
#define _G_REWEIGHT_INTRANUKE_HELPER_H_

class TLorentzVector;

namespace genie {
namespace rew   {

class GReWeightINukeHelper {

public :
  GReWeightINukeHelper();
 ~GReWeightINukeHelper();


private:

  void Initialize(void);
};

}      // rew   namespace
}      // genie namespace

#endif // _G_REWEIGHT_INTRANUKE_HELPER_H_
