//____________________________________________________________________________
/*!

\class   genie::rew::GReWeightINukeHelper

\brief   Helper class for cross section model reweighting

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created October 22, 2005

*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_INTRANUKE_FATEPARAMS_H_
#define _G_REWEIGHT_INTRANUKE_FATEPARAMS_H_

#include <map>

#include "GSyst.h"
#include "HadronTransport/INukeHadroFates.h"

using std::map;

class TLorentzVector;

namespace genie {
namespace rew   {

//namespace inukehelper{
  typedef enum EHadronType {
    kNull = 0,
    kPionInNucleus,
    kNuclInNucleus
  } EHadronType_t;
//}

class GReWeightINukeFateParams {

public :
  GReWeightINukeFateParams(int=5, int=1);
 ~GReWeightINukeFateParams();

  void   Reconfigure();
  void   SetSystematic(EGSyst syst, double val);
  double TwkDialValue(EGSyst, double kinE);
  double GetAverageChi2(){ return fAverageChi2;}

private:

  void   Initialize(int, int);
  void   Reset();
  void   SetFateSyst(EGSyst syst, bool is_cushion, double val);
  double GetFateXSecTwk(EGSyst syst, double kinE, double twk_dial_val);
  bool   CheckFatesUnity(int n_points = 400);
  double GetFatesChi2(double kinE);  
  void   AddCushionTerms();
  bool   CheckSystLists();
  bool   CheckSystType(EGSyst syst);
  EHadronType_t GetHadronType(EGSyst syst);
  EGSyst GetRelevantSyst(EINukeFateHA_t hadron_fate);

  map<EGSyst, double>           fSystListMap; ///< List of systematics included. 
  map<EGSyst, double>::iterator fSystListMap_iter; ///< List of systematics included. 
  map<EGSyst, bool>             fIsCushionMap;
  map<EGSyst, bool>::iterator   fIsCushionMap_iter;
//  
  int    fNSysts;
  int    fNCushionTerms;
  bool   fArePionSysts;
  bool   fAreNuclSysts;
  double fAverageChi2;
  EHadronType_t fHadronType;


};

}      // rew   namespace
}      // genie namespace

#endif // _G_REWEIGHT_INTRANUKE_FATEPARAMS_H_
