//_____________________________________________________________________________
/*!

\class    genie::nuvld::NeuGenFitParams

\brief    Holds the NeuGEN parameters that can be fitted and it is set by the
          NeuGenFitParamsDialog. It has a rather flat structure since this
          is more convenient for use with ROOT's fitting machinery. There is
          a correspondence with NeuGenConfig class - see this for more info
          on the fitted parameters.

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 28, 2005
*/
//_____________________________________________________________________________

#ifndef _NEUGEN_FIT_PARAM_H_
#define _NEUGEN_FIT_PARAM_H_

#include <string>

using std::string;

namespace genie {
namespace nuvld {

const int kNNGFitParams = 22;

typedef enum ENeuGenFitParam {

  kNgfMaQel = 0,
  kNgfMaRes,
  kNgfMaCoh,
  kNgfQelFa0,
  kNgfQelEta,
  kNgfResOmega,
  kNgfResZ,
  kNgfCohR0,
  kNgfCohREI,
  kNgfKnoB,
  kNgfKnoAvp,
  kNgfKnoAvn,
  kNgfKnoAvbp,
  kNgfKnoAvbn,
  kNgfDisResM2vp,
  kNgfDisResM3vp,
  kNgfDisResM2vn,
  kNgfDisResM3vn,
  kNgfDisResM2vbp,
  kNgfDisResM3vbp,
  kNgfDisResM2vbn,
  kNgfDisResM3vbn

} NeuGenFitParam_t;


class NeuGenFitParams {

public:

  friend class NeuGenFitParamsDialog;
  
  NeuGenFitParams();
  NeuGenFitParams(const NeuGenFitParams & ngfp);
  ~NeuGenFitParams();

  int    NFittedParams (void)       const;
  bool   IsFitted      (int iparam) const;
  double RangeMin      (int iparam) const;
  double RangeMax      (int iparam) const;
  double Step          (int iparam) const;
  string ParamAsString (int iparam) const;
  
private:

  void   Init    (void);
  bool   InRange (int iparam) const;
  
  bool   * fIsFitted;
  double * fRangeMin;
  double * fRangeMax;
  double * fStep;

};

} // nuvld namespace
} // genie namespace

#endif 

