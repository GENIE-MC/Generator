//____________________________________________________________________________
/*!

\class    genie::XSecSplineList

\brief    List of cross section vs energy splines

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 12, 2005
*/
//____________________________________________________________________________

#ifndef _XSEC_SPLINE_LIST_H_
#define _XSEC_SPLINE_LIST_H_

#include <map>
#include <string>

using std::map;
using std::pair;
using std::string;

namespace genie {

class XSecAlgorithmI;
class Spline;

class XSecSplineList {

public:

  static XSecSplineList * Instance();

  bool           SplineExists (const XSecAlgorithmI * alg, const Interaction * i) const;
  const Spline * GetSpline    (const XSecAlgorithmI * alg, const Interaction * i) const;
  void           CreateSpline (const XSecAlgorithmI * alg, const Interaction * i, double Emin, double Emax);
  
private:

  XSecSplineList(); 
  XSecSplineList(const XSecSplineList & spline_list);
  virtual ~XSecSplineList();

  string BuildSplineKey(const XSecAlgorithmI * alg, const Interaction * i) const;

  static XSecSplineList * fInstance;

  map<string, Spline *> fSplineMap; ///< xsec_alg_name/param_set/interaction -> Spline
  
  struct Cleaner {
      void DummyMethodAndSilentCompiler() { }
      ~Cleaner() {
         if (XSecSplineList::fInstance !=0) {
            delete XSecSplineList::fInstance;
            XSecSplineList::fInstance = 0;
         }
      }
  };

  friend struct Cleaner;
};

}      // genie namespace

#endif // _XSEC_SPLINE_LIST_H_
