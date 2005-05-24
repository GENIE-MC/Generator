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

#include <ostream>
#include <map>
#include <string>

using std::map;
using std::pair;
using std::string;
using std::ostream;

namespace genie {

class XSecAlgorithmI;
class Spline;

class XSecSplineList {

public:

  static XSecSplineList * Instance();

  //-- query the existence, access or create a spline

  bool           SplineExists (const XSecAlgorithmI * alg, const Interaction * i) const;
  const Spline * GetSpline    (const XSecAlgorithmI * alg, const Interaction * i) const;
  void           CreateSpline (const XSecAlgorithmI * alg, const Interaction * i,
                                   int nknots = -1, double Emin = -1, double Emax = -1);
  //-- set XSecSplineList options

  void   SetLogE   (bool   on); ///< set opt to build splines as f(E) or as f(logE)
  void   SetNKnots (int    nk); ///< set default number of knots for building the spline
  void   SetMinE   (double Ev); ///< set default minimum energy for xsec splines
  void   SetMaxE   (double Ev); ///< set default maximum energy for xsec splines
  void   SetExtrap (double Ev); 
  
  //-- read XSecSplineList options

  bool   UseLogE     (void) const { return fUseLogE;     }
  bool   Extrapolate (void) const { return fExtrapolate; }
  int    NKnots      (void) const { return fNKnots;      }
  double Emin        (void) const { return fEmin;        }
  double Emax        (void) const { return fEmax;        }
  double EExtrap     (void) const { return fEExtrap;     }

  //-- save to / load from file

  void   SaveSplineList (string filename  );
  void   LoadSplineList (bool keep = false);
  
  //-- print available splines
  
  void   Print (ostream & stream) const;
  friend ostream & operator << (ostream & stream, const XSecSplineList & xsl);
  
private:

  XSecSplineList(); 
  XSecSplineList(const XSecSplineList & spline_list);
  virtual ~XSecSplineList();

  string BuildSplineKey(const XSecAlgorithmI * alg, const Interaction * i) const;

  static XSecSplineList * fInstance;

  bool   fUseLogE;
  bool   fExtrapolate;
  int    fNKnots;
  double fEmin;
  double fEmax;
  double fEExtrap;
  
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
