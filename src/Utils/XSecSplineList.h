//____________________________________________________________________________
/*!

\class    genie::XSecSplineList

\brief    List of cross section vs energy splines

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 12, 2005

\cpright  Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _XSEC_SPLINE_LIST_H_
#define _XSEC_SPLINE_LIST_H_

#include <ostream>
#include <map>
#include <vector>
#include <string>

#include "Conventions/XmlParserStatus.h"

using std::map;
using std::pair;
using std::vector;
using std::string;
using std::ostream;

namespace genie {

class XSecAlgorithmI;
class Interaction;
class Spline;

class XSecSplineList {

public:

  static XSecSplineList * Instance();

  // Query the existence, access or create a spline
  bool           SplineExists (const XSecAlgorithmI * alg, const Interaction * i) const;
  bool           SplineExists (string spline_key) const;
  const Spline * GetSpline    (const XSecAlgorithmI * alg, const Interaction * i) const;
  const Spline * GetSpline    (string spline_key) const;
  void           CreateSpline (const XSecAlgorithmI * alg, const Interaction * i,
                                   int nknots = -1, double Emin = -1, double Emax = -1);

  int  NSplines (void) const { return fSplineMap.size();        }
  bool IsEmpty  (void) const { return (fSplineMap.size() == 0); }

  // Set XSecSplineList options
  void   SetLogE   (bool   on); ///< set opt to build splines as f(E) or as f(logE)
  void   SetNKnots (int    nk); ///< set default number of knots for building the spline
  void   SetMinE   (double Ev); ///< set default minimum energy for xsec splines
  void   SetMaxE   (double Ev); ///< set default maximum energy for xsec splines

  // Read XSecSplineList options
  bool   UseLogE     (void) const { return fUseLogE;     }
  int    NKnots      (void) const { return fNKnots;      }
  double Emin        (void) const { return fEmin;        }
  double Emax        (void) const { return fEmax;        }

  // Save/load to/from XML file
  void               SaveAsXml   (string filename) const;
  XmlParserStatus_t  LoadFromXml (string filename, bool keep = false);

  // Autosave/autoload
  bool AutoLoad (void);
  void AutoSave (void);

  // Methods for building / getting keys
  string BuildSplineKey(const XSecAlgorithmI * alg, const Interaction * i) const;
  const vector<string> * GetSplineKeys(void) const;

  // Print available splines
  void   Print (ostream & stream) const;
  friend ostream & operator << (ostream & stream, const XSecSplineList & xsl);


private:

  XSecSplineList();
  XSecSplineList(const XSecSplineList & spline_list);
  virtual ~XSecSplineList();

  static XSecSplineList * fInstance;

  bool   fUseLogE;
  int    fNKnots;
  double fEmin;
  double fEmax;

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
