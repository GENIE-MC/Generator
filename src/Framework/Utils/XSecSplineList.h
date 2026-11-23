//____________________________________________________________________________
/*!

\class    genie::XSecSplineList

\brief    List of cross section vs energy splines

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  May 12, 2005

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _XSEC_SPLINE_LIST_H_
#define _XSEC_SPLINE_LIST_H_

#include <ostream>
#include <map>
#include <memory>
#include <vector>
#include <string>

#include "Framework/Conventions/XmlParserStatus.h"

using std::map;
using std::pair;
using std::vector;
using std::string;
using std::ostream;

namespace genie {

class XSecAlgorithmI;
class Interaction;
class Spline;

class XSecSplineList;
ostream & operator << (ostream & stream, const XSecSplineList & xsl);

class XSecSplineList {

public:

  static XSecSplineList * Instance();

  // Save/load to/from XML file
  void               SaveAsXml   (const string & filename, bool save_init = true) const;
  XmlParserStatus_t  LoadFromXml (const string & filename, bool keep = false);
  XmlParserStatus_t  LoadFromXml (const string & filename, bool keep,
                                   bool require_consistent_em_cuts);

  // Print available splines
  void   Print (ostream & stream) const;
  friend ostream & operator << (ostream & stream, const XSecSplineList & xsl);

  // Set and query current tune.
  // An XSecSplineList can keep splines for numerous tunes and pick the appropriate
  // one for each process, as instructed.
  void   SetCurrentTune (const string & tune) { fCurrentTune = tune; }
  string CurrentTune    (void) const  { return fCurrentTune; }
  bool   HasSplineFromTune( const string & tune ) const { return fSplineMap.count(tune) > 0 ; }
  // A zero value means no configured cut. A missing stamp returns false.
  bool   GetSplineQ2MinKinematics(const string & key, double & Q2min) const;
  void   EnableQ2CutValidation(void) { fValidateQ2Cuts = true; }
  bool   Q2CutValidationEnabled(void) const { return fValidateQ2Cuts; }

  // Query the existence, access or create a spline
  // The results of the following methods depend on the current tune setting
  bool           SplineExists (const XSecAlgorithmI * alg, const Interaction * i) const;
  bool           SplineExists (string spline_key) const;
  const Spline * GetSpline    (const XSecAlgorithmI * alg, const Interaction * i) const;
  const Spline * GetSpline    (string spline_key) const;
  void           CreateSpline (const XSecAlgorithmI * alg, const Interaction * i,
                               int nknots = -1, double e_min = -1, double e_max = -1);
  int  NSplines (void) const;
  bool IsEmpty  (void) const;

  // Methods for building / getting keys
  // The results of the following methods depend on the current tune setting
  string BuildSplineKey(const XSecAlgorithmI * alg, const Interaction * i) const;
  const vector<string> * GetSplineKeys(void) const;


  // XSecSplineList options
  void   SetLogE   (bool   on); ///< set opt to build splines as f(E) or as f(logE)
  void   SetNKnots (int    nk); ///< set default number of knots for building the spline
  void   SetMinE   (double Ev); ///< set default minimum energy for xsec splines
  void   SetMaxE   (double Ev); ///< set default maximum energy for xsec splines
  bool   UseLogE   (void) const { return fUseLogE;  }
  int    NKnots    (void) const { return fNKnots;   }
  double Emin      (void) const { return fEmin;     }
  double Emax      (void) const { return fEmax;     }

private:

  XSecSplineList();
  XSecSplineList(const XSecSplineList & spline_list);
  virtual ~XSecSplineList();

  static XSecSplineList * fInstance;

  bool   fUseLogE;
  int    fNKnots;
  double fEmin;
  double fEmax;

  string fCurrentTune; ///< The `active' tune, out the many that can co-exist

  struct KinematicsMetadata {
    bool   known = false;
    double Q2_min = 0.;
    string Q2_source;
  };

  struct SplineEntry {
    std::unique_ptr<Spline> spline;
    KinematicsMetadata cuts;
    bool loaded = false;
    mutable bool validated = false;
  };

  struct EMCutOrigin {
    double Q2_min;
    string filename;
  };

  using SplineMap = map<string, SplineEntry>;
  map<string, SplineMap> fSplineMap; ///< tune -> spline key -> spline and metadata
  map<string, EMCutOrigin> fMergedEMCuts; ///< First known EM cut per tune in a checked merge
  bool fValidateQ2Cuts = false;

  const SplineEntry * FindSpline(const string & key) const;
  void ValidateQ2Cut(const string & key, const Interaction * interaction,
                     const SplineEntry & entry) const;
  bool CheckEMCutForMerge(const string & tune, const string & key,
                         const KinematicsMetadata & cuts, const string & filename);

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
