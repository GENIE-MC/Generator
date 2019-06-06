//____________________________________________________________________________
/*!

\class    genie::Spline

\brief    A numeric analysis tool class for interpolating 1-D functions.

          Uses ROOT's TSpline3 for the actual interpolation and can retrieve
          function (x,y(x)) pairs from an XML file, a flat ascii file, a
          TNtuple, a TTree or an SQL database.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 04, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _SPLINE_H_
#define _SPLINE_H_

#include <string>
#include <fstream>
#include <ostream>

#include <TObject.h>
#include <TSpline.h>

class TNtupleD;
class TTree;
class TSQLServer;
class TGraph;

using std::string;
using std::ostream;
using std::ofstream;

namespace genie {

class Spline;
ostream & operator << (ostream & stream, const Spline & spl);

class Spline : public TObject {

public:
  using TObject::Print; // suppress clang 'hides overloaded virtual function [-Woverloaded-virtual]' warnings

  // ctors & dtor
  Spline();
  Spline(string filename, string xtag="", string ytag="", bool is_xml = false);
  Spline(TNtupleD * ntuple, string xy, string cut="");
  Spline(TTree * tree,     string xy, string cut="");
  Spline(TSQLServer * db, string query);
  Spline(int nentries, double x[], double y[]);
  Spline(int nentries, float  x[], float  y[]);
  Spline(const Spline & spline);
  Spline(const TSpline3 & spline, int nknots);
  virtual ~Spline();

  // Load the Spline from XML, flat ASCII, ROOT ntuple/tree/tspline3, or SQL DB
  bool   LoadFromXmlFile    (string filename, string xtag, string ytag);
  bool   LoadFromAsciiFile  (string filename);
  bool   LoadFromNtuple     (TNtupleD * nt, string xy, string cut = "");
  bool   LoadFromTree       (TTree *    tr, string xy, string cut = "");
  bool   LoadFromDBase      (TSQLServer * db,  string query);
  bool   LoadFromTSpline3   (const TSpline3 & spline, int nknots);

  // Get xmin,xmax,nknots, check x variable against valid range and evaluate spline
  int    NKnots             (void) const {return fNKnots;}
  void   GetKnot            (int iknot, double & x, double & y) const;
  double GetKnotX           (int iknot) const;
  double GetKnotY           (int iknot) const;
  double XMin               (void) const {return fXMin;  }
  double XMax               (void) const {return fXMax;  }
  double YMax               (void) const {return fYMax;  }
  double Evaluate           (double x) const;
  bool   IsWithinValidRange (double x) const;

  void   SetName (string name) { fName = name; }
  string Name (void) const     { return fName; }

  void   YCanBeNegative(bool tf) { fYCanBeNegative = tf; }

  // Save the Spline in XML, flat ASCII or ROOT format
  void   SaveAsXml (string filename, string xtag, string ytag, string name="") const;
  void   SaveAsXml (ofstream & str,  string xtag, string ytag, string name="") const;
  void   SaveAsText(string filename, string format="%10.6f\t%10.6f") const;
  void   SaveAsROOT(string filename, string name="", bool recreate=false) const;

  // Export Spline as TGraph or TSpline3
  TGraph *   GetAsTGraph  (int np = 500, bool xscaling = false,
                           bool inlog=false, double fx=1., double fy=1.) const;
  TSpline3 * GetAsTSpline (void) const { return fInterpolator; }

  // Knot manipulation methods in additions to the TSpline3 ones
  void FindClosestKnot(double x, double & xknot, double & yknot, Option_t * opt="-+") const;
  bool ClosestKnotValueIsZero(double x, Option_t * opt="-+") const;

  // Common mathematical operations applied simultaneously on all spline knots
  void Add      (const Spline & spl, double c=1);
  void Multiply (const Spline & spl, double c=1);
  void Divide   (const Spline & spl, double c=1);
  void Add      (double a);
  void Multiply (double a);
  void Divide   (double a);

  // Print knots
  void Print(ostream & stream) const;

  // Overloaded operators
  friend ostream & operator << (ostream & stream, const Spline & spl);

private:

  // Initialize and build spline
  void InitSpline  (void);
  void ResetSpline (void);
  void BuildSpline (int nentries, double x[], double y[]);

  // Private data members
  string     fName;
  int        fNKnots;
  double     fXMin;
  double     fXMax;
  double     fYMax;
  TSpline3 * fInterpolator;
  bool       fYCanBeNegative;

ClassDef(Spline,1)
};

}

#endif
