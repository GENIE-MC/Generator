//____________________________________________________________________________
/*!

\class    genie::Spline

\brief    A numeric analysis tool class for interpolating 1-D functions.

          Uses ROOT's TSpline3 for the actual interpolation and can retrieve
          function (x,y(x)) pairs from an XML file, a flat ascii file, a
          TNtuple, a TTree or an SQL database.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 04, 2004

*/
//____________________________________________________________________________

#ifndef _SPLINE_H_
#define _SPLINE_H_

#include <string>
#include <fstream>

#include <TObject.h>
#include <TSpline.h>

class TNtuple;
class TTree;
class TSQLServer;
class TGraph;

using std::string;
using std::ofstream;

namespace genie {

class Spline : public TObject {

public:

  Spline();
  Spline(string filename, string xtag="", string ytag="", bool is_xml = false);
  Spline(TNtuple * ntuple, string xy, string cut="");
  Spline(TTree * tree,     string xy, string cut="");
  Spline(TSQLServer * db, string query);
  Spline(int nentries, double x[], double y[]);
  Spline(int nentries, float  x[], float  y[]);
  Spline(const Spline & spline);
  Spline(const TSpline3 & spline, int nknots);
  virtual ~Spline();

  // load the Spline from XML, flat ASCII, ROOT ntuple/tree/tspline3, or SQL DB
  bool   LoadFromXmlFile    (string filename, string xtag, string ytag);
  bool   LoadFromAsciiFile  (string filename);
  bool   LoadFromNtuple     (TNtuple * nt, string xy, string cut = "");
  bool   LoadFromTree       (TTree *   tr, string xy, string cut = "");
  bool   LoadFromDBase      (TSQLServer * db,  string query);
  bool   LoadFromTSpline3   (const TSpline3 & spline, int nknots);
  
  // get xmin,xmax,nknots, check x variable against valid range and evaluate spline
  int    NKnots(void) const {return fNKnots;}
  double XMin  (void) const {return fXMin;  }
  double XMax  (void) const {return fXMax;  }
  double Evaluate           (double x) const;
  bool   IsWithinValidRange (double x) const;

  // save the Spline in XML, flat ASCII or ROOT format
  void   SaveAsXml (string filename, string xtag, string ytag, string name) const;
  void   SaveAsXml (ofstream & str,  string xtag, string ytag,
                                          string name, bool insert = false) const;
  void   SaveAsText(string filename, string format="%10.6f\t%10.6f") const;
  void   SaveAsROOT(string filename, string name, bool recreate=false) const;

  // export Spline as TGraph or TSpline3
  TGraph *   GetAsTGraph  (int npoints = 100, bool scale_with_x = false) const;
  TSpline3 * GetAsTSpline (void) const { return fInterpolator; }

private:

  // initialize and build spline
  void InitSpline  (void);
  void BuildSpline (int nentries, double x[], double y[]);

  int        fNKnots;  
  double     fXMin;
  double     fXMax;
  TSpline3 * fInterpolator;

ClassDef(Spline,1)
};

}

#endif
