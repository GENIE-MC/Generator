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

#include <TObject.h>
#include <TSpline.h>

class TNtuple;
class TTree;
class TSQLServer;
class TGraph;
//class TSpline3;

namespace genie {

class Spline : public TObject {

public:

  Spline();
  Spline(const char * filename);
  Spline(TNtuple * ntuple, const char * xy, const char * constraint = 0);
  Spline(TTree * tree,     const char * xy, const char * constraint = 0);
  Spline(TSQLServer * db, const char * query);
  Spline(int nentries, double x[], double y[]);
  Spline(int nentries, float  x[], float  y[]);
  Spline(const Spline & spline);
  Spline(const TSpline3 & spline);
  virtual ~Spline();
  
  void LoadFromFile   (const char * filename);
  void LoadFromNtuple (TNtuple * nt, const char * xy, const char * constraint = 0);
  void LoadFromTree   (TTree *   tr, const char * xy, const char * constraint = 0);
  void LoadFromDBase  (TSQLServer * db,  const char * query);

  bool   IsWithinValidRange (double x) const;
  double Evaluate           (double x) const;
  
  TGraph *   GetAsTGraph  (int npoints = 100, bool scale_with_x = false) const;
  TSpline3 * GetAsTSpline (void) const { return fInterpolator; }

private:

  void InitSpline       (void);
  void BuildSpline      (int nentries, double x[], double y[]);
  
  double     fXMin;
  double     fXMax;
  TSpline3 * fInterpolator;

ClassDef(Spline,1)
};

}

#endif
