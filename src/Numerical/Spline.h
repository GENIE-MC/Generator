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

#include <map>

#include <TSpline.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TSQLServer.h>

using std::pair;

namespace genie {

class Spline {

public:

  Spline();
  Spline(const char * filename);
  Spline(TNtuple & ntuple, const char * xy, const char * constraint = 0);
  Spline(TTree & tree,     const char * xy, const char * constraint = 0);
  Spline(TSQLServer * db, const char * query);
  Spline(int nentries, double x[], double y[]);
  Spline(const Spline & spline);
  Spline(const TSpline3 & spline);
  virtual ~Spline();

  void LoadFromFile     (const char * filename);
  void LoadFromNtuple   (TNtuple & ntuple, const char * xy, const char * constraint = 0);
  void LoadFromTTree    (TTree & tree,     const char * xy, const char * constraint = 0);
  void LoadFromSQLDBase (TSQLServer * db,  const char * query);
    
  bool   IsWithinValidRange (double x) const;
  double Evaluate           (double x) const;

private:

  void InitSpline(void);
  void BuildSpline(int nentries, double x[], double y[]);
  
  TSpline3 *            fInterpolator;
  pair<double, double>  fValidityRange; // (x_min, x_max)
  
};

}

#endif
