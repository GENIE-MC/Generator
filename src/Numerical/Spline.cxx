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

#include "Spline.h"

using namespace genie;

//___________________________________________________________________________
Spline::Spline()
{

}
//___________________________________________________________________________
Spline::Spline(const char * filename)
{
  LoadFromFile(filename);
}
//___________________________________________________________________________
Spline::Spline(
          TNtuple & ntuple, const char * var_string, const char * constraint)
{
  LoadFromNtuple(ntuple, var_string, constraint);
}
//___________________________________________________________________________
Spline::Spline(
              TTree & tree, const char * var_string, const char * constraint)
{
  LoadFromTTree(tree, var_string, constraint);
}
//___________________________________________________________________________
Spline::Spline(TSQLServer * db, const char * query)
{
  LoadFromSQLDBase(db, query);
}
//___________________________________________________________________________
Spline::Spline(int n_entries, double x[], double y[])
{
  BuildSpline(n_entries, x, y);
}
//___________________________________________________________________________
Spline::Spline(const Spline & spline)
{

}
//___________________________________________________________________________
Spline::Spline(const TSpline3 & spline)
{

}
//___________________________________________________________________________
Spline::~Spline()
{
  if(fInterpolator) delete fInterpolator;
}
//___________________________________________________________________________
void Spline::LoadFromFile(const char * filename)
{

}
//___________________________________________________________________________
void Spline::LoadFromNtuple(
          TNtuple & ntuple, const char * var_string, const char * constraint)
{
  ntuple.Draw(var_string,"","GOFF"); // select all entries - graphics off

  int n_entries     = ntuple.GetSelectedRows();
  double * x_values = ntuple.GetV1();
  double * y_values = ntuple.GetV2();

  BuildSpline(n_entries, x_values, y_values);
}
//___________________________________________________________________________
void Spline::LoadFromTTree(
              TTree & tree, const char * var_string, const char * constraint)
{
  tree.Draw(var_string,"","GOFF"); // select all entries - graphics off

  int n_entries     = tree.GetSelectedRows();
  double * x_values = tree.GetV1();
  double * y_values = tree.GetV2();

  BuildSpline(n_entries, x_values, y_values);
}
//___________________________________________________________________________
void Spline::LoadFromSQLDBase(TSQLServer * db,  const char * query)
{

}
//___________________________________________________________________________
bool Spline::IsWithinValidRange(double x) const
{
  bool is_in_range = (fValidityRange.first < x && x < fValidityRange.second);

  return is_in_range;
}
//___________________________________________________________________________
double Spline::Evaluate(double x) const
{
  double y = -1;
  
  if( IsWithinValidRange(x) ) return y = fInterpolator->Eval(x);
  
  return y;
}
//___________________________________________________________________________
void Spline::InitSpline(void)
{
  fValidityRange.first  = 0.0;
  fValidityRange.second = 0.0;
  
  if(fInterpolator) delete fInterpolator;
  
  fInterpolator         = 0;
}
//___________________________________________________________________________
void Spline::BuildSpline(int n_entries, double x[], double y[])
{
  double x_min = x[ TMath::LocMin(n_entries, x) ]; // minimum x in spline
  double x_max = x[ TMath::LocMax(n_entries, x) ]; // maximum x in spline

  fValidityRange.first  = x_min;
  fValidityRange.second = x_max;

  if(fInterpolator) delete fInterpolator;

  fInterpolator = new TSpline3("spl3", x, y, n_entries, "0", x_min, x_max);
}
//___________________________________________________________________________

  
