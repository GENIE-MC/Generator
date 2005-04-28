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

#include <TNtuple.h>
#include <TTree.h>
#include <TSQLServer.h>
#include <TGraph.h>
#include <TSpline.h>

#include "Messenger/Messenger.h"
#include "Numerical/Spline.h"

using namespace genie;

//___________________________________________________________________________
Spline::Spline()
{

}
//___________________________________________________________________________
Spline::Spline(const char * filename)
{
  LOG("Spline", pDEBUG) << "Constructing spline from data in text file";

  this->InitSpline();
  this->LoadFromFile(filename);
}
//___________________________________________________________________________
Spline::Spline(TNtuple * ntuple, const char * var, const char * constraint)
{
  LOG("Spline", pDEBUG) << "Constructing spline from data in a TNtuple";

  this->InitSpline();
  this->LoadFromNtuple(ntuple, var, constraint);
}
//___________________________________________________________________________
Spline::Spline(TTree * tree, const char * var, const char * constraint)
{
  LOG("Spline", pDEBUG) << "Constructing spline from data in a TTree";

  this->InitSpline();
  this->LoadFromTree(tree, var, constraint);
}
//___________________________________________________________________________
Spline::Spline(TSQLServer * db, const char * query)
{
  LOG("Spline", pDEBUG) << "Constructing spline from data in a MySQL server";

  this->InitSpline();
  this->LoadFromDBase(db, query);
}
//___________________________________________________________________________
Spline::Spline(int nentries, double x[], double y[])
{
  LOG("Spline", pDEBUG)
                 << "Constructing spline from the arrays passed to the ctor";

  this->InitSpline();
  this->BuildSpline(nentries, x, y);
}
//___________________________________________________________________________
Spline::Spline(int nentries, float x[], float y[])
{
  LOG("Spline", pDEBUG)
                 << "Constructing spline from the arrays passed to the ctor";

  double * dblx = new double[nentries];
  double * dbly = new double[nentries];

  for(int i = 0; i < nentries; i++) {
     dblx[i] = double( x[i] );
     dbly[i] = double( y[i] );
  }
    
  this->InitSpline();
  this->BuildSpline(nentries, dblx, dbly);

  delete [] x;
  delete [] y;
}
//___________________________________________________________________________
Spline::Spline(const Spline & spline)
{
  LOG("Spline", pDEBUG) << "Copy constructor";
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
  LOG("Spline", pDEBUG) << "Retrieving data from file: " << filename;

}
//___________________________________________________________________________
void Spline::LoadFromNtuple(
                     TNtuple * nt, const char * var, const char * constraint)
{
  LOG("Spline", pDEBUG) << "Retrieving data from ntuple: " << nt->GetName();

  nt->Draw(var,"","GOFF"); // select all entries - graphics off

  int      n = nt->GetSelectedRows();
  double * x = nt->GetV1();
  double * y = nt->GetV2();

  this->BuildSpline(n, x, y);
}
//___________________________________________________________________________
void Spline::LoadFromTree(
              TTree * tree, const char * var_string, const char * constraint)
{
  LOG("Spline", pDEBUG) << "Retrieving data from tree: " << tree->GetName();

  tree->Draw(var_string,"","GOFF"); // select all entries - graphics off

  int      n = tree->GetSelectedRows();
  double * x = tree->GetV1();
  double * y = tree->GetV2();

  this->BuildSpline(n, x, y);
}
//___________________________________________________________________________
void Spline::LoadFromDBase(TSQLServer * db,  const char * query)
{
  LOG("Spline", pDEBUG) << "Retrieving data from data-base: ";
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
  
  if( this->IsWithinValidRange(x) ) return y = fInterpolator->Eval(x);
  
  return y;
}
//___________________________________________________________________________
TGraph * Spline::GetAsTGraph(int np, bool scale_with_x) const
{
  double xmin = fValidityRange.first;
  double xmax = fValidityRange.second;

  if(np < 2) np = 2;

  double dx = (xmax - xmin) / (np-1);

  double * x = new double[np];
  double * y = new double[np];

  for(int i=0; i<np; i++) {

      x[i]  = xmin + i*dx;
      y[i] = this->Evaluate( x[i] );

      if (scale_with_x) y[i] /= x[i];
  }

  TGraph * graph = new TGraph(np, x, y);

  return graph;  
}
//___________________________________________________________________________
void Spline::InitSpline(void)
{
  LOG("Spline", pDEBUG) << "Initializing spline...";
  
  fValidityRange.first  = 0.0;
  fValidityRange.second = 0.0;
  
  fInterpolator = 0;
  
  LOG("Spline", pDEBUG) << "...done initializing spline";
}
//___________________________________________________________________________
void Spline::BuildSpline(int nentries, double x[], double y[])
{
  LOG("Spline", pDEBUG) << "Building spline...";

  double xmin = x[ TMath::LocMin(nentries, x) ]; // minimum x in spline
  double xmax = x[ TMath::LocMax(nentries, x) ]; // maximum x in spline

  fValidityRange.first  = xmin;
  fValidityRange.second = xmax;

  if(fInterpolator) delete fInterpolator;

  fInterpolator = new TSpline3("spl3", x, y, nentries, "0", xmin, xmax);

  LOG("Spline", pDEBUG) << "...done building spline";  
}
//___________________________________________________________________________

  
