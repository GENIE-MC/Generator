//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - November 26, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Jun 25, 2008 - CA
   Partial re-write to fix a serious memory leak. Holding x,y values in a map
   rather than a circular ntuple.

*/
//____________________________________________________________________________

#include "Framework/Utils/CacheBranchFx.h"

using namespace genie;

ClassImp(CacheBranchFx);

//____________________________________________________________________________
namespace genie
{
  ostream & operator << (ostream & stream, const CacheBranchFx & cbntp)
  {
     cbntp.Print(stream);
     return stream;
  }
}
//____________________________________________________________________________
CacheBranchFx::CacheBranchFx(void) :
CacheBranchI()
{
  this->Init();
}
//____________________________________________________________________________
CacheBranchFx::CacheBranchFx(string name) :
CacheBranchI()
{
  this->Init();
  fName = name;
}
//____________________________________________________________________________
CacheBranchFx::~CacheBranchFx()
{
  this->CleanUp();
}
//____________________________________________________________________________
void CacheBranchFx::Init(void)
{
  fName   = "";
  fSpline = 0;
}
//____________________________________________________________________________
void CacheBranchFx::CleanUp(void)
{
  if(fSpline) delete fSpline;
  fFx.clear();
}
//____________________________________________________________________________
void CacheBranchFx::Reset(void)
{
  this->CleanUp();
  this->Init();
}
//____________________________________________________________________________
void CacheBranchFx::AddValues(double x, double y)
{
  fFx.insert(map<double,double>::value_type(x,y));
}
//____________________________________________________________________________
void CacheBranchFx::CreateSpline(void)
{
  int n = fFx.size();
  double * x = new double[n];
  double * y = new double[n];

  int i=0;
  map<double,double>::const_iterator iter = fFx.begin();
  for( ; iter !=fFx.end(); ++iter) {
    x[i] = iter->first;
    y[i] = iter->second;
    i++;
  }

  if(fSpline) delete fSpline;
  fSpline = new Spline(n,x,y);

  delete [] x;
  delete [] y;
}
//____________________________________________________________________________
void CacheBranchFx::Print(ostream & stream) const
{
  stream << "type: [CacheBranchFx]  - nentries: " << fFx.size()
           << " / spline: " << ((fSpline) ? "built" : "null");
}
//____________________________________________________________________________
double CacheBranchFx::operator () (double x) const
{
  if(!fSpline) return 0;
  else         return fSpline->Evaluate(x);
}
//____________________________________________________________________________
