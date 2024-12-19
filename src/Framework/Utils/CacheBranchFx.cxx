//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool 
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
void CacheBranchFx::CreateSpline(string type)
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
  fSpline->SetType(type);

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
