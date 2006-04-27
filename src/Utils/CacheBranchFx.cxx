//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - November 26, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TNtupleD.h>

#include "Utils/CacheBranchFx.h"

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
  this->CreateNtuple(name);
}
//____________________________________________________________________________
CacheBranchFx::~CacheBranchFx()
{
  this->CleanUp();
}
//____________________________________________________________________________
void CacheBranchFx::Init(void)
{
  fNtp    = 0;
  fSpline = 0;
}
//____________________________________________________________________________
void CacheBranchFx::CleanUp(void)
{
  if(fNtp)    delete fNtp;
  if(fSpline) delete fSpline;
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
  fNtp->Fill(x,y);
}
//____________________________________________________________________________
void CacheBranchFx::CreateNtuple(string name)
{
  fNtp = new TNtupleD(name.c_str(), "f(x)", "x:y");
  fNtp->SetDirectory(0);
  fNtp->SetCircular(1600000);
}
//____________________________________________________________________________
void CacheBranchFx::CreateSpline(void)
{
  if(fSpline) delete fSpline;
  fSpline = new Spline(fNtp,"x:y");
}
//____________________________________________________________________________
void CacheBranchFx::Print(ostream & stream) const
{
  if(fNtp) {
    stream << "type: [CacheBranchFx]  - nentries: " << fNtp->GetEntries() 
           << " / spline: " << ((fSpline) ? "built" : "null");
  } else {
    stream << " *** NULL ***";
  }
}
//____________________________________________________________________________
double CacheBranchFx::operator () (double x) const
{
  if(!fSpline) return 0;
  else return fSpline->Evaluate(x);
}
//____________________________________________________________________________
