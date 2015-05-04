//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - July 01, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <cassert>

#include <TMath.h>

#include "Messenger/Messenger.h"
#include "Numerical/UnifGridDimension.h"

using std::endl;
using namespace genie;

//____________________________________________________________________________
namespace genie
{
  ostream & operator << (ostream & stream, const UnifGridDimension & gdim)
  {
     gdim.Print(stream);
     return stream;
  }
}
//____________________________________________________________________________
UnifGridDimension::UnifGridDimension(
              GridSpacing_t sp, double min, double max, unsigned int n) :
GridDimension()
{
  this->Set(sp, min, max, n);
}
//____________________________________________________________________________
UnifGridDimension::UnifGridDimension(const UnifGridDimension & gd) :
GridDimension()
{
  this->Set(gd.fSpacing, gd.fMin, gd.fMax, gd.fNPoints);
}
//____________________________________________________________________________
UnifGridDimension::~UnifGridDimension()
{

}
//____________________________________________________________________________
double UnifGridDimension::Coord(unsigned int ipoint) const
{
  assert(ipoint < fNPoints);

  switch(fSpacing) {
  case (kGSpLinear):
          return fLocMin + fStep * ipoint;
          break;
  case (kGSpLoge):
          return TMath::Exp(fLocMin + fStep * ipoint);
          break;
  case (kGSpLog10):
          return TMath::Power(10, fLocMin + fStep * ipoint);
          break;
  default:
     LOG("GridDim", pERROR)
          << "Can't handle spacing = " << GridSpacing::AsString(fSpacing);
  }
  return 0;
}
//____________________________________________________________________________
unsigned int UnifGridDimension::ClosestPoint(double coord) const
{
  double       loc_coord = 0;
  unsigned int ipoint    = 0;

  switch(fSpacing) {
  case (kGSpLinear):
          loc_coord = coord;
          break;
  case (kGSpLoge):
          if(coord<=0) {
            LOG("GridDim", pERROR)
              << "Non-positive coordinate (x=" << coord
                                    << ") with uniform spacing in ln(x)";
            return ipoint;
          }
          loc_coord = TMath::Log(coord);
          break;
  case (kGSpLog10):
          if(coord<=0) {
            LOG("GridDim", pERROR)
              << "Non-positive coordinate (x=" << coord
                                 << ") with uniform spacing in log10(x)";
            return ipoint;
          }
          loc_coord = TMath::Log10(coord);
          break;
  default:
     LOG("GridDim", pERROR)
          << "Can't handle spacing = " << GridSpacing::AsString(fSpacing);
  }

  if(!fStep>0) {
     LOG("GridDim", pERROR) << "Non-positive step (=" << fStep << ")";
     return ipoint;
  }

  ipoint = (unsigned int) TMath::Nint( (loc_coord-fLocMin)/fStep );
  return ipoint;
}
//____________________________________________________________________________
void UnifGridDimension::Set(
            GridSpacing_t sp, double min, double max, unsigned int n)
{
  GridDimension::Set(sp, min, max, n);

  assert(fNPoints>1);
  assert(fSpacing==kGSpLinear || fSpacing==kGSpLoge || fSpacing==kGSpLog10);

  fStep = 0.; // init;

  switch(fSpacing) {
  case (kGSpLinear):
          fLocMin = fMin;
          fLocMax = fMax;
          break;
  case (kGSpLoge):
          if(fMin<=0 || fMax<=0) {
            LOG("GridDim", pERROR)
              << "Non-positive coordinate with uniform spacing in ln(x)";
            return;
          }
          fLocMin = TMath::Log(fMin);
          fLocMax = TMath::Log(fMax);
          break;
  case (kGSpLog10):
          if(fMin<=0 || fMax<=0) {
            LOG("GridDim", pERROR)
              << "Non-positive coordinate with uniform spacing in log10(x)";
            return;
          }
          fLocMin = TMath::Log10(fMin);
          fLocMax = TMath::Log10(fMax);
          break;
  default:
     LOG("GridDim", pERROR)
          << "Can't handle spacing = " << GridSpacing::AsString(fSpacing);
  }
  fStep = (fLocMax - fLocMin) / (fNPoints-1);
}
//____________________________________________________________________________
void UnifGridDimension::ReBin(unsigned int n)
{
  this->Set(fSpacing, fMin, fMax, n);
}
//____________________________________________________________________________
void UnifGridDimension::ReBin(GridSpacing_t sp)
{
  this->Set(sp, fMin, fMax, fNPoints);
}
//____________________________________________________________________________
void UnifGridDimension::Print(ostream & stream) const
{
  stream << "[" << fMin << ", " << fMax << "] - np: "
       << fNPoints << " - spacing: " << GridSpacing::AsString(fSpacing)
                                            << " - step = " << fStep << endl;
}
//____________________________________________________________________________
