//____________________________________________________________________________
/*!

\class    genie::GridDimension

\brief    Grid dimension ABC

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  July 01, 2004

*/
//____________________________________________________________________________

#include <cassert>

#include "Numerical/GridDimension.h"

using namespace genie;

//____________________________________________________________________________
GridDimension::GridDimension()
{

}
//____________________________________________________________________________
GridDimension::GridDimension(
                  GridSpacing_t sp, double min, double max, unsigned int np) :
fSpacing(sp),
fMin(min),
fMax(max),
fNPoints(np)
{

}
//____________________________________________________________________________
GridDimension::~GridDimension()
{
}
//____________________________________________________________________________
void GridDimension::Set(
                 GridSpacing_t sp, double min, double max, unsigned int np)
{
  fSpacing = sp,
  fMin     = min;
  fMax     = max;
  fNPoints = np;
}
//____________________________________________________________________________

