//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - July 01, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

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

