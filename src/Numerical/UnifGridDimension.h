//____________________________________________________________________________
/*!

\class    genie::UnifGridDimension

\brief    A dimension of a multi-dimensional grid with uniform spacing is some
          metric.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  July 01, 2004

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _UNIFORM_GRID_DIMENSION_H_
#define _UNIFORM_GRID_DIMENSION_H_

#include <iostream>

#include "Numerical/GridDimension.h"

using std::ostream;

namespace genie {

class UnifGridDimension : public GridDimension
{
public:
  UnifGridDimension(GridSpacing_t sp, double min, double max, unsigned int n=2);
  UnifGridDimension(const UnifGridDimension & grdim);
  ~UnifGridDimension();

  double       Coord        (unsigned int ipoint) const;
  unsigned int ClosestPoint (double coord)        const;

  void Set (GridSpacing_t sp, double min, double max, unsigned int n);

  double Step  (void) const { return fStep; }
  void   ReBin (unsigned int n);
  void   ReBin (GridSpacing_t sp);

  friend ostream &
           operator << (ostream & stream, const UnifGridDimension & gd);

  void Print (ostream & stream) const;

private:
  double fStep;
  double fLocMin;
  double fLocMax;
};

}        // genie namespace
#endif   // _UNIFORM_GRID_DIMENSION_H_
