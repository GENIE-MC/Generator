//____________________________________________________________________________
/*!

\class    genie::GridDimension

\brief    Grid dimension ABC

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  July 01, 2004

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GRID_DIMENSION_H_
#define _GRID_DIMENSION_H_

#include "Numerical/GridSpacing.h"

namespace genie {

class GridDimension
{
public:

  virtual ~GridDimension();

  virtual double        Coord        (unsigned int ipoint) const = 0;
  virtual unsigned int  ClosestPoint (double coord)        const = 0;

  virtual unsigned int  NPoints (void) const { return fNPoints; }
  virtual double        Min     (void) const { return fMin;     }
  virtual double        Max     (void) const { return fMax;     }
  virtual GridSpacing_t Spacing (void) const { return fSpacing; }

  virtual void Set(GridSpacing_t sp, double min, double max, unsigned int n);

protected:

  GridDimension();
  GridDimension(GridSpacing_t sp, double min, double max, unsigned int n=2);

  GridSpacing_t fSpacing;
  double        fMin;
  double        fMax;
  unsigned int  fNPoints;
};

}        // namespace
#endif   // _GRID_DIMENSION_H_
