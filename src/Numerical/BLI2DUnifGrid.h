//____________________________________________________________________________
/*!

\class    genie::BLI2DUnifGrid

\brief    Bilinear interpolation of 2D functions on a regular grid.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  May 30, 2009

*/
//____________________________________________________________________________

#ifndef _BILLINEAR_INTERPOLATION_2D_UNIFORM_GRID_H_
#define _BILLINEAR_INTERPOLATION_2D_UNIFORM_GRID_H_

#include <TObject.h>

namespace genie {

class BLI2DUnifGrid : public TObject {

public:
  //-- ctors & dtor
  BLI2DUnifGrid();
  BLI2DUnifGrid(int nx, double xmin, double xmax, int ny, double ymin, double ymax);
  BLI2DUnifGrid(int nx, int ny, double *x, double *y, double **z);
  virtual ~BLI2DUnifGrid();

  //-- add another point in the grid
  bool AddPoint(double x, double y, double z);

  //-- evaluate the function at the input position
  double Evaluate (double x, double y) const;

  // report min/max values
  double XMin (void) const { return fXmin; }
  double XMax (void) const { return fXmax; }
  double YMin (void) const { return fYmin; }
  double YMax (void) const { return fYmax; }
  double ZMin (void) const { return fZmin; }
  double ZMax (void) const { return fZmax; }

private:

  void Init (int nx=0, double xmin=0, double xmax=0, int ny=0, double ymin=0, double ymax=0);
  int  IdxZ (int ix, int iy) const;

  //-- private data members
  int      fNX;
  int      fNY;
  int      fNZ;
  double * fX;  //[fNX] 
  double * fY;  //[fNY]  
  double * fZ;  //[fNZ] 
  double   fDX;
  double   fDY;
  double   fXmin;
  double   fXmax;
  double   fYmin;
  double   fYmax;
  double   fZmin;
  double   fZmax;

ClassDef(BLI2DUnifGrid, 1)
};

}
#endif
