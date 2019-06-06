//____________________________________________________________________________
/*!

\class    genie::BLI2DUnifGrid

\brief    Bilinear interpolation of 2D functions on a regular grid.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 30, 2009

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _BILLINEAR_INTERPOLATION_2D_GRID_H_
#define _BILLINEAR_INTERPOLATION_2D_GRID_H_

#include <TObject.h>

namespace genie {

class BLI2DGrid : public TObject {

public:
  //-- ctors & dtor
  BLI2DGrid();
  //BLI2DGrid(int nx, double xmin, double xmax, int ny, double ymin, double ymax);
  //BLI2DGrid(int nx, int ny, double *x, double *y, double **z);
  virtual ~BLI2DGrid();

  //-- add another point in the grid
  virtual bool AddPoint(double x, double y, double z) =0;

  //-- evaluate the function at the input position
  virtual double Evaluate (double x, double y) const =0;

  // report min/max values
  double XMin (void) const { return fXmin; }
  double XMax (void) const { return fXmax; }
  double YMin (void) const { return fYmin; }
  double YMax (void) const { return fYmax; }
  double ZMin (void) const { return fZmin; }
  double ZMax (void) const { return fZmax; }

protected:

  virtual void Init (int nx, double xmin, double xmax, int ny, double ymin, double ymax) =0;
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

  ClassDef(BLI2DGrid, 1)
  };

//____________________________________________________________________________
//____________________________________________________________________________

class BLI2DUnifGrid : public BLI2DGrid {

public:
  //-- ctors & dtor
  BLI2DUnifGrid();
  BLI2DUnifGrid(int nx, double xmin, double xmax, int ny, double ymin, double ymax);
  BLI2DUnifGrid(int nx, int ny, double *x, double *y, double *z);

  //-- add another point in the grid
  bool AddPoint(double x, double y, double z);

  //-- evaluate the function at the input position
  double Evaluate (double x, double y) const;

private:

  void Init (int nx=0, double xmin=0, double xmax=0, int ny=0, double ymin=0, double ymax=0);

  ClassDef(BLI2DUnifGrid, 1)
  };

//____________________________________________________________________________
//____________________________________________________________________________

class BLI2DNonUnifGrid : public BLI2DGrid {

public:
  //-- ctors & dtor
  BLI2DNonUnifGrid();
  BLI2DNonUnifGrid(int nx, double xmin, double xmax, int ny, double ymin, double ymax);
  BLI2DNonUnifGrid(int nx, int ny, double *x, double *y, double *z);

  //-- add another point in the grid
  bool AddPoint(double x, double y, double z);

  //-- evaluate the function at the input position
  double Evaluate (double x, double y) const;

private:

  void Init (int nx=0, double xmin=0, double xmax=0, int ny=0, double ymin=0, double ymax=0);
  int      fNFillX;
  int      fNFillY;

  ClassDef(BLI2DNonUnifGrid, 1)
  };

}
#endif
