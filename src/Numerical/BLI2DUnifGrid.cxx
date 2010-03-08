//____________________________________________________________________________
/*
 Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - May 30, 2009

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ May 30, 2009 - CA
   Was first added in v2.5.1.
*/
//____________________________________________________________________________

#include <TMath.h>

#include <cassert>

#include "Messenger/Messenger.h"
#include "Numerical/BLI2DUnifGrid.h"

using namespace genie;

ClassImp(BLI2DUnifGrid)

//___________________________________________________________________________
BLI2DUnifGrid::BLI2DUnifGrid()
{
  this->Init();
}
//___________________________________________________________________________
BLI2DUnifGrid::BLI2DUnifGrid(
  int nx, double xmin, double xmax, int ny, double ymin, double ymax)
{
  this->Init(nx, xmin, xmax, ny, ymin, ymax);
}
//___________________________________________________________________________
BLI2DUnifGrid::BLI2DUnifGrid(
    int nx, int ny, double *x, double *y, double **z) 
{
  double xmin = x[0];
  double xmax = x[nx-1];
  double ymin = y[0];
  double ymax = y[nx-1];

  this->Init(nx, xmin, xmax, ny, ymin, ymax);

  for(int ix=0; ix<nx; ix++) {
    for(int iy=0; iy<ny; iy++) {
      this->AddPoint(x[ix], y[iy], z[ix][iy]);
    }
  }
}
//___________________________________________________________________________
BLI2DUnifGrid::~BLI2DUnifGrid()
{
  if (fX) { delete [] fX; }
  if (fY) { delete [] fY; }
  if (fZ) { delete [] fZ; }
}
//___________________________________________________________________________
bool BLI2DUnifGrid::AddPoint(double x, double y, double z)
{
  int ix = (x - fXmin + fDX/2) / fDX;
  int iy = (y - fYmin + fDY/2) / fDY;
  int iz = this->IdxZ(ix,iy);

  fZ[iz] = z;

  fZmin = TMath::Min(z, fZmin);
  fZmax = TMath::Min(z, fZmax);

  LOG("BLI2DUnifGrid", pDEBUG) 
    << "Added x = " << x << " (ix = " << ix << ")"
         << " y = " << y << " (iy = " << iy << ") -> "
         << " z = " << z << " (iz = " << iz << ")";

  return true;
}
//___________________________________________________________________________
double BLI2DUnifGrid::Evaluate(double x, double y) const
{
  if(x < fXmin || x > fXmax) return 0.;
  if(y < fYmin || y > fYmax) return 0.;

  int ix_lo  = (x - fXmin) / fDX; 
  int iy_lo  = (y - fYmin) / fDY; 
  int ix_hi  = ix_lo + 1;
  int iy_hi  = iy_lo + 1;

  double x1  = fX[ix_lo];
  double x2  = fX[ix_hi];
  double y1  = fY[iy_lo];
  double y2  = fY[iy_hi];

  double z11 = fZ[ this->IdxZ(ix_lo,iy_lo) ];
  double z21 = fZ[ this->IdxZ(ix_hi,iy_lo) ];
  double z12 = fZ[ this->IdxZ(ix_lo,iy_hi) ];
  double z22 = fZ[ this->IdxZ(ix_hi,iy_hi) ];

  double z1  = z11 * (x2-x)/(x2-x1) + z21 * (x-x1)/(x2-x1);
  double z2  = z12 * (x2-x)/(x2-x1) + z22 * (x-x1)/(x2-x1);
  double z   = z1  * (y2-y)/(y2-y1) + z2  * (y-y1)/(y2-y1);

/*
  LOG("BLI2DUnifGrid", pDEBUG) << "x = " << x << " -> nearby nodes: " << x1 << ", " << x2;
  LOG("BLI2DUnifGrid", pDEBUG) << "y = " << x << " -> nearby nodes: " << y1 << ", " << y2;
  LOG("BLI2DUnifGrid", pDEBUG) << "z11 := z(" << this->IdxZ(ix_lo,iy_lo) << ")  = " << z11;
  LOG("BLI2DUnifGrid", pDEBUG) << "z21 := z(" << this->IdxZ(ix_hi,iy_lo) << ")  = " << z21;
  LOG("BLI2DUnifGrid", pDEBUG) << "z12 := z(" << this->IdxZ(ix_lo,iy_hi) << ")  = " << z12;
  LOG("BLI2DUnifGrid", pDEBUG) << "z22 := z(" << this->IdxZ(ix_hi,iy_hi) << ")  = " << z22;
  LOG("BLI2DUnifGrid", pDEBUG) << "z1 = " << z1 << ", z2 = " << z2;
  LOG("BLI2DUnifGrid", pDEBUG) << "interpolated z(x,y) = " << z;
*/

  return z;
}
//___________________________________________________________________________
void BLI2DUnifGrid::Init(
  int nx, double xmin, double xmax, int ny, double ymin, double ymax)
{
  fNX   = 0;
  fNY   = 0;
  fNZ   = 0;
  fXmin = 0.;
  fXmax = 0.;
  fYmin = 0.;
  fYmax = 0.;
  fZmin = 0.;
  fZmax = 0.;
  fDX   = 0.;
  fDY   = 0.;
  fX    = 0;
  fY    = 0;
  fZ    = 0;

  if(nx>1 && ny>1) {
    fNX = nx;
    fNY = ny;
    fNZ = nx * ny;

    fXmin = xmin;
    fXmax = xmax;
    fYmin = ymin;
    fYmax = ymax;
    fZmin = 0.;
    fZmax = 0.;

    fDX = (xmax-xmin)/(nx-1);
    fDY = (ymax-ymin)/(ny-1);

    fX = new double[fNX];
    fY = new double[fNY];
    fZ = new double[fNZ];

    for(int i=0; i<fNX; i++) { fX[i] = xmin + i*fDX; }
    for(int i=0; i<fNY; i++) { fY[i] = ymin + i*fDY; }
    for(int i=0; i<fNZ; i++) { fZ[i] = 0.; }
  }
}
//___________________________________________________________________________
int BLI2DUnifGrid::IdxZ(int ix, int iy) const
{
  return ix*fNY+iy;
}
//___________________________________________________________________________

