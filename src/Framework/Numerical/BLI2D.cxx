//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - May 30, 2009

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ May 30, 2009 - CA
   Was first added in v2.5.1.
 @ Mar 08, 2010 - CA
   Fix `converting to `int' from `double'' warnings by using TMath::FloorInt
   in all lines similar to `int ix = (xmax-xmin)/dx'.
 @ July 29, 2011 - AM
   Added BLI2DNonUnifGrid.

*/
//____________________________________________________________________________

#include <TMath.h>

#include <cassert>
#include <limits>

#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/BLI2D.h"

using namespace genie;

ClassImp(BLI2DGrid)

//___________________________________________________________________________
BLI2DGrid::BLI2DGrid()
{

}
//___________________________________________________________________________
BLI2DGrid::~BLI2DGrid()
{
  if (fX) { delete [] fX; }
  if (fY) { delete [] fY; }
  if (fZ) { delete [] fZ; }
}
//___________________________________________________________________________
int BLI2DGrid::IdxZ(int ix, int iy) const
{
  return ix*fNY+iy;
}
//___________________________________________________________________________
//___________________________________________________________________________
//___________________________________________________________________________
ClassImp(BLI2DUnifGrid)

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
  int nx, int ny, double *x, double *y, double *z)
{
  double xmin = x[0];
  double xmax = x[nx-1];
  double ymin = y[0];
  double ymax = y[nx-1];

  this->Init(nx, xmin, xmax, ny, ymin, ymax);

  for(int ix=0; ix<nx; ix++) {
    for(int iy=0; iy<ny; iy++) {
      this->AddPoint(x[ix], y[iy], z[this->IdxZ(ix,iy)]);
    }
  }
}
//___________________________________________________________________________
bool BLI2DUnifGrid::AddPoint(double x, double y, double z)
{
  int ix = TMath::FloorNint( (x - fXmin + fDX/2) / fDX );
  int iy = TMath::FloorNint( (y - fYmin + fDY/2) / fDY );
  int iz = this->IdxZ(ix,iy);

  fZ[iz] = z;

  fZmin = TMath::Min(z, fZmin);
  fZmax = TMath::Max(z, fZmax);

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

  int ix_lo  = TMath::FloorNint( (x - fXmin) / fDX ); 
  int iy_lo  = TMath::FloorNint( (y - fYmin) / fDY ); 
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
  fZmin = std::numeric_limits<double>::max();
  fZmax = std::numeric_limits<double>::min();
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
    fZmin = std::numeric_limits<double>::max();
    fZmax = std::numeric_limits<double>::min();

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
//___________________________________________________________________________
//___________________________________________________________________________
ClassImp(BLI2DNonUnifGrid)

BLI2DNonUnifGrid::BLI2DNonUnifGrid()
{
  this->Init();
}
//___________________________________________________________________________
BLI2DNonUnifGrid::BLI2DNonUnifGrid(
  int nx, double xmin, double xmax, int ny, double ymin, double ymax)
{
  this->Init(nx, xmin, xmax, ny, ymin, ymax);
}
//___________________________________________________________________________
BLI2DNonUnifGrid::BLI2DNonUnifGrid(
  int nx, int ny, double *x, double *y, double *z)
{
  double xmin = x[0];
  double xmax = x[nx-1];
  double ymin = y[0];
  double ymax = y[ny-1];

  this->Init(nx, xmin, xmax, ny, ymin, ymax);

  for(int ix=0; ix<nx; ix++) {
    for(int iy=0; iy<ny; iy++) {
      this->AddPoint(x[ix], y[iy], z[this->IdxZ(ix,iy)]);
    }
  }
}
//___________________________________________________________________________
bool BLI2DNonUnifGrid::AddPoint(double x, double y, double z)
{

  // check the x,y values' existence before moving anything
  //   if they do, use them
  //   if they don't, add them in the appropriate place
  //
  int xidx = -1;
  int yidx = -1;
  bool changex = true;
  bool changey = true;
  // first, check and see if the x,y values exist
  // NOTE: == should not be used with double values:
  //   instead, a tolerance of ~5 significant digits is checked
  for(int i=0;i<fNFillX;i++)
    {
      if (!(TMath::Abs(x-fX[i])<=.5*.0000001*(TMath::Abs(x)+TMath::Abs(fX[i]))) && x<fX[i])
	{
	  // missed one, make sure there aren't too many x and move everything up
	  if (fNFillX==fNX)
	    {
	      LOG("BLI2DNonUnifGrid", pWARN) << "Warning: too many x values, cannot add x = "<<x;
	      return false;
	    }
	  else
	    {
	      xidx=i;
	      changex=true;
	      break;
	    }
	}
      else
	{
	  if (TMath::Abs(x-fX[i])<=.5*.0000001*(TMath::Abs(x)+TMath::Abs(fX[i]))) {
	    xidx=i;
	    LOG("BLI2DNonUnifGrid", pDEBUG) << "x value found at index "<<i;
	    changex=false; break;}
	  changex=true;
	}
    }
  if (changex && xidx<0) xidx=fNFillX;
  if (changex) LOG("BLI2DNonUnifGrid", pDEBUG) << "Adding x value at index "<<xidx;

  for(int i=0;i<fNFillY;i++)
    {
      if (!(TMath::Abs(y-fY[i])<=.5*.0000001*(TMath::Abs(y)+TMath::Abs(fY[i]))) && y<fY[i])
	{
	  if (fNFillY==fNY)
	    {
	      LOG("BLI2DNonUnifGrid", pWARN) << "Warning: too many y values, cannot add y = "<<y;
	      return false;
	    }
	  else
	    {
	      yidx=i;
	      changey=true;
	      break;
	    }
	}
      else
	{
	  if (TMath::Abs(y-fY[i])<=.5*.0000001*(TMath::Abs(y)+TMath::Abs(fY[i]))) {
	    yidx=i;
	    LOG("BLI2DNonUnifGrid", pDEBUG) << "y value found at index "<<i;
	    changey=false; break;}
	  changey=true;
	}
    }
  if (changey && yidx<0) yidx=fNFillY;
  if (changey) LOG("BLI2DNonUnifGrid", pDEBUG) << "Adding y value at index "<<yidx;

  // make new entries if needed
  if (changex && xidx>=0)
    {
      for (int j=0;j<fNFillX-xidx;j++)
	{
	  fX[fNFillX-j]=fX[fNFillX-j-1];
	  for (int k=0;k<fNFillY;k++)
	    {
	      fZ[ this->IdxZ(fNFillX-j,k) ] 
		= fZ[ this->IdxZ(fNFillX-j-1,k) ];
	    }
	}
      fX[xidx]=x;
      fNFillX++;
      fXmin = TMath::Min(x,fXmin);
      fXmax = TMath::Max(x,fXmax);
    }
  if (changey && yidx>=0)
    {
      for (int j=0;j<fNFillY-yidx;j++)
	{
	  fY[fNFillY-j]=fY[fNFillY-j-1];
	  for (int k=0;k<fNFillX;k++)
	    {
	      fZ[ this->IdxZ(k,fNFillY-j) ] 
		= fZ[ this->IdxZ(k,fNFillY-j-1) ];
	    }
	}
      fY[yidx]=y;
      fNFillY++;
      fYmin = TMath::Min(y,fYmin);
      fYmax = TMath::Max(y,fYmax);
    }

  int iz = this->IdxZ(xidx,yidx);

  fZ[iz] = z;

  fZmin = TMath::Min(z, fZmin);
  fZmax = TMath::Max(z, fZmax);

  LOG("BLI2DNonUnifGrid", pDEBUG) 
    << "Added x = " << x << " (ix = " << xidx << ")"
         << " y = " << y << " (iy = " << yidx << ") -> "
         << " z = " << z << " (iz = " << iz   << ")";

  return true;
}
//___________________________________________________________________________
double BLI2DNonUnifGrid::Evaluate(double x, double y) const
{
  double evalx=TMath::Min(x,fXmax);
  evalx=TMath::Max(evalx,fXmin);
  double evaly=TMath::Min(y,fYmax);
  evaly=TMath::Max(evaly,fYmin);

  int ix_lo  = -2; 
  int iy_lo  = -2;
  for (int i=0;i<fNFillX;i++)
    {
      if (evalx<=fX[fNFillX-1-i]) ix_lo=fNFillX-2-i;
      else break;
    }
  for (int i=0;i<fNFillY;i++)
    {
      if (evaly<=fY[fNFillY-1-i]) iy_lo=fNFillY-2-i;
      else break;
    }
  int ix_hi  = ix_lo + 1;
  int iy_hi  = iy_lo + 1;

  // in case x = xmin
  if (ix_lo==-1) {ix_lo++; ix_hi++;}
  if (iy_lo==-1) {iy_lo++; iy_hi++;}
  // in case x = xmax
  if (ix_lo==-2) {ix_lo=fNFillX-2; ix_hi=fNFillX-1;}
  if (iy_lo==-2) {iy_lo=fNFillY-2; iy_hi=fNFillY-1;}
  // if an error occurs
  if (ix_lo<0       || iy_lo<0      ) return 0.;
  if (ix_hi>fNFillX || iy_hi>fNFillY) return 0.;

  double x1  = fX[ix_lo];
  double x2  = fX[ix_hi];
  double y1  = fY[iy_lo];
  double y2  = fY[iy_hi];

  double z11 = fZ[ this->IdxZ(ix_lo,iy_lo) ];
  double z21 = fZ[ this->IdxZ(ix_hi,iy_lo) ];
  double z12 = fZ[ this->IdxZ(ix_lo,iy_hi) ];
  double z22 = fZ[ this->IdxZ(ix_hi,iy_hi) ];

  double z1  = z11 * (x2-evalx)/(x2-x1) + z21 * (evalx-x1)/(x2-x1);
  double z2  = z12 * (x2-evalx)/(x2-x1) + z22 * (evalx-x1)/(x2-x1);
  double z   = z1  * (y2-evaly)/(y2-y1) + z2  * (evaly-y1)/(y2-y1);

  /*  
  LOG("BLI2DNonUnifGrid", pINFO) << "x = " << evalx << " -> nearby nodes: " << x1 << ", " << x2;
  LOG("BLI2DNonUnifGrid", pINFO) << "y = " << evaly << " -> nearby nodes: " << y1 << ", " << y2;
  LOG("BLI2DNonUnifGrid", pINFO) << "xmin = " << fXmin << ", xmax = " << fXmax;
  LOG("BLI2DNonUnifGrid", pINFO) << "ymin = " << fYmin << ", ymax = " << fYmax;
  LOG("BLI2DNonUnifGrid", pINFO) << "z11 := z(" << this->IdxZ(ix_lo,iy_lo) << ")  = " << z11;
  LOG("BLI2DNonUnifGrid", pINFO) << "z21 := z(" << this->IdxZ(ix_hi,iy_lo) << ")  = " << z21;
  LOG("BLI2DNonUnifGrid", pINFO) << "z12 := z(" << this->IdxZ(ix_lo,iy_hi) << ")  = " << z12;
  LOG("BLI2DNonUnifGrid", pINFO) << "z22 := z(" << this->IdxZ(ix_hi,iy_hi) << ")  = " << z22;
  LOG("BLI2DNonUnifGrid", pINFO) << "z1 = " << z1 << ", z2 = " << z2;
  LOG("BLI2DNonUnifGrid", pINFO) << "interpolated z(x,y) = " << z;
  */

  return z;
}
//___________________________________________________________________________
void BLI2DNonUnifGrid::Init(
  int nx, double xmin, double xmax, int ny, double ymin, double ymax)
{
  fNX    = 0;
  fNY    = 0;
  fNZ    = 0;
  fNFillX= 0;
  fNFillY= 0;
  fXmin  = 0.;
  fXmax  = 0.;
  fYmin  = 0.;
  fYmax  = 0.;
  fZmin  = std::numeric_limits<double>::max();
  fZmax  = std::numeric_limits<double>::min(); 
  fX     = 0;
  fY     = 0;
  fZ     = 0;

  if(nx>1 && ny>1) {
    fNX = nx;
    fNY = ny;
    fNZ = nx * ny;

    fXmin = xmin;
    fXmax = xmax;
    fYmin = ymin;
    fYmax = ymax;
    fZmin = std::numeric_limits<double>::max();
    fZmax = std::numeric_limits<double>::min();

    fX = new double[fNX];
    fY = new double[fNY];
    fZ = new double[fNZ];

    for(int i=0; i<fNX; i++) { fX[i] = 0.; }
    for(int i=0; i<fNY; i++) { fY[i] = 0.; }
    for(int i=0; i<fNZ; i++) { fZ[i] = 0.; }
  }
}
//___________________________________________________________________________
