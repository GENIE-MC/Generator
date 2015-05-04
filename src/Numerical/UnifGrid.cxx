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

#include "Numerical/UnifGridDimension.h"
#include "Numerical/UnifGrid.h"
#include "Numerical/GSFunc.h"
#include "Messenger/Messenger.h"
#include "Registry/Registry.h"
#include "Utils/Range1.h"

using std::endl;
using namespace genie;

//____________________________________________________________________________
namespace genie
{
  ostream & operator << (ostream & stream, const UnifGrid & grid)
  {
     grid.Print(stream);
     return stream;
  }
}
//____________________________________________________________________________
UnifGrid::UnifGrid()
{
  this->Init();
}
//____________________________________________________________________________
UnifGrid::UnifGrid(const UnifGrid & grid)
{
  this->Init();
  this->Copy(grid);
}
//____________________________________________________________________________
UnifGrid::UnifGrid(const GSFunc & gsfunc, GridSpacing_t sp, unsigned int np)
{
  this->Init();

  unsigned int ndim = gsfunc.NParams();
  UnifGridDimension * dim = 0;
  for(unsigned int idim = 0; idim < ndim; idim++) {
    Range1D_t l = gsfunc.ParamLimits(idim);
    dim = new UnifGridDimension(sp, l.min, l.max, np);
    fGrid.push_back(dim);
  }
}
//____________________________________________________________________________
UnifGrid::~UnifGrid()
{
  this->CleanUp();
}
//____________________________________________________________________________
void UnifGrid::AddDimension(const UnifGridDimension & gd)
{
  UnifGridDimension * dim = new UnifGridDimension(gd);
  fGrid.push_back(dim);
}
//____________________________________________________________________________
unsigned int UnifGrid::GetNDimensions(void) const
{
  return fGrid.size();
}
//____________________________________________________________________________
UnifGridDimension * UnifGrid::GridDimension(unsigned int idim)
{
  if(idim < this->GetNDimensions()) return fGrid[idim];
  else return 0;
}
//____________________________________________________________________________
const UnifGridDimension * UnifGrid::operator[] (unsigned int idim) const
{
  if(idim < this->GetNDimensions()) return fGrid[idim];
  else return 0;
}
//____________________________________________________________________________
double UnifGrid::operator() (unsigned int idim, unsigned int ipoint) const
{
  return this->Coord(idim, ipoint);
}
//____________________________________________________________________________
UnifGrid & UnifGrid::operator = (const UnifGrid & grid)
{
  this->Copy(grid);
  return (*this);
}
//____________________________________________________________________________
Registry & UnifGrid::GridAttributes(void)
{
  if(!fGridAttr) fGridAttr = new Registry("grid attributes", false);
  return *fGridAttr;
}
//____________________________________________________________________________
double UnifGrid::Coord(unsigned int idim, unsigned int ipoint) const
{
  assert(idim < this->GetNDimensions());

  double coord = fGrid[idim]->Coord(ipoint);
  return coord;
}
//____________________________________________________________________________
unsigned int UnifGrid::NPoints(void) const
{
// total number of points in the multi-dimensional grid

  const UnifGrid & grid = *this;

  unsigned int nd = grid.GetNDimensions();
  unsigned int np = 1;

  for(unsigned int id = 0; id < nd; id++) {
    np *= grid[id]->NPoints();
  }
  return np;
}
//____________________________________________________________________________
void UnifGrid::Point(unsigned int ipnt, vector<unsigned int> & grid_pos) const
{
// the ipnt'th grid point in the std UnifGrid tranversal

  const UnifGrid & grid = *this;

  unsigned int nd = grid.GetNDimensions();
  assert(grid_pos.size() == nd);

  unsigned int np = grid.NPoints();
  assert(ipnt<np);

  bool set=false;
  for(int id=nd-1; id>=0; id--) {

    if(set) continue;
    if(grid_pos[id]+1 < grid[id]->NPoints()) {
         set=true;
         grid_pos[id]++;
    } else grid_pos[id]=0;
  }
}
//____________________________________________________________________________
long int UnifGrid::UId(const vector<unsigned int> & grid_pos) const
{
// uid = ... + [4] * N1*N2*N3 + [3] * N1*N2 + [2] * N1 + [1];

  const UnifGrid & grid = *this;

  unsigned int ndim = grid.GetNDimensions();
  assert(grid_pos.size() == ndim);

  long int uid  = 0;
  long int mult = 1;
  for(unsigned int idim = 0; idim < ndim; idim++) {
     uid  += mult * grid_pos[idim];
     mult *= grid[idim]->NPoints();
  }
  return uid;
}
//____________________________________________________________________________
long int UnifGrid::UId(const vector<double> & x_pos) const
{
  const UnifGrid & grid = *this;

  unsigned int ndim = grid.GetNDimensions();
  assert(x_pos.size() == ndim);

  vector<unsigned int> grid_pos(ndim);

  for(unsigned int idim = 0; idim < ndim; idim++) {
    grid_pos[idim] = grid[idim]->ClosestPoint(x_pos[idim]);
  }
  return this->UId(grid_pos);
}
//____________________________________________________________________________
void UnifGrid::Copy(const UnifGrid & grid)
{
  this->Reset();

  unsigned int ndim = grid.GetNDimensions();
  for(unsigned int i = 0; i < ndim; i++) {
    const UnifGridDimension * dim = grid[i];
    this->AddDimension(*dim);
  }
  if(grid.fGridAttr) {
    this->fGridAttr = new Registry(*grid.fGridAttr);
  }
}
//____________________________________________________________________________
void UnifGrid::Reset(void)
{
  this->CleanUp();
  this->Init();
}
//____________________________________________________________________________
void UnifGrid::Print(ostream & stream) const
{
  stream << endl;
  stream << "ndimensions: " << this->GetNDimensions() << endl;

  vector<UnifGridDimension *>::size_type i;
  for(i = 0; i < fGrid.size(); i++) stream << *fGrid[i];
}
//____________________________________________________________________________
void UnifGrid::Init(void)
{
  fGrid.clear();
  fGridAttr = 0;
}
//____________________________________________________________________________
void UnifGrid::CleanUp(void)
{
  vector<UnifGridDimension *>::size_type i;
  for(i = 0; i < fGrid.size(); i++) delete fGrid[i];
  fGrid.clear();

  if(fGridAttr) delete fGridAttr;
}
//____________________________________________________________________________

