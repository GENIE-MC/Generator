//____________________________________________________________________________
/*!

\class    genie::UnifGrid

\brief

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  July 01, 2004
 
*/
//____________________________________________________________________________

#include <cassert>

#include "Numerical/UnifGrid.h"
#include "Messenger/Messenger.h"

using namespace genie;

//____________________________________________________________________________
UnifGrid::UnifGrid()
{

}
//____________________________________________________________________________
UnifGrid::UnifGrid(const UnifGrid & grid)
{
  unsigned int ndim = grid.GetNDimensions();

  for(unsigned int idim = 0; idim < ndim; idim++)         
           this->AddDimension( 
                      grid[idim]->npoints, grid[idim]->min, grid[idim]->max );
}
//____________________________________________________________________________
UnifGrid::~UnifGrid()
{
  vector<UnifGridDimension *>::size_type i;
  
  for(i = 0; i < fGrid.size(); i++) delete fGrid[i];

  fGrid.clear();
}
//____________________________________________________________________________
void UnifGrid::AddDimension(int npoints, double min, double max)
{
  UnifGridDimension * dimension = new UnifGridDimension;

  dimension->npoints = npoints;
  dimension->min     = min;
  dimension->max     = max;
  dimension->step    = (max-min)/(npoints-1);

  fGrid.push_back(dimension);
}
//____________________________________________________________________________
unsigned int UnifGrid::GetNDimensions(void) const
{
  return fGrid.size();
}
//____________________________________________________________________________
const UnifGridDimension * UnifGrid::operator[] (unsigned int idim) const
{
  if( idim >= 0 && idim < GetNDimensions() ) return fGrid[idim];
  else return 0;
}
//____________________________________________________________________________
long int UnifGrid::GridPoint2UId(const vector<int> & point) const
{
  // uid = ... + [4] * N1*N2*N3 + [3] * N1*N2 + [2] * N1 + [1];

  assert(point.size() == GetNDimensions());

  long int uid  = 0;
  long int mult = 1;
  
  for(unsigned int idim = 0; idim < GetNDimensions(); idim++) {

     uid  += mult * point[idim];

     mult *= (*this)[idim]->npoints;
  }

  return uid;
}
//____________________________________________________________________________
