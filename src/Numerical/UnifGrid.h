//____________________________________________________________________________
/*!

\class    genie::UnifGrid

\brief

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  July 01, 2004

*/
//____________________________________________________________________________

#ifndef _UNIFORM_GRID_H_
#define _UNIFORM_GRID_H_

#include <vector>

using std::vector;

namespace genie {

class UnifGridDimension
{
public:
  UnifGridDimension()  { }
  ~UnifGridDimension() { }

  unsigned int npoints;
  double       min;
  double       max;
  double       step;
};

class UnifGrid
{
public:

  UnifGrid();
  UnifGrid(const UnifGrid & grid);
  virtual ~UnifGrid();

  void AddDimension   (int npoints, double min, double max);
  unsigned int GetNDimensions (void) const;

  const UnifGridDimension * operator[] (unsigned int idimension) const;

  long int GridPoint2UId(const vector<int> & point) const;
  
private:

  vector<UnifGridDimension *> fGrid;
};

}        // namespace

#endif   // _UNIFORM_GRID_H_
