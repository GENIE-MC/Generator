//____________________________________________________________________________
/*!

\class    genie::UnifGrid

\brief    A multidimensional grid with uniform axis spacing (in some metric)

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  July 01, 2004

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _UNIFORM_GRID_H_
#define _UNIFORM_GRID_H_

#include <iostream>
#include <vector>

#include "Numerical/GridSpacing.h"

using std::ostream;
using std::vector;

namespace genie {

class GSFunc;
class UnifGridDimension;
class Registry;

class UnifGrid
{
public:
  UnifGrid();
  UnifGrid(const GSFunc & gsfunc, GridSpacing_t sp = kGSpLinear, unsigned int np=2);
  UnifGrid(const UnifGrid & grid);
  ~UnifGrid();

  double       Coord   (unsigned int idim, unsigned int ipoint) const;
  unsigned int NPoints (void) const;
  void         Point   (unsigned int ipnt, vector<unsigned int> & grid_pos) const;

  long int     UId              (const vector<unsigned int> & grid_pos) const;
  long int     UId              (const vector<double> & x_pos) const;
  Registry &   GridAttributes   (void);

  unsigned int GetNDimensions   (void) const;
  void         AddDimension     (const UnifGridDimension & gd);
  UnifGridDimension *
               GridDimension    (unsigned int idim);

  //! operators
  UnifGrid & operator  =  (const UnifGrid & grid);
  double     operator ( ) (unsigned int idim, unsigned int ipoint) const;
  const UnifGridDimension *
             operator [ ] (unsigned int idim) const;
  friend ostream &
             operator << (ostream & stream, const UnifGrid & grid);

  //! Copy, Reset and Print
  void Copy  (const UnifGrid & grid);
  void Reset (void);
  void Print (ostream & stream) const;

private:

  //! Initialize and CleanUp
  void Init    (void);
  void CleanUp (void);

  //! the n-dimensional grid
  vector<UnifGridDimension *> fGrid;

  //! a set of associated named grid attributes (Registry is a type-safe
  //! "name"->"value" associative container)
  Registry * fGridAttr;
};

}        // namespace

#endif   // _UNIFORM_GRID_H_
