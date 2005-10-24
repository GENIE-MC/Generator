//____________________________________________________________________________
/*!

\class    genie::Simpson2D

\brief

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#include "Numerical/Simpson2D.h"
#include "Numerical/Simpson1D.h"
#include "Messenger/Messenger.h"

using namespace genie;

//____________________________________________________________________________
Simpson2D::Simpson2D():
IntegratorI("genie::Simpson2D")
{

}
//____________________________________________________________________________
Simpson2D::~Simpson2D()
{

}
//____________________________________________________________________________
double Simpson2D::Integrate(FunctionMap & func_map) const
{
  const UnifGrid & grid = func_map.GetGrid();

  assert(grid.GetNDimensions() == 2);
  
  int    Nx     = grid[0]->npoints;
  double xmin   = grid[0]->min;
  double xmax   = grid[0]->max;

  int    Ny     = grid[1]->npoints;
  double ymin   = grid[1]->min;
  double ymax   = grid[1]->max;

  Simpson1D simpson_1d;
      
  UnifGrid subset_grid_x; // 1-D

  subset_grid_x.AddDimension( Nx, xmin, xmax );

  FunctionMap func_map_x( subset_grid_x );

  for(int ix = 0; ix < Nx; ix++) {

      LOG("Simpson2D", pDEBUG)  << "1D integration along ix = " << ix;

      UnifGrid subset_grid_y; // 1-D

      subset_grid_y.AddDimension( Ny, ymin, ymax );

      FunctionMap func_map_y( subset_grid_y );

      for(int iy = 0; iy < Ny; iy++)
                        func_map_y.AddPoint( func_map.Func(ix,iy), iy );

      double integral = simpson_1d.Integrate(func_map_y);

      LOG("Simpson2D", pDEBUG)  
                     << "Integral_y { f(x,y)|x=xo * dy} = " << integral;

      func_map_x.AddPoint( simpson_1d.Integrate(func_map_y), ix );
  }  

  double sum = simpson_1d.Integrate(func_map_x);
  
  return sum;
}
//____________________________________________________________________________

