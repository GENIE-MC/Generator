//____________________________________________________________________________
/*!

\program testNumerical

\brief   

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created May 4, 2004
*/
//____________________________________________________________________________

#include "Numerical/UnifGrid.h"
#include "Numerical/FunctionMap.h"
#include "Messenger/Messenger.h"

using namespace genie;

int main(int argc, char ** argv)
{
  LOG("Main",pINFO) << "Defining a UnifGrid";

  UnifGrid grid2d;

  grid2d.AddDimension(20,  0, 10);
  grid2d.AddDimension(10, -5,  5);

  vector<int> some_position(2);

  some_position[0] = 1;
  some_position[1] = 3;
  
  LOG("Main",pINFO) << "Creating a FunctionMap";

  FunctionMap func_map(grid2d);

  LOG("Main",pINFO) << "Adding points";

  for(int p0 = 0; p0 < 20; p0++) {
    for(int p1 = 0; p1 < 10; p1++) {

         LOG("Main",pINFO) << "(" << p0 << ", " << p1 << ") = " << p1+0.5;

         func_map.AddPoint( p1+0.5, p0, p1);
    }
  }

  for(int i=0; i<20; i++) 
       LOG("Main",pINFO) << "......" << func_map.Func(i,9);

  return 0;
}
