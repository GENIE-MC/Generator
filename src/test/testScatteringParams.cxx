//____________________________________________________________________________
/*!

\program testScatteringParams

\brief   test program used for testing / debugging ScatteringParams

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created May 4, 2004
*/
//____________________________________________________________________________

#include "Interaction/ScatteringParams.h"
#include "Messenger/Messenger.h"

using namespace genie;

int main(int argc, char ** argv)
{
 //------------ Build a Registry, unlock it, add some vars, lock it, print it.

 ScatteringParams scp;

 scp.UnLock();

 scp.Set("x",  0.1781);
 scp.Set("y",  0.6892);
 scp.Set("Q2", 3.2218);
 scp.Set("W",  1.9219);

 scp.Lock();

 LOG("Main", pINFO) << scp;

 double x  = scp.x();
 double y  = scp.y();
 double Q2 = scp.Q2();
 double W  = scp.W();

 LOG("Main", pINFO) << "x  = " << x;
 LOG("Main", pINFO) << "y  = " << y;
 LOG("Main", pINFO) << "Q2 = " << Q2;
 LOG("Main", pINFO) << "W  = " << W;
}

