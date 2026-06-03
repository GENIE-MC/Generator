//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Daniel Scully ( d.i.scully \at warwick.ac.uk)
 University of Warwick
*/
//____________________________________________________________________________

#include <iostream>
#include <sstream>
#include <string>
#include <complex>

#include "Physics/Coherent/XSection/ARWFSolution.h"

namespace genie
{
namespace alvarezruso
{

ARWFSolution::ARWFSolution(bool debug): debug_(debug)
{
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  if(debug_) std::cerr << "WFS@ constructor" << std::endl;
#endif
}

ARWFSolution::~ARWFSolution()
{
}

} //namespace alvarezruso
} //namespace genie
