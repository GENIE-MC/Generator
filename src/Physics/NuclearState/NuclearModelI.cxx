 //____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

 @ Mar 18, 2016- Joe Johnston (SD)
   Update GenerateNucleon() and Prob() to accept a radius as the argument,
   and call the corresponding methods in the nuclear model with a radius.

*/
//____________________________________________________________________________


#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/NuclearState/NuclearModelI.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Numerical/RandomGen.h"

using std::ostringstream;
using namespace genie;
using namespace genie::constants;
using namespace genie::controls;

//____________________________________________________________________________

bool NuclearModelI::GenerateNucleon(const Target & tgt,
                                    double /*hitNucleonRadius*/) const
  {
    return GenerateNucleon(tgt);
  }

double NuclearModelI::Prob(double p, double w, const Target & tgt,
                           double /*hitNucleonRadius*/) const
  {
    return Prob(p,w,tgt);
  }

//____________________________________________________________________________
