//____________________________________________________________________________
/*
 Copyright (c) 2023-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
*/
//____________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/Resonance/XSection/MAIDHelicityAmplModelEMn.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
MAIDHelicityAmplModelEMn::MAIDHelicityAmplModelEMn() :
MAIDHelicityAmplModelI("genie::MAIDHelicityAmplModelEMn")
{

}
//____________________________________________________________________________
MAIDHelicityAmplModelEMn::MAIDHelicityAmplModelEMn(string config) :
MAIDHelicityAmplModelI("genie::MAIDHelicityAmplModelEMn", config)
{

}
//____________________________________________________________________________
MAIDHelicityAmplModelEMn::~MAIDHelicityAmplModelEMn()
{

}
//____________________________________________________________________________
const MAIDHelicityAmpl & MAIDHelicityAmplModelEMn::Compute( const Interaction * interaction ) const {

  Resonance_t res = interaction->ExclTag().Resonance();

  switch(res) {

   case (kP33_1232) :
   {

     break;
   }
   case (kS11_1535) :
   {
     break;
   }
   case (kD13_1520) :
   {
     break;
   }
   case (kS11_1650) :
   {
     break;
   }
   case (kD13_1700) :
   {
    
     break;
   }
   case (kD15_1675) :
   {
    
     break;
   }
   case (kS31_1620) :
   {
    
     break;
   }
   case (kD33_1700) :
   {
    
     break;
   }
   case (kP11_1440) :
   {
     break;
   }
   case (kP33_1600) :
   {

     break;
   }
   case (kP13_1720) :
   {

     break;
   }
   case (kF15_1680) :
   {

     break;
   }
   case (kP31_1910) :
   {

     break;
   }
   case (kP33_1920) :
   {

     break;
   }
   case (kF35_1905) :
   {

     break;
   }
   case (kF37_1950) :
   {

     break;
   }
   case (kP11_1710) :
   {

     break;
   }
   case (kF17_1970) :
   {

     break;
   }
   default:
   {
     LOG("MAIDHAmpl", pWARN) << "*** UNRECOGNIZED RESONANCE!";
     break;
   }

  }//switch

  return fAmpl;
}
//____________________________________________________________________________
