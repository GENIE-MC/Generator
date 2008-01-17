//____________________________________________________________________________
/*!

\program testBaryonResDataSet

\brief   test program used for testing the Baryon Resonance Data Sets

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created May 12, 2004

\cpright Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include "Algorithm/Algorithm.h"
#include "Algorithm/AlgFactory.h"
#include "BaryonResonance/BaryonResonance.h"
#include "BaryonResonance/BaryonResParams.h"
#include "BaryonResonance/BaryonResDataSetI.h"
#include "BaryonResonance/BaryonResDataPDG.h"
#include "Messenger/Messenger.h"

using namespace genie;

int main(int /*argc*/, char ** /*argv*/)
{
 //-- Get a baryon resonance
 Resonance_t resonance = kP33_1232;

 //-- Get a concrete implementation of the BaryonResDataSetI interface.
 AlgFactory * algf = AlgFactory::Instance();

 const BaryonResDataSetI * dataset =
       dynamic_cast<const BaryonResDataSetI *> (
             algf->GetAlgorithm("genie::BaryonResDataPDG","Default"));

 //-- Instantiate a BaryonResParams object 
 BaryonResParams res_params;

 //-- Set a baryon resonance data set
 res_params.SetDataSet(dataset);

 //-- Retrieve data for the input resonance
 res_params.RetrieveData(resonance);

 //-- Print the data
 LOG("Main", pINFO) << "Baryon resonance data:\n" << res_params;
}

