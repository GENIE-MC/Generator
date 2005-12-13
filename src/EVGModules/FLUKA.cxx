//____________________________________________________________________________
/*!

\class   genie::FLUKA

\brief   A GENIE interface to the FLUKA hadron transport code.

         Is a concerete implementation of the EventRecordVisitorI interface.

         Note: the FLUKA code is *not included* in your GENIE installation.
         You need to obtain FLUKA from its official distribution point.

\ref     More information at: http://www.fluka.org

         FLUKA Authors:
         G.Battistoni, A.Ferrari, P.R.Sala (INFN & Univ. Milano, CERN)

         The FLUKA code is maintained and developed under INFN-CERN agreement
         and copyright 1989-2005.

         Please cite FLUKA separately if you include this event generation
         module in your event generation threads.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created December 13, 2005

*/
//____________________________________________________________________________

#include <TLorentzVector.h>
#include <TVector3.h>

#include "EVGCore/EVGThreadException.h"
#include "EVGModules/FLUKA.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepParticle.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"

using namespace genie;

//___________________________________________________________________________
FLUKA::FLUKA() :
EventRecordVisitorI("genie::FLUKA")
{

}
//___________________________________________________________________________
FLUKA::FLUKA(string config) :
EventRecordVisitorI("genie::FLUKA",  config)
{

}
//___________________________________________________________________________
FLUKA::~FLUKA()
{

}
//___________________________________________________________________________
void FLUKA::ProcessEventRecord(GHepRecord * event) const
{
  //-- Check that we have an interaction with a nuclear target.
  //   If not, do not call FLUKA
  Interaction * interaction = event->GetInteraction();
  const InitialState & init_state = interaction->GetInitialState();

  if( ! init_state.GetTarget().IsNucleus() ) {
    LOG("FLUKA", pDEBUG)
        << "Not an interaction with a nuclear target. Skipping call FLUKA";
    return;
  }

  //-- Access this module's configuration options and pass them to
  //   the actual FLUKA code

  //
  // ...

  //-- Translate GENIE GHepRecord to whatever FLUKA needs
  LOG("FLUKA", pDEBUG) << "Translating: GENIE GHepRecord ---> FLUKA input";

  //
  // ...

  //-- Run FLUKA
  LOG("FLUKA", pINFO) << "************ Running FLUKA ************";

  //
  // ...

  //-- Get FLUKA output & add it to GENIE GHepRecord
  LOG("FLUKA", pDEBUG) << "Copying: FLUKA output ---> GENIE GHepRecord";

  //
  // ...

}
//___________________________________________________________________________
