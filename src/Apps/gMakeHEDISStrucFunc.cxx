//____________________________________________________________________________
/*!
\program gmkspl
\brief   GENIE utility program building Structure Functions needed for HEDIS
         package.
         Syntax :
           gmksf [-h]
                  --tune genie_tune
                 [--message-thresholds xml_file]
         Note :
           [] marks optional arguments.
           <> marks a list of arguments out of which only one can be
              selected at any given time.
         Options :
           --tune
              Specifies a GENIE comprehensive neutrino interaction model tune.
           --message-thresholds
              Allows users to customize the message stream thresholds.
              The thresholds are specified using an XML file.
              See $GENIE/config/Messenger.xml for the XML schema.
        ***  See the User Manual for more details and examples. ***
\author  Alfonso Garcia <alfonsog \at nikhef.nl>
         NIKHEF (Amsterdam)
\created July 8, 2019
\cpright Copyright (c) 2003-2025, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/EventGen/InteractionList.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/GEVGDriver.h"
#include "Framework/Utils/AppInit.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Interaction/InitialState.h"
#include "Framework/Messenger/Messenger.h"

using namespace genie;

void   PrintSyntax        (void);

//____________________________________________________________________________
int main(int argc, char ** argv)
{

  RunOpt::Instance()->ReadFromCommandLine(argc,argv);

  if ( ! RunOpt::Instance()->Tune() ) {
    LOG("gmkhedissf", pFATAL) << " No TuneId in RunOption";
    exit(-1);
  }
  RunOpt::Instance()->BuildTune();

  utils::app_init::MesgThresholds(RunOpt::Instance()->MesgThresholdFiles());

  GEVGDriver evg_driver;
  InitialState init_state(1000010020, 14);
  evg_driver.SetEventGeneratorList("CCHEDIS");
  evg_driver.Configure(init_state);

  const InteractionList * intlst   = evg_driver.Interactions();
  InteractionList::const_iterator intliter = intlst->begin();
  Interaction * interaction = *intliter;
  const XSecAlgorithmI * xsec_alg = evg_driver.FindGenerator(interaction)->CrossSectionAlg();
  interaction->SetBit(kISkipKinematicChk);
  xsec_alg->XSec(interaction, kPSxQ2fE);

}//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gmkhedissf", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "\n      gmkhedissf [-h]"
    << "\n                  --tune genie_tune"
    << "\n                  --message-thresholds xml_file"
    << "\n";
}
