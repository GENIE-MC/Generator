//____________________________________________________________________________
/*!

\program gtestInteraction

\brief   Program used for testing / debugging the Interaction and its aggregate 
         objects (InitialState, ProcessInfo, Kinematics, XclsTag)

\author  Costas Andreopoulos <C.V.Andreopoulos@@rl.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created May 4, 2004

\cpright Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include "Conventions/Constants.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"

using namespace genie;
using namespace genie::constants;

int main(int /*argc*/, char ** /*argv*/)
{
  //-- build an initial state
  Target         Fe(26,56);
  int            numu(kPdgNuMu); 
  TLorentzVector pnu(0,0,8,8);
  InitialState   init_state(Fe, numu);

  init_state.SetProbeP4(pnu);

  //-- build process info
  ProcessInfo proc(kScDeepInelastic, kIntWeakCC);

  //-- create an interaction & print it
  Interaction interaction(init_state, proc);

  LOG("test", pINFO) << "Printing an interaction object";
  LOG("test", pINFO) << interaction;

  //-- set struck nucleon & quark info in the initial state's Target
  //   note: using methods ending in Ptr -> they return a 'writable'
  //         object that can be modified

  TLorentzVector pnucl(0,0,0,kNucleonMass);

  interaction.InitStatePtr()->TgtPtr()->SetHitNucP4(pnucl);
  interaction.InitStatePtr()->TgtPtr()->SetHitNucPdg(kPdgProton);
  interaction.InitStatePtr()->TgtPtr()->SetHitQrkPdg(kPdgUQuark);

  //-- get a 'read-only' InitialState and print it (check that struck nucleon
  //   and quark were set) 
  //   note: using the methods not ending in Ptr to get a 'read-only' object

  const InitialState & cinit = interaction.InitState();

  LOG("test", pINFO) << "Printing initial state after setting struck nucl/qrk";
  LOG("test", pINFO) << "\n" << cinit;

  //-- take just the Target from the initial state and print it

  LOG("test", pINFO) << "Printing target after setting struck nucl/qrk";
  const Target & ctgt = interaction.InitState().Tgt();
  LOG("test", pINFO) << "\n" <<ctgt;

  //-- change the struck nucleon
  //-- instead of using the long syntax above, get a writable Target object first

  Target * wtgt = interaction.InitStatePtr()->TgtPtr();
  wtgt->SetHitNucPdg(kPdgProton);

  LOG("test", pINFO) << "Printing target after changing struck nucl";
  LOG("test", pINFO) << "\n" << *wtgt;

  //-- take the Kinematics object and set some

  Kinematics * wkine = interaction.KinePtr();
 
  wkine->Setx(0.1781);
  wkine->Sety(0.6892);
  wkine->SetQ2(3.2218);

  LOG("test", pINFO) << "Printing kinematics after setting x,y,Q2";
  LOG("test", pINFO) << "\n" <<*wkine;

  //-- modify some & add a new
  wkine->Setx(0.2);
  wkine->Sety(0.2);
  wkine->SetW(1.9219);

  LOG("test", pINFO) << "Printing kinematics after modifying x,y and adding W";
  LOG("test", pINFO) << "\n" << *wkine;

  //-- now set "selected" kinematics
  wkine->Setx(0.25, true);
  wkine->Sety(0.21, true);
  wkine->SetW(2.89, true);

  LOG("test", pINFO) << "Printing kinematics after setting 'selected'";
  LOG("test", pINFO) << "\n" << *wkine;

  //-- now delete 'running' kinematics
  wkine->ClearRunningValues();

  LOG("test", pINFO) << "Printing kinematics after deleting 'running'";
  LOG("test", pINFO) << "\n" << *wkine;

  //-- copy the 'selected' kinematics to the 'running' ones
  wkine->UseSelectedKinematics();

  LOG("test", pINFO) << "Printing kinematics after copying 'selected'";
  LOG("test", pINFO) << "\n" << *wkine;

  //-- see that the interaction was updated

  LOG("test", pINFO) << "Printing the interaction object after all changes";
  LOG("test", pINFO) << interaction;
}

