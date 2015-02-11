//____________________________________________________________________________
/*!

\program gtestMuELoss

\brief   Program for the MuELoss utility package. The program saves the
         computed data in an output ROOT ntuple. 
        
         Syntax :
           gtestMuELoss -m materials

         Options :
           -m specifies a comma seperated list of materials 
              (the material ids correspond to the enumeration that can be seen
               in $GENIE/src/MuELoss/MuELMaterial.h)

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created March 10, 2006

\cpright Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <cassert>
#include <cstdlib>
#include <string>
#include <vector>

#include <TFile.h>
#include <TNtuple.h>

#include "Algorithm/AlgFactory.h"
#include "Conventions/Units.h"
#include "MuELoss/MuELossI.h"
#include "MuELoss/MuELMaterial.h"
#include "MuELoss/MuELProcess.h"
#include "Messenger/Messenger.h"
#include "Utils/StringUtils.h"
#include "Utils/CmdLnArgParser.h"

using std::string;
using std::vector;
using namespace genie;
using namespace genie::utils;
using namespace genie::mueloss;

void GetCommandLineArgs (int argc, char ** argv);

string gOptMaterials;

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs(argc, argv);

  const int N = 14;
  double E[N] = {1,5,10,15,20,30,50,100,200,500, 1000,2000, 5000, 9000}; //GeV

  // split the comma separated list of materials
  vector<string> mtv = utils::str::Split(gOptMaterials,  ",");
  
  // get muon energy loss algorithms
  AlgFactory * algf = AlgFactory::Instance();

  const MuELossI * betheBloch = 
         dynamic_cast<const MuELossI *> (algf->GetAlgorithm(
                         "genie::mueloss::BetheBlochModel","Default"));

  const MuELossI * petrukhinShestakov = 
         dynamic_cast<const MuELossI *> (algf->GetAlgorithm(
                 "genie::mueloss::PetrukhinShestakovModel","Default"));

  const MuELossI * kokoulinPetroukhin = 
         dynamic_cast<const MuELossI *> (algf->GetAlgorithm(
                  "genie::mueloss::KokoulinPetrukhinModel","Default"));

  const MuELossI * bezroukovBugaev = 
         dynamic_cast<const MuELossI *> (algf->GetAlgorithm(
                     "genie::mueloss::BezrukovBugaevModel","Default"));

  assert ( betheBloch         );
  assert ( petrukhinShestakov );
  assert ( kokoulinPetroukhin );
  assert ( bezroukovBugaev    );

  double myunits_conversion = units::GeV/(units::g/units::cm2); 
  string myunits_name       = " GeV/(gr/cm^2)";

  // open a ROOT file and define the output ntuple.
  TFile froot("./genie-mueloss.root", "RECREATE");
  TNtuple muntp("muntp","muon dE/dx", "material:E:ion:brem:pair:pnucl");
  
  //loop over materials
  vector<string>::iterator iter;
  for(iter = mtv.begin(); iter != mtv.end(); ++iter) {

    MuELMaterial_t mt = (MuELMaterial_t) atoi(iter->c_str());

     LOG("test", pINFO) 
        << "---------- Computing/Printing muon energy losses in "
                       << MuELMaterial::AsString(mt) << " ----------";

     // loop over energies
     for(int i=0; i<N; i++)  {
      
       // ionization 
       double ion = betheBloch->dE_dx(E[i],mt) / myunits_conversion;

       LOG("test", pINFO) 
         << "Process: " << MuELProcess::AsString(betheBloch->Process())
         << ", Model: " << betheBloch->Id().Key() 
         << " : \n -dE/dx(E=" << E[i] << ") = " << ion << myunits_name;

       // bremsstrahlung
       double brem = petrukhinShestakov->dE_dx(E[i],mt) / myunits_conversion;

       LOG("test", pINFO) 
         << "Process: " << MuELProcess::AsString(petrukhinShestakov->Process())
         << ", Model: " << petrukhinShestakov->Id().Key() 
         << " : \n -dE/dx(E=" << E[i] << ") = " << brem << myunits_name;

       // e-e+ pair production
       double pair = kokoulinPetroukhin->dE_dx(E[i],mt) / myunits_conversion;

       LOG("test", pINFO) 
         << "Process: " << MuELProcess::AsString(kokoulinPetroukhin->Process())
         << ", Model: " << kokoulinPetroukhin->Id().Key() 
         << " : \n -dE/dx(E=" << E[i] << ") = " << pair << myunits_name;

       // photonuclear interactions
       double pnucl = bezroukovBugaev->dE_dx(E[i],mt) / myunits_conversion;

       LOG("test", pINFO) 
         << "Process: " << MuELProcess::AsString(bezroukovBugaev->Process())
         << ", Model: " << bezroukovBugaev->Id().Key() 
         << " : \n -dE/dx(E=" << E[i] << ") = " << pnucl << myunits_name
         << "\n\n";

       muntp.Fill( (int)mt,E[i],ion,brem,pair,pnucl);
    }//e
  }//m

  muntp.Write();
  froot.Close();

  return 0;
}
//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("test", pNOTICE) << "Parsing command line arguments";

  CmdLnArgParser parser(argc,argv);

  if ( parser.OptionExists('m') ) {
    LOG("test", pINFO) << "Reading material ids";
    gOptMaterials = parser.ArgAsString('m');
  } else {
    LOG("test", pINFO) << "Unspecified material ids - Exiting";
    exit(1);
  }
}
//____________________________________________________________________________

