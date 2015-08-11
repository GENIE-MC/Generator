//____________________________________________________________________________
/*!

\program gtestFermiP

\brief   Program used for testing / debugging the Fermi momentum distribution
         models

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created June 20, 2004

\cpright Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <TFile.h>
#include <TNtuple.h>
#include <TMath.h>
#include <TVector3.h>

#include "Algorithm/AlgFactory.h"
#include "Interaction/Target.h"
#include "Messenger/Messenger.h"
#include "Nuclear/NuclearModelI.h"
#include "PDG/PDGCodes.h"
#include "Utils/PrintUtils.h"

using namespace genie;

int main(int /*argc*/, char ** /*argv*/)
{
  const unsigned int kNTargets = 2;
  const unsigned int kNModels  = 2;  // 4 if Benhar models are available
  const unsigned int kNEvents  = 3000;

  //-- Get nuclear models
  AlgFactory * algf = AlgFactory::Instance();

  const NuclearModelI * bodritch =
       dynamic_cast<const NuclearModelI *> (algf->GetAlgorithm(
			     "genie::FGMBodekRitchie","Default"));
  // const NuclearModelI * benhsf1d = 
  //      dynamic_cast<const NuclearModelI *> (
  //             algf->GetAlgorithm("genie::BenharSpectralFunc1D","Default"));
  // const NuclearModelI * benhsf2d = 
  //      dynamic_cast<const NuclearModelI *> (
  //                algf->GetAlgorithm("genie::BenharSpectralFunc","Default"));
  const NuclearModelI * effsf = 
       dynamic_cast<const NuclearModelI *> (
                 algf->GetAlgorithm("genie::EffectiveSF","Default"));

  const NuclearModelI * nuclmodel[kNModels] = { bodritch, /* benhsf1d, benhsf2d, */ effsf };

  //-- Create nuclear targets
  Target * nucltgt[kNTargets];
  nucltgt[0] = new Target ( 6,  12); // C12
  nucltgt[1] = new Target (26,  56); // Fe56

  //-- Output ntuple
  
  TNtuple * nuclnt = new TNtuple("nuclnt","","target:model:p:px:py:pz:w");

  //-- Loop over targets/models & generate 'target nucleons'
  for(unsigned int it = 0; it < kNTargets; it++) {
     nucltgt[it]->SetHitNucPdg(kPdgProton);
     const Target & target = *nucltgt[it];
     LOG("test", pNOTICE)  << "** Using target : " << target;;
     for(unsigned int im = 0; im < kNModels; im++) {
        LOG("test", pNOTICE) << "Running model : " << nuclmodel[im]->Id();
        for(unsigned int iev = 0; iev < kNEvents; iev++) {
            nuclmodel[im]->GenerateNucleon(target);
            double   w  = nuclmodel[im]->RemovalEnergy();
            TVector3 p3 = nuclmodel[im]->Momentum3();
            double   px = p3.Px();
            double   py = p3.Py();
            double   pz = p3.Pz();
            double   p  = p3.Mag();
            LOG("test", pDEBUG)
                << "Nucleon 4-P = " << utils::print::Vec3AsString(&p3);
            nuclnt->Fill(it,im,p,px,py,pz,w);
        }//ievents
     }//immodels
  }//itargets

  //-- Save ntuple
  TFile f("./fermip.root","recreate");
  nuclnt->Write();
  f.Close();

  //-- Clean-up
  delete nuclnt;

  return 0;
}
//____________________________________________________________________________
