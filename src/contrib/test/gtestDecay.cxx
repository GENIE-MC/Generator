//____________________________________________________________________________
/*!

\program gtestDecay

\brief   test program used for testing / debugging the GENIE particle decayers

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created June 20, 2004

\cpright Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <ostream>
#include <iomanip>

#include <RVersion.h>
#include <TClonesArray.h>
#include <TParticlePDG.h>
#include <TIterator.h>

#include "Algorithm/Algorithm.h"
#include "Algorithm/AlgFactory.h"
#include "Conventions/Units.h"
#include "Conventions/Constants.h"
#include "Decay/Decayer.h"
#include "Decay/PythiaDecayer.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGLibrary.h"

using namespace genie;

using std::ostream;
using std::endl;
using std::setw;
using std::setprecision;
using std::setfill;
using std::ios;

ostream & operator<< (ostream & stream, const TClonesArray * particle_list);
ostream & operator<< (ostream & stream, const GHepParticle * particle);

void TestPythiaTauDecays(void);
void Decay(const Decayer * decayer, int pdgc, double E,  int ndecays);

//__________________________________________________________________________
int main(int /*argc*/, char ** /*argv*/)
{
  TestPythiaTauDecays();

  return 0;
}
//__________________________________________________________________________
void TestPythiaTauDecays(void)
{
  // Get the pythia decayer
  LOG("test",pINFO)
     << "Asking the AlgFactory for a genie::PythiaDecayer\\Default instance";
  AlgFactory * algf = AlgFactory::Instance();
  const Decayer * pdecayer =
     dynamic_cast<const Decayer *> (
         algf->GetAlgorithm("genie::PythiaDecayer","Default"));

  // Decayer config print-out
  LOG("test",pINFO) << "Algorithm name = " << pdecayer->Id().Name();
  LOG("test",pINFO) << "Parameter set  = " << pdecayer->Id().Config();

  const Registry & conf_registry = pdecayer->GetConfig();
  LOG("test", pINFO) << conf_registry;

  double E = 10;
  int ndec = 3;

  // inhibit pi0 decay
  pdecayer->InhibitDecay(kPdgPi0);

  // decay tau-
  Decay(pdecayer, kPdgTau, E, ndec);

  // now inhibit all but the tau- --> nu_mu_bar + mu- + nu_tau decay channel
  // (see $GENIE/src/contrib/misc/print_decay_channels.C)
  // and perform some more decays
  LOG("test",pINFO) 
    << "\n\n"
    << " **** Inhibiting all but the `tau- --> nu_mu_bar + mu- + nu_tau' decay channel"
    << "\n\n";

  PDGLibrary * pdglib = PDGLibrary::Instance();

  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(0));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(2));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(3));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(4));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(5));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(6));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(7));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(8));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(9));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(10));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(11));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(12));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(13));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(14));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(15));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(16));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(17));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(18));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(19));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(20));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(21));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(22));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(23));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(24));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(25));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(26));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(27));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(28));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(29));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(30));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(31));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(32));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(33));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(34));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(35));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(36));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(37));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(38));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(39));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(40));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(41));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(42));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(43));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(44));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(45));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(46));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(47));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(48));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(49));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(50));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(51));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(52));
  pdecayer->InhibitDecay(kPdgTau, pdglib->Find(kPdgTau)->DecayChannel(53));

  // a few more decays
  Decay(pdecayer, kPdgTau, E, ndec);

  // restore all decay channels
  LOG("test",pINFO) 
    << "\n\n"
    << " **** Restoring all tau- decay channels"
    << "\n\n";

  pdecayer->UnInhibitDecay(kPdgTau);

  // a few more decays
  Decay(pdecayer, kPdgTau, E, ndec);

  // now inhibit all decay channels and try to decay!
  LOG("test",pINFO) 
    << "\n\n"
    << " **** Inhibit all tau- decay channels"
    << "\n\n";

  pdecayer->InhibitDecay(kPdgTau);

  // a few more decays
  Decay(pdecayer, kPdgTau, E, ndec);
}
//__________________________________________________________________________
void Decay(const Decayer * decayer, int pdgc, double E, int ndecays)
{
  DecayerInputs_t dinp;

  TLorentzVector p4;
  p4.SetE(E);
  p4.SetTheta(0.);
  p4.SetPhi(0.);

  dinp.PdgCode = pdgc;
  dinp.P4      = &p4;

  PDGLibrary * pdglib = PDGLibrary::Instance();
  TParticlePDG * pp = pdglib->Find(pdgc);
  if(!pp) return;

  LOG("test",pINFO) 
       << "Decaying a " << pp->GetName() << " with E = " << p4.Energy();

  // Perform the decay a few times & print-out the decay products
  for(int idec = 0; idec < ndecays; idec++) {

        // Decay
  	TClonesArray * particle_list = decayer->Decay(dinp); 

        if(!particle_list) {
          LOG("test",pWARN) 
               << "\n ** Decay nu.: " << idec << " ==> NULL particle list";
           continue;
        }

        // Print-out       
        LOG("test",pINFO) 
             << "\n ** Decay nu.: " << idec 
             << " (weight = " << decayer->Weight() << ") : \n "
             << particle_list;

        // Clean-up
        particle_list->Delete();
        delete particle_list;
  }//ndecays

}
//__________________________________________________________________________
ostream & operator<< (ostream & stream, const TClonesArray * particle_list)
{
  GHepParticle * p = 0;
  TObjArrayIter particle_iter(particle_list);


  stream 
    << setfill(' ') << setw(10) << "name "
    << setfill(' ') << setw(10)  << "PDG"
    << setfill(' ') << setw(10)  << "Status"
    << setfill(' ') << setw(15) << "E (GeV)"
    << setfill(' ') << setw(15) << "Px (GeV/c)"
    << setfill(' ') << setw(15) << "Py (GeV/c)"
    << setfill(' ') << setw(15) << "Pz (GeV/c)"
    << setfill(' ') << setw(15) << "t (mm/c)"
    << setfill(' ') << setw(15) << "x (mm)"
    << setfill(' ') << setw(15) << "y (mm)"
    << setfill(' ') << setw(15) << "z (mm)"
    << endl;

  while( (p = (GHepParticle *) particle_iter.Next()) ) stream << p;

  stream << setfill('-') << setw(100) << "|";

  return stream; 
}
//__________________________________________________________________________
ostream & operator<< (ostream & stream, const GHepParticle * p)
{
  stream 
   << std::scientific << setprecision(6);

  stream 
    << setfill(' ') << setw(10) << p->Name()
    << setfill(' ') << setw(10) << p->Pdg() 
    << setfill(' ') << setw(10) << p->Status()
    << setfill(' ') << setw(15) << p->Energy()
    << setfill(' ') << setw(15) << p->Px() 
    << setfill(' ') << setw(15) << p->Py() 
    << setfill(' ') << setw(15) << p->Pz() 
    << setfill(' ') << setw(15) << p->Vt()   /(units::mm) 
    << setfill(' ') << setw(15) << p->Vx()   /(units::mm) 
    << setfill(' ') << setw(15) << p->Vy()   /(units::mm) 
    << setfill(' ') << setw(15) << p->Vz()   /(units::mm) 
    << endl;

  return stream;
}
//__________________________________________________________________________

