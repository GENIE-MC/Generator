#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_INCL_ENABLED__

#ifndef INCLConvertParticle_hh
#define INCLConvertParticle_hh 1

// INCL++
#include "G4INCLParticleSpecies.hh" // includes G4INCLParticleType
#include "G4INCLParticleTable.hh"

// GENIE
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"

using namespace genie;

namespace G4INCL {

  int INCLpartycleSpecietoPDGCODE(G4INCL::ParticleSpecies theSpecies) {
    if (theSpecies.theType != Composite) {
      if      (ParticleTable::getName(theSpecies.theType) == "pi0") return  111;
      else if (ParticleTable::getName(theSpecies.theType) == "pi+") return  211;
      else if (ParticleTable::getName(theSpecies.theType) == "pi-") return -211;
    }
    if      (theSpecies.theA == 1 && theSpecies.theZ == 1) return 2212;
    else if (theSpecies.theA == 1 && theSpecies.theZ == 0) return 2112;
    return 0; // return _something_
  }

  int INCLtopdgcode(int A, int Z) {
    if      (A == 1 && Z == 1) return 2212;
    else if (A == 1 && Z == 0) return 2112;
    else if (A == 0 && Z == 0) return  111;
    else if (A == 0 && Z == 1) return  211;
    else if (A == 0 && Z ==-1) return -211;
    else {
      return genie::pdg::IonPdgCode( A , Z );
    }
  }

  int FindlargestFragment(G4INCL::EventInfo result){
    int maxA(0), rem_index(0);
    for (size_t ll = 0; ll < result.nParticles; ll++) {
      if (result.A[ll] > maxA) {
        maxA = result.A[ll];
        rem_index = ll;
      }
    }
      return rem_index;
    }

    G4INCL::ParticleType toINCLparticletype(int pdgc) {

      if      (pdgc ==       2212) return G4INCL::Proton;
      else if (pdgc ==       2112) return G4INCL::Neutron;
      else if (pdgc ==        211) return G4INCL::PiPlus;
      else if (pdgc ==       -211) return G4INCL::PiMinus;
      else if (pdgc ==        111) return G4INCL::PiZero;
      else if (pdgc == 1000010020) return G4INCL::Composite;
      else if (pdgc == 1000010030) return G4INCL::Composite;
      else if (pdgc == 1000020030) return G4INCL::Composite;
      else if (pdgc == 1000020040) return G4INCL::Composite;
      else                         return G4INCL::UnknownParticle;

    }

    GHepParticle *INCLtoGenieParticle(G4INCL::EventInfo result,
      int nP, GHepStatus_t ist, int mom1, int mom2) {
      int pdg_codeP(0);
      double E_pnP(0), EKinP = result.EKin[nP];
      TVector3 p3M(result.px[nP]/1000,result.py[nP]/1000,result.pz[nP]/1000);
      TLorentzVector x4null(0.,0.,0.,0.);
      double  m_pnP = 0.5*((result.px[nP])*(result.px[nP]) +
       (result.py[nP])*(result.py[nP]) +
       (result.pz[nP])*(result.pz[nP]) -
       EKinP*EKinP) / (EKinP);
    if (m_pnP < 10 && result.A[nP] == 0 && result.Z[nP] == 0) { // photon
      pdg_codeP = 22;
      E_pnP = TMath::Sqrt((result.px[nP])*(result.px[nP]) +
        (result.py[nP])*(result.py[nP]) +
        (result.pz[nP])*(result.pz[nP])   );
    } else {
      pdg_codeP = INCLtopdgcode(result.A[nP],result.Z[nP]);
      TParticlePDG * pdgRemn=PDGLibrary::Instance()->Find(pdg_codeP,false); 
      if(!pdgRemn)
      {
        LOG("HINCLCascadeIntranuke", pINFO)
        << "NO Particle with pdg = " << pdg_codeP << " in PDGLibrary!";
        pdg_codeP=kPdgHadronicBlob;
        ist=kIStFinalStateNuclearRemnant;
      }
      E_pnP = EKinP + m_pnP;
    }
    TLorentzVector p4tgtf(p3M,E_pnP/1000);
    return new GHepParticle(pdg_codeP,ist,mom1,mom2,-1,-1,p4tgtf,x4null);
  }

}

#endif

#endif // __GENIE_INCL_ENABLED__
