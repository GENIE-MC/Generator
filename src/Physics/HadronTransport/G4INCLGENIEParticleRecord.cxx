#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_INCL_ENABLED__

#include "Physics/HadronTransport/G4INCLGENIEParticleRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"


namespace G4INCL{
  GENIEParticleRecord::GENIEParticleRecord(genie::GHepParticle *p, int scType, GENIERecordCode recordCode){
    fPdgCode        = p->Pdg();
    fCharge         = p->Charge()/3;
    // if(p->Status() == genie::kIStHadronInTheNucleus) fStatus = 1;
    // else fStatus = 0;
    fStatus         = p->Status();
    fFirstMother    = p->FirstMother();
    fLastMother     = p->LastMother();
    fFirstDaughter  = p->FirstDaughter();
    fLastDaughter   = p->LastDaughter();
    fScatteringType = scType;
    fPType          = this->PDG_to_INCLType(p->Pdg());
    fMass           = p->P4()->M() * GeV / MeV;
    fP3.set(p->P4()->Px() * GeV / MeV, p->P4()->Py() * GeV / MeV, p->P4()->Pz() * GeV / MeV);
    fX3.set(p->X4()->X(), p->X4()->Y(), p->X4()->Z());
    fID = -1;
    fGenieRecordCode = recordCode;


  }

  GENIEParticleRecord::~GENIEParticleRecord(){

  }
  void GENIEParticleRecord::setID(int id){
    fID = id;
  }

  int GENIEINCLUtil::INCLPDG_to_GHEPPDG(int pdg, int A, int Z, int S){
    // Convert the INCL pdg id to GHep pdg id
    TDatabasePDG * fDatabasePDG = TDatabasePDG::Instance();
    TParticlePDG * p = fDatabasePDG->GetParticle(pdg);
    // pdg code is not in the data base
    if(!p){     
      // It is a nucleus, using A, Z, S to construct the nucleus ID.
      if(A != 0){  
        int ion_pdg =  genie::pdg::IonPdgCode( A , Z, std::abs(S), 0 );
        TParticlePDG *ion = fDatabasePDG->GetParticle(ion_pdg);
        // the hyper-nuclei is not in the database, add a temp one 
        if(S != 0 && !ion){ 
          genie::PDGLibrary *pdg_library = genie::PDGLibrary::Instance();
          pdg_library->AddHypernucleus(ion_pdg);
        }
        return ion_pdg;
      }
      else{
        exit(1);
      }
    }
    else{ // find the pdg code in database, return it.
      return pdg;
    }
  }
}
#endif // __GENIE_INCL_ENABLED__
