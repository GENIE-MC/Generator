#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_INCL_ENABLED__

//---------------------
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

// GENIE
#include "HINCLCascadeIntranuke.h"

// INCL++
#include "G4INCLConfig.hh"
#include "G4INCLCascade.hh"
#include "G4INCLConfigEnums.hh"
#include "G4INCLParticle.hh"
// signal handler (for Linux and GCC)
#include "G4INCLSignalHandling.hh"

// Generic de-excitation interface
#include "G4INCLIDeExcitation.hh"

// ABLA v3p de-excitation
#ifdef INCL_DEEXCITATION_ABLAXX
#include "G4INCLAblaInterface.hh"
#endif

// ABLA07 de-excitation
#ifdef INCL_DEEXCITATION_ABLA07
#include "G4INCLAbla07Interface.hh"
#endif

// SMM de-excitation
#ifdef INCL_DEEXCITATION_SMM
#include "G4INCLSMMInterface.hh"
#endif

// GEMINIXX de-excitation
#ifdef INCL_DEEXCITATION_GEMINIXX
#include "G4INCLGEMINIXXInterface.hh"
#endif

//#ifdef HAS_BOOST_DATE_TIME
//#include <boost/date_time/posix_time/posix_time.hpp>
//namespace bpt = boost::posix_time;
//#endif

#ifdef HAS_BOOST_TIMER
#include <boost/timer/timer.hpp>
namespace bt = boost::timer;
#endif


// --------------------------------------Include for GENIE---------------------
// GENIE
#include "INCLConvertParticle.hh"
#include "INCLConfigParser.h"
// INCL++
//#include "ConfigParser.hh"

#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/SystemUtils.h"

// ROOT
#include "TSystem.h"

using namespace genie;
using namespace genie::utils;
using namespace G4INCL;
using std::ostringstream;
using namespace std;

HINCLCascadeIntranuke::HINCLCascadeIntranuke() :
EventRecordVisitorI("genie::HINCLCascadeIntranuke"),
theINCLConfig(0), theINCLModel(0), theDeExcitation(0)
{
  LOG("HINCLCascadeIntranuke", pDEBUG)
  << "default ctor";
}

//______________________________________________________________________________
HINCLCascadeIntranuke::HINCLCascadeIntranuke(string config) :
EventRecordVisitorI("genie::HINCLCascadeIntranuke", config),
theINCLConfig(0), theINCLModel(0), theDeExcitation(0)
{
  LOG("HINCLCascadeIntranuke", pDEBUG)
  << "ctor from config " << config;
}

//______________________________________________________________________________
HINCLCascadeIntranuke::~HINCLCascadeIntranuke()
{

  // Config is owned by model once handed over
  if ( theINCLConfig   ) { theINCLConfig=0;   }
  if ( theINCLModel    ) { delete theINCLModel;    theINCLModel=0;    }
  if ( theDeExcitation ) { delete theDeExcitation; theDeExcitation=0; }

}

//______________________________________________________________________________
void HINCLCascadeIntranuke::LoadConfig(void)
{
  LOG("HINCLCascadeIntranuke", pINFO) << "Settings for INCL++ mode: " ;

#ifdef INCL_SIGNAL_HANDLING
  enableSignalHandling();
#endif

  // Config is owned by model once handed over
  if ( theINCLConfig   ) { theINCLConfig=0;   }
  if ( theINCLModel    ) { delete theINCLModel;    theINCLModel=0;    }
  if ( theDeExcitation ) { delete theDeExcitation; theDeExcitation=0; }

  INCLConfigParser theParser;

  size_t maxFlags = 200;
  size_t nflags   = 0;
  char * flags[maxFlags];  // note non-const'ness desired by ConfigParser

  std::string infile;
  GetParamDef( "INCL-infile", infile, std::string("NULL"));
  flags[nflags] = strdup(infile.c_str()); ++nflags;

  std::string pflag;
  GetParamDef( "INCL-pflag",  pflag,  std::string("-pp"));
  flags[nflags] = strdup(pflag.c_str()); ++nflags;

  std::string tflag;
  GetParamDef( "INCL-tflag",  tflag,  std::string("-tFe56"));
  flags[nflags] = strdup(tflag.c_str()); ++nflags;

  std::string Nflag;
  GetParamDef( "INCL-Nflag",  Nflag,  std::string("-N1"));
  flags[nflags] = strdup(Nflag.c_str()); ++nflags;

  std::string Eflag;
  GetParamDef( "INCL-Eflag",  Eflag,  std::string("-E1"));
  flags[nflags] = strdup(Eflag.c_str()); ++nflags;

  std::string dflag;
  GetParamDef( "INCL-dflag",  dflag,  std::string("-dABLA07"));
  flags[nflags] = strdup(dflag.c_str()); ++nflags;

  // arbitary extra flags, need to be tokenized
  std::string extra_incl_flags;
  GetParamDef( "INCL-extra-flags",  extra_incl_flags,  std::string(""));
  std::vector<std::string> extraTokens = genie::utils::str::Split(extra_incl_flags," \t\n");
  for (size_t j=0; j < extraTokens.size(); ++j) {
    std::string& token = extraTokens[j];
    if ( token != "" ) {
      // don't add empty strings
      flags[nflags] = strdup(token.c_str()); ++nflags;
    }
  }

  LOG("HINCLCascadeIntranuke", pDEBUG)
  << "LoadConfig() create theINCLConfig";
  theINCLConfig = theParser.parse(nflags,flags);

  // there's currently no way to update *all* of the data paths in the Config
  // after it's been generated, so check whether default file paths "work",
  // add to the flags if necessary, and then regenerate
  bool updateNeeded = AddDataPathFlags(nflags,flags);
  if (updateNeeded) {
    delete theINCLConfig;
    theINCLConfig = theParser.parse(nflags,flags);
  }

  LOG("HINCLCascadeIntranuke", pDEBUG)
  << "doCascade new G4INCL::INCL";
  theINCLModel = new G4INCL::INCL(theINCLConfig);

  // code copied fomr INCL's main/src/INCLCascade.cc
  // with additions for GENIE messenger system
  switch(theINCLConfig->getDeExcitationType()) {
#ifdef INCL_DEEXCITATION_ABLAXX
  case G4INCL::DeExcitationABLAv3p:
    theDeExcitation = new G4INCLAblaInterface(theINCLConfig);
    LOG("HINCLCascadeIntranuke", pINFO)
    << "using ABLAv3p for DeExcitation";
    break;
#endif
#ifdef INCL_DEEXCITATION_ABLA07
  case G4INCL::DeExcitationABLA07:
    theDeExcitation = new ABLA07CXX::Abla07Interface(theINCLConfig);
    LOG("HINCLCascadeIntranuke", pINFO)
    << "using ABLA07CXX for DeExcitation";
    break;
#endif
#ifdef INCL_DEEXCITATION_SMM
  case G4INCL::DeExcitationSMM:
    theDeExcitation = new SMMCXX::SMMInterface(theINCLConfig);
    LOG("HINCLCascadeIntranuke", pINFO)
    << "using SMMCXX for DeExcitation";
    break;
#endif
#ifdef INCL_DEEXCITATION_GEMINIXX
  case G4INCL::DeExcitationGEMINIXX:
    theDeExcitation = new G4INCLGEMINIXXInterface(theINCLConfig);
    LOG("HINCLCascadeIntranuke", pINFO)
    << "using GEMINIXX for DeExcitation";
    break;
#endif
  default:
    std::stringstream ss;
    ss << "########################################################\n"
    << "###              WARNING WARNING WARNING             ###\n"
    << "###                                                  ###\n"
    << "### You are running the code without any coupling to ###\n"
    << "###              a de-excitation model!              ###\n"
    << "###    Results will be INCOMPLETE and UNPHYSICAL!    ###\n"
    << "###    Are you sure this is what you want to do?     ###\n"
    << "########################################################\n";
    INCL_INFO(ss.str());
    LOG("HINCLCascadeIntranuke", pWARN)
    << '\n' << ss.str();
    //std::cout << ss.str();
    break;
  }

  // try not to leak strings
  for (size_t i=0; i < nflags; ++i) {
    char * p = flags[i];
    free(p);
    flags[i] = 0;
  }

}

//______________________________________________________________________________
bool HINCLCascadeIntranuke::AddDataPathFlags(size_t& nflags, char** flags) {

  // check if the default INCL path works
  bool needed_update = false;
  size_t nflags_in = nflags;

  std::string validpath;
  std::vector<std::string> datapaths;

  LOG("HINCLCascadeIntranuke", pINFO)
  << "check for data location of INCL";

  datapaths.clear();
  datapaths.push_back(theINCLConfig->getINCLXXDataFilePath());
  datapaths.push_back("!INCL-incl-data-alt1");
  datapaths.push_back("!INCL-incl-data-alt2");
  needed_update |=
  LookForAndAddValidPath(datapaths,0,"--inclxx-datafile-path",nflags,flags);

  switch(theINCLConfig->getDeExcitationType()) {
#ifdef INCL_DEEXCITATION_ABLAXX
  case G4INCL::DeExcitationABLAv3p:
    LOG("HINCLCascadeIntranuke", pINFO)
    << "using ABLAv3p for DeExcitation -- check for data location";
    datapaths.clear();
    datapaths.push_back(theINCLConfig->getABLAv3pCxxDataFilePath());
    datapaths.push_back("!INCL-ablav3p-data-alt1");
    datapaths.push_back("!INCL-ablav3p-data-alt2");
    needed_update |=
    LookForAndAddValidPath(datapaths,0,"--ablav3p-cxx-datafile-path",nflags,flags);
    break;
#endif
#ifdef INCL_DEEXCITATION_ABLA07
  case G4INCL::DeExcitationABLA07:
    LOG("HINCLCascadeIntranuke", pINFO)
    << "using ABLA07 for DeExcitation -- check for data location";
    datapaths.clear();
    datapaths.push_back(theINCLConfig->getABLA07DataFilePath());
    datapaths.push_back("!INCL-abla07-data-alt1");
    datapaths.push_back("!INCL-abla07-data-alt2");
    needed_update |=
    LookForAndAddValidPath(datapaths,0,"--abla07-datafile-path",nflags,flags);
    break;
#endif
#ifdef INCL_DEEXCITATION_SMM
  case G4INCL::DeExcitationSMM:
    LOG("HINCLCascadeIntranuke", pINFO)
    << "using SMMCXX for DeExcitation -- no data files to check for";
    break;
#endif
#ifdef INCL_DEEXCITATION_GEMINIXX
  case G4INCL::DeExcitationGEMINIXX:
    LOG("HINCLCascadeIntranuke", pINFO)
    << "using GEMINIXX for DeExcitation -- check for data location";
    datapaths.clear();
    datapaths.push_back(theINCLConfig->getGEMINIXXDataFilePath());
    datapaths.push_back("!INCL-gemini-data-alt1");
    datapaths.push_back("!INCL-gemini-data-alt2");
    needed_update |=
    LookForAndAddValidPath(datapaths,0,"--geminixx-datafile-path",nflags,flags);
    break;
#endif
  default:
    LOG("HINCLCascadeIntranuke", pINFO)
    << "using no DeExcitation -- no data files to check for";
    break;
  }

  // print info
  std::stringstream ss;
  for (size_t i=0; i < nflags; ++i) {
    if ( i == nflags_in ) ss << "---- list prior to AddDataPathFlags()\n";
    ss << "[" << setw(3) << i << "] " << flags[i] << '\n';
  }
  LOG("HINCLCascadeIntranuke", pNOTICE)
  << "Flags passed to INCLConfigParser"
  << '\n' << ss.str();

  return needed_update;
}

//______________________________________________________________________________
bool HINCLCascadeIntranuke::LookForAndAddValidPath(std::vector<std::string>& datapaths,
 size_t defaultIndx,
 const char* optString,
 size_t& nflags, char** flags) {

  // okay, we have a series of paths _OR_ parameter names ("!" as first char)
  // loop over possibilities
  //    convert parameter names to their actual values
  //    expand string to evaluate ${} env variables
  //    check if the path exists
  // first instance that exists is our choice
  // if it isn't the defaultIndx entry, then we must add to the
  // flags passed to ConfigParser
  // add the optString, then the path

  bool needed_update = false;

  LOG("HINCLCascadeIntranuke", pINFO)
  << "looking for a valid path for " << optString
  << " (default [" << defaultIndx << "]";

  size_t foundIndx = SIZE_MAX; // flag as unfound
  size_t npaths = datapaths.size();
  for (size_t ipath = 0; ipath < npaths; ++ipath) {
    std::string& apath = datapaths[ipath];
    LOG("HINCLCascadeIntranuke", pINFO)
    << "looking at " << apath;
    if ( apath.at(0) == '!' ) {
      apath.erase(0,1); // remove the "!"
      // parameter name, returned value, default value
      // default if not found is parameter name, just to make clear what failed
      std::string newpath = "";
      std::string notfoundvalue = std::string("no-such-param-") + apath;
      GetParamDef(apath,newpath,notfoundvalue);
      LOG("HINCLCascadeIntranuke", pINFO)
      << "fetch param "<< "[" << ipath << "] "<< apath << " got " << newpath;
      apath = newpath;
    }
    const char* expandedPath =  gSystem->ExpandPathName(apath.c_str());
    if ( ! expandedPath ) {
      LOG("HINCLCascadeIntranuke", pINFO)
      << "expandedPath was NULL";
      continue;
    }
    bool valid = utils::system::DirectoryExists(expandedPath);
    LOG("HINCLCascadeIntranuke", pINFO)
    << "expandedPath " << expandedPath << " "
    << ((valid)?"valid":"doesn't exist");
    if ( valid ) {
      apath = std::string(expandedPath);
      foundIndx = ipath;
      break;
    }
  }
  if ( foundIndx == defaultIndx ) {
    // nothing to do ... the default works
  } else if ( foundIndx > npaths-1 ) {
    // nothing valid found      // yikes
    std::stringstream ss;
    for (size_t ipath = 0; ipath < npaths; ++ipath) {
      std::string& apath = datapaths[ipath];
      ss << "[" << ipath << "] " << apath << "\n";
    }
    LOG("HINCLCascadeIntranuke", pWARN)
    << "no valid path found for " << optString
    << ", tried: \n" << ss.str();
  } else {
    // add the flag with the valid path
    flags[nflags] = strdup(optString); ++nflags;
    flags[nflags] = strdup(datapaths[foundIndx].c_str()); ++nflags;
    needed_update = true;
  }

  return needed_update;
}

//______________________________________________________________________________
int HINCLCascadeIntranuke::doCascade(GHepRecord * evrec) const {

  if ( ! theINCLConfig || ! theINCLModel ) return 0;

  int tpos = evrec->TargetNucleusPosition();
  GHepParticle * target = evrec->Particle(tpos);
  GHepParticle * pprobe = evrec->Probe();

  const ParticleType theType = toINCLparticletype(pprobe->Pdg());

  double  E    = (pprobe->E())*1000;
  double massp = G4INCL::ParticleTable::getRealMass(theType);
  double EKin = E - massp;

  G4INCL::ParticleSpecies  theSpecies;
  theSpecies.theType = theType;
  theSpecies.theA    = pdgcpiontoA(pprobe->Pdg());
  theSpecies.theZ    = pdgcpiontoZ(pprobe->Pdg());

  G4INCL::Random::SeedVector const theInitialSeeds = G4INCL::Random::getSeeds();
  GHepStatus_t    ist1  = kIStStableFinalState;
  int pdg_codeProbe = 0;
  pdg_codeProbe =  INCLpartycleSpecietoPDGCODE(theSpecies);

  G4INCL::EventInfo result;
  result = theINCLModel->processEvent(theSpecies,EKin,target->A(),target->Z());

  double m_target = ParticleTable::getTableMass(result.At, result.Zt);
  GHepParticle * fsProbe = evrec->Probe();

  TLorentzVector p4h   (0.,0.,fsProbe->Pz(),fsProbe->E());
  TLorentzVector x4null(0.,0.,0.,0.);
  TLorentzVector p4tgt (0.,0.,0.,m_target/1000);
  int pdg_codeTarget= genie::pdg::IonPdgCode(target->A(), target->Z());

  if ( result.transparent ) {
    evrec->AddParticle(pdg_codeProbe, ist1, 0,-1,-1,-1, p4h,x4null);
    evrec->AddParticle(pdg_codeTarget,kIStStableFinalState,1,-1,-1,-1,p4tgt,x4null);
    INCL_DEBUG("Transparent event" << std::endl);
  } else {
    INCL_DEBUG("Number of produced particles: " << result.nParticles << "\n");
    if ( theDeExcitation != 0 ) {
      theDeExcitation->deExcite(&result);
    }

    for (int nP = 0; nP < result.nParticles; nP++){
      if ( nP == result.nParticles-1 ) {
        int pdgRem=INCLtopdgcode(result.A[nP],result.Z[nP]);
        TParticlePDG * pdgRemn=PDGLibrary::Instance()->Find(pdgRem,false); 
        if(!pdgRemn)
        {
          LOG("HINCLCascadeIntranuke", pINFO)
          << "NO Particle with pdg = " << pdgRem << " in PDGLibrary!";
              // Add the particle with status id=15 and change it to HadroBlob
          TVector3 p3M(result.px[nP]/1000,result.py[nP]/1000,result.pz[nP]/1000);

          double MassRem=0.5*((result.px[nP])*(result.px[nP]) +
           (result.py[nP])*(result.py[nP]) +
           (result.pz[nP])*(result.pz[nP]) -
           result.EKin[nP]*result.EKin[nP]) / (result.EKin[nP]);
          float ERem=result.EKin[nP]+MassRem;
          TLorentzVector p4tgtf(p3M,ERem/1000);
          GHepParticle p_outR(kPdgHadronicBlob, kIStFinalStateNuclearRemnant,
            1,-1,-1,-1, p4tgtf, x4null);;
          evrec->AddParticle(p_outR);   
        }
        else
        {
          GHepParticle *p_outR =
          INCLtoGenieParticle(result,nP,kIStStableFinalState,1,-1);
          evrec->AddParticle(*p_outR);
          delete p_outR;
        }
      } else {
        GHepParticle *p_outR =
        INCLtoGenieParticle(result,nP,kIStStableFinalState,0,-1);
        evrec->AddParticle(*p_outR);
        delete p_outR;
      }

    }
  }
  return 0;
}

void HINCLCascadeIntranuke::ProcessEventRecord(GHepRecord * evrec) const {

  LOG("HINCLCascadeIntranuke", pNOTICE)
  << "************ Running HINCLCascadeIntranuke MODE INTRANUKE ************";

  fGMode = evrec->EventGenerationMode();
  if ( fGMode == kGMdHadronNucleus ||
   fGMode == kGMdPhotonNucleus    ) {
    HINCLCascadeIntranuke::doCascade(evrec);
} else if ( fGMode == kGMdLeptonNucleus ) {
  HINCLCascadeIntranuke::TransportHadrons(evrec);
}

LOG("HINCLCascadeIntranuke", pINFO) << "Done with this event";
}

bool HINCLCascadeIntranuke::CanRescatter(const GHepParticle * p) const {

  // checks whether a particle that needs to be rescattered, can in fact be
  // rescattered by this cascade MC
  assert(p);
  return  ( p->Pdg() == kPdgPiP     ||
    p->Pdg() == kPdgPiM     ||
    p->Pdg() == kPdgPi0     ||
    p->Pdg() == kPdgProton  ||
    p->Pdg() == kPdgNeutron
    );
}

void HINCLCascadeIntranuke::TransportHadrons(GHepRecord * evrec) const {

  int inucl = -1;
  if ( fGMode == kGMdHadronNucleus ||
   fGMode == kGMdPhotonNucleus    ) {
    inucl = evrec->TargetNucleusPosition();
} else {
  if ( fGMode == kGMdLeptonNucleus ||
   fGMode == kGMdNucleonDecay  ||
   fGMode == kGMdNeutronOsc       ) {
    inucl = evrec->RemnantNucleusPosition();
}
}

LOG("HINCLCascadeIntranuke", pNOTICE)
<< "Propagating hadrons within nucleus found in position = " << inucl;
int tpos = evrec->TargetNucleusPosition();
GHepParticle * target = evrec->Particle(tpos);
GHepParticle * nucl = evrec->Particle(inucl);
if ( ! nucl ) {
  LOG("HINCLCascadeIntranuke", pERROR)
  << "No nucleus found in position = " << inucl;
  LOG("HINCLCascadeIntranuke", pERROR)
  << *evrec;
  return;
}
fRemnA = nucl->A();
fRemnZ = nucl->Z();
GHepParticle * Incident_particle = evrec->Particle(0);

int A_f(0), Z_f(0), Aft(0), A_i(target->A()),Z_i(0), Charge_probe(Incident_particle->Charge());
if(Charge_probe==0) Z_i=target->Z();
else if(Charge_probe<0) Z_i=target->Z()-1;
else if(Charge_probe>0)Z_i=target->Z()+1;


LOG("HINCLCascadeIntranuke", pNOTICE)
<< "Nucleus (A,Z) = (" << fRemnA << ", " << fRemnZ << ")";

const TLorentzVector & p4nucl = *(nucl->P4());
TLorentzVector x4null(0.,0.,0.,0.);
fRemnP4 = p4nucl;

TObjArrayIter piter(evrec);
GHepParticle * p = 0;

int icurr = -1;

bool is_QE = evrec->Summary()->ProcInfo().IsQuasiElastic();

TLorentzVector * p_4 = nucl->P4();
  // momentum of the remnant nucleus.
double pxRemn = p_4->Px();
double pyRemn = p_4->Py();
double pzRemn = p_4->Pz();
int pdg_codeTargetRemn= genie::pdg::IonPdgCode(nucl->A(),nucl->Z());
TLorentzVector p4tgf(p_4->Px(),p_4->Py(),p_4->Pz(),0.0);

  // Loop over GHEP and run intranuclear rescattering on handled particles
std::vector<G4INCL::EventInfo>ListeOfINCLresult;
std::vector<int> Postion_evrec,num_of_AZexception;
  GHepParticle * fsl = evrec->FinalStatePrimaryLepton();  // primary Lepton
  double ExcitaionE(0), the_pxRemn(0), the_pyRemn(0), the_pzRemn(0);
  int Zl(0), Aresult(0), Zresult(0), Aexception(0), Zexception(0),Pos(0),
  mother1(0),mother2(0),theA_Remn(0), theZ_Remn(0);

  if ( fsl->Charge() == 0. ) Zl =  0;
  else if ( fsl->Charge()  < 0. ) Zl = -1;
  else if ( fsl->Charge()  > 0. ) Zl =  1;
  bool has_remnant = false;
  while ( (p = (GHepParticle *) piter.Next() ) ) {

    icurr++;

    // Check whether the particle needs rescattering, otherwise skip it
    if( ! this->NeedsRescattering(p) ) continue;
    GHepParticle * sp = new GHepParticle(*p);

    // Set clone's mom to be the hadron that was cloned
    sp->SetFirstMother(icurr);

    // Check whether the particle can be rescattered
    if ( ! this->CanRescatter(sp) ) {

      // if I can't rescatter it, I will just take it out of the nucleus
      LOG("HINCLCascadeIntranuke", pNOTICE)
      << "... Current version can't rescatter a " << sp->Name();
      sp->SetFirstMother(icurr);
      sp->SetStatus(kIStStableFinalState);
      if ( sp->Charge() == 0. ) {
       Zl+=0;
       Aft+=1;
     }
     else if ( sp->Charge()  < 0. ) {
       Zl-=1;
       Aft+=1;
     }
     else if ( sp->Charge()  > 0. ) {
       Zl+=1;
       Aft+=1;
     }
     
     evrec->AddParticle(*sp);
     evrec->Particle(sp->FirstMother())->SetRescatterCode(1);
     delete sp;
      continue; // <-- skip to next GHEP entry
    }

    TLorentzVector *v4= sp->GetX4();

    ThreeVector thePosition(0.,0.,0.);
    ThreeVector momentum (0.,0.,0.);
    thePosition.setX(v4->X());
    thePosition.setY(v4->Y());
    thePosition.setZ(v4->Z());
    TLorentzVector * p4 = sp->P4();

    momentum.setX(p4->Px()*1000);
    momentum.setY(p4->Py()*1000);
    momentum.setZ(p4->Pz()*1000);

    int pdgc = sp->Pdg();

    const ParticleType theType = toINCLparticletype(pdgc);

    if ( theType == G4INCL::UnknownParticle) {
      // INCL++ can't handle the particle
      sp->SetFirstMother(icurr);
      sp->SetStatus(kIStStableFinalState);
      evrec->AddParticle(*sp);
      evrec->Particle(sp->FirstMother())->SetRescatterCode(1);
      delete sp;
      continue;
    }

    double  E    = (p4->Energy())*1000;
    double massp = G4INCL::ParticleTable::getRealMass(theType);

    double EKin = E - massp;

    G4INCL::ParticleSpecies  theSpecies;
    theSpecies.theType=theType;
    theSpecies.theA=pdgcpiontoA(sp->Pdg());
    theSpecies.theZ=pdgcpiontoZ(sp->Pdg());

    G4INCL::Particle *pincl =
    new G4INCL::Particle( theType , E , momentum, thePosition);

    G4INCL::Random::SeedVector const theInitialSeeds =
    G4INCL::Random::getSeeds();

    G4INCL::Random::saveSeeds();
    G4INCL::EventInfo result;


    result=theINCLModel->processEvent(theSpecies,pincl,EKin,fRemnA,fRemnZ,"FSI");

    // Exception get remnant with Z and A <0
    Aresult += (fRemnA + pdgcpiontoA(sp->Pdg())- result.ARem[0]); // Aresult & Zresult are particles from cascade
    Zresult += (fRemnZ + pdgcpiontoZ(sp->Pdg())- result.ZRem[0]);
    Aexception = A_i - Aresult; // remaining particles in the nucleus
    Zexception = Z_i - Zresult;
    bool AZexception(false);
    if ( Zexception <= 0 || Aexception <= 0 || Aexception<=Zexception) {
      // Make sure that ARemn and Zremn >0
      Zl+=pdgcpiontoZ(sp->Pdg());
      Aft+=pdgcpiontoA(sp->Pdg());
      sp->SetFirstMother(icurr);
      sp->SetStatus(kIStStableFinalState);
      evrec->AddParticle(*sp);
      evrec->Particle(sp->FirstMother())->SetRescatterCode(1);
      delete sp;
      for (size_t it=0; it<ListeOfINCLresult.size();it++){
     Pos=Postion_evrec.at(it);// Position of the mother in evrec

     GHepParticle * pinthenucleus = evrec->Particle(Pos);
     theA_Remn+= (fRemnA + pdgcpiontoA(pinthenucleus->Pdg())- ListeOfINCLresult.at(it).ARem[0]);
     theZ_Remn+= (fRemnZ + pdgcpiontoZ(pinthenucleus->Pdg())- ListeOfINCLresult.at(it).ZRem[0]);
     if ( (A_i-theA_Remn-Aft) < (Z_i-theZ_Remn-Zl) ) {
      theA_Remn-= (fRemnA + pdgcpiontoA(pinthenucleus->Pdg())- ListeOfINCLresult.at(it).ARem[0]);
      theZ_Remn-= (fRemnZ + pdgcpiontoZ(pinthenucleus->Pdg())- ListeOfINCLresult.at(it).ZRem[0]);
      Zl+=pdgcpiontoZ(pinthenucleus->Pdg());
      Aft+=pdgcpiontoA(pinthenucleus->Pdg());
      
      pinthenucleus->SetFirstMother(Pos);
      pinthenucleus->SetStatus(kIStStableFinalState);
      evrec->AddParticle(*pinthenucleus);
      evrec->Particle(pinthenucleus->FirstMother())->SetRescatterCode(1);
      AZexception=true;
      num_of_AZexception.push_back(it);
    } else {
      the_pxRemn+=ListeOfINCLresult.at(it).pxRem[0];
      the_pyRemn+=ListeOfINCLresult.at(it).pyRem[0];
      the_pzRemn+=ListeOfINCLresult.at(it).pzRem[0];
      ExcitaionE+=ListeOfINCLresult.at(it).EStarRem[0];
    }
  }
  if (AZexception) {
    for (size_t it=0;it<num_of_AZexception.size();it++) {
      ListeOfINCLresult.pop_back();
    }}

    if(ListeOfINCLresult.size() != 0) {
      int Number_of_Sec=ListeOfINCLresult.size();
      int mom2(0),mom1(0);
      mom1=Postion_evrec.at(0);
      if(Number_of_Sec==1) mom2=-1;
      if(Number_of_Sec>1) mom2=Postion_evrec.at(Number_of_Sec-1);

      for (size_t it=0; it<ListeOfINCLresult.size();it++){

        if ( it < ListeOfINCLresult.size()-1 ) {
          for (int nP=0; nP < ListeOfINCLresult.at(it).nParticles; nP++ ) {
            GHepParticle *p_outD = INCLtoGenieParticle(ListeOfINCLresult.at(it),
             nP,kIStStableFinalState,mom1,mom2);
            evrec->AddParticle(*p_outD);
            delete p_outD;
          } //Add result without the remnant nucleus
        } else {
          ListeOfINCLresult.at(it).ARem[0]=A_i-theA_Remn- Aft;
          ListeOfINCLresult.at(it).ZRem[0]=Z_i-theZ_Remn- Zl;

          ListeOfINCLresult.at(it).pxRem[0]= the_pxRemn + (pxRemn*1000);
          ListeOfINCLresult.at(it).pyRem[0]= the_pyRemn + (pyRemn*1000);
          ListeOfINCLresult.at(it).pzRem[0]= the_pzRemn + (1000*pzRemn);
          ListeOfINCLresult.at(it).EStarRem[0]=ExcitaionE;
          theDeExcitation->deExcite(&ListeOfINCLresult.at(it));
          for (int nP=0;nP<ListeOfINCLresult.at(it).nParticles;nP++ ) {
            int rem_index=FindlargestFragment(ListeOfINCLresult.at(it));
            if(nP==rem_index||ListeOfINCLresult.at(it).A[nP]>1) {
             mom1=inucl;
             mom2=-1;
           }
           else{
            mom1=Postion_evrec.at(0);
            if(Number_of_Sec==1) mom2=-1;
            if(Number_of_Sec>1) mom2=Postion_evrec.at(Number_of_Sec-1);
          }
          GHepParticle *p_outFinal = INCLtoGenieParticle(ListeOfINCLresult.at(it),
           nP,kIStStableFinalState,mom1,mom2);
          evrec->AddParticle(*p_outFinal);
          delete p_outFinal;
          has_remnant=true;
          }//Add all the result with the correct remnant nucleus
        }
      }
      ListeOfINCLresult.clear();
      num_of_AZexception.clear();
    }//  
  } else {
      //if result  give a transparent event FSI=1
      // Store *sp
    if ( result.transparent ) {
      Zl+=pdgcpiontoZ(sp->Pdg());
      Aft+=pdgcpiontoA(sp->Pdg());
        //sp->SetFirstMother(icurr);
      sp->SetStatus(kIStStableFinalState);
      evrec->AddParticle(*sp);
      evrec->Particle(sp->FirstMother())->SetRescatterCode(1);
      delete sp;
    } else {
      Postion_evrec.push_back(icurr);
      ListeOfINCLresult.push_back(result);
      delete sp;
    }

  }

  } //Ghep-entry

  if ( ListeOfINCLresult.size() != 0 ) {
    bool AZexception=false;
    for (size_t it=0; it<ListeOfINCLresult.size();it++){
      Pos = Postion_evrec.at(it);  // Position of the mother in evrec


      GHepParticle * pinthenucleus = evrec->Particle(Pos);
      theA_Remn+= (fRemnA + pdgcpiontoA(pinthenucleus->Pdg())- ListeOfINCLresult.at(it).ARem[0]);
      theZ_Remn+= (fRemnZ + pdgcpiontoZ(pinthenucleus->Pdg())- ListeOfINCLresult.at(it).ZRem[0]);
      if ( (A_i-theA_Remn-Aft) < (Z_i-theZ_Remn-Zl) ) {
        theA_Remn-= (fRemnA + pdgcpiontoA(pinthenucleus->Pdg())- ListeOfINCLresult.at(it).ARem[0]);
        theZ_Remn-= (fRemnZ + pdgcpiontoZ(pinthenucleus->Pdg())- ListeOfINCLresult.at(it).ZRem[0]);
        Zl+=pdgcpiontoZ(pinthenucleus->Pdg());
        Aft+=pdgcpiontoA(pinthenucleus->Pdg());

        pinthenucleus->SetFirstMother(Pos);
        pinthenucleus->SetStatus(kIStStableFinalState);
        evrec->AddParticle(*pinthenucleus);
        AZexception=true;
        num_of_AZexception.push_back(it);
      } else {
        the_pxRemn+=ListeOfINCLresult.at(it).pxRem[0];
        the_pyRemn+=ListeOfINCLresult.at(it).pyRem[0];
        the_pzRemn+=ListeOfINCLresult.at(it).pzRem[0];
        ExcitaionE+=ListeOfINCLresult.at(it).EStarRem[0];
      }
    }
    if (AZexception) {
      for (size_t it=0;it<num_of_AZexception.size();it++) {
        ListeOfINCLresult.pop_back();
      }}
      for (size_t it=0; it < ListeOfINCLresult.size(); it++) {
        if ( is_QE) {
        // QE - event
          ListeOfINCLresult.at(it).pxRem[0] += pxRemn*1000;
          ListeOfINCLresult.at(it).pyRem[0] += pyRemn*1000;
          ListeOfINCLresult.at(it).pzRem[0] += 1000*pzRemn;
          if ( theDeExcitation != 0 ) theDeExcitation->deExcite(&ListeOfINCLresult.at(it));
          int  rem_index=FindlargestFragment(ListeOfINCLresult.at(it));

          for (int nP=0; nP < ListeOfINCLresult.at(it).nParticles; nP++ ) {
           // if(nP==ListeOfINCLresult.at(it).nParticles-1)            Pos=inucl;
            if(nP==rem_index) Pos=inucl;
            else if(nP!=rem_index) Pos=Postion_evrec.at(it);
            GHepParticle *p_out = INCLtoGenieParticle(ListeOfINCLresult.at(it),
              nP,kIStStableFinalState,Pos,-1);
            evrec->AddParticle(*p_out);
            delete p_out;
        }// Add to evrec the result
      } else {
        if ( it < ListeOfINCLresult.size()-1 ) {
          for (int nP=0; nP < ListeOfINCLresult.at(it).nParticles; nP++ ) {
            A_f+=ListeOfINCLresult.at(it).A[nP];
            Z_f+=ListeOfINCLresult.at(it).Z[nP];
          if(ListeOfINCLresult.at(it).A[nP]>1) { // if the event emits cluster during cascade
            mother1=inucl;
            mother2=-1;
          }
          else{
            mother1=Postion_evrec.at(0);
            if(ListeOfINCLresult.size()==1) mother2=-1;
            else if(ListeOfINCLresult.size()>1){
              mother2=Postion_evrec.at(ListeOfINCLresult.size()-1);
            }
          }
          GHepParticle *p_outD =
          INCLtoGenieParticle(ListeOfINCLresult.at(it),nP,kIStStableFinalState,mother1,mother2);
          evrec->AddParticle(*p_outD);
          delete p_outD;
            } //Add result without the remnant nucleus
          } else {
            ListeOfINCLresult.at(it).ARem[0]=A_i-theA_Remn-Aft;
            ListeOfINCLresult.at(it).ZRem[0]=Z_i-theZ_Remn-Zl;
            ListeOfINCLresult.at(it).pxRem[0]= the_pxRemn + (pxRemn*1000);
            ListeOfINCLresult.at(it).pyRem[0]= the_pyRemn + (pyRemn*1000);
            ListeOfINCLresult.at(it).pzRem[0]= the_pzRemn + (1000*pzRemn);
            ListeOfINCLresult.at(it).EStarRem[0]=ExcitaionE;
            if ( theDeExcitation != 0 ) theDeExcitation->deExcite(&ListeOfINCLresult.at(it));
            int rem_index=FindlargestFragment(ListeOfINCLresult.at(it));
            for (int nP=0; nP < ListeOfINCLresult.at(it).nParticles; nP++ ) {
              A_f+=ListeOfINCLresult.at(it).A[nP];
              Z_f+=ListeOfINCLresult.at(it).Z[nP];
              if(nP==rem_index||ListeOfINCLresult.at(it).A[nP]>1) {
                mother1=inucl;
                mother2=-1;
              }
              else{
                mother1=Postion_evrec.at(0);
                if(ListeOfINCLresult.size()==1) mother2=-1;
                else if(ListeOfINCLresult.size()>1){
                  mother2=Postion_evrec.at(ListeOfINCLresult.size()-1);
                }
              }
              GHepParticle *p_outR = INCLtoGenieParticle(ListeOfINCLresult.at(it),
               nP,kIStStableFinalState,mother1,mother2);
              evrec->AddParticle(*p_outR);
              has_remnant=true;
              delete p_outR;
          } //Add all the result with the correct remnant nucleus
        }
      }
    }// Loop over all the result from INCLCascade
  }

  if ( ListeOfINCLresult.size() == 0 && !has_remnant) {
   TParticlePDG * pdgRemn=PDGLibrary::Instance()->Find(pdg_codeTargetRemn,false); 
   if(!pdgRemn)
   {
    LOG("HINCLCascadeIntranuke", pINFO)
    << "NO Particle with pdg = " << pdg_codeTargetRemn << " in PDGLibrary!";
    GHepParticle remnant(kPdgHadronicBlob, kIStFinalStateNuclearRemnant, inucl,-1,-1,-1, fRemnP4, x4null);
    evrec->AddParticle(remnant);
  }else{
    GHepParticle remnant(pdg_codeTargetRemn, kIStStableFinalState, inucl,-1,-1,-1, fRemnP4, x4null);
    evrec->AddParticle(remnant);
  }
}
evrec->Particle(inucl)->SetStatus(kIStIntermediateState);


int dau1(0), dau2(0),Nsec=ListeOfINCLresult.size();
if(Nsec>1){
  GHepParticle * pinN = evrec->Particle(Postion_evrec.at(0));
  dau1=pinN->FirstDaughter();
  dau2=pinN->LastDaughter(); 
  for(int ii=1;ii<Nsec;ii++){
    evrec->Particle(Postion_evrec.at(ii))->SetFirstDaughter(dau1);
    evrec->Particle(Postion_evrec.at(ii))->SetLastDaughter(dau2);
  }
}
}

//______________________________________________________________________________
int HINCLCascadeIntranuke::pdgcpiontoA(int pdgc) const {

  if      ( pdgc == 2212 || pdgc == 2112 ) return 1;
  else if ( pdgc ==  211 || pdgc == -211 || pdgc == 111 ) return 0;
  return 0;  // return something

}

//______________________________________________________________________________
int HINCLCascadeIntranuke::pdgcpiontoZ(int pdgc) const {

  if      ( pdgc == 2212 || pdgc == 211 ) return 1;
  else if ( pdgc == 2112 || pdgc == 111 ) return 0;
  else if ( pdgc == -211 ) return -1;
  return 0; // return something

}

//______________________________________________________________________________
bool HINCLCascadeIntranuke::NeedsRescattering(const GHepParticle * p) const {

  // checks whether the particle should be rescattered
  assert(p);
  // attempt to rescatter anything marked as 'hadron in the nucleus'
  return ( p->Status() == kIStHadronInTheNucleus );

}

//______________________________________________________________________________
void HINCLCascadeIntranuke::Configure(const Registry & config) {

  LOG("HINCLCascadeIntranuke", pDEBUG)
  << "Configure from Registry: '" << config.Name() << "'\n"
  << config;

  Algorithm::Configure(config);
  this->LoadConfig();

}

//___________________________________________________________________________
void HINCLCascadeIntranuke::Configure(string param_set) {

  LOG("HINCLCascadeIntranuke", pDEBUG)
  << "Configure from param_set name: " << param_set;

  Algorithm::Configure(param_set);
  this->LoadConfig();

}

#endif // __GENIE_INCL_ENABLED__
