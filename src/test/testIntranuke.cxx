//____________________________________________________________________________
/*!

\program testIntranuke

\brief   testIntranuke from hadtest.F

         Syntax :
           testIntranuke [-n nev] [-e E] [-p nupdg] [-t tgtpdg] [-r run#] [-a R0] -i inputpdg -k KE

         Options :
           [] denotes an optional argument
           -n specifies the number of events to generate
           -e specifies the neutrino energy
              (if what follows the -e option is a comma separated pair of
               values it will be interpreted as an energy range and neutrino
               energies will be generated uniformly in this range).
              Note that this is the energy of an "interacting" neutrino, not
              the energy of a "flux" neutrino -- the  neutrino is forced to
              interact anyway (see GMCJDriver for how to input fluxes)
           -p specifies the neutrino PDG code
           -t specifies the target PDG code (std format: 1aaazzz000)
           -r specifies the MC run number
           -i specifies the rescattering particle PDG code as a input
           -k specifies the rescattering particle kinetic energy
           -a specifies an effective nucleus size given in fm (R=Ro*A^1/3)

         Example:
           gtestIntranuke -n 300 -e 6.5 -p 14 -t 1056026000 -a 1.3 -i 211 -k 0.165

           will generate 300 events of muon neutrinos (pdg = 14) on Iron
           (A=56,Z=26) at E = 6.5 GeV REGARDLESS OF cross section splines
           and will only generate a rescattering pion with input of KE,
           as well as vertex at around nuclear boundry using an effective nucleus size
           and will show cross section.

\author  Minsuk Kim and Steve Dytman
         University of Pittsburgh

\version 1.0

\created Mar 1, 2007

\cpright Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
         All rights reserved.
         For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#include <cassert>
#include <sstream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

#include <iomanip>

#include "EVGCore/EventRecord.h"
#include "EVGDrivers/GEVGDriver.h"
#include "EVGDrivers/GMCJMonitor.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
//#include "Ntuple/NtpWriter.h"
//#include "Ntuple/NtpMCFormat.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/XSecSplineList.h"
#include "Utils/StringUtils.h"
#include "Utils/CmdLineArgParserUtils.h"
#include "Utils/CmdLineArgParserException.h"

#include "Algorithm/Algorithm.h"
#include "Algorithm/AlgFactory.h"
#include "Algorithm/AlgConfigPool.h"
#include "Conventions/Controls.h"
#include "Conventions/Units.h"
#include "Conventions/Constants.h"
#include "Numerical/Spline.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodeList.h"
#include "Utils/PrintUtils.h"
#include "Utils/NuclearUtils.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepStatus.h"
#include "HadronTransport/Intranuke.h"
#include "HadronTransport/INukeHadroData.h"
#include "HadronTransport/INukeHadroFates.h"
#include "EVGCore/EventRecordVisitorI.h"

using std::string;
using std::vector;
using std::ostringstream;

using namespace genie;

void GetCommandLineArgs (int argc, char ** argv);
void PrintSyntax        (void);

using namespace genie::utils;
using namespace genie::constants;
using namespace genie::controls;

using std::cout;
using std::endl;
using std::ios;
using std::setw;

//Default options (override them using the command line arguments):
int           kDefOptNevents   = 100;            // n-events to generate
Long_t        kDefOptRunNu     = 0;              // default run number

double        kDefOptMinNuEnergy  = 5.2;
int           kDefOptNuPdgCode    = 14;
int           kDefOptTgtPdgCode   = 1056026000;
double        kDefOptR0           = 1.2;

//User-specified options:
int           gOptNevents;           // n-events to generate
bool          gOptBuildSplines;      // spline building option
double        gOptMinNuEnergy;       // min neutrino energy
double        gOptNuEnergyRange;     // max-min neutrino energy
int           gOptNuPdgCode;         // neutrino PDG code
int           gOptTgtPdgCode;        // target PDG code
Long_t        gOptRunNu;             // run number

int           gOptInputPdgCode;      // rescattering particle PDG code as a input
double        gOptInputKE;           // This is KE = E - M. So E = M + KE
double        gOptR0;                // R0


namespace genie {

//class IntranukeTester;
//__________________________________________________
class IntranukeTester
{
public:
  IntranukeTester() {
    this->Init();
  }
  ~IntranukeTester() {
    this->Cleanup();
  }
  void   TransportHadrons   (GHepRecord* ev) {
    fINTRANUKE->TransportHadrons(ev);
  }
  void   GenerateVertex     (GHepRecord * ev) {
    fINTRANUKE->GenerateVertex(ev);
  }
  void Init() {
    AlgFactory * algf = AlgFactory::Instance();
    fINTRANUKE = dynamic_cast<Intranuke *> (algf->AdoptAlgorithm("genie::Intranuke","Default"));
    assert(fINTRANUKE);
  }
  void Cleanup() {
    if(fINTRANUKE) delete fINTRANUKE;
  }

private:

  Intranuke * fINTRANUKE;

};

} // genie namespace

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  //-- parse command line arguments
  GetCommandLineArgs(argc,argv);

  //-- print the options you got from command line arguments
  LOG("testIntranuke", pINFO) << "Number of events requested = " << gOptNevents;
  //LOG("testIntranuke", pINFO) << "Building splines at init.  = " << gOptBuildSplines; 
  LOG("testIntranuke", pINFO) << "Neutrino PDG code          = " << gOptNuPdgCode;
  LOG("testIntranuke", pINFO) << "Target PDG code            = " << gOptTgtPdgCode;
  LOG("testIntranuke", pINFO) << "MC Run Number              = " << gOptRunNu;
  if(gOptNuEnergyRange>0) {
    LOG("testIntranuke", pINFO) << "Neutrino energy            = ["
        << gOptMinNuEnergy << ", " << gOptMinNuEnergy+gOptNuEnergyRange << "]";
  } else {
    LOG("testIntranuke", pINFO) << "Neutrino energy            = " << gOptMinNuEnergy;
  }

  GHepStatus_t status = kIStHadronInTheNucleus;

  string gHadAlg    = "genie::Intranuke";
  string gHadConfig = "Default";

  //-- get the algorithm factory & config pool
  AlgConfigPool * cnfp = AlgConfigPool::Instance();
  AlgFactory *    algf = AlgFactory::Instance();
  //cout << *algf;
  AlgId id("genie::Intranuke","Default");
  //const Algorithm * alg = algf->GetAlgorithm(id);
  //string config = alg->Id().Config();
  Algorithm * alg = algf->AdoptAlgorithm(id);
  //cout << *alg;
  Registry * reg = cnfp->FindRegistry(alg);
  if(gOptR0!=1.2) {
    reg->UnLock();

    ////reg->Set("Kpt2",1);
    reg->Set("R0",gOptR0);
    ////reg->Set("ct0","2");
    //reg->Set("mode","hN");
    //reg->Set("nuc-removal-energy",0.007); //reg->Set("nucleon-removal-energy",0.007); 

    //reg->Set("INUKE-KPt2",1.1);
    //reg->Set("INUKE-Ro",gOptR0);         [fm]
    //reg->Set("INUKE-FormationZone",2.1); [fm]
    //reg->Set("INUKE-Mode","hN");
    //reg->Set("INUKE-NucRemovalE",0.007); [GeV]

    //-- force reconfiguration
    ////algf->ForceReconfiguration();
    alg->Configure(*reg);
  }
  // frac = Pi Absorption Scale factor

  int NNUCLN = genie::pdg::IonPdgCodeToA(gOptTgtPdgCode); // from $GENIE/PDG/PDGUtils.h
  TParticlePDG * input = PDGLibrary::Instance()->Find(gOptInputPdgCode);
  TParticlePDG * tgt = PDGLibrary::Instance()->Find(gOptTgtPdgCode);
  //double nuclear_density = GV_RHONUC;
  ////gOptR0 = gOptR0 + 0.00000005; (in neugen)
  double nuclear_radius = gOptR0*TMath::Power(NNUCLN,1./3.); //exactly rsiz*TMath::Power(NNUCLN,1/3.);

  double area = TMath::Pi()*TMath::Power(nuclear_radius,2);
  double nuclear_vertex[3];
  for(int k=0; k<3; k++) nuclear_vertex[k] = 0;
  const int nfates = 12;
  int countfate[nfates]; string FateType[nfates];
  double sigma[nfates], err[nfates];
  for(int k=0; k<nfates; k++) { countfate[k] = 0; sigma[k] = 0; err[k] = 0; }

  //-- this driver produces events for monoenergetic neutrinos
  TLorentzVector nu_p4(0.,0.,gOptMinNuEnergy,gOptMinNuEnergy); // px,py,pz,E (GeV)

  //-- create an MC Job Monitor
  GMCJMonitor mcjmonitor(gOptRunNu);

  //-- get a random number generator
  RandomGen * r = RandomGen::Instance();
  time_t t; (void) time(&t);
  r->SetSeed(t);

  //-- generate events / print the GHEP record / add it to the ntuple
  int ievent = 0;
  while ( ievent < gOptNevents) {

    // generate neutrino energy (if an energy range was defined)
    if(gOptNuEnergyRange>0) {
      double Ev = gOptMinNuEnergy + gOptNuEnergyRange * r->RndEvg().Rndm();
      nu_p4.SetPxPyPzE(0,0,Ev,Ev);
    }

    EventRecord * evrec = new EventRecord();

    //-- Generate a vertex within the nucleus (from $NEUGEN3PATH/gen/nuclear_rescattering.F)
    double random[2];
    double xpysq = 10000, rs = 0;
    while ( xpysq > TMath::Power((nuclear_radius-0.01),2)) {

      for(int k=0; k<2; k++) random[k] = r->RndGen().Rndm();
      random[1] = 2*random[1] - 1;
      nuclear_vertex[0] = random[0]*(nuclear_radius-0.01);
      nuclear_vertex[1] = random[1]*(nuclear_radius-0.01);	
      xpysq = TMath::Power(nuclear_vertex[0],2) + TMath::Power(nuclear_vertex[1],2);

    }

    nuclear_vertex[2] = -1.*TMath::Sqrt(TMath::Power((nuclear_radius-0.01),2) - xpysq);
    rs = TMath::Sqrt(TMath::Power(nuclear_vertex[0],2) + TMath::Power(nuclear_vertex[1],2) + TMath::Power(nuclear_vertex[2],2));

    double mass = input->Mass();
    double energy = mass + gOptInputKE;
    double pz = TMath::Sqrt(energy*energy - mass*mass);
    
    TLorentzVector newp4(0.,0.,pz,energy);
    TLorentzVector newv4(nuclear_vertex[0],nuclear_vertex[1],nuclear_vertex[2],0.);
  
    TLorentzVector v4(0,0,0,0);
    TLorentzVector tgt_p4(0,0,0,tgt->Mass());
    evrec->AddParticle(gOptNuPdgCode,kIStInitialState,-1,-1,0,0,nu_p4,v4);
    evrec->AddParticle(gOptTgtPdgCode,kIStInitialState,-1,-1,0,0,tgt_p4,v4);

    int tgt_n = gOptTgtPdgCode - 1000000; 
    int tgt_p = gOptTgtPdgCode - 1001000;
    TParticlePDG * tgtn = PDGLibrary::Instance()->Find(tgt_n);
    TParticlePDG * tgtp = PDGLibrary::Instance()->Find(tgt_p);

    TLorentzVector p4(0.036,-0.143,0.009,0.94);            // neutron(0.94) or proton(0.938)
    TLorentzVector p4_n(-0.036,0.143,-0.009,tgtn->Mass()); // Fe55(51.172) or O15(13.945)
    TLorentzVector p4_p(-0.036,0.143,-0.009,tgtp->Mass()); // Mn55(51.175) or N15(13.973)
    TLorentzVector p4_f(0.415,-0.069,0.425,0.598);         // nu_mu or mu-
    TLorentzVector p4_h(-0.379,-0.074,4.583,5.334);        // HardSyst

    if(r->RndGen().Rndm()<0.65) {
      evrec->AddParticle(2112,kIStNucleonTarget,1,-1,5,5,p4,v4); // neutron
      evrec->AddParticle(tgt_n,kIStStableFinalState,1,-1,-1,-1,p4_n,v4);
    } else {
      evrec->AddParticle(2212,kIStNucleonTarget,1,-1,5,5,p4,v4); // proton
      evrec->AddParticle(tgt_p,kIStStableFinalState,1,-1,-1,-1,p4_p,v4);
    }
    if(r->RndGen().Rndm()<0.35) evrec->AddParticle(gOptNuPdgCode,kIStStableFinalState,0,-1,-1,-1,p4_f,v4); // NC
    else { // CC
      if(gOptNuPdgCode>0) evrec->AddParticle(13,kIStStableFinalState,0,-1,-1,-1,p4_f,v4);
      else evrec->AddParticle(-13,kIStStableFinalState,0,-1,-1,-1,p4_f,v4);
    }
    evrec->AddParticle(1111111002,kIStDISPreFragmHadronicState,2,-1,6,6,p4_h,v4);
    evrec->AddParticle(gOptInputPdgCode,status,5,-1,7,7,newp4,newv4);

    int idx = evrec->ParticlePosition(gOptInputPdgCode,status,0);

    IntranukeTester tester;
    //-- Generate and set a vertex in the nucleus coordinate system
    tester.GenerateVertex(evrec);
    GHepParticle * input14 = evrec->Particle(idx);
    input14->SetPosition(newv4);
    //-- Transport all particles outside the nucleus and exit
    tester.TransportHadrons(evrec);

    bool interacts = false;
    INukeFateHA_t fate = kIHAFtUndefined;
    int nn = 0, np = 0, npip = 0, npi0 = 0, npim = 0;
    TObjArrayIter piter(evrec);
    GHepParticle * p = 0;
    int icurr = -1;
    while( (p = (GHepParticle *) piter.Next()) ) {
      icurr++;
      if(icurr<8) continue;
      if(icurr==8) {
	//cout << " -> Event No. " << setw(4) << ievent << setw(12) << p->Pdg() << setw(5) << p->Name() << setw(3) << p->Status() << endl;
	//cout << " -> KinE = " << p->KinE() << " " << print::P4AsShortString(p->P4()) << endl;
	//cout << " -> rs = " << p->X4()->Vect().Mag() << " " << print::X4AsString(p->X4()) << ", NuclearRadius=" << nuclear_radius << endl;
	if(p->X4()->Vect().Mag() < nuclear_radius) interacts = true;
        if(p->Status() == kIStStableFinalState) {
          if(p->Pdg() != gOptInputPdgCode) fate = kIHAFtCEx;
          else {
	    double E = TMath::Sqrt(2)*p->Mass(); //TMath::Sqrt(2*p->P4()->Vect().Mag2());
            //printf(" ==> %f %f %f %f %f\n",p->P4()->Vect().Mag(),p->Mass(),p->E(),E,p->E()-E);
            fate = (TMath::Abs(p->E()-E)<1E-10) ? kIHAFtElas : kIHAFtInelas;
          }
          break;
        }
      } else if(p->Status() == kIStStableFinalState) {
        if(p->Pdg() == kPdgProton) np++;
        if(p->Pdg() == kPdgNeutron) nn++;
        if(p->Pdg() == kPdgPiP) npip++;
        if(p->Pdg() == kPdgPiM) npim++;
        if(p->Pdg() == kPdgPi0) npi0++;
      }
    }
    int npi = npip + npi0 + npim;
    if(npi==0) { 
      if(nn==1 && np==1) fate = kIHAFtAbsNP;
      if(np==2 && nn==0) fate = kIHAFtAbsPP;
      if(nn==1 && np==2) fate = kIHAFtAbsNPP;
      if(nn==2 && np==1) fate = kIHAFtAbsNNP;
      if(nn==2 && np==2) fate = kIHAFtAbs2N2P; // only for rescattering pion
      if(nn==2 && np==3) fate = kIHAFtAbs2N3P; // either rescattering neutron or proton
    } else {
      if(nn==1 && npip==1 && npi0==0) fate = kIHAFtNPip; // either rescattering neutron or proton
      if(nn==1 && npip==1 && npi0==1) fate = kIHAFtNPipPi0;
    }

    if(interacts) {

      ////GHepParticle * genp = evrec->Particle(idx+2); if(genp) fate = tester.HadronFateHA(sp);
      //cout << " => Event No. " << setw(7) << ievent << ": " << INukeHadroFates::AsString(fate) << endl;
      //cout << " (nn=" << nn << ", np=" << np << ", npip=" << npip << ", npi0=" << npi0 << ", npim=" << npim << ")" << endl;
      countfate[fate]++;
      FateType[fate] = INukeHadroFates::AsString(fate).c_str();

    }
    //cout << " ==> Generated Event GHEP Record: " << *evrec << endl;

    //-- refresh the mc job monitor
    mcjmonitor.Update(ievent,evrec);
    
    ievent++;
    delete evrec;
    
  } // end loop events

  cout << endl << endl;
  cout << " Total cross section results for " << input->GetName() << " (KE=" << gOptInputKE << ") + " << tgt->GetName() << " interaction " << endl;
  cout << " A = " << NNUCLN << " nuclear radius, area are: " << nuclear_radius << " fm, " << area << " fm**2 " << endl;

  double fm2tomb = (units::fm2 / units::mb);

  int cnttot = 0;
  //int nullint = 0;
  double sigtot = 0, sigtoterr = 0;
  double sigtotScat = 0, sigtotAbs = 0, sigtotProd = 0;
  for(int k=0; k<nfates; k++) {

    //if(k==0) nullint += countfate[k];
    if(k!=0) {
      cnttot += countfate[k];
      double ratio = countfate[k]/(double)ievent;
      sigma[k] = fm2tomb * area * ratio;
      err[k]   = fm2tomb * area * TMath::Sqrt(ratio*(1-ratio)/(double)ievent);
      if(err[k]==0) err[k] = fm2tomb * area * TMath::Sqrt(countfate[k])/(double)ievent;
    } 
    if(countfate[k]>0) {
      cout << " --> " << setw(26) << FateType[k] << ": " << setw(7) << countfate[k] << " events -> " << setw(7) << sigma[k] << " +- " << err[k] << " (mb)" << endl;
    }
    if(k>=1 && k<=3) sigtotScat += sigma[k];
    if(k>=4 && k<=9) sigtotAbs += sigma[k];
    if(k>=10 && k<=11) sigtotProd += sigma[k];

  }
  sigtot    = fm2tomb * area * cnttot/(double)ievent;
  sigtoterr = fm2tomb * area * TMath::Sqrt(cnttot)/(double)ievent;
  double ratioAS = (sigtotScat==0) ? 0 : sigtotAbs/(double)sigtotScat;

  cout << " -------------------------------------------------------------------------- " << endl;
  cout << " ==> " << setw(28) << " Total: " << setw(7) << cnttot << " events -> " << setw(7) << sigtot << " +- " << sigtoterr << " (mb)" << endl;
  cout << " ==> " << setw(28) << " Ratio (abs/scat) = " << setw(7) << ratioAS << endl;
  cout << " ==> " << setw(28) << " avg. num of int. = " << setw(7) << cnttot/(double)ievent << endl;
  cout << " ==> " << setw(28) << " no interaction   = " << setw(7) << (ievent-cnttot)/(double)ievent << endl; 
  cout << endl;

  cout << *alg << endl; delete alg;

  return 0;
}
//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("testIntranuke", pNOTICE) << "Parsing command line arguments";

  //-- Optional arguments

  //rescattering particle PDG code:
  try {
    LOG("testIntranuke", pINFO) << "Reading rescattering particle PDG code";
    gOptInputPdgCode = genie::utils::clap::CmdLineArgAsInt(argc,argv,'i');
    if(gOptInputPdgCode != 111 && TMath::Abs(gOptInputPdgCode) != 211 && gOptInputPdgCode != 2112 && gOptInputPdgCode != 2212) {
      LOG("testIntranuke", pFATAL) << "Incorrect PDG code - Exiting";
      PrintSyntax();
      exit(1);
    }
  } catch(exceptions::CmdLineArgParserException e) {
      if(!e.ArgumentFound()) {
      LOG("testIntranuke", pFATAL) << "Unspecified PDG code - Exiting";
      PrintSyntax();
      exit(1);
    }
  }

  //rescattering particle kinetic energy
  try {
    LOG("testIntranuke", pINFO) << "Reading rescattering particle KE energy";
    gOptInputKE = genie::utils::clap::CmdLineArgAsDouble(argc,argv,'k');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("testIntranuke", pFATAL) << "Unspecified KE - Exiting";
      PrintSyntax();
      exit(1);
    }
  }

  try {
    LOG("testIntranuke", pINFO) << "Reading RO";
    gOptR0 = genie::utils::clap::CmdLineArgAsDouble(argc,argv,'a');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      gOptR0 = kDefOptR0;
    }
  }

  //number of events:
  try {
    LOG("testIntranuke", pINFO) << "Reading number of events to generate";
    gOptNevents = genie::utils::clap::CmdLineArgAsInt(argc,argv,'n');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("testIntranuke", pINFO)
            << "Unspecified number of events to generate - Using default";
      gOptNevents = kDefOptNevents;
    }
  }

  //run number:
  try {
    LOG("testIntranuke", pINFO) << "Reading MC run number";
    gOptRunNu = genie::utils::clap::CmdLineArgAsInt(argc,argv,'r');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("testIntranuke", pINFO) << "Unspecified run number - Using default";
      gOptRunNu = kDefOptRunNu;
    }
  }

  //spline building option
  gOptBuildSplines = genie::utils::clap::CmdLineArgAsBool(argc,argv,'s');

  //-- Required arguments

  //neutrino energy:
  try {
    LOG("testIntranuke", pINFO) << "Reading neutrino energy";
    string nue = genie::utils::clap::CmdLineArgAsString(argc,argv,'e');

    // is it just a value or a range (comma separated set of values)
    if(nue.find(",") != string::npos) {
       // split the comma separated list
       vector<string> nurange = utils::str::Split(nue, ",");
       assert(nurange.size() == 2);
       double emin = atof(nurange[0].c_str());
       double emax = atof(nurange[1].c_str());
       assert(emax>emin && emin>0);
       gOptMinNuEnergy   = emin;
       gOptNuEnergyRange = emax-emin;
    } else {
       gOptMinNuEnergy   = atof(nue.c_str());
       gOptNuEnergyRange = -1;
    }
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      gOptMinNuEnergy = kDefOptMinNuEnergy;
      //LOG("testIntranuke", pFATAL) << "Unspecified neutrino energy - Exiting";
      //PrintSyntax();
      //exit(1);
    }
  }

  //neutrino PDG code:
  try {
    LOG("testIntranuke", pINFO) << "Reading neutrino PDG code";
    gOptNuPdgCode = genie::utils::clap::CmdLineArgAsInt(argc,argv,'p');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      gOptNuPdgCode = kDefOptNuPdgCode;
      //LOG("testIntranuke", pFATAL) << "Unspecified neutrino PDG code - Exiting";
      //PrintSyntax();
      //exit(1);
    }
  }

  //target PDG code:
  try {
    LOG("testIntranuke", pINFO) << "Reading target PDG code";
    gOptTgtPdgCode = genie::utils::clap::CmdLineArgAsInt(argc,argv,'t');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      gOptTgtPdgCode = kDefOptTgtPdgCode;
      //LOG("testIntranuke", pFATAL) << "Unspecified target PDG code - Exiting";
      //PrintSyntax();
      //exit(1);
    }
  }

}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("testIntranuke", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "   gtestIntranuke [-n nev] [-e energy] [-p nupdg] [-t tgtpdg] [-r run] [-a R0] -i inputpdg -k KE\n"
    << "    for inputpdg : piplus  =  211\n"
    << "                   pizero  =  111\n"
    << "                   piminus = -211\n"
    << "                   neutron = 2112\n"
    << "                   proton  = 2212\n\n"
    << "    to generate 1k events with intranuclear rescattering in the default target (pi+,Fe56)\n"
    << "    gtestIntranuke -n 1000 -i 211 -k .165\n\n"
    << "    to generate 1k events with intranuclear rescattering in the oxegen target (pi+,O16)\n"
    << "    gtestIntranuke -n 1000 -i 211 -k .165 -t 1016008000\n\n";
}
//____________________________________________________________________________
