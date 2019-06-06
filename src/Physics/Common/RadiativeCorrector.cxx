//____________________________________________________________________________
/*
 Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 19, 2007 - CA
   Added code to generate the interaction according to a realistic nuclear
   density and made that the default setting.
 @ Dec 01, 2007 - CA
   For COH and ve- interactions setting the vertex on the nuclear boundary
   rather than inside the nucleus.
 @ Sep 15, 2009 - CA
   IsFake() and IsNucleus() are no longer available in GHepParticle. 
   Use pdg::IsPseudoParticle() and pdg::IsIon().
 @ Feb 12, 2013 - CA (code from Rosen Matev)
   Handle the IMD annihilation channel.
 @ Jun 03, 2016 - JJ (SD)
   Move code to generate the position to a public method, so this code
   can easily be reused by other classes. (Specifically LwlynSmithQELCCPXSec
   and NievesQELCCPXSec to generate a position before calculating the xsec
   when making splines).
*/
//____________________________________________________________________________
#include <algorithm>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <TMath.h>

#include <RVersion.h>
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,6)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif
#include <TParticlePDG.h>
#include <TF1.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Physics/Common/RadiativeCorrector.h"
#include "Physics/Decay/DecayModelI.h"
#include "Physics/Decay/UnstableParticleDecayer.h"
#include "Physics/Decay/RadiativeDecayer.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Registry/Registry.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/RunOpt.h"

using std::count;
using std::ostringstream;

using namespace genie;
using namespace genie::utils;
using namespace genie::controls;
using namespace genie::constants;

//___________________________________________________________________________
RadiativeCorrector::RadiativeCorrector() :
EventRecordVisitorI("genie::RadiativeCorrector")
{
  fDecayers = 0;
  fInitState = 0;
  //rad_gOptRunNu = 0;
  //rad_kDefOptNtpFormat = kNFGHEP;
}
//___________________________________________________________________________
RadiativeCorrector::RadiativeCorrector(string config) :
EventRecordVisitorI("genie::RadiativeCorrector", config)
{
  fDecayers = 0;
}
//___________________________________________________________________________
RadiativeCorrector::~RadiativeCorrector()
{
  if(fDecayers) {
   fDecayers->clear();
   delete fDecayers;
  }
  if (fInitState)        delete fInitState;
  //rad_gOptRunNu = 0;
  //rad_kDefOptNtpFormat = kNFGHEP;
}
//___________________________________________________________________________
void RadiativeCorrector::BuildInitialState(const InitialState & init_state)
{
  LOG("RadiativeCorrector", pINFO) << "Setting the initial state";

  if(fInitState) delete fInitState;
  fInitState = new InitialState(init_state);

  //this->AssertIsValidInitState();
}
//___________________________________________________________________________
void RadiativeCorrector::ProcessEventRecord(GHepRecord * evrec) const 
{
  // decay a photon with dE from the electron 
  evrec->Print();
  bool radDone = false;
 
  bool vanderhagen = true;
  bool simc = false;
 
  Interaction * interaction = evrec->Summary();
  Kinematics * kine = interaction->KinePtr();
  InitialState * init_state_ptr = interaction -> InitStatePtr();
  Target * target_ptr = init_state_ptr -> TgtPtr();
  int Z = target_ptr->Z();
  //const InitialState & const_init_state = interaction -> InitState();
  //InitialState init_state(const_init_state);
  TObjArrayIter piter(evrec);
  GHepParticle * p = 0;
  unsigned int ipos = 0;
  //unsigned int hitNucleonPosition; 
  while( (p = (GHepParticle *) piter.Next()) )
  {
     if( this->ToBeDecayed(p)) {
      if( init_state_ptr->Probe()->PdgCode() == 11 || init_state_ptr->Probe()->PdgCode() == -11) {
        TLorentzVector p4 = *(p->P4());
        TLorentzVector x4 = *(p->X4());
        DecayerInputs_t dinp;
        dinp.PdgCode = p->Pdg();
        dinp.P4      = &p4;
	double e = dinp.P4->E();
	double energyLoss = 0.;

	//double targetThickness = init_state.TgtThickness();
	//if (fabs(p4.E() - init_state_ptr->ProbeE(kRfLab)) > 0.1) {
	//	LOG("RadiativeCorrector", pNOTICE) << "Decaying ISR particle is not the probe, not performing ISR correction";
	//	continue;
	//}	
	/////double energyLoss = fCurrDecayer->findEnergyLoss(dinp);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// find energy loss
        
	if (vanderhagen) {
           //double energyLossLimit = init_state_ptr->ProbeE(kRfLab) - kine->FSLeptonP4().E();
           double energyLossLimit = evrec->Probe()->GetP4()->E() - kine->FSLeptonP4().E();
	   LOG("RadiativeCorrector", pINFO) << "probE  "<<init_state_ptr->ProbeE(kRfLab)<<" evrec->Probe()->GetP4()->E() "<<evrec->Probe()->GetP4()->E()<<" final state lepton e "<< kine->FSLeptonP4().E()<<" energy loss limit "<<energyLossLimit;
           double a = (kAem/kPi)*(TMath::Log(kine->Q2(true)/pow(kElectronMass,2)) - 1.);
           if (energyLossLimit<0) energyLossLimit = 0.;
           TF1 *f = new TF1("f","([0]/x)*TMath::Power(x/[1],[0])",1E-10,energyLossLimit);
           f->SetParameter(0,a);
           f->SetParameter(1,e);
           energyLoss = f->GetRandom();
	   delete f;
	   LOG("RadiativeCorrector", pINFO) << "Energy loss is "<<energyLoss;
	}
	double L1,L2,b,thickness,lambda_e,g,e_gamma_max,e_gamma_min,power_hi,power_lo;
	if (simc) { 
           L1 = log(184.15) - (1./3)*log(Z);
	   L2 = log(1194.) - (2./3)*log(Z);
           b = (1./9)*(12 + float(Z+1)/(Z*(Z*L1 + L2)));
	   lambda_e = (kAem/kPi)*( TMath::Log(4*pow(p->P4()->P(),2)/pow(kElectronMass,2)-1) + 2*TMath::Log(init_state_ptr->GetProbeP4(kRfLab)->P()/kine->FSLeptonP4().P()) + TMath::Log(0.5*(1-kine->FSLeptonP4().CosTheta()) ) );
	   thickness =  0.005176; //thickness in radiation length (0.1 cm carbon)
           g = b*thickness + lambda_e;
           std::cout<<" lambda e  is "<<lambda_e<<" g "<<g<<std::endl;

	   e_gamma_max = evrec->Probe()->GetP4()->E() - kine->FSLeptonP4().E();
	   e_gamma_min = 0.2;
	   std::cout<<"e_gamma_max "<<e_gamma_max<<std::endl;
	   power_hi = pow(e_gamma_max,g);
	   power_lo  = pow(e_gamma_min,g);
	   std::cout<<"power_hi "<<power_hi<<" power_lo "<<power_lo<<std::endl;
	   double ymin = power_lo/power_hi;
	   std::cout<<"ymin "<<ymin<<std::endl; 
	   TRandom3 rnd;
	   rnd.SetSeed(0);
	   double y = ymin + rnd.Rndm()*(1.-ymin);
	   std::cout<<"y "<<y<<std::endl;
	   double x = pow(y,1./g);
	   energyLoss = x*e_gamma_max;
	   //double c_int = lambda_e/pow(init_state_ptr->GetProbeP4(kRfLab)->E()*kine->FSLeptonP4().E(),lambda_e/2);
	   //TF1 *f = new TF1("f","([0]*pow(x,[0]-1)/[1]"
	}

	double momentumLoss = energyLoss;
	double ptLoss;
	double pzLoss;
	if (p4.Pz()==p4.E()) // for the z direction going probe theta = -nan 
	{
	   ptLoss = 0.;
	   pzLoss = momentumLoss;
	}
	else
	{
	   ptLoss = momentumLoss*sin(p4.Theta()); 
	   pzLoss = momentumLoss*cos(p4.Theta());
	}
	TLorentzVector p4RadGamma;
	p4RadGamma.SetPxPyPzE(ptLoss*cos(p4.Phi()),ptLoss*sin(p4.Phi()),pzLoss,energyLoss);
	//LOG("RadiativeCorrector", pINFO) << "old momentum e "<<p4.E()<<" pz "<<p4.Pz()<<" theta "<<p4.Theta()<<" eta "<<p4.Eta();
	//LOG("RadiativeCorrector", pINFO) << "gamma        e "<<p4RadGamma.E()<<" pz "<<p4RadGamma.Pz();
	p4 = p4 - p4RadGamma;
	//LOG("RadiativeCorrector", pINFO) << "new momentum e "<<p4.E()<<" pz "<<p4.Pz();

	if (fISR && !radDone) {
	  LOG("RadiativeCorrector", pINFO) << "performing ISR correction for: " << p->Name() << " reduced energy is : "<<p4.E();
	  // changing the probe energy instead of decaying it 
	  init_state_ptr->SetProbeP4(p4);
	  //p->SetMomentum(p4);

	  //-- Mark it as a 'decayed state' & add its daughter links
	  p->SetStatus(kIStDecayedState);
	  ////-- Add the mom & daughters to the event record
	  LOG("RadiativeCorrector", pINFO) << "Adding daughter... PDG=" << p->Pdg();
          evrec->AddParticle(p->Pdg(), kIStIntermediateState, ipos,-1,-1,-1, p4, x4);
          LOG("RadiativeCorrector", pINFO) << "Adding daughter... PDG= 22";
          evrec->AddParticle(22, kIStStableFinalState, ipos,-1,-1,-1, p4RadGamma, x4);
	  radDone = true;
          //LOG("RadiativeCorrector", pINFO) <<"TESTING PROBE "<<init_state_ptr->ProbeE(kRfLab)<<" from event rec mother "<<evrec->Probe()->FirstMother();
          //evrec->Print();
	}
	
        if (!fISR && !radDone) {
	   LOG("RadiativeCorrector", pINFO) << "performing FSR correction for: " << p->Name();
           //-- Mark it as a 'decayed state' & add its daughter links
       	   p->SetStatus(kIStDecayedState);

           ////-- Add the mom & daughters to the event record
           LOG("RadiativeCorrector", pINFO) << "Adding daughter... PDG=" << p->Pdg();
	   evrec->AddParticle(p->Pdg(), kIStStableFinalState, ipos,-1,-1,-1, p4, x4);
           LOG("RadiativeCorrector", pINFO) << "Adding daughter... PDG= 22";
	   evrec->AddParticle(22, kIStStableFinalState, ipos,-1,-1,-1, p4RadGamma, x4);
           
 	   //-- Update the event weight for each weighted particle decay
	   float radcor_weight = 1.;
	   if(vanderhagen) 
	   {
	   	TF1 *fsp = new TF1("fsp","1./x * TMath::Log(1-x)");
	   	double SP = -1*fsp->Integral(0.03,pow(cos(p4.Theta())/2,2));
		delete fsp;
	   	//radcor_weight = 1. + Z*(kAem/kPi)*((13./6)*TMath::Log(kine->Q2(true)/pow(kElectronMass,2))-(28./9)-0.5*pow(TMath::Log(init_state_ptr->ProbeE(kRfLab)/p4.E()),2)-pow(kPi,2)/6 + SP);
	   	radcor_weight = 1. + Z*(kAem/kPi)*( (TMath::Log(kine->Q2(true)/pow(kElectronMass,2))-1)*TMath::Log(pow(init_state_ptr->ProbeE(kRfLab),2)/(init_state_ptr->ProbeE(kRfLab)*kine->FSLeptonP4().E())) 
					       + (13./6)*TMath::Log(kine->Q2(true)/pow(kElectronMass,2))
					       - (28./9)
					       - 0.5*pow(TMath::Log(init_state_ptr->ProbeE(kRfLab)/kine->FSLeptonP4().E()),2)
					       - pow(kPi,2)/6 
					       + SP );
	   }
           if(simc)
	   {
		double c = lambda_e/pow(evrec->Probe()->GetP4()->E() * kine->FSLeptonP4().E(),(lambda_e/2.));
		radcor_weight = c/(g * (power_hi-power_lo));
	   }
           evrec->SetWeight(evrec->Weight() * radcor_weight);
	   LOG("RadiativeCorrector", pINFO) << "Applying radiative correction weight "<<evrec->Weight() * radcor_weight;
	   //---------------------------------
	   evrec->Print();
	   radDone = true;
	}
	if (!fISR && radDone && (p->Status() == kIStStableFinalState) && (p->Pdg()==11)) { 
           LOG("RadiativeCorrector", pINFO) <<"FSR stable electron in pos "<<ipos; 
	   //-- Correct the final state lepton 
           evrec->Probe()->SetFirstDaughter(ipos);
           evrec->Probe()->SetLastDaughter(-1);
           LOG("RadiativeCorrector", pINFO) <<"After setting the lepton first and last mother to -1  "<<evrec->FinalStatePrimaryLepton()->P4()->E();
           evrec->Print(); 
	}
      }
     }
     //if (fISR && (p->Status() ==1)) p->SetStatus(kIStIntermediateState);
     if (fISR && (p->Status() ==1) && (p->Pdg()==11)) {
        //if (p->Pdg() == 2212 || p->Pdg() == 2112 || p->Pdg() == 1000000010 || p->Pdg() == 1000010000 ) {
	//	p->SetMomentum(0.,0.,0.,p->Mass());
	//	//p->SetEnergy(p->Mass());
	//}
	//-- Correct the final state lepton 
        //evrec->Probe()->SetFirstDaughter(-1);
        //evrec->Probe()->SetLastDaughter(-1);
        
	p->SetStatus(kIStDecayedState);
	LOG("RadiativeCorrector", pINFO) <<"After setting the particles status to decyed and momentum of the nucleon to be zero, pdg "<<p->Pdg()<<" mass "<<p->Mass();
        evrec->Print();
        //kine->ClearRunningValues();
        //kine->Reset();
     }
     ipos++;
  } // loop over particles
  //LOG("RadiativeCorrector", pINFO) << "The final state lepton energy is "<<evrec->FinalStatePrimaryLepton()->P4()->E();

  /*if (fISR)
  {
     evrec->Print();
     //evrec->RemoveIntermediateParticles();
     int i=0;
     while( (p = (GHepParticle *)piter.Next()) ) {
       if(!p) continue;
       GHepStatus_t ist = p->Status();
       bool keep = (ist==kIStInitialState) ||
                  (ist==kIStStableFinalState) || (ist==kIStNucleonTarget);
       if (!keep) {
         evrec->RemoveAt(i);
       }
       i++;
     }
     evrec->Print();
  }*/
}
//___________________________________________________________________________
bool RadiativeCorrector::ToBeDecayed(GHepParticle * particle) const
{
   if(particle->Pdg() != 0) {
     if (fISR) {
     	if ((particle->Status() == 0) && (pdg::IsElectron(particle->Pdg()) || pdg::IsPositron(particle->Pdg()))) return true;
     	else return false;	
     }
     else {
     	if ((particle->Status() == 1) && (pdg::IsElectron(particle->Pdg()) || pdg::IsPositron(particle->Pdg()))) return true;
     	else return false;	
     }
   }
   return false;
}
//___________________________________________________________________________
//bool RadiativeCorrector::findEnergyLoss(const DecayerInputs_t & inp) const
//{
//        double e = inp.P4->E();
//        double Q = 1.; // typical momentum scale,
//        double a = (kAem/kPi)*(2*TMath::Log(Q/kElectronMass) - 1.);
//        TF1 *f = new TF1("f","([0]/x)*TMath::Power(x/[1],[0])",0,e);
//        f->SetParameter(0,a);
//        f->SetParameter(1,e);
//        double energyloss = f->GetRandom();
//        LOG("RadiativeDecay", pNOTICE) << "energy loss " << energyloss;
//        return energyloss;
//}
//___________________________________________________________________________
void RadiativeCorrector::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void RadiativeCorrector::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void RadiativeCorrector::LoadConfig(void)
{
  //AlgConfigPool * confp = AlgConfigPool::Instance();
  //const Registry * gc = confp->GlobalParameterList();

  //int ndec = fConfig->GetIntDef("NDecayers",0);
  //assert(ndec>0);
  //GetParam( "NDecayers", fNDecayers);
  //assert(fNDecayers>0);

  //if(fDecayers) {
  //  fDecayers->clear();
  //  delete fDecayers;
  //}
  //fDecayers = new vector<const RadiativeDecayer *>(ndec);

  //for(int idec = 0; idec < ndec; idec++) {
  //   ostringstream alg_key;
  //   alg_key     << "Decayer-" << idec;
  //   const RadiativeDecayer * decayer =
  //            dynamic_cast<const RadiativeDecayer *> (this->SubAlg(alg_key.str()));
  //   (*fDecayers)[idec] = decayer;
  //}
  
  //fISR = fConfig->GetBool("ISR");
  GetParam( "ISR",fISR);

  //Target target(1, 1);     //ADITODO change 
  //InitialState init_state(target, 11); // ADITODO change
  
  //ostringstream mesg;
 // mesg << "Configuring for initial state: `"
   //    << init_state.AsString();

  //LOG("RadiativeCorrector", pNOTICE)
   //     << utils::print::PrintFramedMesg(mesg.str(), 0, '*');

  //this -> BuildInitialState            (init_state);
  //LOG("RadiativeCorrector", pINFO) << "Done configuring. \n"; 

}
//____________________________________________________________________________
void RadiativeCorrector::Configure(const InitialState & is)
{
  /*InitialState init_state(is.TgtPdg(), is.ProbePdg()); // filter any other init state info

  ostringstream mesg;
  mesg << "Configuring for initial state: `"
       << init_state.AsString();

  LOG("RadiativeCorrector", pNOTICE)
        << utils::print::PrintFramedMesg(mesg.str(), 0, '*');

  this -> BuildInitialState            (init_state);
  //this -> BuildGeneratorList           ();
  //this -> BuildInteractionGeneratorMap ();
  //this -> BuildInteractionSelector     ();

  LOG("RadiativeCorrector", pINFO) << "Done configuring. \n";*/
}
//____________________________________________________________________________
/*void RadiativeCorrector::SetE(double e) 
{
	r_e = e;
}*/
//____________________________________________________________________________
