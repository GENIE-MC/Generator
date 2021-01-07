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
#include <TCanvas.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Physics/Common/RadiativeCorrector.h"
#include "Physics/Decay/UnstableParticleDecayer.h"
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
  fInitState = 0;
  //rad_gOptRunNu = 0;
  //rad_kDefOptNtpFormat = kNFGHEP;
}
//___________________________________________________________________________
RadiativeCorrector::RadiativeCorrector(string config) :
EventRecordVisitorI("genie::RadiativeCorrector", config)
{
}
//___________________________________________________________________________
RadiativeCorrector::~RadiativeCorrector()
{
  if (fInitState)        delete fInitState;
  //rad_gOptRunNu = 0;
  //rad_kDefOptNtpFormat = kNFGHEP;
}
//___________________________________________________________________________
void RadiativeCorrector::BuildInitialState(const InitialState & init_state)
{

  if(fInitState) delete fInitState;
  fInitState = new InitialState(init_state);

  //this->AssertIsValidInitState();
}
//___________________________________________________________________________
void RadiativeCorrector::ProcessEventRecord(GHepRecord * evrec) const 
{
  LOG("RadiativeCorrector", pDEBUG) << "Welcome to the Radiative corrector using the "<<fModel<<" model";
  if (fISR) LOG("RadiativeCorrector", pDEBUG) <<"This is ISR";
  else LOG("RadiativeCorrector", pDEBUG) <<"This is FSR";
  // decay a photon with dE from the electron
  evrec->SetPrintLevel(13); 
  //evrec->Print();
  //std::cout<<"\n";
  bool radDone = false;
 
  Interaction * interaction = evrec->Summary();
  Kinematics * kine = interaction->KinePtr();
  InitialState * init_state_ptr = interaction -> InitStatePtr();
  Target * target_ptr = init_state_ptr -> TgtPtr();
  int Z = target_ptr->Z();
  //double thickness = 0.1;
  //map<int,double>::const_iterator it = fThicknesses.find(Z);
  //if(it != fThicknesses.end()) thickness = it->second;
  //const InitialState & const_init_state = interaction -> InitState();
  //InitialState init_state(const_init_state);
  TObjArrayIter piter(evrec);
  GHepParticle * p = 0;
  unsigned int ipos = 0;
  //unsigned int hitNucleonPosition; 
  while( (p = (GHepParticle *) piter.Next()) )
  {
     if( this->ToBeDecayed(p)) {
        TLorentzVector p4 = *(p->P4());
        TLorentzVector x4 = *(p->X4());
	double e = p4.E();
	double energyLoss = 0.;
        double e_gamma_max;
	double e_gamma_min = 1E-25;
        //if (fISR) e_gamma_max = init_state_ptr->ProbeE(kRfLab) - fP4l.E();
        //else 
	//{
	//  e_gamma_max = init_state_ptr->ProbeE(kRfLab) - kine->FSLeptonP4().E();
	//  if (e_gamma_max > kine->FSLeptonP4().E()) e_gamma_max = kine->FSLeptonP4().E();
        //}
        e_gamma_max = 0.2*p4.E();
        LOG("RadiativeCorrector", pDEBUG) << " particle for decay "<<p->Pdg()<<" e_gamma_max "<<e_gamma_max;
	if (fModel == "vanderhagen") {
           if (e_gamma_max<0) e_gamma_max = 1E-10;
           double a;
 	   if (fISR) a = (kAem/kPi)*(TMath::Log(fQ2/pow(kElectronMass,2)) - 1.);
	   else a = (kAem/kPi)*(TMath::Log(kine->Q2(true)/pow(kElectronMass,2)) - 1.);
	   LOG("RadiativeCorrector", pDEBUG) << "\n ipos "<<ipos<<" probE  "<<init_state_ptr->ProbeE(kRfLab)<<" evrec->CorrectProbe()->GetP4()->E() "<<evrec->CorrectProbe()->GetP4()->E()<<" energy loss limit "<<e_gamma_max<<" a parameter "<<a<<" the radiated lepton e "<<e;
           TF1 *f = new TF1("f","([0]/x)*TMath::Power(x/[1],[0])",e_gamma_min,e_gamma_max);
           f->SetParameter(0,a);
           f->SetParameter(1,e);
           energyLoss = f->GetRandom();
	   delete f;
	   LOG("RadiativeCorrector", pINFO) << "Vanderhagen Energy loss is "<<energyLoss;
	}
	double L1,L2,b,lambda_e,g,power_hi,power_lo;
	if (fModel =="simc" || fModel == "simple") {
 
	   if (Z==1) {
	      L1 = 5.31;
	      L2 = 6.144; 
	   }
	   else {
              L1 = TMath::Log(184.15) - (1./3)*TMath::Log(Z);
	      L2 = TMath::Log(1194.) - (2./3)*TMath::Log(Z);
	   }
           b = (1./9)*(12 + float(Z+1)/(Z*L1 + L2));
	}
	if (fModel == "simple") {
	   double lambda = (kAem/kPi)*(2*TMath::Log(2*init_state_ptr->ProbeE(kRfLab)/kElectronMass) - 1) + b*fThickness;
	   TF1 *f = new TF1("f","[0]*pow(x,[0]-1)/[1]",e_gamma_min,e_gamma_max);
           f->SetParameter(0,lambda);
           f->SetParameter(1,pow(init_state_ptr->ProbeE(kRfLab),-1.*lambda));
           energyLoss = f->GetRandom();
           delete f;
           LOG("RadiativeCorrector", pDEBUG) << "Simple Energy loss is "<<energyLoss;    
	}
	if (fModel == "simc") {
	   if (fISR) lambda_e = (kAem/kPi)*( 2*TMath::Log(2*p->P4()->P()/kElectronMass) -1 + TMath::Log(0.5*(1-fP4l.CosTheta())) );//+ 2*TMath::Log(init_state_ptr->GetProbeP4(kRfLab)->P()/fP4l.E()) + TMath::Log(0.5*(1-fP4l.CosTheta() ) ) );
	   else lambda_e =      (kAem/kPi)*( 2*TMath::Log(2*p->P4()->P()/kElectronMass) -1 + TMath::Log(0.5*(1-kine->FSLeptonP4().CosTheta())) );//+ 2*TMath::Log(init_state_ptr->GetProbeP4(kRfLab)->P()/kine->FSLeptonP4().E()) + TMath::Log(0.5*(1-kine->FSLeptonP4().CosTheta() ) ) );
           g = b*fThickness + lambda_e;
	   LOG("RadiativeCorrector", pINFO) << "SIMC chosen Thickness "<<fThickness<<" lambda e  is "<<lambda_e<<" g "<<g<<" b "<<b;

	   power_hi = pow(e_gamma_max,g);
	   //power_lo  = pow(e_gamma_min,g);
	   power_lo  = pow(fCutoff,g);
	   //TRandom3 rnd;
	   //rnd.SetSeed(0);
	   //double ymin = power_lo/power_hi;
	   //std::cout<<"ymin "<<ymin<<std::endl; 
	   //double y = ymin + rnd.Rndm()*(1.-ymin);
	   //std::cout<<"y "<<y<<std::endl;
	   //double x = pow(y,1./g);
	   //energyLoss = x*e_gamma_max;
	   //double c_int = lambda_e/pow(init_state_ptr->GetProbeP4(kRfLab)->E()*kine->FSLeptonP4().E(),lambda_e/2);
	   TF1 *f = new TF1("f","[0]*pow(x,[0]-1)/[1]",e_gamma_min,e_gamma_max);
           f->SetParameter(0,g);
           f->SetParameter(1,power_hi - power_lo);
           energyLoss = f->GetRandom();
           delete f;
           LOG("RadiativeCorrector", pINFO) << "SIMC Energy loss is "<<energyLoss;

	}
	//std::cout<<"CHECK energyLoss "<<energyLoss<<" fCutoff "<<fCutoff<<" original energy is "<<p4.E()<<" the loss is "<<energyLoss*100/p4.E()<<"%"<<std::endl;
	if (energyLoss<fCutoff) energyLoss = 0.; 

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
	TLorentzVector p4tag = p4 - p4RadGamma;


	if (fISR && !radDone && energyLoss>0) {
	  if(fModel =="simc") {
            double C = g/(TMath::Gamma(1+b*fThickness)*pow(p4.P(),b*fThickness)*pow(p4.P()*p4tag.P(),lambda_e/2)); 
            double W_rad_e = (C/g)*(power_hi-power_lo);
            //double Phi_ext_e = 1. - ( (b*fThickness/2)/(b*fThickness/2 + lambda_e) )*( p4.E()/p4.P()  );
            double Phi_ext_e = 1. - ( (b*fThickness) / p4.E() / g * p4RadGamma.E()) ;
            double radcor_weight = W_rad_e*Phi_ext_e;
	    if( fDoInternal) evrec->SetWeight(evrec->Weight() * radcor_weight);
            //std::cout<<"weights C "<<C<<" g "<<g<<" W_rad_e "<<W_rad_e<<" Phi_ext_e "<<Phi_ext_e<<std::endl; 
	    //std::cout<<"weights b "<<b<<" fThickness "<<fThickness<<" p4.E() "<<p4.E()<<" p4RadGamma.E() "<<p4RadGamma.E()<<std::endl;
            LOG("RadiativeCorrector", pINFO) << "Applying ISR part of the radiative correction weight "<<evrec->Weight() * radcor_weight;
          }
	  
	  if (energyLoss<fCutoff) std::cout<<"CHECK BAD "<<std::endl;
	  LOG("RadiativeCorrector", pDEBUG) << "performing ISR correction for: " << p->Name() << " reduced energy is : "<<p4tag.E();
	  // changing the probe energy for the initial state 
	  init_state_ptr->SetProbeP4(p4tag);
	  //-- Mark it as a 'decayed state' & add its daughter links
	  //-- Add the mom & daughters to the event record
	  LOG("RadiativeCorrector", pDEBUG) << "Adding daughter... PDG=" << p->Pdg();
          evrec->AddParticle(p->Pdg(), kIStCorrectedProbe, ipos,-1,-1,-1, p4tag, x4);
          LOG("RadiativeCorrector", pDEBUG) << "Adding daughter... PDG= 22";
          evrec->AddParticle(22, kIStStableFinalState, ipos,-1,-1,-1, p4RadGamma, x4);
	  //evrec->SetWeight(evrec->Weight() * 1.000000001); //changing the weight to identify events with radiated photon
	  radDone = true;
          //evrec->Print();
	}
	
        if (!fISR && !radDone) {
 	   //-- Update the event weight for each weighted particle decay
	   float radcor_weight = 1.;
	   if(fModel == "vanderhagen" ) 
	   {
	   	TF1 *fsp = new TF1("fsp","TMath::Log(1-x)* (1./x)");
	   	double SP = -1*fsp->Integral(e_gamma_min,pow(cos(p4tag.Theta())/2,2));
		delete fsp;
	   	radcor_weight = 1. + Z*(kAem/kPi)*( (TMath::Log(kine->Q2(true)/pow(kElectronMass,2))-1)*TMath::Log(pow(init_state_ptr->ProbeE(kRfLab),2)/(init_state_ptr->ProbeE(kRfLab)*kine->FSLeptonP4().E())) 
					             + (13./6.)*TMath::Log(kine->Q2(true)/pow(kElectronMass,2))
					             - (28./9.)
					             - 0.5*pow(TMath::Log(init_state_ptr->ProbeE(kRfLab)/kine->FSLeptonP4().E()),2)
					             - pow(kPi,2)/6. 
					             + SP );
   		LOG("RadiativeCorrector", pINFO) << "radcor_weight "<<radcor_weight;
	   }

	   if(fModel=="simple") 
	   {
		//double QSq = 2.*init_state_ptr->ProbeE(kRfLab)* kine->FSLeptonP4().E() * (1.-TMath::Cos(p4tag.Theta()));
  		//double eta = init_state_ptr->ProbeE(kRfLab)/kine->FSLeptonP4().E();

		//radcor_weight = - (kAem /kPi) * ( 13.0/6.*TMath::Logkine->Q2(true)/pow(kElectronMass,2)) 
		//			        - 28.0/9. 
		//				- 0.5*TMath::Power(TMath::Log(eta),2.)
		//	 			+ TMath::DiLog(TMath::Power(TMath::Cos(0.5*theta_e),2.)) - M_PI*M_PI/6.);
	
		radcor_weight = 1 + (2*kAem /kPi) * ( (13./12.)* (TMath::Log(kine->Q2(true)/pow(kElectronMass,2)) - 1)
						    - (17./36.)
						    - (1./4.) * pow(TMath::Log(init_state_ptr->ProbeE(kRfLab)*kine->FSLeptonP4().E()),2)
						    - (1./2.) * ( (pow(kPi,2)/6) -  TMath::DiLog(TMath::Power(TMath::Cos(0.5*p4tag.Theta()),2.)) ) );


		std::cout<<"(13./12.)* (TMath::Log(kine->Q2(true)/pow(kElectronMass,2)) - 1) "<<(13./12.)* (TMath::Log(kine->Q2(true)/pow(kElectronMass,2)) - 1)<<std::endl;
                std::cout<<"(1./4.) * pow(TMath::Log(init_state_ptr->ProbeE(kRfLab)*kine->FSLeptonP4().E()),2) "<< (1./4.) * pow(TMath::Log(init_state_ptr->ProbeE(kRfLab)*kine->FSLeptonP4().E()),2)<<std::endl;
                std::cout<<"TMath::DiLog(TMath::Power(TMath::Cos(0.5*p4tag.Theta()),2.)) " <<TMath::DiLog(TMath::Power(TMath::Cos(0.5*p4tag.Theta()),2.))<<std::endl;
		std::cout<<"radcor_weight "<<radcor_weight<<std::endl;
	   }

           if(fModel=="simc")
	   {
                //double tot_delta    = (1./(3*kPi))*( -5./3 + TMath::Log(kine->Q2(true)/pow(kElectronMass,2)) );
		//double delta_hard   = 2*kAem*( (-3./(4*kPi))*TMath::Log(kine->Q2(true)/pow(kElectronMass,2)) + 1./kPi - 2*tot_delta );
		double delta_hard   = -1.*(kAem/kPi)*(-28./9.+(13./6.)*TMath::Log(kine->Q2(true)/pow(kElectronMass,2)));
		//double c          = lambda_e/pow(evrec->CorrectProbe()->GetP4()->E() * kine->FSLeptonP4().E(),(lambda_e/2.));
		double C            = g/(TMath::Gamma(1+b*fThickness)*pow(p4.P(),b*fThickness)*pow(p4.P()*p4tag.P(),lambda_e/2));
                double W_rad_el     = (C/g)*(power_hi-power_lo);
                //double Phi_ext_el = 1. - ( (b*fThickness/2)/(b*fThickness/2 + lambda_e) )*( p4.E()/p4.P() );
                double Phi_ext_el   = 1. - ( (b*fThickness) / p4.E() / g * p4RadGamma.E()) ;
		double W_rad_pl     = 1.;
                double Phi_ext_pl   = 1.;
                if (fDoInternal) radcor_weight       = W_rad_el*W_rad_pl*Phi_ext_el*Phi_ext_pl*(1-delta_hard);
		else radcor_weight       = Phi_ext_el*Phi_ext_pl;
                //radcor_weight       = W_rad_el*W_rad_pl*Phi_ext_el*Phi_ext_pl;
		//std::cout<<"Weight: log "<<TMath::Log(kine->Q2(true)/pow(kElectronMass,2))<<" () "<<(-28./9.+(13./6.)*TMath::Log(kine->Q2(true)/pow(kElectronMass,2)))<<" E2 "<<pow(p4.E(),2)<<" delta_hard "<<delta_hard<<std::endl;
	   	//std::cout<<"C "<<C<<" g "<<g<<" W_rad_el "<<W_rad_el<<" W_rad_pl "<<W_rad_pl<<" Phi_ext_el"<<Phi_ext_el<<" delta hard "<<delta_hard<<std::endl;

		// crisis attempt 
		double lambda_e_ISR   = (kAem/kPi)*( 2*TMath::Log(2*init_state_ptr->GetProbeP4(kRfLab)->P()/kElectronMass) -1 + TMath::Log(0.5*(1-fP4l.CosTheta())) );
                double lambda_e_FSR   =      (kAem/kPi)*( 2*TMath::Log(2*p->P4()->P()/kElectronMass) -1 + TMath::Log(0.5*(1-kine->FSLeptonP4().CosTheta())) );
                double g_ISR          = b*fThickness + lambda_e_ISR;
                double g_FSR          = b*fThickness + lambda_e_FSR;
		double extrad_phi_ISR = 1. - (b*fThickness/init_state_ptr->ProbeE(kRfLab) + b*fThickness/kine->FSLeptonP4().E() ) / (g_ISR+g_FSR) * (init_state_ptr->ProbeE(kRfLab) - evrec->CorrectProbe()->GetP4()->P() );
	        double extrad_phi_FSR = 1. - b*fThickness/kine->FSLeptonP4().E() / g_FSR * p4RadGamma.E();
		double C_ISR = b*fThickness/pow(init_state_ptr->ProbeE(kRfLab),b*fThickness)/TMath::Gamma(1+b*fThickness);
		double C_FSR = b*fThickness/pow(kine->FSLeptonP4().E()        ,b*fThickness)/TMath::Gamma(1+b*fThickness);
		double g_ext = 2*b*fThickness; 
		double c_ext = C_ISR*C_FSR * g_ext; //c_ext(1)*c_ext(2) * g_ext / bt(1)/bt(2)
		c_ext = c_ext * TMath::Gamma(1+b*fThickness) * TMath::Gamma(1+b*fThickness) / TMath::Gamma(1 + g_ext);
		//std::cout<<"Weight try extrad_phi_ISR "<<extrad_phi_ISR<<" extrad_phi_FSR "<<extrad_phi_FSR<<" c_ext "<<c_ext<<std::endl;
		//std::cout<<" C_ISR "<<C_ISR<<" C_FSR "<<C_FSR<<" g_ext "<<g_ext<<std::endl;
		//bei= (-1./(2.*pi))*log(ein/egamma)
		//bef= (-1./(2.*pi))*log(eout/egamma)
		//adot= ein*eout*(1.-cos(eang))
		//alpha= 2.*ame**2 - 2.*adot
		//ar1= 0.5+sqrt(adot**2 - ame**4)/alpha
		//ar2= 0.5-sqrt(adot**2 - ame**4)/alpha
		//bee= -1.*adot*inter(calculate_spence,alpha,ar1,ar2,ak,akp,de)
		//b= 2.*e2*(bei+bef+bee)
	        //double peaked_rad_weight = c_ext/g_ext * (exp(-dsoft_intmax)*emax**g_ext - exp(-dsoft_intmin)*emin**g_ext)
		//peaked_rad_weight = peaked_rad_weight * exp(-eul*g(4))/gamma(1.+g(4))
    		//	* gamma(1.+g(4)-bt(1)-bt(2))*gamma(1.+bt(1))
    		//	* gamma(1.+bt(2))/gamma(1.+g(4))

	   }

	   if (energyLoss>0) {
	       LOG("RadiativeCorrector", pDEBUG) << "performing FSR correction for: " << p->Name();
               //-- Mark it as a 'decayed state' & add its daughter links
       	       p->SetStatus(kIStDecayedState);
               ////-- Add the mom & daughters to the event record
               LOG("RadiativeCorrector", pDEBUG) << "Adding daughter... PDG=" << p->Pdg();
	       evrec->AddParticle(p->Pdg(), kIStStableFinalState, ipos,-1,-1,-1, p4tag, x4);
               LOG("RadiativeCorrector", pDEBUG) << "Adding daughter... PDG= 22";
	       evrec->AddParticle(22, kIStStableFinalState, ipos,-1,-1,-1, p4RadGamma, x4);
               ////-- Printouts
	       //evrec->Print();
	       //std::cout<<"\n";
	       radDone = true;
	   }

	   //bool radiatedPhoton = false;
	   //GHepParticle * p1 = 0;
           //while( (p1 = (GHepParticle *) piter.Next()) ) { if (ParticleWasRadiated(p1)) radiatedPhoton = true;}

 	   //if ( energyLoss > 0 || evrec->Weight() > 1)
	   //{ LOG("RadiativeCorrector", pDEBUG) << "Applying radiative correction weight "<<evrec->Weight()<<" radcor_weight " <<radcor_weight;
	     std::cout << "Applying radiative correction weight "<<evrec->Weight()<<" radcor_weight " <<radcor_weight<<std::endl;
	     evrec->SetWeight(evrec->Weight() * radcor_weight);
	   //}
	}
	/* //code no longer needed now with FSL before and after radiation are keptand PrimaryLeptonPosition method is changed  
           if (!fISR && radDone && (p->Status() == kIStStableFinalState) && (p->Pdg()==11)) {    
           LOG("RadiativeCorrector", pINFO) <<"FSR stable electron in pos "<<ipos; 
	   //-- Correct the final state lepton 
           evrec->CorrectProbe()->SetFirstDaughter(ipos);
           evrec->CorrectProbe()->SetLastDaughter(-1);
           LOG("RadiativeCorrector", pINFO) <<"After setting the lepton first and last mother to -1  "<<evrec->FinalStatePrimaryLepton()->P4()->E();
           evrec->Print(); 
	}*/
     }
     //if (fISR && (p->Status() ==1)) p->SetStatus(kIStIntermediateState);
     if (fISR && (p->Status() ==1) && (p->Pdg()==11)) {
        //if (p->Pdg() == 2212 || p->Pdg() == 2112 || p->Pdg() == 1000000010 || p->Pdg() == 1000010000 ) {
	//	p->SetMomentum(0.,0.,0.,p->Mass());
	//	//p->SetEnergy(p->Mass());
	//}
	//-- Correct the final state lepton 
        //evrec->CorrectProbe()->SetFirstDaughter(-1);
        //evrec->CorrectProbe()->SetLastDaughter(-1);
        
	p->SetStatus(kIStDecayedState);
	LOG("RadiativeCorrector", pDEBUG) <<"After setting the particles status to decyed and momentum of the nucleon to be zero, pdg "<<p->Pdg()<<" mass "<<p->Mass();
        //evrec->Print();
	//std::cout<<"\n";
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
   LOG("RadiativeCorrector", pDEBUG) <<"Particle ToBeDecayed pdg "<<particle->Pdg()<<" status "<<particle->Status()<<" mother "<<particle->FirstMother();
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
bool RadiativeCorrector::ParticleWasRadiated(GHepParticle * particle) const
{
  std::cout<<"particle was radiated pdg "<<particle->Pdg()<<" status "<<particle->Status()<<" mother "<<particle->FirstMother() <<std::endl;
  if ( (particle->Pdg() ==11) && ((particle->Status() == 5) || (particle->Status() ==3))) return true;
  else return false;
}
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
  GetParam( "ISR",fISR);
  GetParam( "RadiativeCorrectionModel" , fModel);
  GetParam( "RadiativeCorrectionThickness_1000010010",fThickness_1000010010);
  GetParam( "RadiativeCorrectionThickness_1000020040",fThickness_1000020040);
  GetParam( "RadiativeCorrectionThickness_1000060120",fThickness_1000060120);
  GetParam( "RadiativeCorrectionThickness_1000260560",fThickness_1000260560);
  GetParam( "RadiativeCorrectionCutoff",fCutoff);
  GetParam( "RadiativeCorrectionDoInternal",fDoInternal);
  
  /*
  const std::string keyStart = "Thickness@Pdg=";

  RgIMap entries = GetConfig().GetItemMap();

  for(RgIMap::const_iterator it = entries.begin(); it != entries.end(); ++it){
    const std::string& key = it->first;
    int pdg = 0;
    int Z = 0;
    if (0 == key.compare(0, keyStart.size(), keyStart.c_str())) {
      pdg = atoi(key.c_str() + keyStart.size());
      Z = pdg::IonPdgCodeToZ(pdg);
    }
    if (0 != pdg && 0 != Z) {
      ostringstream key_ss ;
      key_ss << keyStart << pdg;
      RgKey rgkey   = key_ss.str();
      double thickness ;
      GetParam( rgkey, thickness ) ;
      thickness = TMath::Max(thickness, 0.);
      LOG("RadiativeCorrector", pINFO) << "Nucleus: " << pdg << " -> using Thickness =  " << thickness << " radiation lengths";
      fThicknesses.insert(map<int,double>::value_type(Z,thickness));
    }
  }*/
}
//____________________________________________________________________________
void RadiativeCorrector::SetISR(bool isr)
{
  fISR = isr;
}
//____________________________________________________________________________
void RadiativeCorrector::SetModel(std::string model)
{
  fModel = model;
}
//____________________________________________________________________________
void RadiativeCorrector::SetQ2(double Q2)
{
  fQ2 = Q2;
}
//____________________________________________________________________________
void RadiativeCorrector::SetP4l(TLorentzVector p4l)
{
  fP4l = p4l;
}
//____________________________________________________________________________
void RadiativeCorrector::SetCutoff(double cutoff)
{
  fCutoff = cutoff;
}
//____________________________________________________________________________
void RadiativeCorrector::SetThickness(int tgtpdg, double thickness_1000010010, double thickness_1000020040, double thickness_1000060120, double thickness_1000260560)
{
  if      (tgtpdg == 1000010010) fThickness = thickness_1000010010;
  else if (tgtpdg == 1000020040) fThickness = thickness_1000020040;
  else if (tgtpdg == 1000060120) fThickness = thickness_1000060120;
  else if (tgtpdg == 1000260560) fThickness = thickness_1000260560;
}
//____________________________________________________________________________
void RadiativeCorrector::SetDoInternalRad(bool doInternal)
{
  fDoInternal = doInternal;
}
//____________________________________________________________________________
void RadiativeCorrector::Configure(const InitialState & is)
{
  /*InitialState init_state(is.TgtPdg(), is.CorrectProbePdg()); // filter any other init state info

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
