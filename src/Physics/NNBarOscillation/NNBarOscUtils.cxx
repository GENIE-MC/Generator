//____________________________________________________________________________
/*
	Copyright (c) 2003-2023, The GENIE Collaboration
	For the full text of the license visit http://copyright.genie-mc.org

	Jeremy Hewes, Georgia Karagiorgi
	University of Manchester

	Linyan Wan
	Fermilab
	*/
//____________________________________________________________________________

#include "Framework/Messenger/Messenger.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/PrintUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Physics/NNBarOscillation/NNBarOscUtils.h"

using namespace genie;

//____________________________________________________________________________
string genie::utils::nnbar_osc::AsString(NNBarOscMode_t ndm)
{
  // this just maps the decay mode integers to string descriptors. replaced. -j
  switch(ndm) {
	 case (kNORandom)          :    return "Random mode";
											  break;
	 case (kNOpto1pip1pi0)     :    return "p + nbar --> pi+ pi0";
											  break;
	 case (kNOpto1pip2pi0)     :    return "p + nbar --> pi+ 2pi0";
											  break;
	 case (kNOpto1pip3pi0)     :    return "p + nbar --> pi+ 3pi0";
											  break;
	 case (kNOpto1pip4pi0)     :    return "p + nbar --> pi+ 4pi0";
											  break;
	 case (kNOpto2pip1pim) 		:    return "p + nbar --> 2pi+ pi-";
											  break;
	 case (kNOpto2pip1pim1pi0) :    return "p + nbar --> 2pi+ pi- pi0";
											  break;
	 case (kNOpto2pip1pim2pi0) :    return "p + nbar --> 2pi+ pi- 2pi0";
											  break;
	 case (kNOpto2pip1pim3pi0) :    return "p + nbar --> 2pi+ pi- 3pi0";
											  break;
	 case (kNOpto3pip2pim)		:    return "p + nbar --> 3pi+ 2pi";
											  break;
	 case (kNOpto3pip2pim1pi0) :    return "p + nbar --> 3pi+ 2pi- pi0";
											  break;
	 case (kNOpto1pip1pi01o)   :    return "p + nbar --> pi+ pi0 omega";
											  break;
	 case (kNOpto2pip1pim1o)   :    return "p + nbar --> 2pi+ pi- omega";
											  break;

	 case (kNOnto2pi0)         :    return "n + nbar --> 2pi0";
											  break;
	 case (kNOnto3pi0)         :    return "n + nbar --> 3pi0";
											  break;
	 case (kNOnto4pi0)         :    return "n + nbar --> 4pi0";
											  break;
	 case (kNOnto5pi0)         :    return "n + nbar --> 5pi0";
											  break;
	 case (kNOnto7pi0)         :    return "n + nbar --> 7pi0";
											  break;
	 case (kNOnto1pip1pim)     :    return "n + nbar --> pi+ pi-";
											  break;
	 case (kNOnto1pip1pim1pi0) :    return "n + nbar --> pi+ pi- pi0";
											  break;
	 case (kNOnto1pip1pim2pi0) :    return "n + nbar --> pi+ pi- 2pi0";
											  break;
	 case (kNOnto1pip1pim3pi0) :    return "n + nbar --> pi+ pi- 3pi0";
											  break;
	 case (kNOnto1pip1pim4pi0) :    return "n + nbar --> pi+ pi- 4pi0";
											  break;
	 case (kNOnto1pip1pim5pi0) :    return "n + nbar --> pi+ pi- 5pi0";
											  break;
	 case (kNOnto2pip2pim)     :    return "n + nbar --> 2pi+ 2pi-";
											  break;
	 case (kNOnto2pip2pim1pi0) :    return "n + nbar --> 2pi+ 2pi- pi0";
											  break;
	 case (kNOnto2pip2pim2pi0) :    return "n + nbar --> 2pi+ 2pi- 2pi0";
											  break;
	 case (kNOnto2pip2pim3pi0) :    return "n + nbar --> 2pi+ 2pi- 3pi0";
											  break;
	 case (kNOnto3pip3pim) 		:    return "n + nbar --> 3pi+ 3pi-";
											  break;
	 case (kNOnto3pip3pim1pi0) :    return "n + nbar --> 3pi+ 3pi- 1pi0";
											  break;
	 case (kNOnto1rho1pi0)   	:    return "n + nbar --> rho0 pi0";
											  break;
	 case (kNOnto1rhopm1pimp)  :    return "n + nbar --> rho+- pi-+";
											  break;
	 case (kNOnto2o)   			:	  return "n + nbar --> 2omega";
											  break;
	 case (kNOnto1rho1o)   		:    return "n + nbar --> rho0 omega0";
											  break;
	 case (kNOnto2pi01o)   		:    return "n + nbar --> 2pi0 omega0";
											  break;
	 case (kNOnto1pip1pim1o)   :    return "n + nbar --> pi+ pi- omega0";
											  break;
	 case (kNOnto1eta1o)   		:    return "n + nbar --> eta0 omega0";
											  break;
	 case (kNOnto1pip1pim1eta) :    return "n + nbar --> pi+ pi- eta0";
											  break;
	 case (kNOnto1Kp1Km)   		:    return "n + nbar --> K+ K-";
											  break;
	 case (kNOnto1Kp1Km1o)   	:    return "n + nbar --> K+ K- omega0";
											  break;
	 case (kNOnto2Ks1o)   	:    return "n + nbar --> Ks Ks omega0";
										  break;
	 case (kNOnto1pipm1Kmp1Kl) :    return "n + nbar --> pi+- K-+ Kl";
											  break;
	 case (kNOnto1pipm1Kmp1Ks) :    return "n + nbar --> pi+- K-+ Ks";
											  break;
	 case (kNOnto1pipm1Kmp1pi0K0):  return "n + nbar --> pi+- K-+ pi0 K0";
											  break;
	 case (kNOnto1pip1pim1Ks1Kl):    return "n + nbar --> pi+ pi- Ks Kl";
												break;

	 default                   :    return "?";
											  break;
  }
  return "??";
}
//____________________________________________________________________________
bool genie::utils::nnbar_osc::IsValidMode(NNBarOscMode_t ndm)
{
  // checks if a mode is valid. just straight replaced. -j
  switch(ndm) {
	 case (kNORandom)          :
	 case (kNOpto1pip1pi0)     :
	 case (kNOpto1pip2pi0)     :
	 case (kNOpto1pip3pi0)     :
	 case (kNOpto1pip4pi0)     :
	 case (kNOpto2pip1pim) :
	 case (kNOpto2pip1pim1pi0) :
	 case (kNOpto2pip1pim2pi0) :
	 case (kNOpto2pip1pim3pi0) :
	 case (kNOpto3pip2pim) :
	 case (kNOpto3pip2pim1pi0) :
	 case (kNOpto1pip1pi01o)   :
	 case (kNOpto2pip1pim1o)   :

	 case (kNOnto2pi0)         :    
	 case (kNOnto3pi0)         :    
	 case (kNOnto4pi0)         :    
	 case (kNOnto5pi0)         :    
	 case (kNOnto7pi0)         :    
	 case (kNOnto1pip1pim)     :    
	 case (kNOnto1pip1pim1pi0) :    
	 case (kNOnto1pip1pim2pi0) :    
	 case (kNOnto1pip1pim3pi0) :    
	 case (kNOnto1pip1pim4pi0) :    
	 case (kNOnto1pip1pim5pi0) :    
	 case (kNOnto2pip2pim)     :    
	 case (kNOnto2pip2pim1pi0) :    
	 case (kNOnto2pip2pim2pi0) :    
	 case (kNOnto2pip2pim3pi0) :    
	 case (kNOnto3pip3pim) 		:    
	 case (kNOnto3pip3pim1pi0) :    
	 case (kNOnto1rho1pi0)   	:    
	 case (kNOnto1rhopm1pimp)  :    
	 case (kNOnto2o)   			:	  
	 case (kNOnto1rho1o)   		:    
	 case (kNOnto2pi01o)   		:    
	 case (kNOnto1pip1pim1o)   :    
	 case (kNOnto1eta1o)   		:    
	 case (kNOnto1pip1pim1eta) :    
	 case (kNOnto1Kp1Km)   		:    
	 case (kNOnto1Kp1Km1o)   	:    
	 case (kNOnto2Ks1o)   	:    
	 case (kNOnto1pipm1Kmp1Kl) :    
	 case (kNOnto1pipm1Kmp1Ks) :    
	 case (kNOnto1pipm1Kmp1pi0K0):  
	 case (kNOnto1pip1pim1Ks1Kl):    

		return true;
		break;
	 default :
		return false;
		break;
  }
  return false;
}
//____________________________________________________________________________
int genie::utils::nnbar_osc::AnnihilatingNucleonPdgCode(NNBarOscMode_t ndm)
{
  // name isn't really accurate any more. instead of decayed nucleon, function
  // returns what particle the oscillated neutron annihilated with -j
  switch(ndm) {
	 case (kNOpto1pip1pi0)     : return kPdgProton;  break;
	 case (kNOpto1pip2pi0)     : return kPdgProton;  break;
	 case (kNOpto1pip3pi0)     : return kPdgProton;  break;
	 case (kNOpto1pip4pi0)     : return kPdgProton;  break;
	 case (kNOpto2pip1pim)		: return kPdgProton;  break;
	 case (kNOpto2pip1pim1pi0) : return kPdgProton;  break;
	 case (kNOpto2pip1pim2pi0) : return kPdgProton;  break;
	 case (kNOpto2pip1pim3pi0) : return kPdgProton;  break;
	 case (kNOpto3pip2pim)		: return kPdgProton;  break;
	 case (kNOpto3pip2pim1pi0) : return kPdgProton;  break;
	 case (kNOpto1pip1pi01o)   : return kPdgProton;  break;
	 case (kNOpto2pip1pim1o)   : return kPdgProton;  break;

	 case (kNOnto2pi0)         :  return kPdgNeutron; break;  
	 case (kNOnto3pi0)         :  return kPdgNeutron; break;  
	 case (kNOnto4pi0)         :  return kPdgNeutron; break;  
	 case (kNOnto5pi0)         :  return kPdgNeutron; break;  
	 case (kNOnto7pi0)         :  return kPdgNeutron; break;  
	 case (kNOnto1pip1pim)     :  return kPdgNeutron; break;  
	 case (kNOnto1pip1pim1pi0) :  return kPdgNeutron; break;  
	 case (kNOnto1pip1pim2pi0) :  return kPdgNeutron; break;  
	 case (kNOnto1pip1pim3pi0) :  return kPdgNeutron; break;  
	 case (kNOnto1pip1pim4pi0) :  return kPdgNeutron; break;  
	 case (kNOnto1pip1pim5pi0) :  return kPdgNeutron; break;  
	 case (kNOnto2pip2pim)     :  return kPdgNeutron; break;  
	 case (kNOnto2pip2pim1pi0) :  return kPdgNeutron; break;  
	 case (kNOnto2pip2pim2pi0) :  return kPdgNeutron; break;  
	 case (kNOnto2pip2pim3pi0) :  return kPdgNeutron; break;  
	 case (kNOnto3pip3pim) 		:  return kPdgNeutron; break;  
	 case (kNOnto3pip3pim1pi0) :  return kPdgNeutron; break;  
	 case (kNOnto1rho1pi0)   	:  return kPdgNeutron; break;  
	 case (kNOnto1rhopm1pimp)  :  return kPdgNeutron; break;  
	 case (kNOnto2o)   			:	return kPdgNeutron; break;  
	 case (kNOnto1rho1o)   		:  return kPdgNeutron; break;  
	 case (kNOnto2pi01o)   		:  return kPdgNeutron; break;  
	 case (kNOnto1pip1pim1o)   :  return kPdgNeutron; break;  
	 case (kNOnto1eta1o)   		:  return kPdgNeutron; break;  
	 case (kNOnto1pip1pim1eta) :  return kPdgNeutron; break;  
	 case (kNOnto1Kp1Km)   		:  return kPdgNeutron; break;  
	 case (kNOnto1Kp1Km1o)   	:  return kPdgNeutron; break;  
	 case (kNOnto2Ks1o)   	:  return kPdgNeutron; break;  
	 case (kNOnto1pipm1Kmp1Kl) :  return kPdgNeutron; break;  
	 case (kNOnto1pipm1Kmp1Ks) :  return kPdgNeutron; break;  
	 case (kNOnto1pipm1Kmp1pi0K0):return kPdgNeutron; break;  
	 case (kNOnto1pip1pim1Ks1Kl):  return kPdgNeutron; break;  
	 default                   : return 0;           break;
  }
  return 0;
}
//____________________________________________________________________________
PDGCodeList genie::utils::nnbar_osc::DecayProductList(
	 NNBarOscMode_t ndm)
{
  // ok so i think this is the first function where a straight replacement
  // isn't gonna cut it. all the nucleon decay modes are two-body, but that is
  // painfully untrue for nnbar. i just threw aaaaaaallll of the final state
  // particles into the vector, so let's just hope for the best -j

  // need to implement a lorentz boost into rest frame of two nucleons -j

  bool allow_duplicate = true;
  PDGCodeList decay_products(allow_duplicate);

  switch(ndm) {
	 case (kNOpto1pip1pi0) :
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPi0);
		break;
	 case (kNOpto1pip2pi0) :
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		break;
	 case (kNOpto1pip3pi0) :
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		break;
	 case (kNOpto1pip4pi0) :
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		break;
	 case (kNOpto2pip1pim) :
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiM);
		break;
	 case (kNOpto2pip1pim1pi0) :
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiM);
		decay_products.push_back(kPdgPi0);
		break;
	 case (kNOpto2pip1pim2pi0) :
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiM);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		break;
	 case (kNOpto2pip1pim3pi0) :
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiM);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		break;
	 case (kNOpto3pip2pim) :
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiM);
		decay_products.push_back(kPdgPiM);
		break;
	 case (kNOpto3pip2pim1pi0) :
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiM);
		decay_products.push_back(kPdgPiM);
		decay_products.push_back(kPdgPi0);
		break;
	 case (kNOpto1pip1pi01o) :
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgomega);
		break;
	 case (kNOpto2pip1pim1o) :
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiM);
		decay_products.push_back(kPdgomega);
		break;

	 case (kNOnto2pi0) :
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		break;
	 case (kNOnto3pi0) :
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		break;
	 case (kNOnto4pi0) :
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		break;
	 case (kNOnto5pi0) :
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		break;
	 case (kNOnto7pi0) :
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		break;
	 case (kNOnto1pip1pim) :
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiM);
		break;
	 case (kNOnto1pip1pim1pi0) :
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiM);
		decay_products.push_back(kPdgPi0);
		break;
	 case (kNOnto1pip1pim2pi0) :
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiM);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		break;
	 case (kNOnto1pip1pim3pi0) :
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiM);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		break;
	 case (kNOnto1pip1pim4pi0) :
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiM);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		break;
	 case (kNOnto1pip1pim5pi0) :
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiM);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		break;
	 case (kNOnto2pip2pim) :
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiM);
		decay_products.push_back(kPdgPiM);
		break;
	 case (kNOnto2pip2pim1pi0) :
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiM);
		decay_products.push_back(kPdgPiM);
		decay_products.push_back(kPdgPi0);
		break;
	 case (kNOnto2pip2pim2pi0) :
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiM);
		decay_products.push_back(kPdgPiM);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		break;
	 case (kNOnto2pip2pim3pi0) :
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiM);
		decay_products.push_back(kPdgPiM);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		break;
	 case (kNOnto3pip3pim) :
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiM);
		decay_products.push_back(kPdgPiM);
		decay_products.push_back(kPdgPiM);
		break;
	 case (kNOnto3pip3pim1pi0) :
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiM);
		decay_products.push_back(kPdgPiM);
		decay_products.push_back(kPdgPiM);
		decay_products.push_back(kPdgPi0);
		break;
	 case (kNOnto1rho1pi0) :
		decay_products.push_back(kPdgRho0);
		decay_products.push_back(kPdgPi0);
		break;
	 case (kNOnto1rhopm1pimp) :
		decay_products.push_back(kPdgRhoP);
		decay_products.push_back(kPdgPiM);
		break;
	 case (kNOnto2o) :
		decay_products.push_back(kPdgomega);
		decay_products.push_back(kPdgomega);
		break;
	 case (kNOnto1rho1o) :
		decay_products.push_back(kPdgRho0);
		decay_products.push_back(kPdgomega);
		break;
	 case (kNOnto2pi01o) :
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgomega);
		break;
	 case (kNOnto1pip1pim1o) :
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiM);
		decay_products.push_back(kPdgomega);
		break;
	 case (kNOnto1eta1o) :
		decay_products.push_back(kPdgEta);
		decay_products.push_back(kPdgomega);
		break;
	 case (kNOnto1pip1pim1eta) :
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiM);
		decay_products.push_back(kPdgEta);
		break;
	 case (kNOnto1Kp1Km) :
		decay_products.push_back(kPdgKP);
		decay_products.push_back(kPdgKM);
		break;
	 case (kNOnto1Kp1Km1o) :
		decay_products.push_back(kPdgKP);
		decay_products.push_back(kPdgKM);
		decay_products.push_back(kPdgomega);
		break;
	 case (kNOnto2Ks1o) :
		decay_products.push_back(kPdgK0S);
		decay_products.push_back(kPdgK0S);
		decay_products.push_back(kPdgomega);
		break;
	 case (kNOnto1pipm1Kmp1Kl) :
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgKM);
		decay_products.push_back(kPdgK0L);
		break;
	 case (kNOnto1pipm1Kmp1Ks) :
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgKM);
		decay_products.push_back(kPdgK0S);
		break;
	 case (kNOnto1pipm1Kmp1pi0K0) :
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgKM);
		decay_products.push_back(kPdgPi0);
		decay_products.push_back(kPdgK0);
		break;
	 case (kNOnto1pip1pim1Ks1Kl) :
		decay_products.push_back(kPdgPiP);
		decay_products.push_back(kPdgPiM);
		decay_products.push_back(kPdgK0S);
		decay_products.push_back(kPdgK0L);
		break;
	 default :
		break;
  }
  return decay_products;
}
//____________________________________________________________________________
GHepStatus_t genie::utils::nnbar_osc::DecayProductStatus(
	 bool in_nucleus, int pdgc)
{
  // took out all the irrelevant particles -j
  if(in_nucleus) {
	 if( pdgc == kPdgPi0   ||
		  pdgc == kPdgPiM   ||
		  pdgc == kPdgPiP   ||
		  pdgc == kPdgRho0   ||
		  pdgc == kPdgRhoP   ||
		  pdgc == kPdgEta   ||
		  pdgc == kPdgKP   ||
		  pdgc == kPdgKM   ||
		  pdgc == kPdgK0    ||
		  pdgc == kPdgK0S   ||
		  pdgc == kPdgK0L   ||
		  pdgc == kPdgomega)
	 {
		return kIStHadronInTheNucleus;
	 }
  }

  return kIStStableFinalState;
}
//____________________________________________________________________________
