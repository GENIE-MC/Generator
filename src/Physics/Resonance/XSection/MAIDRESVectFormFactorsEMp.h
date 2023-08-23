//____________________________________________________________________________
/*!

  \class    genie::MAIDRESVectFormFactorsEMp

  \brief    The Helicity Amplitudes, for all baryon resonances, for Electro-
  Magnetic (EM) interactions on free neutrons, as computed in MAID
  parameterization.

  \author   Julia Tena Vidal <jtenavidal \at tauex.tau.ac.il>
  Tel Aviv Universtiy

  \created  January 2023

  \cpright  Copyright (c) 2023-2025, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org        
*/
//____________________________________________________________________________

#ifndef _MAID_RES_VFF_EM_P_H_
#define _MAID_RES_VFF_EM_P_H_

#include "Physics/Resonance/XSection/RESVectFormFactorsI.h"

namespace genie {

  class MAIDRESVectFormFactorsEMp : public RESVectFormFactorsI {

  public:
    MAIDRESVectFormFactorsEMp();
    MAIDRESVectFormFactorsEMp(string config);
    virtual ~MAIDRESVectFormFactorsEMp();

    void Configure( const Registry & config );
    void Configure( string param_set ) ; 
  
    RESVectFFAmplitude Compute( const Interaction interaction ) ;

  private:
    void LoadConfig(void) ; 

    double fAM0_P33_1232 ; 
    double fAE0_P33_1232 ; 
    double fAC0_P33_1232 ; 
    double fBetaM_P33_1232 ; 
    double fBetaE_P33_1232 ; 
    double fBetaC_P33_1232 ; 
    double fGammaM_P33_1232 ; 
    double fGammaE_P33_1232 ; 
    double fGammaC_P33_1232 ; 
    double fDC_P33_1232 ; 
    double fDipoleMass ;

    std::map<Resonance_t,double> fA120P ;
    std::map<Resonance_t,double> fA12Alpha1P ;
    std::map<Resonance_t,double> fA12Alpha2P ;
    std::map<Resonance_t,double> fA12Alpha4P ;
    std::map<Resonance_t,double> fA12BetaP ;

    std::map<Resonance_t,double> fS120P ;
    std::map<Resonance_t,double> fS12Alpha1P ;
    std::map<Resonance_t,double> fS12Alpha2P ;
    std::map<Resonance_t,double> fS12Alpha4P ;
    std::map<Resonance_t,double> fS12BetaP ;
 
    std::map<Resonance_t,double> fA320P ;
    std::map<Resonance_t,double> fA32Alpha1P ;
    std::map<Resonance_t,double> fA32Alpha2P ;
    std::map<Resonance_t,double> fA32Alpha4P ;
    std::map<Resonance_t,double> fA32BetaP ;
  };

}        // genie namespace
#endif
