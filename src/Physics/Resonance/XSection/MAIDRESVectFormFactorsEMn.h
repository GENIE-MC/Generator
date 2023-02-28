//____________________________________________________________________________
/*!

  \class    genie::MAIDRESVectFormFactorsEMn

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

#ifndef _MAID_RESVECTFF_EM_N_H_
#define _MAID_RESVECTFF_EM_N_H_

#include "Physics/Resonance/XSection/RESVectFormFactorsI.h"

namespace genie {

  class MAIDRESVectFormFactorsEMn : public RESVectFormFactorsI {

  public:
    MAIDRESVectFormFactorsEMn();
    MAIDRESVectFormFactorsEMn(string config);
    virtual ~MAIDRESVectFormFactorsEMn();

    void Configure(const Registry & config);
    void Configure( string param_set ) ; 

    RESVectFFAmplitude Compute( const Interaction interaction ) ;

  private:
    void LoadConfig(void) ; 

    // Defining constants from fit 
    double fAM0_P33_1232 ; 
    double fAE0_P33_1232 ; 
    double fAS0_P33_1232 ; 
    double fBetaM_P33_1232 ; 
    double fBetaE_P33_1232 ; 
    double fBetaS_P33_1232 ; 
    double fGammaM_P33_1232 ; 
    double fGammaE_P33_1232 ; 
    double fGammaS_P33_1232 ; 
    //    double fNM0_P33_1232 ; 
    // double fNE0_P33_1232 ; ??


    std::map<Resonance_t,double> fA120N ;
    std::map<Resonance_t,double> fA12AlphaN ;
    std::map<Resonance_t,double> fA12BetaN ;

    std::map<Resonance_t,double> fS120N ;
    std::map<Resonance_t,double> fS12AlphaN ;
    std::map<Resonance_t,double> fS12BetaN ;

    std::map<Resonance_t,double> fA320N ;
    std::map<Resonance_t,double> fA32AlphaN ;
    std::map<Resonance_t,double> fA32BetaN ;

  };

}        // genie namespace
#endif
