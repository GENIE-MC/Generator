/**
 * @brief Correction to free NN xsec in nuclear matter
 * 
 * @author Kyle Bachinski, Tomasz Golan
 * @date 2015
 * @remarks V.R. Pandharipande and S. C. Pieper, Phys. Rev. C45 (1992) 791
 * 
*/

#ifndef INUKE_NUCLEON_CORR_H
#define INUKE_NUCLEON_CORR_H

#include <iostream>

#include <TGenPhaseSpace.h>
#include "PDG/PDGCodes.h"

class INukeNucleonCorr
{
  public:
    
    //! get single instance of INukeNucleonCorr; create if necessary
    static INukeNucleonCorr* getInstance() {return fInstance ? fInstance : (fInstance = new INukeNucleonCorr);}
    
    //! get the correction for given four-momentum and density
    double AvgCorrection (const double rho, const int A, const int Z, const int pdg, const double Ek);
    
    double getAvgCorrection (double rho, double A, double ke);

    void OutputFiles(int A, int Z);

          
  private:
  
    static INukeNucleonCorr *fInstance; //!< single instance of INukeNucleonCorr
  
    // ----- MODEL PARAMETERS ----- //
    
    static const unsigned int fRepeat; //!< number of repetition to get average correction

    // ----- POTENTIAL PARAMETERS (from paper) ----- //
    
    static const double fRho0; //!< equilibrium density

    static const double fAlpha1;  //!<  alpha coefficient as defined by Eq. 2.17
    static const double fAlpha2;  //!<  alpha coefficient as defined by Eq. 2.17
    static const double fBeta1;   //!<   beta coefficient as defined by Eq. 2.18
    static const double fLambda0; //!< lambda coefficient as defined by Eq. 2.19
    static const double fLambda1; //!< lambda coefficient as defined by Eq. 2.19
    


    // ----- CALC VARIABLES ----- //
    
    double fFermiMomProton;  // local Fermi momentum for protons
    double fFermiMomNeutron; // local Fermi momentum for neutrons
    
    // ----- SINGLETON "BLOCKADES"----- //
        
    INukeNucleonCorr () {}                                 //!< private constructor (called only by getInstance())
    INukeNucleonCorr (const INukeNucleonCorr&);            //!< block copy constructor
    INukeNucleonCorr& operator= (const INukeNucleonCorr&); //!< block assignment operator

    // ----- CALCULATIONS ----- //

    inline double   beta (const double rho) {return fBeta1 * rho;}                //!< potential component (Eq. 2.18)
    inline double lambda (const double rho) {return (fLambda0 + fLambda1 * rho);} //!< potential component (Eq. 2.19)
    
    inline void setFermiLevel (const double rho, const int A, const int Z) //!< set up Fermi momenta
    {
      fFermiMomProton  = localFermiMom (rho, A, Z, genie::kPdgProton);  // local Fermi momentum for protons
      fFermiMomNeutron = localFermiMom (rho, A, Z, genie::kPdgNeutron); // local Fermi momentum for neutrons
    }
    
    //! return proper Fermi momentum based on nucleon PDG
    inline double fermiMomentum (const int pdg) {return pdg == genie::kPdgProton ? fFermiMomProton : fFermiMomNeutron;}
    
    double mstar (const double rho, const double k2); //!< m* calculated based on Eqs. 2.6 and 2.16
    
    double localFermiMom (const double rho, const int A, const int Z, const int pdg); //!< calculate local Fermi momentum 
    
    TLorentzVector generateTargetNucleon (const double mass, const double fermiMomentum); //!< generate target nucleon
    
    double getCorrection (const double mass, const double rho,
                          const TVector3 &k1, const TVector3 &k2,
                          const TVector3 &k3, const TVector3 &k4); //!< calculate xsec correction
};

#endif // INUKE_NUCLEON_CORR_H
