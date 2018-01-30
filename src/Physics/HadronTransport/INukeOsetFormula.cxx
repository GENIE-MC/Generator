#include "INukeOsetFormula.h"

// constants 
const double INukeOsetFormula :: fCouplingConstant = 0.36 * 4.0 * M_PI;

const double INukeOsetFormula :: fNormalDensity = 0.17; // fm-3

const double INukeOsetFormula :: fNormFactor = 197.327 * 197.327 * 10.0;

// particles mass
                                                            
const double INukeOsetFormula :: fNucleonMass  = kNucleonMass * 1000.0; // [MeV]
const double INukeOsetFormula :: fNucleonMass2 = fNucleonMass * fNucleonMass;

const double INukeOsetFormula :: fDeltaMass = 1232.0; // MeV

// s-wave parametrization (see sec. 3.1)

const double INukeOsetFormula :: fCoefSigma[3] = {-0.01334, 0.06889, 0.19753};
const double INukeOsetFormula :: fCoefB[3] = {-0.01866, 0.06602, 0.21972};
const double INukeOsetFormula :: fCoefD[3] = {-0.08229, 0.37062,-0.03130};
const double INukeOsetFormula :: ImB0 = 0.035;

//! delta parametrization coefficients (NuclPhys A468 (1987) 631-652)
const double INukeOsetFormula ::  fCoefCQ[3] = { -5.19, 15.35,   2.06};
//! delta parametrization coefficients (NuclPhys A468 (1987) 631-652)
const double INukeOsetFormula :: fCoefCA2[3] = {  1.06, -6.64,  22.66};
//! delta parametrization coefficients (NuclPhys A468 (1987) 631-652)
const double INukeOsetFormula :: fCoefCA3[3] = {-13.46, 46.17, -20.34};
//! delta parametrization coefficients (NuclPhys A468 (1987) 631-652)
const double INukeOsetFormula :: fCoefAlpha[3] = {0.382, -1.322, 1.466};
//! delta parametrization coefficients (NuclPhys A468 (1987) 631-652)
const double INukeOsetFormula :: fCoefBeta[3]  = {-0.038,  0.204, 0.613};

using namespace osetUtils;

//! @f$ p_{F} = \left(\frac{3}{2}\pi\rho\right)^{1/3} @f$
void INukeOsetFormula :: setNucleus (const double &density)
{
  static const double constFactor = 3.0 / 2.0 * M_PI * M_PI;

  fNuclearDensity = density;
  fFermiMomentum  = pow (constFactor * fNuclearDensity, 1.0 / 3.0) * 197.327;
  fFermiEnergy    = sqrt (fFermiMomentum*fFermiMomentum + fNucleonMass2);
}

/*! <ul>
 *  <li> set pointer to proper mass (charged vs neutral)
 *  <li> calculate energy, momentum (in LAB and CMS) and invariant mass
 *  </ul>
 */ 
void INukeOsetFormula :: setKinematics (const double &pionTk, const bool &isPi0)
{
  if (isPi0)
  {
    fPionMass  = PDGLibrary::Instance()->Find(kPdgPi0)->Mass() * 1000.0; // [MeV]
    fPionMass2 = fPionMass * fPionMass;
  }
  else
  {
    fPionMass  = PDGLibrary::Instance()->Find(kPdgPiP)->Mass() * 1000.0; // [MeV]
    fPionMass2 = fPionMass * fPionMass;
  }

  fPionKineticEnergy = pionTk;
  fPionEnergy        = fPionKineticEnergy + fPionMass;
  fPionMomentum      = sqrt (fPionEnergy * fPionEnergy - fPionMass2);

  // assuming nucleon at rest
  fInvariantMass = sqrt (fPionMass2 + 2.0 * fNucleonMass * fPionEnergy +
                        fNucleonMass2);

  fMomentumCMS  = fPionMomentum * fNucleonMass / fInvariantMass;
  fMomentumCMS2 = fMomentumCMS * fMomentumCMS;
}

/*! <ul>
 *  <li> calculale Delta half width (in nuclear matter)
 *  <li> calculalte Delta self energy
 *  <li> calculalte Delta propagator
 *  </ul>
 */ 
void INukeOsetFormula :: setDelta ()
{
  static const double constFactor = 1.0 / 12.0 / M_PI;

  fCouplingFactor = fCouplingConstant / fPionMass2; // (f*/m_pi)^2, e.g. eq. 2.6

  // see eq. 2.7 and sec. 2.3
  fReducedHalfWidth = constFactor * fCouplingFactor * fNucleonMass * fMomentumCMS *
                     fMomentumCMS2 / fInvariantMass * deltaReduction ();

  setSelfEnergy (); // calculate qel and absorption part of delta self-energy

  // real and imaginary part of delta denominator (see eq. 2.23)
  const double ReDelta = fInvariantMass - fDeltaMass;
  const double ImDelta = fReducedHalfWidth + fSelfEnergyTotal;

  fDeltaPropagator2 = 1.0 / (ReDelta * ReDelta + ImDelta * ImDelta);
}

/*! <ul>
 *  <li> calculalte p- and s-wave absorption cross section
 *  <li> calculalte total qel (el+cex) p- and s-wave cross section
 *  <li> calculalte fraction of cex in total qel (based on isospin relation)
 *  </ul>
 */
void INukeOsetFormula :: setCrossSections ()
{

  // common part for all p-wave cross sections
  const double pXsecCommon = fNormFactor * fCouplingFactor * fDeltaPropagator2 *
                             fMomentumCMS2 / fPionMomentum;

  // ----- ABSORPTION ----- //

  // constant factor for p-wave absorption cross section
  static const double pAborptionFactor = 4.0 / 9.0;

  // absorption p-wave cross section (see eq. 2.24 and eq. 2.21)
  const double pXsecAbsorption = pAborptionFactor * pXsecCommon *
                                 fSelfEnergyAbsorption;

  // constant factor for s-wave absorption cross section 
  static const double sAbsorptionFactor = 4.0 * M_PI * 197.327 * 10.0 * ImB0;

  // absorption s-wave cross section (see sec. 3.3)
  const double sXsecAbsorption = sAbsorptionFactor / fPionMomentum * fNuclearDensity *
                                 (1.0 + fPionEnergy / 2.0 / fNucleonMass) / pow (fPionMass / 197.327, 4.0);

  // total absorption cross section coming from both s- and p-waves
  fAbsorptionCrossSection = pXsecAbsorption + sXsecAbsorption;

  // ----- TOTAL ----- //  

  // el+cex p-wave cross section (will be multipled by proper factor later)
  const double pXsecTotalQel = pXsecCommon * fReducedHalfWidth;

  // see eq. 3.7 
  const double ksi = (fInvariantMass - fNucleonMass - fPionMass) / fPionMass;

  // el+cex s-wave cross section
  const double sXsecTotalQel = quadraticFunction (ksi, fCoefSigma) /
                               fPionMass2 * fNormFactor;

  // see eq. 3.8
  const double B = quadraticFunction (ksi, fCoefB);
  const double D = quadraticFunction (ksi, fCoefD);
  const double A = 0.5 + 0.5*D;
  const double C = 1.0 - A;

  // see 3.4 (chi = 1 for neutron, chi = -1 for proton, chi = 0 for N=Z)
  const double sTotalQelFactor[fNChannels] = {C - B, A + B, 1.0};

  // see eq. 2.18 (chi = 1 for neutron, chi = -1 for proton, chi = 0 for N=Z)
  static const double pTotalQelFactor[fNChannels] = {2.0 / 9.0, 2.0 / 3.0,
                                                    5.0 / 9.0};

  // total cross section for each channel = qel_s + qel_p + absorption
  for (unsigned int i = 0; i < fNChannels; i++)
    fQelCrossSections[i] = pTotalQelFactor[i] * pXsecTotalQel + sTotalQelFactor[i] * sXsecTotalQel;

  // ----- CEX ----- //

  // see 2.18 (chi = 1 for neutron, chi = -1 for proton, chi = 0 for N=Z)
  static const double pCexFactor[fNChannels] = {4.0 / 27.0, 0.0, 5.0 / 27.0};

  // see eq. 3.4 (chi = 1 for neutron, chi = -1 for proton, chi = 0 for N=Z)
  const double sCexFactor[fNChannels] = {2.0 * C, 0.0, 2.0 * C};

  // total cex cross section coming from both s- and p-waves
  for (unsigned int i = 0; i < fNChannels; i++)
    fCexCrossSections[i] = pCexFactor[i] * pXsecTotalQel + sCexFactor[i] * sXsecTotalQel;
}
                                  
//! related to Pauli blocking, see sec. 2.3
double INukeOsetFormula :: deltaReduction () const
{
  // assuming nucleon at rest 
  const double deltaEnergy = fPionEnergy + fNucleonMass;
  // nucleon energy in CMS
  const double energyCMS   = sqrt(fMomentumCMS * fMomentumCMS + fNucleonMass2);
  // fermiEnergy calculated from density
  const double mu0 = (deltaEnergy * energyCMS - fFermiEnergy * fInvariantMass) /
                     fPionMomentum / fMomentumCMS;

  // see eq. 2.13
  if (mu0 < -1.0) return 0.0;
  if (mu0 >  1.0) return 1.0;
  
  return (mu0*mu0*mu0 + mu0 + 2) / 4.0;
}

//! based on eq. 2.21
void INukeOsetFormula :: setSelfEnergy ()
{
  const double x = fPionKineticEnergy / fPionMass;
  const double densityFraction = fNuclearDensity / fNormalDensity;

  const double alpha = quadraticFunction (x, fCoefAlpha);
  const double  beta = quadraticFunction (x, fCoefBeta);

  double  absNN = quadraticFunction (x, fCoefCA2);
  double absNNN = quadraticFunction (x, fCoefCA3);

  absNN *= pow (densityFraction, beta);

  if (absNNN < 0.0) absNNN = 0.0; // it may happen for Tk < 50 MeV
  else absNNN *= pow (densityFraction, 2.0 * beta);

  // absNN and absNNN are multiplied by number of nucleons to get proper
  // cross section normalization
  //fSelfEnergyAbsorption = 2.0 * absNN + 3.0 * absNNN;
  fSelfEnergyAbsorption = absNN + absNNN;

  // this one goes to the delta propagator
  fSelfEnergyTotal = absNN + absNNN + quadraticFunction (x, fCoefCQ) * pow (densityFraction, alpha);
}

/*! <ul>
 *  <li> set up nucleus
 *  <li> set up kinematics
 *  <li> set up Delta (width, propagator)
 *  <li> calculate cross sections
 *  </ul> 
 */
void INukeOsetFormula :: setupOset (const double &density, const double &pionTk, const int &pionPDG, 
                                    const double &protonFraction)
{
  setNucleus (density);
  setKinematics (pionTk, pionPDG == kPdgPi0);
  setDelta();
  setCrossSections ();
  INukeOset::setCrossSections (pionPDG, protonFraction);  
}

