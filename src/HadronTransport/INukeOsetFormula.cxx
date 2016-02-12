#include "INukeOsetFormula.h"

// constants 
const double INukeOsetFormula :: couplingConstant = 0.36 * 4.0 * M_PI;

const double INukeOsetFormula :: normalDensity = 0.17; // fm-3

const double INukeOsetFormula :: normFactor = 197.327 * 197.327 * 10.0;

// particles mass
                                                            
const double INukeOsetFormula :: nucleonMass  = kNucleonMass * 1000.0; // [MeV]
const double INukeOsetFormula :: nucleonMass2 = nucleonMass * nucleonMass;

const double INukeOsetFormula :: deltaMass = 1232.0; // MeV

// s-wave parametrization (see sec. 3.1)

const double INukeOsetFormula :: coefSigma[3] = {-0.01334, 0.06889, 0.19753};
const double INukeOsetFormula :: coefB[3] = {-0.01866, 0.06602, 0.21972};
const double INukeOsetFormula :: coefD[3] = {-0.08229, 0.37062,-0.03130};
const double INukeOsetFormula :: ImB0 = 0.035;

//! delta parametrization coefficients (NuclPhys A468 (1987) 631-652)
const double INukeOsetFormula ::  coefCQ[3] = { -5.19, 15.35,   2.06};
//! delta parametrization coefficients (NuclPhys A468 (1987) 631-652)
const double INukeOsetFormula :: coefCA2[3] = {  1.06, -6.64,  22.66};
//! delta parametrization coefficients (NuclPhys A468 (1987) 631-652)
const double INukeOsetFormula :: coefCA3[3] = {-13.46, 46.17, -20.34};
//! delta parametrization coefficients (NuclPhys A468 (1987) 631-652)
const double INukeOsetFormula :: coefAlpha[3] = {0.382, -1.322, 1.466};
//! delta parametrization coefficients (NuclPhys A468 (1987) 631-652)
const double INukeOsetFormula :: coefBeta[3]  = {-0.038,  0.204, 0.613};

using namespace osetUtils;

//! @f$ p_{F} = \left(\frac{3}{2}\pi\rho\right)^{1/3} @f$
void INukeOsetFormula :: setNucleus (const double &density)
{
  static const double constFactor = 3.0 / 2.0 * M_PI * M_PI;

  nuclearDensity = density;
  fermiMomentum  = pow (constFactor * nuclearDensity, 1.0 / 3.0) * 197.327;
  fermiEnergy    = sqrt (fermiMomentum*fermiMomentum + nucleonMass2);
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
    pionMass  = PDGLibrary::Instance()->Find(kPdgPi0)->Mass() * 1000.0; // [MeV]
    pionMass2 = pionMass * pionMass;
  }
  else
  {
    pionMass  = PDGLibrary::Instance()->Find(kPdgPiP)->Mass() * 1000.0; // [MeV]
    pionMass2 = pionMass * pionMass;
  }

  pionKineticEnergy = pionTk;
  pionEnergy        = pionKineticEnergy + pionMass;
  pionMomentum      = sqrt (pionEnergy * pionEnergy - pionMass2);

  // assuming nucleon at rest
  invariantMass = sqrt (pionMass2 + 2.0 * nucleonMass * pionEnergy +
                        nucleonMass2);

  momentumCMS  = pionMomentum * nucleonMass / invariantMass;
  momentumCMS2 = momentumCMS * momentumCMS;
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

  couplingFactor = couplingConstant / pionMass2; // (f*/m_pi)^2, e.g. eq. 2.6

  // see eq. 2.7 and sec. 2.3
  reducedHalfWidth = constFactor * couplingFactor * nucleonMass * momentumCMS *
                     momentumCMS2 / invariantMass * deltaReduction ();

  setSelfEnergy (); // calculate qel and absorption part of delta self-energy

  // real and imaginary part of delta denominator (see eq. 2.23)
  const double ReDelta = invariantMass - deltaMass;
  const double ImDelta = reducedHalfWidth + selfEnergyTotal;

  deltaPropagator2 = 1.0 / (ReDelta * ReDelta + ImDelta * ImDelta);
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
  const double pXsecCommon = normFactor * couplingFactor * deltaPropagator2 *
                             momentumCMS2 / pionMomentum;

  // ----- ABSORPTION ----- //

  // constant factor for p-wave absorption cross section
  static const double pAborptionFactor = 4.0 / 9.0;

  // absorption p-wave cross section (see eq. 2.24 and eq. 2.21)
  const double pXsecAbsorption = pAborptionFactor * pXsecCommon *
                                 selfEnergyAbsorption;

  // constant factor for s-wave absorption cross section 
  static const double sAbsorptionFactor = 4.0 * M_PI * 197.327 * 10.0 * ImB0;

  // absorption s-wave cross section (see sec. 3.3)
  const double sXsecAbsorption = sAbsorptionFactor / pionMomentum *
                                 nuclearDensity *
                                 (1.0 + pionEnergy / 2.0 / nucleonMass) /
                                 pow (pionMass / 197.327, 4.0);

  // total absorption cross section coming from both s- and p-waves
  absorptionCrossSection = pXsecAbsorption + sXsecAbsorption;

  // ----- TOTAL ----- //  

  // el+cex p-wave cross section (will be multipled by proper factor later)
  const double pXsecTotalQel = pXsecCommon * reducedHalfWidth;

  // see eq. 3.7 
  const double ksi = (invariantMass - nucleonMass - pionMass) / pionMass;

  // el+cex s-wave cross section
  const double sXsecTotalQel = quadraticFunction (ksi, coefSigma) /
                               pionMass2 * normFactor;

  // see eq. 3.8
  const double B = quadraticFunction (ksi, coefB);
  const double D = quadraticFunction (ksi, coefD);
  const double A = 0.5 + 0.5*D;
  const double C = 1.0 - A;

  // see 3.4 (chi = 1 for neutron, chi = -1 for proton, chi = 0 for N=Z)
  const double sTotalQelFactor[nChannels] = {C - B, A + B, 1.0};

  // see eq. 2.18 (chi = 1 for neutron, chi = -1 for proton, chi = 0 for N=Z)
  static const double pTotalQelFactor[nChannels] = {2.0 / 9.0, 2.0 / 3.0,
                                                    5.0 / 9.0};

  // total cross section for each channel = qel_s + qel_p + absorption
  for (unsigned int i = 0; i < nChannels; i++)
    qelCrossSections[i] = pTotalQelFactor[i] * pXsecTotalQel +
                          sTotalQelFactor[i] * sXsecTotalQel;

  // ----- CEX ----- //

  // see 2.18 (chi = 1 for neutron, chi = -1 for proton, chi = 0 for N=Z)
  static const double pCexFactor[nChannels] = {4.0 / 27.0, 0.0, 5.0 / 27.0};

  // see eq. 3.4 (chi = 1 for neutron, chi = -1 for proton, chi = 0 for N=Z)
  const double sCexFactor[nChannels] = {2.0 * C, 0.0, 2.0 * C};

  // total cex cross section coming from both s- and p-waves
  for (unsigned int i = 0; i < nChannels; i++)
    cexCrossSections[i] = pCexFactor[i] * pXsecTotalQel +
                          sCexFactor[i] * sXsecTotalQel;
}
                                  
//! related to Pauli blocking, see sec. 2.3
double INukeOsetFormula :: deltaReduction () const
{
  // assuming nucleon at rest 
  const double deltaEnergy = pionEnergy + nucleonMass;
  // nucleon energy in CMS
  const double energyCMS   = sqrt(momentumCMS * momentumCMS + nucleonMass2);
  // fermiEnergy calculated from density
  const double mu0 = (deltaEnergy * energyCMS - fermiEnergy * invariantMass) /
                     pionMomentum / momentumCMS;

  // see eq. 2.13
  if (mu0 < -1.0) return 0.0;
  if (mu0 >  1.0) return 1.0;
  
  return (mu0*mu0*mu0 + mu0 + 2) / 4.0;
}

//! based on eq. 2.21
void INukeOsetFormula :: setSelfEnergy ()
{
  const double x = pionKineticEnergy / pionMass;
  const double densityFraction = nuclearDensity / normalDensity;

  const double alpha = quadraticFunction (x, coefAlpha);
  const double  beta = quadraticFunction (x, coefBeta);

  double  absNN = quadraticFunction (x, coefCA2);
  double absNNN = quadraticFunction (x, coefCA3);

  absNN *= pow (densityFraction, beta);

  if (absNNN < 0.0) absNNN = 0.0; // it may happen for Tk < 50 MeV
  else absNNN * pow (densityFraction, 2.0 * beta);

  // absNN and absNNN are multiplied by number of nucleons to get proper
  // cross section normalization
  //selfEnergyAbsorption = 2.0 * absNN + 3.0 * absNNN;
  selfEnergyAbsorption = absNN + absNNN;

  // this one goes to the delta propagator
  selfEnergyTotal = absNN + absNNN +
                    quadraticFunction (x, coefCQ) *
                    pow (densityFraction, alpha);
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

