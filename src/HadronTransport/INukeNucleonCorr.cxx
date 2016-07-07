#include "HadronTransport/INukeNucleonCorr.h"
#include "HadronTransport/INukeUtils2015.h"
#include "HadronTransport/INukeHadroData2015.h"
#include "PDG/PDGLibrary.h"
#include "Conventions/Units.h"
#include "Numerical/RandomGen.h"
#include "Messenger/Messenger.h"
using namespace genie;

INukeNucleonCorr* INukeNucleonCorr::fInstance = NULL; // initialize instance with NULL

// ----- STATIC CONSTANTS ----- //

const unsigned int INukeNucleonCorr::fRepeat = 1000; // how many times kinematics should be generated to get avg corr

const double INukeNucleonCorr::fRho0  = 0.16; // fm^-3

const double INukeNucleonCorr::fBeta1   = -116.00 / fRho0 / 1000.0;          // converted to GeV
const double INukeNucleonCorr::fLambda0 =    3.29 / (units::fermi);          // converted to GeV
const double INukeNucleonCorr::fLambda1 =  -0.373 / (units::fermi) / fRho0;  // converted to GeV

//Three Cache arrays give fine binning at low energy and coarse binning at high energy//
const int INukeNucleonCorr::fNDensityBins = 16;  // number of bins for density cache array
const double INukeNucleonCorr::fDensityStep = fRho0 / fNDensityBins;      // [fm^{-3}]

const int INukeNucleonCorr::fNEnergyBins1  = 20;  // number of bins for energy cache1 array
const double INukeNucleonCorr::fMaxEnergy1    = 0.1; // [GeV]
const double INukeNucleonCorr::fEnergyStep1  = fMaxEnergy1 / fNEnergyBins1; // [GeV]

const int INukeNucleonCorr::fNEnergyBins2  = 8;  // number of bins for energy cache2 array
const double INukeNucleonCorr::fMaxEnergy2    = 0.4; // [GeV]
const double INukeNucleonCorr::fEnergyStep2  = fMaxEnergy2 / fNEnergyBins2; // [GeV]

const int INukeNucleonCorr::fNEnergyBins3  = 1;  // number of bins for energy cache3 array
const double INukeNucleonCorr::fMaxEnergy3    = 1.0; // [GeV]
const double INukeNucleonCorr::fEnergyStep3  = fMaxEnergy3 / fNEnergyBins3; // [GeV]


// ----- CALCULATIONS ----- //

//! \f$m^* (k,\rhp) = m \frac{(\Lambda^2 + k^2)^2}{\Lambda^2 + k^2)^2 - 2\Lambda^2\beta m}\f$
double INukeNucleonCorr::mstar (const double rho, const double k2)
{
  // density [fm^-3], momentum square [GeV^2]
  
  static const double m = (PDGLibrary::Instance()->Find(kPdgProton)->Mass() + 
                           PDGLibrary::Instance()->Find(kPdgNeutron)->Mass()) / 2.0;
  
  const double L = lambda (rho); // potential coefficient lambda
  const double B =   beta (rho); // potential coefficient beta
    
  const double L2 = L * L; // lambda^2 used twice in equation
  
  const double num = (L2 + k2) * (L2 + k2); // numerator
    
  return m * num / (num - 2.0 * L2 * B * m);
}

//! \f$k_F = (\frac{3}{2}\pi^2\rho)^{1/3}\f$
double INukeNucleonCorr :: localFermiMom (const double rho, const int A, const int Z, const int pdg)
{
  static double factor = 3.0 * M_PI * M_PI / 2.0;
    
  return pdg == kPdgProton ? pow (factor * rho * Z / A, 1.0 / 3.0) / (units::fermi) :
                             pow (factor * rho * (A - Z) / A, 1.0 / 3.0) / (units::fermi);
}

//! generate random momentum direction and return 4-momentum of target nucleon
TLorentzVector INukeNucleonCorr :: generateTargetNucleon (const double mass, const double fermiMom)
{
  RandomGen * rnd = RandomGen::Instance();
    
  // get random momentum direction
  const double costheta = 2.0 * rnd->RndGen().Rndm() - 1.0;  // random cos (theta)
  const double sintheta = sqrt (1.0 - costheta * costheta);  // sin (theta)  
  const double      phi = 2.0 * M_PI * rnd->RndGen().Rndm(); // random phi

  // set nucleon 4-momentum
  const double p = rnd->RndGen().Rndm() * fermiMom; // random nucleon momentum up to Fermi level

  const TVector3   p3 = TVector3 (p * sintheta * cos (phi), p * sintheta * sin (phi), p * costheta); // 3-momentum
  const double energy = sqrt (p3.Mag2() + mass * mass); // energy

  return TLorentzVector (p3, energy);
}

//! calculate correction given by eq. 2.9
double INukeNucleonCorr :: getCorrection (const double mass, const double rho,
                                          const TVector3 &k1, const TVector3 &k2,
                                          const TVector3 &k3, const TVector3 &k4)
{
  const double num = (k1 - k2).Mag() * mstar (rho, (k3.Mag2() + k4.Mag2()) / 2.0) / mass / mass;
  const double den = (k1 * (1.0 / mstar (rho, k1.Mag2())) - k2 * (1.0 / mstar (rho, k2.Mag2()))).Mag();
  
  return num / den;
}

//! generate kinematics fRepeat times to calculate average correction
double INukeNucleonCorr :: getAvgCorrection (const double rho, const int A, const int Z, const int pdg, const double Ek)
{
  static double cache1[fNDensityBins][fNEnergyBins1] = {{-1}}; // corrections cache
  static double cache2[fNDensityBins][fNEnergyBins2] = {{-1}}; // corrections cache
  static double cache3[fNDensityBins][fNEnergyBins3] = {{-1}}; // corrections cache
  bool useCache1 = false; 
  bool useCache2 = false;


  const int densityBin = rho / fDensityStep < fNDensityBins ? rho / fDensityStep : fNDensityBins - 1;
  int energyBin;

  if(Ek < fMaxEnergy1) {energyBin = Ek / fEnergyStep1; useCache1 = true;}  // Determine which cache to use(0-100MeV, 100-400MeV, or 400-1000MeV). Then determine which bin to use
  else if(Ek < fMaxEnergy2) {energyBin = Ek / fEnergyStep2; useCache2 = true;} 
  else {energyBin = (fNEnergyBins3 - 1);}
  
  if (Ek < fMaxEnergy1 && cache1[densityBin][energyBin] > 0.0) return cache1[densityBin][energyBin]; // correction already calculated
  if (Ek < fMaxEnergy2 && cache2[densityBin][energyBin] > 0.0) return cache2[densityBin][energyBin]; // correction already calculated
  if (Ek < fMaxEnergy3 && cache3[densityBin][energyBin] > 0.0) return cache3[densityBin][energyBin]; // correction already calculated

  RandomGen * rnd = RandomGen::Instance();

  setFermiLevel (rho, A, Z); // set Fermi momenta for protons and neutrons
  
  const double mass   = PDGLibrary::Instance()->Find(pdg)->Mass(); // mass of incoming nucleon
  const double energy = Ek + mass;

  TLorentzVector p (0.0, 0.0, sqrt (energy * energy - mass * mass), energy); // incoming particle 4-momentum
  GHepParticle incomingParticle (pdg, kIStUndefined, -1,-1,-1,-1, p, TLorentzVector ()); // for IntBounce
  
  double corrPauliBlocking = 0.0; // correction coming from Pauli blocking
  double     corrPotential = 0.0; // correction coming from potential
  
  for (unsigned int i = 0; i < fRepeat; i++) // generate kinematics fRepeat times to get avg corrections
  {
    // get proton vs neutron randomly based on Z/A
    const int targetPdg = rnd->RndGen().Rndm() < (double) Z / A ? kPdgProton : kPdgNeutron;
    
    const double targetMass = PDGLibrary::Instance()->Find(targetPdg)->Mass(); // set nucleon mass
    
    const TLorentzVector target = generateTargetNucleon (targetMass, fermiMomentum (targetPdg)); // generate target nucl
    
    TLorentzVector outNucl1, outNucl2, RemnP4; // final 4-momenta
    
    // random scattering angle
    double C3CM = INukeHadroData2015::Instance()->IntBounce (&incomingParticle, targetPdg, pdg, kIHNFtElas);
    
    // generate kinematics
    utils::intranuke2015::TwoBodyKinematics (mass, targetMass, p, target, outNucl1, outNucl2, C3CM, RemnP4);

    // update Pauli blocking correction
    corrPauliBlocking += (outNucl1.Vect().Mag() > fermiMomentum (pdg) and outNucl2.Vect().Mag() > fermiMomentum (targetPdg));
    
    // update potential-based correction
    corrPotential += getCorrection (mass, rho, p.Vect(), target.Vect(), outNucl1.Vect(), outNucl2.Vect());
  }
  
  corrPauliBlocking /= fRepeat;
      corrPotential /= fRepeat;
      
      if(useCache1 == true){
	cache1[densityBin][energyBin] = corrPauliBlocking * corrPotential; // save this results for future calls
	return cache1[densityBin][energyBin]; } 
      else if(useCache2 == true){
	cache2[densityBin][energyBin] = corrPauliBlocking * corrPotential; // save this results for future calls
	return cache2[densityBin][energyBin]; } 
      else{
	cache3[densityBin][energyBin] = corrPauliBlocking * corrPotential; // save this results for future calls
	return cache3[densityBin][energyBin]; } 

}
