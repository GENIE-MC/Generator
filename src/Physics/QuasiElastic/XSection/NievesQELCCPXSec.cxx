//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Joe Johnston, University of Pittsburgh
         Steven Dytman, University of Pittsburgh

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <Math/IFunction.h>
#include <Math/Integrator.h>
#include <complex>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/QuasiElastic/XSection/QELFormFactors.h"
#include "Physics/QuasiElastic/XSection/QELFormFactorsModelI.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/Units.h"
#include "Physics/Common/VertexGenerator.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/NuclearState/FermiMomentumTablePool.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Physics/QuasiElastic/XSection/NievesQELCCPXSec.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Numerical/GSLUtils.h"

#include <iostream> // Used for testing code
#include <fstream> // Used for testing code
#include "Physics/NuclearState/NuclearModelI.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;
using namespace genie::utils;

//____________________________________________________________________________
NievesQELCCPXSec::NievesQELCCPXSec() :
XSecAlgorithmI("genie::NievesQELCCPXSec")
{

}
//____________________________________________________________________________
NievesQELCCPXSec::NievesQELCCPXSec(string config) :
XSecAlgorithmI("genie::NievesQELCCPXSec", config)
{

}
//____________________________________________________________________________
NievesQELCCPXSec::~NievesQELCCPXSec()
{

 }
//____________________________________________________________________________
double NievesQELCCPXSec::XSec(const Interaction * interaction,
  KinePhaseSpace_t kps) const
{
  /*// TESTING CODE:
  // The first time this method is called, output tensor elements and other
  // kinmeatics variables for various kinematics. This can the be compared
  // to Nieves' fortran code for validation purposes
  if(fCompareNievesTensors){
    LOG("Nieves",pNOTICE) << "Printing tensor elements for specific "
                          << "kinematics for testing purposes";
    CompareNievesTensors(interaction);
    fCompareNievesTensors = false;
    exit(0);
  }
  // END TESTING CODE*/


  if ( !this->ValidProcess   (interaction) ) return 0.;
  if ( !this->ValidKinematics(interaction) ) return 0.;

  // Get kinematics and init-state parameters
  const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();
  const Target & target = init_state.Tgt();

  // HitNucMass() looks up the PDGLibrary (on-shell) value for the initial
  // struck nucleon
  double mNi = target.HitNucMass();

  // Hadronic matrix element for CC neutrino interactions should really use
  // the "nucleon mass," i.e., the mean of the proton and neutrino masses.
  // This expression would also work for NC and EM scattering (since the
  // initial and final on-shell nucleon masses would be the same)
  double mNucleon = ( mNi + interaction->RecoilNucleon()->Mass() ) / 2.;

  // Create a copy of the struck nucleon 4-momentum that is forced
  // to be on-shell (this will be needed later for the tensor contraction,
  // in which the nucleon is treated in this way)
  double inNucleonOnShellEnergy = std::sqrt( std::pow(mNi, 2)
    + std::pow(target.HitNucP4().P(), 2) );

  // The Nieves CCQE model follows the de Forest prescription: free nucleon
  // (i.e., on-shell) form factors and spinors are used, but an effective
  // value of the 4-momentum transfer "qTilde" is used when computing the
  // contraction of the hadronic tensor. See comments in the
  // FullDifferentialXSec() method of LwlynSmithQELCCPXSec for more details.
  TLorentzVector inNucleonMomOnShell( target.HitNucP4().Vect(),
    inNucleonOnShellEnergy );

  // Get the four kinematic vectors and caluclate GFactor
  // Create copies of all kinematics, so they can be rotated
  // and boosted to the nucleon rest frame (because the tensor
  // constraction below only applies for the initial nucleon
  // at rest and q in the z direction)

  // All 4-momenta should already be stored, with the hit nucleon off-shell
  // as appropriate
  TLorentzVector* tempNeutrino = init_state.GetProbeP4(kRfLab);
  TLorentzVector neutrinoMom = *tempNeutrino;
  delete tempNeutrino;
  TLorentzVector inNucleonMom = target.HitNucP4();
  TLorentzVector leptonMom = kinematics.FSLeptonP4();
  TLorentzVector outNucleonMom = kinematics.HadSystP4();

  // Apply Pauli blocking if enabled
  if ( fDoPauliBlocking && target.IsNucleus() && !interaction->TestBit(kIAssumeFreeNucleon) ) {
    int final_nucleon_pdg = interaction->RecoilNucleonPdg();
    double kF = fPauliBlocker->GetFermiMomentum(target, final_nucleon_pdg,
      target.HitNucPosition());
    double pNf = outNucleonMom.P();
    if ( pNf < kF ) return 0.;
  }

  // Use the lab kinematics to calculate the Gfactor, in order to make
  // the XSec differential in initial nucleon momentum and energy
  // Divide by 4.0 because Nieves' conventions for the leptonic and hadronic
  // tensor contraction differ from LwlynSmith by a factor of 4
  double Gfactor = kGF2*fCos8c2 / (8.0*kPi*kPi*inNucleonOnShellEnergy
    *neutrinoMom.E()*outNucleonMom.E()*leptonMom.E()) / 4.0;

  // Calculate Coulomb corrections
  double ml = interaction->FSPrimLepton()->Mass();
  double ml2 = TMath::Power(ml, 2);
  double coulombFactor = 1.0;
  double plLocal = leptonMom.P();

  bool is_neutrino = pdg::IsNeutrino(init_state.ProbePdg());
  double r = target.HitNucPosition();

  if ( fCoulomb ) {
    // Coulomb potential
    double Vc = vcr(& target, r);

    // Outgoing lepton energy and momentum including Coulomb potential
    int sign = is_neutrino ? 1 : -1;
    double El = leptonMom.E();
    double pl = leptonMom.P();
    double ElLocal = El - sign*Vc;

    if ( ElLocal - ml <= 0. ) {
      LOG("Nieves", pDEBUG) << "Event should be rejected. Coulomb effects"
        << " push kinematics below threshold. Returning xsec = 0.0";
      return 0.0;
    }

    // The Coulomb correction factor blows up as pl -> 0. To guard against
    // unphysically huge corrections here, require that the lepton kinetic energy
    // (at infinity) is larger than the magnitude of the Coulomb potential
    // (should be around a few MeV)
    double KEl = El - ml;
    if ( KEl <= std::abs(Vc) ) {
      LOG("Nieves", pDEBUG) << "Outgoing lepton has a very small kinetic energy."
        << " Protecting against near-singularities in the Coulomb correction"
        << " factor by returning xsec = 0.0";
      return 0.0;
    }

    // Local value of the lepton 3-momentum magnitude for the Coulomb
    // correction
    plLocal = TMath::Sqrt( ElLocal*ElLocal - ml2 );

    // Correction factor
    coulombFactor= (plLocal * ElLocal) / (pl * El);

  }

  // When computing the contraction of the leptonic and hadronic tensors,
  // we need to use an effective value of the 4-momentum transfer q.
  // The energy transfer (q0) needs to be modified to account for the binding
  // energy of the struck nucleon, while the 3-momentum transfer needs to
  // be corrected for Coulomb effects.
  //
  // See the original Valencia model paper:
  // https://journals.aps.org/prc/abstract/10.1103/PhysRevC.70.055503

  double q0Tilde = outNucleonMom.E() - inNucleonMomOnShell.E();

  // If binding energy effects pull us into an unphysical region, return
  // zero for the differential cross section
  if ( q0Tilde <= 0. && target.IsNucleus() && !interaction->TestBit(kIAssumeFreeNucleon) ) return 0.;

  // Note that we're working in the lab frame (i.e., the rest frame
  // of the target nucleus). We can therefore use Nieves' explicit
  // form of the Amunu tensor if we rotate the 3-momenta so that
  // qTilde is in the +z direction
  TVector3 neutrinoMom3 = neutrinoMom.Vect();
  TVector3 leptonMom3 = leptonMom.Vect();

  TVector3 inNucleonMom3 = inNucleonMom.Vect();
  TVector3 outNucleonMom3 = outNucleonMom.Vect();

  // If Coulomb corrections are being used, adjust the lepton 3-momentum used
  // to get q3VecTilde so that its magnitude matches the local
  // Coulomb-corrected value calculated earlier. Note that, although the
  // treatment of Coulomb corrections by Nieves et al. doesn't change the
  // direction of the lepton 3-momentum, it *does* change the direction of the
  // 3-momentum transfer, and so the correction should be applied *before*
  // rotating coordinates into a frame where q3VecTilde lies along the positive
  // z axis.
  TVector3 leptonMomCoulomb3 = (! fCoulomb ) ? leptonMom3
    : plLocal * leptonMom3 * (1. / leptonMom3.Mag());
  TVector3 q3VecTilde = neutrinoMom3 - leptonMomCoulomb3;

  // Find the rotation angle needed to put q3VecTilde along z
  TVector3 zvec(0.0, 0.0, 1.0);
  TVector3 rot = ( q3VecTilde.Cross(zvec) ).Unit(); // Vector to rotate about
  // Angle between the z direction and q
  double angle = zvec.Angle( q3VecTilde );

  // Handle the edge case where q3VecTilde is along -z, so the
  // cross product above vanishes
  if ( q3VecTilde.Perp() == 0. && q3VecTilde.Z() < 0. ) {
    rot = TVector3(0., 1., 0.);
    angle = kPi;
  }

  // Rotate if the rotation vector is not 0
  if ( rot.Mag() >= kASmallNum ) {

    neutrinoMom3.Rotate(angle,rot);
    neutrinoMom.SetVect(neutrinoMom3);

    leptonMom3.Rotate(angle,rot);
    leptonMom.SetVect(leptonMom3);

    inNucleonMom3.Rotate(angle,rot);
    inNucleonMom.SetVect(inNucleonMom3);
    inNucleonMomOnShell.SetVect(inNucleonMom3);

    outNucleonMom3.Rotate(angle,rot);
    outNucleonMom.SetVect(outNucleonMom3);

  }

  // Calculate q and qTilde
  TLorentzVector qP4 = neutrinoMom - leptonMom;
  TLorentzVector qTildeP4(0., 0., q3VecTilde.Mag(), q0Tilde);

  double Q2 = -1. * qP4.Mag2();
  double Q2tilde = -1. * qTildeP4.Mag2();

  // Store Q2tilde in the interaction so that we get the correct
  // values of the form factors (according to the de Forest prescription)
  interaction->KinePtr()->SetQ2(Q2tilde);

  double q2 = -Q2tilde;

  // Check that q2 < 0 (accounting for rounding errors)
  if ( q2 >= kASmallNum ) {
    LOG("Nieves", pWARN) << "q2 >= 0, returning xsec = 0.0";
    return 0.0;
  }

  // Calculate form factors
  fFormFactors.Calculate( interaction );

  // Now that the form factors have been calculated, store Q2
  // in the event instead of Q2tilde
  interaction->KinePtr()->SetQ2( Q2 );

  // Do the contraction of the leptonic and hadronic tensors. See the
  // RPA-corrected expressions for the hadronic tensor elements in appendix A
  // of Phys. Rev. C 70, 055503 (2004). Note that the on-shell 4-momentum of
  // the initial struck nucleon should be used in the calculation, as well as
  // the effective 4-momentum transfer q tilde (corrected for the nucleon
  // binding energy and Coulomb effects)
  double LmunuAnumuResult = LmunuAnumu(neutrinoMom, inNucleonMomOnShell,
    leptonMom, qTildeP4, mNucleon, is_neutrino, target,
    interaction->TestBit( kIAssumeFreeNucleon ));

  // Calculate xsec
  double xsec = Gfactor*coulombFactor*LmunuAnumuResult;

  // Apply the factor that arises from elimination of the energy-conserving
  // delta function
  xsec *= genie::utils::EnergyDeltaFunctionSolutionQEL( *interaction );

  // Apply given scaling factor
  xsec *= fXSecScale;

  LOG("Nieves",pDEBUG) << "TESTING: RPA=" << fRPA
                       << ", Coulomb=" << fCoulomb
                       << ", q2 = " << q2 << ", xsec = " << xsec;

  //----- The algorithm computes dxsec/dQ2 or kPSQELEvGen
  //      Check whether variable tranformation is needed
  if ( kps != kPSQELEvGen ) {

    // Compute the appropriate Jacobian for transformation to the requested
    // phase space
    double J = utils::kinematics::Jacobian(interaction, kPSQELEvGen, kps);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("Nieves", pDEBUG)
     << "Jacobian for transformation to: "
                  << KinePhaseSpace::AsString(kps) << ", J = " << J;
#endif
    xsec *= J;
  }

  // Number of scattering centers in the target
  int nucpdgc = target.HitNucPdg();
  int NNucl = (pdg::IsProton(nucpdgc)) ? target.Z() : target.N();

  xsec *= NNucl; // nuclear xsec

  return xsec;
}
//____________________________________________________________________________
double NievesQELCCPXSec::Integral(const Interaction * in) const
{
  // If we're using the new spline generation method (which integrates
  // over the kPSQELEvGen phase space used by QELEventGenerator) then
  // let the cross section integrator do all of the work. It's smart
  // enough to handle free nucleon vs. nuclear targets, different
  // nuclear models (including the local Fermi gas model), etc.
  if ( fXSecIntegrator->Id().Name() == "genie::NewQELXSec" ) {
    return fXSecIntegrator->Integrate(this, in);
  }
  else {
    LOG("Nieves", pFATAL) << "Splines for the Nieves CCQE model must be"
      << " generated using genie::NewQELXSec";
    std::exit(1);
  }
}
//____________________________________________________________________________
bool NievesQELCCPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();

  if(!proc_info.IsQuasiElastic()) return false;

  int  nuc = init_state.Tgt().HitNucPdg();
  int  nu  = init_state.ProbePdg();

  bool isP   = pdg::IsProton(nuc);
  bool isN   = pdg::IsNeutron(nuc);
  bool isnu  = pdg::IsNeutrino(nu);
  bool isnub = pdg::IsAntiNeutrino(nu);

  bool prcok = proc_info.IsWeakCC() && ((isP&&isnub) || (isN&&isnu));
  if(!prcok) return false;

  return true;
}
//____________________________________________________________________________
void NievesQELCCPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NievesQELCCPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NievesQELCCPXSec::LoadConfig(void)
{
  double thc;
  GetParam( "CabibboAngle", thc ) ;
  fCos8c2 = TMath::Power(TMath::Cos(thc), 2);

  // Cross section scaling factor
  GetParam( "QEL-CC-XSecScale", fXSecScale ) ;

  // hbarc for unit conversion, GeV*fm
  fhbarc = kLightSpeed*kPlankConstant/genie::units::fermi;

  // load QEL form factors model
  fFormFactorsModel = dynamic_cast<const QELFormFactorsModelI *> (
    this->SubAlg("FormFactorsAlg") );
  assert(fFormFactorsModel);
  fFormFactors.SetModel( fFormFactorsModel ); // <-- attach algorithm

  // load XSec Integrator
  fXSecIntegrator = dynamic_cast<const XSecIntegratorI*>(
    this->SubAlg("XSec-Integrator") );
  assert(fXSecIntegrator);

  // Load settings for RPA and Coulomb effects

  // RPA corrections will not affect a free nucleon
  GetParamDef("RPA", fRPA, true ) ;

  // Coulomb Correction- adds a correction factor, and alters outgoing lepton
  // 3-momentum magnitude (but not direction)
  // Correction only becomes sizeable near threshold and/or for heavy nuclei
  GetParamDef( "Coulomb", fCoulomb, true ) ;

  LOG("Nieves", pNOTICE) << "RPA=" << fRPA << ", useCoulomb=" << fCoulomb;

  // Get nuclear model for use in Integral()
  RgKey nuclkey = "IntegralNuclearModel";
  fNuclModel = dynamic_cast<const NuclearModelI *> (this->SubAlg(nuclkey));
  assert(fNuclModel);

  // Check if the model is a local Fermi gas
  fLFG = fNuclModel->ModelType(Target()) == kNucmLocalFermiGas;

  if(!fLFG){
    // get the Fermi momentum table for relativistic Fermi gas
    GetParam( "FermiMomentumTable", fKFTableName ) ;

    fKFTable = 0;
    FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
    fKFTable = kftp->GetTable( fKFTableName );
    assert( fKFTable );
  }

  // TESTING CODE
  GetParamDef( "PrintDebugData", fCompareNievesTensors, false ) ;
  // END TESTING CODE

  // Nuclear radius parameter (R = R0*A^(1/3)) to use when computing
  // the maximum radius to use to integrate the Coulomb potential
  GetParam("NUCL-R0", fR0) ; // fm

  std::string temp_mode;
  GetParamDef( "RmaxMode", temp_mode, std::string("VertexGenerator") ) ;

  // Translate the string setting the Rmax mode to the appropriate
  // enum value, or complain if one couldn't be found
  if ( temp_mode == "VertexGenerator" ) {
    fCoulombRmaxMode = kMatchVertexGeneratorRmax;
  }
  else if ( temp_mode == "Nieves" ) {
    fCoulombRmaxMode = kMatchNieves;
  }
  else {
    LOG("Nieves", pFATAL) << "Unrecognized setting \"" << temp_mode
      << "\" requested for the RmaxMode parameter in the"
      << " configuration for NievesQELCCPXSec";
    gAbortingInErr = true;
    std::exit(1);
  }

  // Method to use to calculate the binding energy of the initial hit nucleon when
  // generating splines
  std::string temp_binding_mode;
  GetParamDef( "IntegralNucleonBindingMode", temp_binding_mode, std::string("UseNuclearModel") );
  fIntegralNucleonBindingMode = genie::utils::StringToQELBindingMode( temp_binding_mode );

  // Cutoff energy for integrating over nucleon momentum distribution (above this
  // lab-frame probe energy, the effects of Fermi motion and binding energy
  // are taken to be negligible for computing the total cross section)
  GetParamDef("IntegralNuclearInfluenceCutoffEnergy", fEnergyCutOff, 2.5 ) ;

  // Get PauliBlocker for possible use in XSec()
  fPauliBlocker = dynamic_cast<const PauliBlocker*>( this->SubAlg("PauliBlockerAlg") );
  assert( fPauliBlocker );

  // Decide whether or not it should be used in XSec()
  GetParamDef( "DoPauliBlocking", fDoPauliBlocking, true );
}
//___________________________________________________________________________
void NievesQELCCPXSec::CNCTCLimUcalc(TLorentzVector qTildeP4,
  double M, double r, bool is_neutrino, bool tgtIsNucleus, int tgt_pdgc,
  int A, int Z, int N, bool hitNucIsProton, double & CN, double & CT, double & CL,
  double & imaginaryU, double & t0, double & r00, bool assumeFreeNucleon) const
{
  if ( tgtIsNucleus && !assumeFreeNucleon ) {
    double dq = qTildeP4.Vect().Mag();
    double dq2 = TMath::Power(dq,2);
    double q2 = 1 * qTildeP4.Mag2();
    //Terms for polarization coefficients CN,CT, and CL
    double hbarc2 = TMath::Power(fhbarc,2);
    double c0 = 0.380/fhbarc;//Constant for CN in natural units

    //Density gives the nuclear density, normalized to 1
    //Input radius r must be in fm
    double rhop = nuclear::Density(r,A)*Z;
    double rhon = nuclear::Density(r,A)*N;
    double rho = rhop + rhon;
    double rho0 = A*nuclear::Density(0,A);

    double fPrime = (0.33*rho/rho0+0.45*(1-rho/rho0))*c0;

    // Get Fermi momenta
    double kF1, kF2;
    if(fLFG){
      if(hitNucIsProton){
        kF1 = TMath::Power(3*kPi2*rhop, 1.0/3.0) *fhbarc;
        kF2 = TMath::Power(3*kPi2*rhon, 1.0/3.0) *fhbarc;
      }else{
        kF1 = TMath::Power(3*kPi2*rhon, 1.0/3.0) *fhbarc;
        kF2 = TMath::Power(3*kPi2*rhop, 1.0/3.0) *fhbarc;
      }
    }else{
      if(hitNucIsProton){
        kF1 = fKFTable->FindClosestKF(tgt_pdgc, kPdgProton);
        kF2 = fKFTable->FindClosestKF(tgt_pdgc, kPdgNeutron);
      }else{
        kF1 = fKFTable->FindClosestKF(tgt_pdgc, kPdgNeutron);
        kF2 = fKFTable->FindClosestKF(tgt_pdgc, kPdgProton);
      }
    }

    double kF = TMath::Power(1.5*kPi2*rho, 1.0/3.0) *fhbarc;

    std::complex<double> imU(relLindhardIm(qTildeP4.E(),dq,kF1,kF2,
                                           M,is_neutrino,t0,r00));

    imaginaryU = imag(imU);

    std::complex<double> relLin(0,0),udel(0,0);

    // By comparison with Nieves' fortran code
    if(imaginaryU < 0.){
      relLin = relLindhard(qTildeP4.E(),dq,kF,M,is_neutrino,imU);
      udel = deltaLindhard(qTildeP4.E(),dq,rho,kF);
    }
    std::complex<double> relLinTot(relLin + udel);

  /* CRho = 2
     DeltaRho = 2500 MeV, (2.5 GeV)^2 = 6.25 GeV^2
     mRho = 770 MeV, (0.770 GeV)^2 = 0.5929 GeV^2
     g' = 0.63 */
    double Vt = 0.08*4*kPi/kPionMass2 *
      (2* TMath::Power((6.25-0.5929)/(6.25-q2),2)*dq2/(q2-0.5929) + 0.63);
  /* f^2/4/Pi = 0.08
     DeltaSubPi = 1200 MeV, (1.2 GeV)^2 = 1.44 GeV^2
     g' = 0.63 */
    double Vl = 0.08*4*kPi/kPionMass2 *
      (TMath::Power((1.44-kPionMass2)/(1.44-q2),2)*dq2/(q2-kPionMass2)+0.63);

    CN = 1.0/TMath::Power(abs(1.0-fPrime*relLin/hbarc2),2);

    CT = 1.0/TMath::Power(abs(1.0-relLinTot*Vt),2);
    CL = 1.0/TMath::Power(abs(1.0-relLinTot*Vl),2);
  }
  else {
    //Polarization Coefficients: all equal to 1.0 for free nucleon
    CN = 1.0;
    CT = 1.0;
    CL = 1.0;
    imaginaryU = 0.0;
  }
}
//____________________________________________________________________________
// Gives the imaginary part of the relativistic lindhard function in GeV^2
// and sets the values of t0 and r00
std::complex<double> NievesQELCCPXSec::relLindhardIm(double q0, double dq,
                                                     double kFn, double kFp,
                                                     double M,
                                                     bool isNeutrino,
                                                     double & t0,
                                                     double & r00) const
{
  double M2 = TMath::Power(M,2);
  double EF1,EF2;
  if(isNeutrino){
    EF1 = TMath::Sqrt(M2+TMath::Power(kFn,2)); //EFn
    EF2 = TMath::Sqrt(M2+TMath::Power(kFp,2)); //EFp
  }else{
    EF1 = TMath::Sqrt(M2+TMath::Power(kFp,2)); //EFp
    EF2 = TMath::Sqrt(M2+TMath::Power(kFn,2)); //EFn
  }

  double q2 = TMath::Power(q0,2) - TMath::Power(dq,2);
  double a = (-q0+dq*TMath::Sqrt(1-4.0*M2/q2))/2.0;
  double epsRP = TMath::Max(TMath::Max(M,EF2-q0),a);

  // Other theta functions for q are handled by nuclear suppression
  // That is, q0>0 and -q2>0 are always handled, and q0>EF2-EF1 is
  // handled if pauli blocking is on, because otherwise the final
  // nucleon would be below the fermi sea
  //if(fNievesSuppression && !interaction->TestBit(kIAssumeFreeNucleon )
  //&& !EF1-epsRP<0){
  //LOG("Nieves", pINFO) << "Average value of E(p) above Fermi sea";
  //return 0;
  //}else{
  t0 = 0.5*(EF1+epsRP);
  r00 = (TMath::Power(EF1,2)+TMath::Power(epsRP,2)+EF1*epsRP)/3.0;
  std::complex<double> result(0.0,-M2/2.0/kPi/dq*(EF1-epsRP));
  return result;
  //}
}
//____________________________________________________________________________
//Following obtained from fortran code by J Nieves, which contained the following comment:
/*
 NUCLEON relativistic Lindhard Function
 Same normalization as ULIN
 Real part
 taken from Eur.Phys.J.A25:299-318,2005 (Barbaro et al)
 Eq. 61

 Im. part: Juan.

 INPUT: Real*8
  q0:Energy   [fm]
  qm: modulus 3mom [fm]
  kf: Fermi mom [fm]

 OUTPUT: Complex*16 [fm]

 USES: ruLinRelX, relLindhardIm
 */
//Takes inputs in GeV (with imU in GeV^2), and gives output in GeV^2
std::complex<double> NievesQELCCPXSec::relLindhard(double q0gev,
                        double dqgev, double kFgev, double M,
                        bool isNeutrino,
                        std::complex<double> /* relLindIm */) const
{
  double q0 = q0gev/fhbarc;
  double qm = dqgev/fhbarc;
  double kf = kFgev/fhbarc;
  double m = M/fhbarc;

  if(q0>qm){
    LOG("Nieves", pWARN) << "relLindhard() failed";
    return 0.0;
  }

  std::complex<double> RealLinRel(ruLinRelX(q0,qm,kf,m)+ruLinRelX(-q0,qm,kf,m));
  double t0,r00;
  std::complex<double> ImLinRel(relLindhardIm(q0gev,dqgev,kFgev,kFgev,M,isNeutrino,t0,r00));
  //Units of GeV^2
  return(RealLinRel*TMath::Power(fhbarc,2) + 2.0*ImLinRel);
}
//____________________________________________________________________________
//Inputs assumed to be in natural units
std::complex<double> NievesQELCCPXSec::ruLinRelX(double q0, double qm,
                                                 double kf, double m) const
{
  double q02 = TMath::Power(q0, 2);
  double qm2 = TMath::Power(qm, 2);
  double kf2 = TMath::Power(kf, 2);
  double m2  = TMath::Power(m,  2);
  double m4  = TMath::Power(m,  4);

  double ef = TMath::Sqrt(m2+kf2);
  double q2 = q02-qm2;
  double q4 = TMath::Power(q2,2);
  double ds = TMath::Sqrt(1.0-4.0*m2/q2);
  double L1 = log((kf+ef)/m);
  std::complex<double> uL2(
       TMath::Log(TMath::Abs(
                    (ef + q0 - TMath::Sqrt(m2+TMath::Power(kf-qm,2)))/
                    (ef + q0 - TMath::Sqrt(m2 + TMath::Power(kf + qm,2))))) +
       TMath::Log(TMath::Abs(
                    (ef + q0 + TMath::Sqrt(m2 + TMath::Power(kf - qm,2)))/
                    (ef + q0 + TMath::Sqrt(m2 + TMath::Power(kf + qm,2))))));

  std::complex<double> uL3(
       TMath::Log(TMath::Abs((TMath::Power(2*kf + q0*ds,2)-qm2)/
                             (TMath::Power(2*kf - q0*ds,2)-qm2))) +
       TMath::Log(TMath::Abs((TMath::Power(kf-ef*ds,2) - (4*m4*qm2)/q4)/
                             (TMath::Power(kf+ef*ds,2) - (4*m4*qm2)/q4))));

  std::complex<double> RlinrelX(-L1/(16.0*kPi2)+
                                uL2*(2.0*ef+q0)/(32.0*kPi2*qm)-
                                uL3*ds/(64.0*kPi2));

  return RlinrelX*16.0*m2;
}
//____________________________________________________________________________
//Following obtained from fortran code by J Nieves, which contained the following comment:
/*
   complex Lindhard function for symmetric nuclear matter:
                    from Appendix of
                    E.Oset et al Phys. Rept. 188:79, 1990
                    formula A.4

   input variables:
     q_zero [fm^-1] : Energy
     q_mod  [fm^-1] : Momentum
     rho    [fm^3]  : Nuclear density
     k_fermi[fm^-1] : Fermi momentum

   All variables are real*8

   output variable:
     delta_lind [fm^-2]

            ATTENTION!!!
 Only works properly for real q_zero,
 if q_zero has an imaginary part calculates the L. function
 assuming Gamma= 0.
 Therefore this subroutine provides two different functions
 depending on whether q_zero is real or not!!!!!!!!!!!
*/
std::complex<double> NievesQELCCPXSec::deltaLindhard(double q0,
                                double dq, double rho, double kF) const
{
  double q_zero = q0/fhbarc;
  double q_mod = dq/fhbarc;
  double k_fermi = kF/fhbarc;
  //Divide by hbarc in order to use natural units (rho is already in the correct units)

  //m = 939/197.3, md = 1232/197.3, mpi = 139/197.3
  double m = 4.7592;
  double md = 6.2433;
  double mpi = 0.7045;

  double fdel_f = 2.13;
  double wr = md-m;
  double gamma = 0;
  double gammap = 0;

  double q_zero2 =  TMath::Power(q_zero,  2);
  double q_mod2 =   TMath::Power(q_mod,   2);
  double k_fermi2 = TMath::Power(k_fermi, 2);

  double m2 =       TMath::Power(m,       2);
  double m4 =       TMath::Power(m,       4);
  double mpi2 =     TMath::Power(mpi,     2);
  double mpi4 =     TMath::Power(mpi,     4);

  double fdel_f2 =  TMath::Power(fdel_f,  2);

  //For the current code q_zero is always real
  //If q_zero can have an imaginary part then only the real part is used
  //until z and zp are calculated

  double s = m2+q_zero2-q_mod2+
    2.0*q_zero *TMath::Sqrt(m2+3.0/5.0*k_fermi2);

  if(s>TMath::Power(m+mpi,2)){
    double srot = TMath::Sqrt(s);
    double qcm = TMath::Sqrt(TMath::Power(s,2)+mpi4+m4-2.0*(s*mpi2+s*m2+
                mpi2*m2)) /(2.0*srot);
    gamma = 1.0/3.0 * 1.0/(4.0*kPi) * fdel_f2*
     TMath::Power(qcm,3)/srot*(m+TMath::Sqrt(m2+TMath::Power(qcm,2)))/mpi2;
  }
  double sp = m2+q_zero2-q_mod2-
    2.0*q_zero *TMath::Sqrt(m2+3.0/5.0*k_fermi2);


  if(sp > TMath::Power(m+mpi,2)){
    double srotp = TMath::Sqrt(sp);
    double qcmp = TMath::Sqrt(TMath::Power(sp,2)+mpi4+m4-2.0*(sp*mpi2+sp*m2+
                 mpi2*m2))/(2.0*srotp);
    gammap = 1.0/3.0 * 1.0/(4.0*kPi) * fdel_f2*
      TMath::Power(qcmp,3)/srotp*(m+TMath::Sqrt(m2+TMath::Power(qcmp,2)))/mpi2;
  }
  //}//End if statement
  const std::complex<double> iNum(0,1.0);

  std::complex<double> z(md/(q_mod*k_fermi)*(q_zero-q_mod2/(2.0*md)
                         -wr +iNum*gamma/2.0));
  std::complex<double> zp(md/(q_mod*k_fermi)*(-q_zero-q_mod2/(2.0*md)
                          -wr +iNum*gammap/2.0));

  std::complex<double> pzeta(0.0);
  if(abs(z) > 50.0){
    pzeta = 2.0/(3.0*z)+2.0/(15.0*z*z*z);
  }else if(abs(z) < TMath::Power(10.0,-2)){
    pzeta = 2.0*z-2.0/3.0*z*z*z-iNum*kPi/2.0*(1.0-z*z);
  }else{
    pzeta = z + (1.0-z*z) * log((z+1.0)/(z-1.0))/2.0;
  }

  std::complex<double> pzetap(0);
  if(abs(zp) > 50.0){
    pzetap = 2.0/(3.0*zp)+2.0/(15.0*zp*zp*zp);
  }else if(abs(zp) < TMath::Power(10.0,-2)){
    pzetap = 2.0*zp-2.0/3.0*zp*zp*zp-iNum*kPi/2.0*(1.0-zp*zp);
  }else{
    pzetap = zp+ (1.0-zp*zp) * log((zp+1.0)/(zp-1.0))/2.0;
  }

  //Multiply by hbarc^2 to give answer in units of GeV^2
  return 2.0/3.0 * rho * md/(q_mod*k_fermi) * (pzeta +pzetap) * fdel_f2 *
    TMath::Power(fhbarc,2);
}

//____________________________________________________________________________
// Gives coulomb potential in units of GeV
double NievesQELCCPXSec::vcr(const Target * target, double Rcurr) const{
  if(target->IsNucleus()){
    int A = target->A();
    int Z = target->Z();
    double Rmax = 0.;

    if ( fCoulombRmaxMode == kMatchNieves ) {
      // Rmax calculated using formula from Nieves' fortran code and default
      // charge and neutron matter density parameters from NuclearUtils.cxx
      if (A > 20) {
        double c = TMath::Power(A,0.35), z = 0.54;
        Rmax = c + 9.25*z;
      }
      else {
        // c = 1.75 for A <= 20
        Rmax = TMath::Sqrt(20.0)*1.75;
      }
    }
    else if ( fCoulombRmaxMode == kMatchVertexGeneratorRmax ) {
      // TODO: This solution is fragile. If the formula used by VertexGenerator
      // changes, then this one will need to change too. Switch to using
      // a common function to get Rmax for both.
      Rmax = 3. * fR0 * std::pow(A, 1./3.);
    }
    else {
      LOG("Nieves", pFATAL) << "Unrecognized setting for fCoulombRmaxMode encountered"
        << " in NievesQELCCPXSec::vcr()";
      gAbortingInErr = true;
      std::exit(1);
    }

    if(Rcurr >= Rmax){
      LOG("Nieves",pNOTICE) << "Radius greater than maximum radius for coulomb corrections."
                          << " Integrating to max radius.";
      Rcurr = Rmax;
    }

    ROOT::Math::IBaseFunctionOneDim * func = new
      utils::gsl::wrap::NievesQELvcrIntegrand(Rcurr,A,Z);
    ROOT::Math::IntegrationOneDim::Type ig_type =
      utils::gsl::Integration1DimTypeFromString("adaptive");

    double abstol = 1; // We mostly care about relative tolerance;
    double reltol = 1E-4;
    int nmaxeval = 100000;
    ROOT::Math::Integrator ig(*func,ig_type,abstol,reltol,nmaxeval);
    double result = ig.Integral(0,Rmax);
    delete func;

    // Multiply by Z to normalize densities to number of protons
    // Multiply by hbarc to put result in GeV instead of fm
    return -kAem*4*kPi*result*fhbarc;
  }else{
    // If target is not a nucleus the potential will be 0
    return 0.0;
  }
}
//____________________________________________________________________________
int NievesQELCCPXSec::leviCivita(int input[]) const{
  int copy[4] = {input[0],input[1],input[2],input[3]};
  int permutations = 0;
  int temp;

  for(int i=0;i<4;i++){
    for(int j=i+1;j<4;j++){
      //If any two elements are equal return 0
      if(input[i] == input[j])
        return 0;
      //If a larger number is before a smaller one, use permutations
      //(exchanges of two adjacent elements) to move the smaller element
      //so it is before the larger element, eg 2341->2314->2134->1234
      if(copy[i]>copy[j]){
        temp = copy[j];
        for(int k=j;k>i;k--){
          copy[k] = copy[k-1];
          permutations++;
        }
        copy[i] = temp;
      }
    }
  }

  if(permutations % 2 == 0){
    return 1;
  }else{
    return -1;
  }
}
//____________________________________________________________________________
// Calculates the constraction of the leptonic and hadronic tensors. The
// expressions used here are valid in a frame in which the
// initial nucleus is at rest, and qTilde must be in the z direction.
double NievesQELCCPXSec::LmunuAnumu(const TLorentzVector neutrinoMom,
const TLorentzVector inNucleonMomOnShell, const TLorentzVector leptonMom,
const TLorentzVector qTildeP4, double M, bool is_neutrino,
const Target& target, bool assumeFreeNucleon) const
{
  double r = target.HitNucPosition();
  bool tgtIsNucleus = target.IsNucleus();
  int tgt_pdgc = target.Pdg();
  int A = target.A();
  int Z = target.Z();
  int N = target.N();
  bool hitNucIsProton = pdg::IsProton( target.HitNucPdg() );

  const double k[4] = {neutrinoMom.E(),neutrinoMom.Px(),neutrinoMom.Py(),neutrinoMom.Pz()};
  const double kPrime[4] = {leptonMom.E(),leptonMom.Px(),
                            leptonMom.Py(),leptonMom.Pz()};

  double q2 = qTildeP4.Mag2();

  const double q[4] = {qTildeP4.E(),qTildeP4.Px(),qTildeP4.Py(),qTildeP4.Pz()};
  double q0 = q[0];
  double dq = TMath::Sqrt(TMath::Power(q[1],2)+
                          TMath::Power(q[2],2)+TMath::Power(q[3],2));

  int sign = (is_neutrino) ? 1 : -1;

  // Get the QEL form factors (were calculated before this method was called)
  double F1V   = 0.5*fFormFactors.F1V();
  double xiF2V = 0.5*fFormFactors.xiF2V();
  double FA    = -fFormFactors.FA();
  // According to Nieves' paper, Fp = 2.0*M*FA/(kPionMass2-q2), but Llewelyn-
  // Smith uses Fp = 2.0*M^2*FA/(kPionMass2-q2), so I divide by M
  // This gives units of GeV^-1
  double Fp    = -1.0/M*fFormFactors.Fp();

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Nieves", pDEBUG) << "\n" << fFormFactors;
#endif

  // Calculate auxiliary parameters
  double M2      = TMath::Power(M,     2);
  double FA2     = TMath::Power(FA,    2);
  double Fp2     = TMath::Power(Fp,    2);
  double F1V2    = TMath::Power(F1V,   2);
  double xiF2V2  = TMath::Power(xiF2V, 2);
  double q02     = TMath::Power(q[0],  2);
  double dq2     = TMath::Power(dq,    2);
  double q4      = TMath::Power(q2,    2);

  double t0,r00;
  double CN=1.,CT=1.,CL=1.,imU=0;
  CNCTCLimUcalc(qTildeP4, M, r, is_neutrino, tgtIsNucleus,
    tgt_pdgc, A, Z, N, hitNucIsProton, CN, CT, CL, imU,
    t0, r00, assumeFreeNucleon);

  if ( imU > kASmallNum )
    return 0.;

  if ( !fRPA || assumeFreeNucleon ) {
    CN = 1.0;
    CT = 1.0;
    CL = 1.0;
  }

  double tulin[4] = {0.,0.,0.,0.};
  double rulin[4][4] = { {0.,0.,0.,0.},
                         {0.,0.,0.,0.},
                         {0.,0.,0.,0.},
                         {0.,0.,0.,0.} };

  // TESTING CODE:
  if(fCompareNievesTensors){
    // Use average values for initial momentum to calculate A, as given
    // in Appendix B of Nieves' paper. T gives average values of components
    // of p, and R gives the average value of two components multiplied
    // together
    double t3 = (0.5*q2 + q0*t0)/dq; // Average pz

    // Vector of p

    tulin[0] = t0;
    tulin[3] = t3;

    // R is a 4x4 matrix, with R[mu][nu] is the average
    // value of p[mu]*p[nu]
    double aR = r00-M2;
    double bR = (q4+4.0*r00*q02+4.0*q2*q0*t0)/(4.0*dq2);
    double gamma = (aR-bR)/2.0;
    double delta = (-aR+3.0*bR)/2.0/dq2;

    double r03 = (0.5*q2*t0 + q0*r00)/dq; // Average E(p)*pz

    rulin[0][0] = r00;
    rulin[0][3] = r03;
    rulin[1][1] = gamma;
    rulin[2][2] = gamma;
    rulin[3][0] = r03;
    rulin[3][3] = gamma+delta*dq2; // END TESTING CODE
  }
  else {
    // For normal code execulation, tulin is the initial nucleon momentum
    tulin[0] = inNucleonMomOnShell.E();
    tulin[1] = inNucleonMomOnShell.Px();
    tulin[2] = inNucleonMomOnShell.Py();
    tulin[3] = inNucleonMomOnShell.Pz();

    for(int i=0; i<4; i++)
      for(int j=0; j<4; j++)
        rulin[i][j] = tulin[i]*tulin[j];
  }

  //Additional constants and variables
  const int g[4][4] = {{1,0,0,0},{0,-1,0,0},{0,0,-1,0},{0,0,0,-1}};
  const std::complex<double> iNum(0,1);
  int leviCivitaIndexArray[4];
  double imaginaryPart = 0;

  std::complex<double> sum(0.0,0.0);

  double kPrimek = k[0]*kPrime[0]-k[1]*kPrime[1]-k[2]*kPrime[2]-k[3]*kPrime[3];

  std::complex<double> Lmunu(0.0,0.0),Lnumu(0.0,0.0);
  std::complex<double> Amunu(0.0,0.0),Anumu(0.0,0.0);

  // Calculate LmunuAnumu by iterating over mu and nu
  // In each iteration, add LmunuAnumu to sum if mu=nu, and add
  // LmunuAnumu + LnumuAmunu if mu != nu, since we start nu at mu
  double axx=0.,azz=0.,a0z=0.,a00=0.,axy=0.;
  for(int mu=0;mu<4;mu++){
    for(int nu=mu;nu<4;nu++){
      imaginaryPart = 0;
      if(mu == nu){
        //if mu==nu then levi-civita = 0, so imaginary part = 0
        Lmunu = g[mu][mu]*kPrime[mu]*g[nu][nu]*k[nu]+g[nu][nu]*kPrime[nu]*g[mu][mu]*k[mu]-g[mu][nu]*kPrimek;
      }else{
        //if mu!=nu, then g[mu][nu] = 0
        //This same leviCivitaIndex array can be used in the else portion when
        //calculating Anumu
        leviCivitaIndexArray[0] = mu;
        leviCivitaIndexArray[1] = nu;
        for(int a=0;a<4;a++){
          for(int b=0;b<4;b++){
            leviCivitaIndexArray[2] = a;
            leviCivitaIndexArray[3] = b;
            imaginaryPart += - leviCivita(leviCivitaIndexArray)*kPrime[a]*k[b];
          }
        }
        //real(Lmunu) is symmetric, and imag(Lmunu) is antisymmetric
        //std::complex<double> num(g[mu][mu]*kPrime[mu]*g[nu][nu]*k[nu]+g[nu][nu]*kPrime[nu]*g[mu][mu]*k[mu],imaginaryPart);
        Lmunu = g[mu][mu]*kPrime[mu]*g[nu][nu]*k[nu]+g[nu][nu]*kPrime[nu]*g[mu][mu]*k[mu] + iNum*imaginaryPart;
        Lnumu = g[nu][nu]*kPrime[nu]*g[mu][mu]*k[mu]+g[mu][mu]*kPrime[mu]*g[nu][nu]*k[nu ]- iNum*imaginaryPart;
      } // End Lmunu calculation

      if(mu ==0 && nu == 0){
        Amunu = 16.0*F1V2*(2.0*rulin[0][0]*CN+2.0*q[0]*tulin[0]+q2/2.0)+
          2.0*q2*xiF2V2*
          (4.0-4.0*rulin[0][0]/M2-4.0*q[0]*tulin[0]/M2-q02*(4.0/q2+1.0/M2)) +
          4.0*FA2*(2.0*rulin[0][0]+2.0*q[0]*tulin[0]+(q2/2.0-2.0*M2))-
          (2.0*CL*Fp2*q2+8.0*FA*Fp*CL*M)*q02-16.0*F1V*xiF2V*(-q2+q02)*CN;
        a00 = real(Amunu); // TESTING CODE
        sum += Lmunu*Amunu;
      }else if(mu == 0 && nu == 3){
        Amunu = 16.0*F1V2*((2.0*rulin[0][3]+tulin[0]*dq)*CN+tulin[3]*q[0])+
          2.0*q2*xiF2V2*
          (-4.0*rulin[0][3]/M2-2.0*(dq*tulin[0]+q[0]*tulin[3])/M2-dq*q[0]*(4.0/q2+1.0/M2))+
          4.0*FA2*((2.0*rulin[0][3]+dq*tulin[0])*CL+q[0]*tulin[3])-
          (2.0*CL*Fp2*q2+8.0*FA*Fp*CL*M)*dq*q[0]-
          16.0*F1V*xiF2V*dq*q[0];
        a0z= real(Amunu); // TESTING CODE
        Anumu = Amunu;
        sum += Lmunu*Anumu + Lnumu*Amunu;
      }else if(mu == 3 && nu == 3){
        Amunu = 16.0*F1V2*(2.0*rulin[3][3]+2.0*dq*tulin[3]-q2/2.0)+
          2.0*q2*xiF2V2*(-4.0-4.0*rulin[3][3]/M2-4.0*dq*tulin[3]/M2-dq2*(4.0/q2+1.0/M2))+
          4.0*FA2*(2.0*rulin[3][3]+2.0*dq*tulin[3]-(q2/2.0-2.0*CL*M2))-
          (2.0*CL*Fp2*q2+8.0*FA*Fp*CL*M)*dq2-
          16.0*F1V*xiF2V*(q2+dq2);
        azz = real(Amunu); // TESTING CODE
        sum += Lmunu*Amunu;
      }else if(mu ==1 && nu == 1){
        Amunu = 16.0*F1V2*(2.0*rulin[1][1]-q2/2.0)+
          2.0*q2*xiF2V2*(-4.0*CT-4.0*rulin[1][1]/M2) +
          4.0*FA2*(2.0*rulin[1][1]-(q2/2.0-2.0*CT*M2))-
          16.0*F1V*xiF2V*CT*q2;
        axx = real(Amunu); // TESTING CODE
        sum += Lmunu*Amunu;
      }else if(mu == 2 && nu == 2){
        // Ayy not explicitly listed in paper. This is included so rotating the
        // coordinates of k and k' about the z-axis does not change the xsec.
        Amunu = 16.0*F1V2*(2.0*rulin[2][2]-q2/2.0)+
          2.0*q2*xiF2V2*(-4.0*CT-4.0*rulin[2][2]/M2) +
          4.0*FA2*(2.0*rulin[2][2]-(q2/2.0-2.0*CT*M2))-
          16.0*F1V*xiF2V*CT*q2;
        sum += Lmunu*Amunu;
      }else if(mu ==1 && nu == 2){
        Amunu = sign*16.0*iNum*FA*(xiF2V+F1V)*(-dq*tulin[0]*CT + q[0]*tulin[3]);
        Anumu = -Amunu; // Im(A) is antisymmetric
        axy = imag(Amunu); // TESTING CODE
        sum += Lmunu*Anumu+Lnumu*Amunu;
      }
      // All other terms will be 0 because the initial nucleus is at rest and
      // qTilde is in the z direction

    } // End loop over nu
  } // End loop over mu

  // TESTING CODE
  if(fCompareNievesTensors){
    // get tmu
    double tmugev = leptonMom.E() - leptonMom.Mag();
    // Print Q2, form factors, and tensor elts
    std::ofstream ffstream;
    ffstream.open(fTensorsOutFile, std::ios_base::app);
    if(q0 > 0){
      ffstream << -q2 << "\t" << q[0] << "\t" << dq
               << "\t" << axx << "\t" << azz << "\t" << a0z
               << "\t" << a00 << "\t" << axy << "\t"
               << CT << "\t" << CL << "\t" << CN << "\t"
               << tmugev << "\t" << imU << "\t"
               << F1V << "\t" << xiF2V << "\t"
               << FA << "\t" << Fp << "\t"
               << tulin[0] << "\t"<< tulin[1] << "\t"
               << tulin[2] << "\t"<< tulin[3] << "\t"
               << rulin[0][0]<< "\t"<< rulin[0][1]<< "\t"
               << rulin[0][2]<< "\t"<< rulin[0][3]<< "\t"
               << rulin[1][0]<< "\t"<< rulin[1][1]<< "\t"
               << rulin[1][2]<< "\t"<< rulin[1][3]<< "\t"
               << rulin[2][0]<< "\t"<< rulin[2][1]<< "\t"
               << rulin[2][2]<< "\t"<< rulin[2][3]<< "\t"
               << rulin[3][0]<< "\t"<< rulin[3][1]<< "\t"
               << rulin[3][2]<< "\t"<< rulin[3][3]<< "\t"
               << fVc << "\t" << fCoulombFactor << "\t";
        ffstream << "\n";
    }
    ffstream.close();
  }
  // END TESTING CODE

  // Since the real parts of A and L are both symmetric and the imaginary
  // parts are antisymmetric, the contraction should be real
  if ( imag(sum) > kASmallNum )
    LOG("Nieves",pWARN) << "Imaginary part of tensor contraction is nonzero "
                        << "in QEL XSec, real(sum) = " << real(sum)
                        << "imag(sum) = " << imag(sum);

  return real(sum);
}

//___________________________________________________________________________
// Auxiliary scalar function for internal integration
//____________________________________________________________________________
utils::gsl::wrap::NievesQELvcrIntegrand::NievesQELvcrIntegrand(
                                               double Rcurr, int A, int Z):
ROOT::Math::IBaseFunctionOneDim()
{
  fRcurr = Rcurr;
  fA = A;
  fZ = Z;
}
//____________________________________________________________________________
utils::gsl::wrap::NievesQELvcrIntegrand::~NievesQELvcrIntegrand()
{

}
//____________________________________________________________________________
unsigned int utils::gsl::wrap::NievesQELvcrIntegrand::NDim(void) const
{
  return 1;
}
//____________________________________________________________________________
double utils::gsl::wrap::NievesQELvcrIntegrand::DoEval(double rin) const
{
  double rhop = fZ*nuclear::Density(rin,fA);
  if(rin<fRcurr){
    return rhop*rin*rin/fRcurr;
  }else{
    return rhop*rin;
  }
}
//____________________________________________________________________________
ROOT::Math::IBaseFunctionOneDim *
  utils::gsl::wrap::NievesQELvcrIntegrand::Clone(void) const
{
  return new utils::gsl::wrap::NievesQELvcrIntegrand(fRcurr,fA,fZ);
}
//____________________________________________________________________________

//____________________________________________________________________________
//
// NOTE: THE REMAINING IS TESTING CODE
//
// This method prints the tensor elements (as well as various inputs) for
// different kinematics. The tensor elements can then be compared to the
// output of Nieves' fortran code.
//
// The results of this code will only agree exactlly with Nieves' fortran
// if Dipole form factors are set (in UserPhysicsOptions).
//
void NievesQELCCPXSec::CompareNievesTensors(const Interaction* in)
  const {
  Interaction * interaction = new Interaction(*in); // copy in

  // Set input values here
  double ein = 0.2;
  double ctl = 0.5;
  double rmaxfrac = 0.25;

  bool carbon = false; // true -> C12, false -> Pb208

  if(fRPA)
    fTensorsOutFile = "gen.RPA";
  else
    fTensorsOutFile = "gen.noRPA";

  // Calculate radius
  bool klave;
  double rp,ap,rn,an;
  if(carbon){
    klave = true;
    rp = 1.692;
    ap = 1.082;
    rn = 1.692;
    an = 1.082;
  }else{
    // Pb208
    klave = false;
    rp = 6.624;
    ap = 0.549;
    rn = 6.890;
    an = 0.549;
  }
  double rmax;
  if(!klave)
    rmax = TMath::Max(rp,rn) + 9.25*TMath::Max(ap,an);
  else
    rmax = TMath::Sqrt(20.0)*TMath::Max(rp,rn);
  double r = rmax *  rmaxfrac;

  // Relevant objects and parameters
  //const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();
  const Target & target = init_state.Tgt();

  // Parameters required for LmunuAnumu
  double M  = target.HitNucMass();
  double ml = interaction->FSPrimLepton()->Mass();
  bool is_neutrino = pdg::IsNeutrino(init_state.ProbePdg());

  // Iterate over lepton energy (which then affects q, which is passed to
  // LmunuAnumu using in and out NucleonMom
  double delta = (ein-0.025)/100.0;
  for(int it=0;it<100;it++){
    double tmu = it*delta;
    double eout = ml + tmu;
    double pout = TMath::Sqrt(eout*eout-ml*ml);

    double pin = TMath::Sqrt(ein*ein); // Assume massless neutrinos

    double q0 = ein-eout;
    double dq = TMath::Sqrt(pin*pin+pout*pout-2.0*ctl*pin*pout);
    double q2 = q0*q0-dq*dq;
    interaction->KinePtr()->SetQ2(-q2);

    // When this method is called, inNucleonMomOnShell is unused.
    // I can thus provide the calculated values using a null vector for
    // inNucleonMomOnShell. I also need to put qTildeP4 in the z direction, as
    // Nieves does in his paper.
    TLorentzVector qTildeP4(0, 0, dq, q0);
    TLorentzVector inNucleonMomOnShell(0,0,0,0);

    // neutrinoMom and leptonMom only directly affect the leptonic tensor, which
    // we are not calculating now. Use them to transfer q.
    TLorentzVector neutrinoMom(0,0,pout+dq,eout+q0);
    TLorentzVector leptonMom(0,0,pout,eout);

    if(fCoulomb){ // Use same steps as in XSec()
      // Coulomb potential
      double Vc = vcr(& target, r);
      fVc = Vc;

      // Outgoing lepton energy and momentum including coulomb potential
      int sign = is_neutrino ? 1 : -1;
      double El = leptonMom.E();
      double ElLocal = El - sign*Vc;
      if(ElLocal - ml <= 0.0){
        LOG("Nieves",pINFO) << "Event should be rejected. Coulomb effects "
                            << "push kinematics below threshold";
        return;
      }
      double plLocal = TMath::Sqrt(ElLocal*ElLocal-ml*ml);

      // Correction factor
      double coulombFactor= plLocal*ElLocal/leptonMom.Vect().Mag()/El;
      fCoulombFactor = coulombFactor; // Store and print
    }

    // TODO: apply Coulomb correction to 3-momentum transfer dq

    fFormFactors.Calculate(interaction);
    LmunuAnumu(neutrinoMom, inNucleonMomOnShell, leptonMom, qTildeP4,
               M, is_neutrino, target, false);
  }
  return;
} // END TESTING CODE
//____________________________________________________________________________
