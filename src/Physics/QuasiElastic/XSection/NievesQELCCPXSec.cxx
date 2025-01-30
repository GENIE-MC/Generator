//____________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Joe Johnston, University of Pittsburgh
 Steven Dytman, University of Pittsburgh
*/
//____________________________________________________________________________

#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <Math/IFunction.h>
#include <Math/Integrator.h>

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
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Numerical/GSLUtils.h"
#include "Physics/Common/PrimaryLeptonUtils.h"

#include <iostream> // Used for testing code
#include <fstream> // Used for testing code
#include "Physics/NuclearState/NuclearModelI.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;
using namespace genie::utils;
using namespace std::complex_literals;

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

  // Shift the q0Tilde if required:
  if( fQvalueShifter ) q0Tilde += q0Tilde * fQvalueShifter->Shift(*interaction) ;

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
  double xsec_scale = 1 ;
  const ProcessInfo& proc_info = interaction->ProcInfo();

  if( proc_info.IsWeakCC() ) xsec_scale = fXSecCCScale;
  else if( proc_info.IsWeakNC() ) xsec_scale = fXSecNCScale;

  xsec *= xsec_scale ;

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
  bool good_config = true ;
  double thc;
  GetParam( "CabibboAngle", thc ) ;
  fCos8c2 = TMath::Power(TMath::Cos(thc), 2);
  
  // Do precise calculation of lepton polarization
  GetParamDef( "PreciseLeptonPol", fIsPreciseLeptonPolarization, false ) ;

  // Cross section scaling factor
  GetParam( "QEL-CC-XSecScale", fXSecCCScale ) ;
  GetParam( "QEL-NC-XSecScale", fXSecNCScale ) ;

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

  // Read optional QvalueShifter:
  fQvalueShifter = nullptr;
  if( GetConfig().Exists("QvalueShifterAlg") ) {
    fQvalueShifter = dynamic_cast<const QvalueShifter *> ( this->SubAlg("QvalueShifterAlg") );
    if( !fQvalueShifter ) {
      good_config = false ;
      LOG("NievesQELCCPXSec", pERROR) << "The required QvalueShifterAlg does not exist. AlgID is : " << SubAlg("QvalueShifterAlg")->Id() ;
    }
  }

  if( ! good_config ) {
    LOG("NievesQELCCPXSec", pERROR) << "Configuration has failed.";
    exit(78) ;
  }

  // Scaling factor for the Coulomb potential
  GetParamDef( "CoulombScale", fCoulombScale, 1.0 );
}
//___________________________________________________________________________
void NievesQELCCPXSec::CNCTCLimUcalc(TLorentzVector qTildeP4,
  double M, double r, bool is_neutrino, bool tgtIsNucleus, int tgt_pdgc,
  int A, int Z, int N, double & CN, double & CT, double & CL,
  double & imaginaryU, double & t0, double & r00, bool assumeFreeNucleon) const
{
  if ( tgtIsNucleus && !assumeFreeNucleon ) 
  {
    double dq  = qTildeP4.Vect().Mag();
    double dq2 = TMath::Sq(dq);
    double q2  = qTildeP4.Mag2();
    //Terms for polarization coefficients CN,CT, and CL
    double hbarc2 = TMath::Sq(fhbarc);
    double c0 = 0.380/fhbarc;//Constant for CN in natural units

    //Density gives the nuclear density, normalized to 1
    //Input radius r must be in fm
    double rhop = nuclear::Density(r,A)*Z;
    double rhon = nuclear::Density(r,A)*N;
    double rho = rhop + rhon;
    double rho0 = A*nuclear::Density(0,A);

    double fPrime = (0.33*rho/rho0 + 0.45*(1 - rho/rho0))*c0;

    // Get Fermi momenta
    double kFn, kFp;
    if(fLFG)
    {
        kFn = TMath::Power(3*kPi2*rhon, 1.0/3.0) *fhbarc;
        kFp = TMath::Power(3*kPi2*rhop, 1.0/3.0) *fhbarc;
    }
    else
    {
        kFn = fKFTable->FindClosestKF(tgt_pdgc, kPdgNeutron);
        kFp = fKFTable->FindClosestKF(tgt_pdgc, kPdgProton);
    }

    double kF = TMath::Power(1.5*kPi2*rho, 1./3.)*fhbarc;

    imaginaryU = relLindhardIm(qTildeP4.E(),dq,kFn,kFp,
                               M,is_neutrino,t0,r00);
 

    std::complex<double> relLin(0, 0), udel(0, 0);

    // By comparison with Nieves' fortran code
    if(imaginaryU < 0)
    {
      relLin = relLindhard(qTildeP4.E(), dq, kF, M);
      udel = deltaLindhard(qTildeP4.E(), dq, rho, kF);
    }
    std::complex<double> relLinTot(relLin + udel);
  /* CRho = 2
     DeltaRho = 2500 MeV, (2.5 GeV)^2 = 6.25 GeV^2
     mRho = 770 MeV, (0.770 GeV)^2 = 0.5929 GeV^2
     g' = 0.63 */
    double Vt = 0.08*4*kPi/kPionMass2*
      (2*TMath::Sq( (6.25 - 0.5929)/(6.25 - q2) )*dq2/(q2 - 0.5929) + 0.63);
  /* f^2/4/Pi = 0.08
     DeltaSubPi = 1200 MeV, (1.2 GeV)^2 = 1.44 GeV^2
     g' = 0.63 */
    double Vl = 0.08*4*kPi/kPionMass2*
      (TMath::Sq( (1.44 - kPionMass2)/(1.44 - q2) )*dq2/(q2 - kPionMass2) + 0.63);

    CN = 1/TMath::Sq(abs(1. - fPrime*relLin/hbarc2));
    CT = 1/TMath::Sq(abs(1. - relLinTot*Vt));
    CL = 1/TMath::Sq(abs(1. - relLinTot*Vl));
  }
  else 
  {
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
double NievesQELCCPXSec::relLindhardIm(double q0, double dq,
                                                     double kFn, double kFp,
                                                     double M,
                                                     bool isNeutrino,
                                                     double & t0,
                                                     double & r00) const
{
  double M2 = TMath::Sq(M);
  double EF1,EF2;
  if(isNeutrino){
    EF1 = TMath::Sqrt(M2 + TMath::Sq(kFn)); //EFn
    EF2 = TMath::Sqrt(M2 + TMath::Sq(kFp)); //EFp
  }else{
    EF1 = TMath::Sqrt(M2 + TMath::Sq(kFp)); //EFp
    EF2 = TMath::Sqrt(M2 + TMath::Sq(kFn)); //EFn
  }

  double q2 = TMath::Sq(q0) - TMath::Sq(dq);
  double a = (-q0 + dq*TMath::Sqrt(1 - 4*M2/q2))/2;
  double epsRP = TMath::Max(TMath::Max(M,EF2 - q0),a);
  int factor = (EF1 - EF2 + q0 >= 0)*(EF1 - epsRP >= 0);
  // theta functions q0>0 and -q2>0 are always handled
  t0  = factor*(EF1 + epsRP)/2;
  r00 = factor*(TMath::Sq(EF1) + TMath::Sq(epsRP) + EF1*epsRP)/3;
  return -factor*M2/2/kPi/dq*(EF1 - epsRP);
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
                        double dqgev, double kFgev, double M) const
{
  double q0 = q0gev/fhbarc;
  double qm = dqgev/fhbarc;
  double kf = kFgev/fhbarc;
  double m = M/fhbarc;

  double dummy;
  double relLindIm = relLindhardIm(q0gev, dqgev, kFgev, kFgev, M, true, dummy, dummy);
  //Units of GeV^2
  std::complex<double> relLind(TMath::Sq(fhbarc)*(ruLinRelX(q0,qm,kf,m) + ruLinRelX(-q0,qm,kf,m)), relLindIm);
  return relLind;
}
//____________________________________________________________________________
//Inputs assumed to be in natural units
double NievesQELCCPXSec::ruLinRelX(double q0, double qm,
                                   double kf, double m) const
{
  double q02 = TMath::Sq(q0);
  double qm2 = TMath::Sq(qm);
  double kf2 = TMath::Sq(kf);
  double m2  = TMath::Sq(m);
  double m4  = TMath::Sq(m2);

  double ef = TMath::Sqrt(m2 + kf2);
  double q2 = q02 - qm2;
  double q4 = TMath::Sq(q2);
  double ds = TMath::Sqrt(1 - 4*m2/q2);
  double L1 = TMath::Log((kf + ef)/m);
  
  double aux1 = TMath::Sq(ef + q0);
  double aux2 = (m2 + TMath::Sq(kf - qm))/aux1;
  double aux3 = (m2 + TMath::Sq(kf + qm))/aux1;
  double L2 = TMath::Log(TMath::Abs((1 - aux2)/(1 - aux3)));
  
  double aux5 = TMath::Sq(2*kf + q0*ds)/qm2;
  double aux6 = TMath::Sq(2*kf - q0*ds)/qm2;
  double aux7 = 4*m4*qm2/q4;
  double aux8 = TMath::Sq(kf - ef*ds)/aux7;
  double aux9 = TMath::Sq(kf + ef*ds)/aux7;
  double L3 = TMath::Log(TMath::Abs((aux5 - 1)/(aux6 - 1)*
                                    (aux8 - 1)/(aux9 - 1)));

  return m2*(-L1 + L2*(2*ef + q0)/2/qm - L3*ds/4)/kPi2;
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
std::complex<double> NievesQELCCPXSec::deltaLindhard(double q0, double dq, 
                                                     double rho, double kF) const
{
  double q_zero  = q0/fhbarc;
  double q_mod   = dq/fhbarc;
  double k_fermi = kF/fhbarc;
  //Divide by hbarc in order to use natural units (rho is already in the correct units)

  //m = 939/197.3, md = 1232/197.3, mpi = 139/197.3
  double m   = 4.7592;
  double md  = 6.2433;
  double mpi = 0.7045;

  double fdel_f2 = 0.45;
  double fdel2   = 4*0.36*kPi;
  double wr = md - m;
  double gamma = 0;
  double gammap = 0;

  double q_zero2 =  TMath::Sq(q_zero);
  double q_mod2  =  TMath::Sq(q_mod);
  double k_fermi2 = TMath::Sq(k_fermi);

  double m2 =       TMath::Sq(m);
  double m4 =       TMath::Sq(m2);
  double mpi2 =     TMath::Sq(mpi);
  double mpi4 =     TMath::Sq(mpi2);


  //For the current code q_zero is always real
  //If q_zero can have an imaginary part then only the real part is used
  //until z and zp are calculated
  double aux1 = m2 + q_zero2 - q_mod2;
  double aux2 = 2*q_zero*TMath::Sqrt(m2 + 3*k_fermi2/5);
  double aux3 = TMath::Sq(m + mpi);
  double s    = aux1 + aux2;
  double sp   = aux1 - aux2;
  
  if(s > aux3)
  {
    double srot = TMath::Sqrt(s);
    double qcm = TMath::Sqrt(TMath::Sq(s) + mpi4 + m4 - 2*(s*mpi2 + s*m2 + mpi2*m2))/2/srot;
    gamma = fdel2*TMath::Power(qcm,3)*(m + TMath::Sqrt(m2 + TMath::Sq(qcm)))/mpi2/12/kPi/srot;
  }
  
  if(sp > aux3)
  {
    double srotp = TMath::Sqrt(sp);
    double qcmp  = TMath::Sqrt(TMath::Sq(sp) + mpi4 + m4 - 2*(sp*mpi2 + sp*m2 + mpi2*m2))/2/srotp;
    gammap = fdel2*TMath::Power(qcmp,3)*(m + TMath::Sqrt(m2 + TMath::Sq(qcmp)))/mpi2/12/kPi/srotp;
  }
  
  std::complex<double> z (md/q_mod/k_fermi*( q_zero - q_mod2/2./md - wr + 1i*gamma/2.));
  std::complex<double> zp(md/q_mod/k_fermi*(-q_zero - q_mod2/2./md - wr + 1i*gammap/2.));
  
  std::complex<double> pzeta(0, 0);
  double one_sqrt10 = 1/TMath::Sq(10);
  std::complex<double> z2(z*z);
  double abs_z = abs(z);
  if(abs_z > 50)
  {
    pzeta = 2.*(1. + 1./5./z2)/3./z;
  }
  else if(abs_z < one_sqrt10)
  {
    pzeta = 2.*z*(1. - z2/3.) - 1i*kPi*(1. - z2)/2.;
  }
  else
  {
    pzeta = z + (1. - z2)*log((z + 1.)/(z - 1.))/2.;
  }

  std::complex<double> pzetap(0,0);
  std::complex<double> zp2(zp*zp);
  double abs_zp = abs(zp);
  if(abs_zp > 50)
  {
    pzetap = 2.*(1. + 1./5./zp2)/3./zp;
  }
  else if(abs_zp < one_sqrt10)
  {
    pzetap = 2.*zp*(1. - zp2/3.) - 1i*kPi*(1. - zp2)/2.;
  }
  else
  {
    pzetap = zp + (1. - zp2)*log((zp + 1.)/(zp - 1.))/2.;
  }

  //Multiply by hbarc^2 to give answer in units of GeV^2
  return 2.*rho*md*(pzeta + pzetap)*fdel_f2*TMath::Sq(fhbarc)/q_mod/k_fermi/3.;
}

//____________________________________________________________________________
// Gives coulomb potential in units of GeV
double NievesQELCCPXSec::vcr(const Target * target, double Rcurr) const
{
    double Rmax = MaximalRadius(target);
    if (Rmax == 0) return 0;
    if(Rcurr >= Rmax)
    {
      LOG("Nieves",pNOTICE) << "Radius greater than maximum radius for coulomb corrections."
                          << " Integrating to max radius.";
      Rcurr = Rmax;
    }
    
    int A = target->A();
    int Z = target->Z();
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
    // Multiply by an extra configurable scaling factor that defaults to unity
    return -kAem*4*kPi*result*fhbarc*fCoulombScale;

}
//____________________________________________________________________________
double NievesQELCCPXSec::MaximalRadius(const Target * target) const
{
  if(target->IsNucleus())
  {
    int A = target->A();
    double Rmax = 0.;

    if ( fCoulombRmaxMode == kMatchNieves ) 
    {
      // Rmax calculated using formula from Nieves' fortran code and default
      // charge and neutron matter density parameters from NuclearUtils.cxx
      if (A > 20) 
      {
        double c = TMath::Power(A,0.35), z = 0.54;
        Rmax = c + 9.25*z;
      }
      else 
      {
        // c = 1.75 for A <= 20
        Rmax = TMath::Sqrt(20.0)*1.75;
      }
    }
    else if ( fCoulombRmaxMode == kMatchVertexGeneratorRmax ) 
    {
      // TODO: This solution is fragile. If the formula used by VertexGenerator
      // changes, then this one will need to change too. Switch to using
      // a common function to get Rmax for both.
      Rmax = 3. * fR0 * std::pow(A, 1./3.);
    }
    else 
    {
      LOG("Nieves", pFATAL) << "Unrecognized setting for fCoulombRmaxMode encountered"
        << " in NievesQELCCPXSec::vcr()";
      gAbortingInErr = true;
      std::exit(1);
    }
    
    return Rmax;

  }
  else
  {
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
  double F1V2    = TMath::Power(F1V,   2);
  double xiF2V2  = TMath::Power(xiF2V, 2);
  double q02     = TMath::Power(q[0],  2);
  double dq2     = TMath::Power(dq,    2);
  double q4      = TMath::Power(q2,    2);

  double t0,r00;
  double CN=1.,CT=1.,CL=1.,imU=0;
  CNCTCLimUcalc(qTildeP4, M, r, is_neutrino, tgtIsNucleus,
    tgt_pdgc, A, Z, N, CN, CT, CL, imU,
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
        Lmunu = g[mu][mu]*kPrime[mu]*g[nu][nu]*k[nu]+g[nu][nu]*kPrime[nu]*g[mu][mu]*k[mu] + 1i*imaginaryPart;
        Lnumu = g[nu][nu]*kPrime[nu]*g[mu][mu]*k[mu]+g[mu][mu]*kPrime[mu]*g[nu][nu]*k[nu ]- 1i*imaginaryPart;
      } // End Lmunu calculation
      double aux1 = 2*CL*Fp*(Fp*q2 + 4*FA*M);
      if(mu ==0 && nu == 0){
        Amunu = 16.0*F1V2*(2.0*rulin[0][0]*CN+2.0*q[0]*tulin[0]+q2/2.0)+
          2.0*q2*xiF2V2*
          (4.0-4.0*rulin[0][0]/M2-4.0*q[0]*tulin[0]/M2-q02*(4.0/q2+1.0/M2)) +
          4.0*FA2*(2.0*rulin[0][0]+2.0*q[0]*tulin[0]+(q2/2.0-2.0*M2))-
          aux1*q02-16.0*F1V*xiF2V*(-q2+q02)*CN;
        a00 = real(Amunu); // TESTING CODE
        sum += Lmunu*Amunu;
      }else if(mu == 0 && nu == 3){
        Amunu = 16.0*F1V2*((2.0*rulin[0][3]+tulin[0]*dq)*CN+tulin[3]*q[0])+
          -4.0*q2*xiF2V2*(2.0*rulin[0][3]/M2+(dq*tulin[0]+q[0]*tulin[3])/M2+dq*q[0]*(2.0/q2+0.5/M2))+
          4.0*FA2*((2.0*rulin[0][3]+dq*tulin[0])*CL+q[0]*tulin[3])-dq*q[0]*(aux1+16.0*F1V*xiF2V);
        a0z= real(Amunu); // TESTING CODE
        Anumu = Amunu;
        sum += Lmunu*Anumu + Lnumu*Amunu;
      }else if(mu == 3 && nu == 3){
        Amunu = 16.0*F1V2*(2.0*rulin[3][3]+2.0*dq*tulin[3]-q2/2.0)+
          2.0*q2*xiF2V2*(-4.0-4.0*rulin[3][3]/M2-4.0*dq*tulin[3]/M2-dq2*(4.0/q2+1.0/M2))+
          4.0*FA2*(2.0*rulin[3][3]+2.0*dq*tulin[3]-(q2/2.0-2.0*CL*M2))-
          aux1*dq2-16.0*F1V*xiF2V*q02;
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
        Amunu = sign*16.0*1i*FA*(xiF2V+F1V)*(-dq*tulin[0]*CT + q[0]*tulin[3]);
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
utils::gsl::wrap::NievesQELSmithMonizIntegrand::NievesQELSmithMonizIntegrand(
              const NievesQELCCPXSec* alg_, const Interaction* interaction_, int mu_, int nu_):
ROOT::Math::IBaseFunctionOneDim(),
alg(alg_),
interaction(interaction_), 
mu(mu_),
nu(nu_)
{

}
//____________________________________________________________________________
utils::gsl::wrap::NievesQELSmithMonizIntegrand::~NievesQELSmithMonizIntegrand()
{

}
//____________________________________________________________________________
unsigned int utils::gsl::wrap::NievesQELSmithMonizIntegrand::NDim(void) const
{
  return 1;
}
//____________________________________________________________________________
double utils::gsl::wrap::NievesQELSmithMonizIntegrand::DoEval(double rin) const
{
   return alg->IntegratedAmunuOverMomentum(interaction, rin, mu, nu); 
}
//____________________________________________________________________________
ROOT::Math::IBaseFunctionOneDim *
  utils::gsl::wrap::NievesQELSmithMonizIntegrand::Clone(void) const
{
  return new utils::gsl::wrap::NievesQELSmithMonizIntegrand(alg, interaction, mu, nu);
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
const TVector3 & NievesQELCCPXSec::FinalLeptonPolarization (const Interaction* interaction) const
{
  if (!fIsPreciseLeptonPolarization) return XSecAlgorithmI::FinalLeptonPolarization(interaction);
  
  // Get kinematics and init-state parameters
  const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();
  const Target & target = init_state.Tgt();
  
  const double M = 1; // polarization doesn't depend on M
  // common factor will be reduced
  //double factor  = fCos8c2*M*Omega/pl/El/kPi/4;
  double Rmax = MaximalRadius(&target);
  double W00, Wxx, Wzz, ReW0z, ImWxy;
  if (!target.IsNucleus() || Rmax <= 0)
  {
      FinalLeptonPolarizationOnFreeNucleon (interaction, W00, Wxx, Wzz, ReW0z, ImWxy);
  }
  else
  {
    ROOT::Math::IntegrationOneDim::Type ig_type = utils::gsl::Integration1DimTypeFromString("adaptive");
    double reltol = 1E-4;
    int nmaxeval = 100000;
    
    ROOT::Math::IBaseFunctionOneDim * func = new utils::gsl::wrap::NievesQELSmithMonizIntegrand(this, interaction, 0, 0);
    ROOT::Math::Integrator ig(*func,ig_type,0,reltol,nmaxeval);
    W00 = ig.Integral(0, Rmax);
    delete func;
    
    func = new utils::gsl::wrap::NievesQELSmithMonizIntegrand(this, interaction, 1, 1);
    ig.SetFunction(*func);
    Wxx = ig.Integral(0, Rmax);
    delete func;
    
    func = new utils::gsl::wrap::NievesQELSmithMonizIntegrand(this, interaction, 3, 3);
    ig.SetFunction(*func);
    Wzz = ig.Integral(0, Rmax);
    delete func;
    
    func = new utils::gsl::wrap::NievesQELSmithMonizIntegrand(this, interaction, 0, 3);
    ig.SetFunction(*func);
    ReW0z = ig.Integral(0, Rmax);
    delete func;
    
    func = new utils::gsl::wrap::NievesQELSmithMonizIntegrand(this, interaction, 1, 2);
    ig.SetFunction(*func);
    ImWxy= ig.Integral(0, Rmax);
    delete func;
  }
    
  bool is_neutrino = pdg::IsNeutrino(init_state.ProbePdg());
  TLorentzVector * tempNeutrino = init_state.GetProbeP4(kRfLab);
  TLorentzVector   neutrinoMom  = *tempNeutrino;
  delete tempNeutrino;
  const TLorentzVector leptonMom = kinematics.FSLeptonP4();
  
  TLorentzVector q4 = neutrinoMom - leptonMom;
  
  double q0  = q4.E();
  double q   = q4.Vect().Mag();
  double q02 = q0*q0;
  double q2  = q*q;
  
  double T1 = Wxx/2/M;
  double T2 = (W00 + Wxx + q02/q2*(Wzz - Wxx) - 2*q0/q*ReW0z)/2/M;
  double T3 = -ImWxy/q;
  double T4 = M/2/q2*(Wzz - Wxx);
  double T5 = (ReW0z - q0/q*(Wzz - Wxx))/q;
  
  CalculatePolarizationVectorInTargetRestFrame(
                        fFinalLeptonPolarization, 
                        neutrinoMom, 
                        leptonMom,
                        is_neutrino, 
                        M, T1,T2,T3,T4,T5,0);
  

  
//  std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
//  std::cout << fFinalLeptonPolarization.Mag() << "\n";
//  std::cout << "NV@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n" << std::endl;
  
  return fFinalLeptonPolarization;

}
//___________________________________________________________________________________
void NievesQELCCPXSec::FinalLeptonPolarizationOnFreeNucleon (const Interaction* interaction, double & A00, double & Axx, double & Azz, double & A0z, double & Axy) const
{
  
  A00 = Axx = Azz = A0z = Axy = 0;
  // Get kinematics and init-state parameters
  const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();
  const Target & target = init_state.Tgt();

  // HitNucMass() looks up the PDGLibrary (on-shell) value for the initial
  // struck nucleon
  double Mi = target.HitNucMass();
  
  // On-shell mass of final nucleon (from PDGLibrary)
  double Mf = interaction->RecoilNucleon()->Mass();

  // Isoscalar mass of nucleon
  double M = (Mi + Mf)/2;
  
  // Note that GetProbeP4 defaults to returning the probe 4-momentum in the
  // struck nucleon rest frame, so we have to explicitly ask for the lab frame
  // here
  TLorentzVector * tempNeutrino = init_state.GetProbeP4(kRfLab);
  TLorentzVector neutrinoMom = *tempNeutrino;
  delete tempNeutrino;
  TLorentzVector inNucleonMom(*init_state.TgtPtr()->HitNucP4Ptr());

  TLorentzVector leptonMom = kinematics.FSLeptonP4();
  TLorentzVector outNucleonMom = kinematics.HadSystP4();
  

  double outNucleonEnergy = TMath::Hypot(M, outNucleonMom.P() );
  
  double q0Tilde = outNucleonEnergy - M;

  // Note that we're working in the lab frame (i.e., the rest frame
  // of the target nucleus). We can therefore use Nieves' explicit
  // form of the Amunu tensor if we rotate the 3-momenta so that
  // qTilde is in the +z direction
  TVector3 neutrinoMom3 = neutrinoMom.Vect();
  TVector3 leptonMom3 = leptonMom.Vect();

  // If Coulomb corrections are being used, adjust the lepton 3-momentum used
  // to get q3VecTilde so that its magnitude matches the local
  // Coulomb-corrected value calculated earlier. Note that, although the
  // treatment of Coulomb corrections by Nieves et al. doesn't change the
  // direction of the lepton 3-momentum, it *does* change the direction of the
  // 3-momentum transfer, and so the correction should be applied *before*
  // rotating coordinates into a frame where q3VecTilde lies along the positive
  // z axis.
  TVector3 q3VecTilde = neutrinoMom3 - leptonMom3;

  // Find the rotation angle needed to put q3VecTilde along z
  TVector3 zvec(0, 0, 1);
  TVector3 rot = ( q3VecTilde.Cross(zvec) ).Unit(); // Vector to rotate about
  // Angle between the z direction and q
  double angle = zvec.Angle( q3VecTilde );

  // Handle the edge case where q3VecTilde is along -z, so the
  // cross product above vanishes
  if ( q3VecTilde.Perp() == 0. && q3VecTilde.Z() < 0. ) 
  {
    rot = TVector3(0., 1., 0.);
    angle = kPi;
  }

  // Rotate if the rotation vector is not 0
  if ( rot.Mag() > 0 ) 
  {
    neutrinoMom3.Rotate(angle,rot);
    neutrinoMom.SetVect(neutrinoMom3);
    leptonMom3.Rotate(angle,rot);
    leptonMom.SetVect(leptonMom3);
  }

  // Calculate qTilde
  TLorentzVector qTildeP4(0., 0., q3VecTilde.Mag(), q0Tilde);
  
  double Q2      = interaction->KinePtr()->Q2(true);
  double Q2tilde = -qTildeP4.Mag2();
  
  // Check that Q2tilde > 0 (accounting for rounding errors)
  if (Q2tilde < 0) 
  {
      return;
  }

  // Store Q2tilde in the kinematic variable representing Q2.
  // This will ensure that the form factors are calculated correctly
  // using the de Forest prescription (Q2tilde instead of Q2).
  interaction->KinePtr()->SetQ2(Q2tilde);
  // Calculate form factors
  fFormFactors.Calculate( interaction );
  // Now that the form factors have been calculated, store Q2
  // in the event instead of Q2tilde
  interaction->KinePtr()->SetQ2(Q2);


  // Get the QEL form factors (were calculated before this method was called)
  double F1V   = 0.5*fFormFactors.F1V();
  double xiF2V = 0.5*fFormFactors.xiF2V();
  double FA    = fFormFactors.FA();
  // According to Nieves' paper, Fp = 2.0*M*FA/(kPionMass2-q2), but Llewelyn-
  // Smith uses Fp = 2.0*M^2*FA/(kPionMass2-q2), so I divide by M
  // This gives units of GeV^-1
  double Fp = fFormFactors.Fp()/M;

  // Calculate auxiliary parameters
  // Off shell mass of initial nucleon
  double M2      = M*M;
  double FA2     = FA*FA;
  double F1V2    = F1V*F1V;
  double xiF2V2  = xiF2V*xiF2V;
  double q0      = qTildeP4.E();
  double dq      = qTildeP4.Pz();
  double q02     = q0*q0;
  double dq2     = dq*dq;
  double q2      = q02 - dq2;

  double aux1 = 2*Fp*(Fp*q2 + 4*FA*M);
  A00 = 32*F1V2*(M2 + q0*M + q2/4) + 8*q2*xiF2V2*(-q0/M - q02*(1/q2 + 1/M2/4)) +
                           8*FA2*(q0*M + q2/4) - aux1*q02 - 16*F1V*xiF2V*dq2;
  Axx = 8*(FA2*M2 - q2*( (F1V+xiF2V)*(F1V+xiF2V) + FA2/4 ) );
  Azz =  -8*F1V2*q2 - 8*q2*xiF2V2*(1 + dq2*(1/q2 + 1/M2/4))+
          8*FA2*(M2 - q2/4) - aux1*dq2 - 16*F1V*xiF2V*q02;
  A0z = 16*F1V2*M*dq - 4*q2*xiF2V2*(dq/M + dq*q0*(2/q2 + 1/M2/2)) 
          + 4*FA2*M*dq - dq*q0*(aux1 + 16*F1V*xiF2V);
  Axy = 16*FA*(xiF2V+F1V)*dq*M;
  
  return;
}
//___________________________________________________________________________________
double NievesQELCCPXSec::IntegratedAmunuOverMomentum (const Interaction* interaction, double r, int mu, int nu) const
{
  // Get kinematics and init-state parameters
  const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();
  const Target & target = init_state.Tgt();

  // HitNucMass() looks up the PDGLibrary (on-shell) value for the initial
  // struck nucleon
  double Mi_onshell = target.HitNucMass();
  
  // On-shell mass of final nucleon (from PDGLibrary)
  double Mf = interaction->RecoilNucleon()->Mass();

  // Isoscalar mass of nucleon
  double M = (Mi_onshell + Mf)/2;
  
  // Note that GetProbeP4 defaults to returning the probe 4-momentum in the
  // struck nucleon rest frame, so we have to explicitly ask for the lab frame
  // here
  TLorentzVector * tempNeutrino = init_state.GetProbeP4(kRfLab);
  TLorentzVector neutrinoMom = *tempNeutrino;
  delete tempNeutrino;
  TLorentzVector leptonMom = kinematics.FSLeptonP4();
 
  // Calculate Coulomb corrections
  double ml = interaction->FSPrimLepton()->Mass();
  double ml2 = ml*ml;
  double plLocal = leptonMom.P();

  bool is_neutrino = pdg::IsNeutrino(init_state.ProbePdg());
  int sign = is_neutrino ? 1 : -1;
  
  double ElLocal = leptonMom.E();
  if ( fCoulomb ) 
  {
    // Coulomb potential
    double Vc = vcr(&target, r);

    // Outgoing lepton energy and momentum including Coulomb potential
    double El      = ElLocal;
    ElLocal = El - sign*Vc;

    if ( ElLocal - ml <= 0 ) return 0;

    // The Coulomb correction factor blows up as pl -> 0. To guard against
    // unphysically huge corrections here, require that the lepton kinetic energy
    // (at infinity) is larger than the magnitude of the Coulomb potential
    // (should be around a few MeV)
    double KEl = El - ml;
    if ( KEl <= TMath::Abs(Vc) ) return 0;

    // Local value of the lepton 3-momentum magnitude for the Coulomb correction
    plLocal = TMath::Sqrt(ElLocal*ElLocal - ml2);

  }

  double q0Tilde = neutrinoMom.E() - leptonMom.E();
  
  int nucl_pdg_ini = target.HitNucPdg();
  
  double kFi     = fPauliBlocker->GetFermiMomentum(target, nucl_pdg_ini, r);
  double kFf     = fPauliBlocker->GetFermiMomentum(target, interaction->RecoilNucleonPdg(), r);
  double EFi     = TMath::Hypot(M, kFi);
  double EFf     = TMath::Hypot(M, kFf);
  
  bool tgtIsNucleus = target.IsNucleus();
  int tgt_pdgc = target.Pdg();
  int A = target.A();
  int Z = target.Z();
  int N = target.N();
  bool hitNucIsProton = pdg::IsProton( nucl_pdg_ini );
  
  // This part of the code is strictly in accordance with the original Nieves' paper
  double Mi = target.Mass();
  int Zf = (hitNucIsProton) ? Z - 1 : Z + 1;
  PDGLibrary * pdglib = PDGLibrary::Instance();
  TParticlePDG * nucl_f = pdglib->Find( pdg::IonPdgCode(A, Zf) );
  double Q = 0;
  if(nucl_f) Q = nucl_f -> Mass() - Mi;
  double Q_LFG   = EFf - EFi;
  q0Tilde -= (Q - Q_LFG);
  
  // Shift the q0Tilde if required:
  if( fQvalueShifter ) q0Tilde += q0Tilde * fQvalueShifter->Shift(*interaction) ;

  // If binding energy effects pull us into an unphysical region, return
  // zero
  if ( q0Tilde <= 0 && tgtIsNucleus && !interaction->TestBit(kIAssumeFreeNucleon) ) return 0;

  // Note that we're working in the lab frame (i.e., the rest frame
  // of the target nucleus). We can therefore use Nieves' explicit
  // form of the Amunu tensor if we rotate the 3-momenta so that
  // qTilde is in the +z direction
  TVector3 neutrinoMom3 = neutrinoMom.Vect();
  TVector3 leptonMom3 = leptonMom.Vect();

  // If Coulomb corrections are being used, adjust the lepton 3-momentum used
  // to get q3VecTilde so that its magnitude matches the local
  // Coulomb-corrected value calculated earlier. Note that, although the
  // treatment of Coulomb corrections by Nieves et al. doesn't change the
  // direction of the lepton 3-momentum, it *does* change the direction of the
  // 3-momentum transfer, and so the correction should be applied *before*
  // rotating coordinates into a frame where q3VecTilde lies along the positive
  // z axis.
  TVector3 leptonMomCoulomb3 = !fCoulomb ? leptonMom3: plLocal*leptonMom3*(1/leptonMom3.Mag());
  TVector3 q3VecTilde = neutrinoMom3 - leptonMomCoulomb3;
  // Calculate qTilde
  TLorentzVector qTildeP4(0., 0., q3VecTilde.Mag(), q0Tilde);

  // Find the rotation angle needed to put q3VecTilde along z
  TVector3 zvec(0, 0, 1);
  TVector3 rot = ( q3VecTilde.Cross(zvec) ).Unit(); // Vector to rotate about
  // Angle between the z direction and q
  double angle = zvec.Angle( q3VecTilde );

  // Handle the edge case where q3VecTilde is along -z, so the
  // cross product above vanishes
  if ( q3VecTilde.Perp() == 0. && q3VecTilde.Z() < 0. ) 
  {
    rot = TVector3(0., 1., 0.);
    angle = kPi;
  }

  // Rotate if the rotation vector is not 0
  if ( rot.Mag() > 0 ) 
  {
    neutrinoMom3.Rotate(angle,rot);
    neutrinoMom.SetVect(neutrinoMom3);
    leptonMom3.Rotate(angle,rot);
    leptonMom.SetVect(leptonMom3);
  }

  double Q2      = interaction->KinePtr()->Q2(true);
  double Q2tilde = -qTildeP4.Mag2();
  
  // Check that Q2tilde > 0 (accounting for rounding errors)
  if (Q2tilde < 0) return 0;

  // Store Q2tilde in the kinematic variable representing Q2.
  // This will ensure that the form factors are calculated correctly
  // using the de Forest prescription (Q2tilde instead of Q2).
  interaction->KinePtr()->SetQ2(Q2tilde);
  // Calculate form factors
  fFormFactors.Calculate( interaction );
  // Now that the form factors have been calculated, store Q2
  // in the event instead of Q2tilde
  interaction->KinePtr()->SetQ2(Q2);


  // Get the QEL form factors (were calculated before this method was called)
  double F1V   = 0.5*fFormFactors.F1V();
  double xiF2V = 0.5*fFormFactors.xiF2V();
  double FA    = fFormFactors.FA();
  // According to Nieves' paper, Fp = 2.0*M*FA/(kPionMass2-q2), but Llewelyn-
  // Smith uses Fp = 2.0*M^2*FA/(kPionMass2-q2), so I divide by M
  // This gives units of GeV^-1
  double Fp = fFormFactors.Fp()/M;

  
  double t0,r00;
  double CN(1), CT(1), CL(1), imU(0);
  CNCTCLimUcalc(qTildeP4, M, r, is_neutrino, tgtIsNucleus,
    tgt_pdgc, A, Z, N, CN, CT, CL, imU,
    t0, r00, interaction->TestBit( kIAssumeFreeNucleon ));
    
  if ( imU > 0 ) return 0;
  
  if ( !fRPA || interaction->TestBit( kIAssumeFreeNucleon ) ) 
  {
    CN = 1.0;
    CT = 1.0;
    CL = 1.0;
  }
  
  // Calculate auxiliary parameters
  // Off shell mass of initial nucleon
  double M2      = M*M;
  double FA2     = FA*FA;
  double F1V2    = F1V*F1V;
  double xiF2V2  = xiF2V*xiF2V;
  double q0      = qTildeP4.E();
  double dq      = qTildeP4.Pz();
  double dq2     = dq*dq;
  double q02     = q0*q0;
  double q2      = q02 - dq2;


  double factor  = plLocal*ElLocal*r*r/dq;

  double c       = q0/dq;
  double d       = q2/2/M/dq;
  
  double Elow    = TMath::Max( EFf - q0, M*(c*d + TMath::Sqrt(1- c*c + d*d))/(1 - c*c) );
  double Eup     = EFi;
  
  if (Elow >= Eup) return 0;
  
  double Elow2   = Elow*Elow;
  double Elow3   = Elow*Elow2;
  double Eup2    = Eup*Eup;
  double Eup3    = Eup*Eup2;
  
  double b0      = factor*(Eup - Elow);
  double b1      = factor/M/2*(Eup2 - Elow2);
  double b2      = factor/M2/3*(Eup3 - Elow3);
  
  double a1      = b0;
  double a2      = b2 - b0;
  double a3      = c*c*b2 + 2*c*d*b1 + d*d*b0;
  double a4      = b2;
  double a5      = c*b2 + d*b1;
  double a6      = c*b1 + d*b0;
  double a7      = b1;
  
  
  double Ep   = M*a7;
  double Ep2  = M2*a4;
  double px2  = 0.5*M2*(a2 - a3); //py2=px2
  double pz   = M*a6;
  double Eppz = M2*a5;
  double pz2  = M2*a3;
  
  double aux1 = 2*CL*Fp*(Fp*q2 + 4*FA*M);

  if (mu == 0 && nu == 0)
  {
      return  32*F1V2*(Ep2*CN + Ep*q0 + a1*q2/4)+
              8*q2*xiF2V2*(a1*(1 - q02*(1/q2 + 1/M2/4)) - Ep2/M2 - Ep*q0/M2) +
              8*FA2*(Ep2 + Ep*q0 + a1*(q2/4 - M2)) - a1*(aux1*q02 + 16*F1V*xiF2V*(q02 - q2)*CN);
  }
  if ( (mu == 1 && nu == 1) || (mu == 2 && nu == 2) )
  {
      return  32*F1V2*(px2 - a1*q2/4)-
              8*q2*xiF2V2*(a1*CT + px2/M2) +
              8*FA2*(px2 + a1*(CT*M2 - q2/4))-
              16*a1*F1V*xiF2V*CT*q2;
  }
  if (mu == 3 && nu == 3)
  {
      return  32*F1V2*(pz2 + pz*dq - a1*q2/4)-
              8*q2*xiF2V2*(a1 + pz2/M2 + pz*dq/M2 + a1*dq2*(1/q2 + 1/M2/4))+
              8*FA2*(pz2 + pz*dq + a1*(CL*M2 - q2/4)) - a1*(aux1*dq2 + 16*F1V*xiF2V*q02);
  }
  if (mu == 0 && nu == 3)
  {
      return  16*F1V2*((2*Eppz + Ep*dq)*CN + pz*q0)
              -4*q2*xiF2V2*(2*Eppz/M2 + (Ep*dq + pz*q0)/M2 + a1*dq*q0*(2/q2 + 1/M2/2))+
              4*FA2*((2*Eppz + Ep*dq)*CL + pz*q0) - a1*dq*q0*(aux1 + 16*F1V*xiF2V);
  }
  if (mu == 1 && nu == 2) // should be multiplied by i
  {
      return -16*FA*(xiF2V+F1V)*(pz*q0 - Ep*dq*CT);
  }
  return 0;
}
//___________________________________________________________________________________
