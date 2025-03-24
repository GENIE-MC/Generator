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
#include "Physics/QuasiElastic/XSection/NievesQELCCXSec.h"
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
    int sign = is_neutrino ? -1 : 1;
    double El = leptonMom.E();
    double pl = leptonMom.P();
    double ElLocal = El + sign*Vc;

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

  LOG("Nieves",pDEBUG) << "RPA=" << fRPA
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
  if ( fXSecIntegrator->Id().Name() == "genie::NievesQELCCXSec" ) 
  {
    Target * tgt = in->InitStatePtr()->TgtPtr();
    tgt->SetHitNucPosition(MaximalRadius(tgt) );
    return fXSecIntegrator->Integrate(this, in);
  }
  else 
  {
    LOG("Nieves", pFATAL) << "Splines for the Nieves CCQE model must be"
      << " generated using genie::NievesQELCCPXSec";
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
  
  // The calculation is only appropriate for complex nuclear targets,
  // not free nucleons.
  if ( !init_state.Tgt().IsNucleus() ) return false;

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

  // Nuclear radius parameter (R = R0*A^(1/3)) to use when computing
  // the maximum radius to use to integrate the Coulomb potential
  GetParam("NUCL-R0", fR0) ; // fm

  std::string temp_mode;
  GetParamDef( "RmaxMode", temp_mode, std::string("VertexGenerator") ) ;
  
  GetParamDef( "LindhardFunction", fLindhardFunction, std::string("OriginalByNieves") ) ;

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
  double M, double r, bool tgtIsNucleus, int A, int Z, int N, 
  double & CN, double & CT, double & CL, bool assumeFreeNucleon) const
{
  if ( tgtIsNucleus && !assumeFreeNucleon ) 
  {
    double q0  = qTildeP4.E();
    double dq  = qTildeP4.Vect().Mag();
    double dq2 = TMath::Sq(dq);
    double q2  = qTildeP4.Mag2();
    //Terms for polarization coefficients CN,CT, and CL
    double hbarc2 = TMath::Sq(fhbarc);
    double c0 = 0.380/fhbarc/hbarc2; //Constant for CN in natural units

    //Density gives the nuclear density, normalized to 1
    //Input radius r must be in fm
    double rhop = nuclear::Density(r,A)*Z;
    double rhon = nuclear::Density(r,A)*N;
    double rho = rhop + rhon;
    double rho0 = A*nuclear::Density(0,A);

    double fPrime = (0.33*rho/rho0 + 0.45*(1 - rho/rho0));

    double kF = TMath::Power(1.5*kPi2*rho, 1./3.)*fhbarc;

    std::complex<double> Unuc( LindhardNuclear(q0, dq, kF, M) );
    std::complex<double> Udel( LindhardDelta(q0, dq, kF, M, rho*hbarc2*fhbarc) );
    std::complex<double> Utot = Unuc + Udel;
 
// CRho = 2, DeltaRho = 2500 MeV, (2.5 GeV)^2 = 6.25 GeV^2, mRho = 770 MeV, (0.770 GeV)^2 = 0.5929 GeV^2, g' = 0.63 
    double aux = 0.08*4*kPi/kPionMass2;
    double Vt = aux*(2*TMath::Sq( (6.25 - 0.5929)/(6.25 - q2) )*dq2/(q2 - 0.5929) + 0.63);
// f^2/4/Pi = 0.08, DeltaSubPi = 1200 MeV, (1.2 GeV)^2 = 1.44 GeV^2, g' = 0.63 
    double Vl = aux*(TMath::Sq( (1.44 - kPionMass2)/(1.44 - q2) )*dq2/(q2 - kPionMass2) + 0.63);

    CN = 1/std::norm(1. - c0*fPrime*Unuc);
    CT = 1/std::norm(1. - Vt*Utot);
    CL = 1/std::norm(1. - Vl*Utot);
  }
  else 
  {
    //Polarization Coefficients: all equal to 1.0 for free nucleon
    CN = 1.0;
    CT = 1.0;
    CL = 1.0;
  }
}
//____________________________________________________________________________
// Gives the imaginary part of the relativistic lindhard function in GeV^2, Ref.1, Eq.B2
double NievesQELCCPXSec::relLindhardIm(double q0, double dq,
                                       double kFn, double kFp,
                                       double M, bool isNeutrino) const
{
  double M2 = TMath::Sq(M);
  double EF1,EF2;
  if(isNeutrino)
  {
    EF1 = TMath::Sqrt(M2 + TMath::Sq(kFn)); //EFn
    EF2 = TMath::Sqrt(M2 + TMath::Sq(kFp)); //EFp
  }
  else
  {
    EF1 = TMath::Sqrt(M2 + TMath::Sq(kFp)); //EFp
    EF2 = TMath::Sqrt(M2 + TMath::Sq(kFn)); //EFn
  }

  double q2 = TMath::Sq(q0) - TMath::Sq(dq);
  double a = (-q0 + dq*TMath::Sqrt(1 - 4*M2/q2))/2;
  double epsRP = TMath::Max(TMath::Max(a,EF2 - q0),M);
  // theta functions q0>0 and -q2>0 are always handled
  if ( (EF2 - q0 >= EF1) || (a >= EF1) ) return 0;

  return -M2/2/kPi/dq*(EF1 - epsRP);
}
//____________________________________________________________________________
//Inputs assumed to be in natural units, Ref.2, Eq.61
double NievesQELCCPXSec::ruLinRelX(double q0, double dq,
                                   double kF, double M) const
{
  double q02 = q0*q0;
  double dq2 = dq*dq;
  double kF2 = kF*kF;
  double M2  = M*M;
  double M4  = M2*M2;

  double EF = TMath::Sqrt(M2 + kF2);
  double q2 = q02 - dq2;
  double q4 = q2*q2;
  double ds = TMath::Sqrt(1 - 4*M2/q2);
  double L1 = TMath::Log((kF + EF)/M);
  
  double aux1 = TMath::Sq(EF + q0);
  double aux2 = (M2 + TMath::Sq(kF - dq))/aux1;
  double aux3 = (M2 + TMath::Sq(kF + dq))/aux1;
  double L2 = TMath::Log(TMath::Abs((1 - aux2)/(1 - aux3)));
  
  double aux5 = TMath::Sq(2*kF + q0*ds)/dq2;
  double aux6 = TMath::Sq(2*kF - q0*ds)/dq2;
  double aux7 = 4*M4*dq2/q4;
  double aux8 = TMath::Sq(kF - EF*ds)/aux7;
  double aux9 = TMath::Sq(kF + EF*ds)/aux7;
  double L3 = TMath::Log(TMath::Abs((aux5 - 1)/(aux6 - 1)*
                                    (aux8 - 1)/(aux9 - 1)));

  return M2*(-L1 + L2*(2*EF + q0)/2/dq - L3*ds/4)/kPi2;
}
//____________________________________________________________________________  
//Takes inputs in GeV and gives output in GeV^2
std::complex<double> NievesQELCCPXSec::LindhardNuclear(double q0, double dq, double kF, double M) const
{
    if (fLindhardFunction == "CJP46")
    {
        // Ref.4, Eqs.27-29 
        double v = q0*M/kF/kF;
        double q = dq/kF;
        double q2 = q*q;
        double auxm = v/q - q/2;
        double auxp = v/q + q/2;
        double auxm2 = auxm*auxm;
        double auxp2 = auxp*auxp;
        double ReUnuc =  M*kF/kPi2*(-1 + ( (1 - auxm2)*TMath::Log( TMath::Abs( (auxm + 1)/(auxm - 1) ) ) -
                                        (1 - auxp2)*TMath::Log( TMath::Abs( (auxp + 1)/(auxp - 1) ) )  )/2/q);
        
        double uplim  = q + q2/2;
        double lowlim = q - q2/2;                                   
        double ImUnuc = 0;
        if ((q > 2 && uplim >= v && v >= -lowlim ) || (q < 2 && uplim >= v && v >= lowlim) ) 
            ImUnuc = -M*kF*(1 - auxm2 )/2/kPi/q;
        else if (q < 2 && 0 <= v && v <= lowlim)
            ImUnuc = -M*kF*v/kPi/q; 
        
        return std::complex(ReUnuc, ImUnuc);
    }
    
    //Following obtained from fortran code by J Nieves, which contained the following comment:
    // NUCLEON relativistic Lindhard Function
    // Same normalization as ULIN
    // Real part
    // taken from Eur.Phys.J.A25:299-318,2005 (Barbaro et al)
    // Eq. 61
    // Im. part: Juan.
    double relLindIm = relLindhardIm(q0, dq, kF, kF, M, true);
    std::complex<double> relLind(ruLinRelX(q0,dq,kF,M) + ruLinRelX(-q0,dq,kF,M), 2*relLindIm);
    return relLind;
}
//____________________________________________________________________________
std::complex<double> NievesQELCCPXSec::LindhardDelta(double q0, double dq, double kF, double M, double rho) const
{
    double MD = 1.232;
    double mpi = kPionMass;
    
    double fs2_f2 = 4.5;
    double wR = MD - M;
    double gamma  = 0;
    double gammap = 0;
    
    double q02  = q0*q0;
    double dq2  = dq*dq;
    double kF2  = kF*kF;
    
    double M2 =       M*M;
    double M4 =       M2*M2;
    double mpi2 =     mpi*mpi;
    double mpi4 =     mpi2*mpi2;
    //For the current code q0 is always real
    //If q0 can have an imaginary part then only the real part is used
    double aux1 = M2 + q02 - dq2;
    double aux2 = 2*q0*TMath::Sqrt(M2 + 3*kF2/5);
    double aux3 = TMath::Sq(M + mpi);
    double s    = aux1 + aux2;
    double sp   = aux1 - aux2;
    
    if (fLindhardFunction == "CJP46")
    {
        // Ref.4, Eq.30-31
        const double fs2 = 4*kPi*0.36;
        if(s > aux3)
        {
            double srot = TMath::Sqrt(s);
            double qcm  = TMath::Sqrt(TMath::Sq(s) + mpi4 + M4 - 2*(s*mpi2 + s*M2 + mpi2*M2))/2/srot;
            double qcm3 = qcm*qcm*qcm;
            gamma       = fs2*M*qcm3/mpi2/12/kPi/srot;
        }
        
        if(sp > aux3)
        {
            double srotp = TMath::Sqrt(sp);
            double qcmp  = TMath::Sqrt(TMath::Sq(sp) + mpi4 + M4 - 2*(sp*mpi2 + sp*M2 + mpi2*M2))/2/srotp;
            double qcmp3 = qcmp*qcmp*qcmp;
            gammap       = fs2*M*qcmp3/mpi2/12/kPi/srotp;
        }
        
        double b = dq/MD;
        double b3 = b*b*b;  
        std::complex<double> a ( q0 - dq2/2/MD - wR, gamma );
        std::complex<double> ap(-q0 - dq2/2/MD - wR, gammap);
        a  /= kF;
        ap /= kF;
        
        std::complex<double> a_plus_b   = a + b;
        std::complex<double> a_minus_b  = a - b;
        std::complex<double> ap_plus_b  = ap + b;
        std::complex<double> ap_minus_b = ap - b;
        
        std::complex<double> L, Lp;
        
        if (gamma <= 0)
            L = log(abs(a_plus_b/a_minus_b));
        else
            L = log(a_plus_b/a_minus_b);
            
        if (gammap <= 0)
            Lp = log(abs(ap_plus_b/ap_minus_b));
        else
            Lp = log(ap_plus_b/ap_minus_b);
            
        double factor = 4*kF2*fs2_f2/9/kPi2;
        return factor*(b*(a + ap) - (a_plus_b*a_minus_b*L + ap_plus_b*ap_minus_b*Lp)/2.)/b3;
    }
    
    //   Following obtained from fortran code by J Nieves, which contained the following comment:
    //   complex Lindhard function for symmetric nuclear matter:
    //                    from Appendix of
    //                    E.Oset et al Phys. Rept. 188:79, 1990
    //                    formula A.4
    //
    //            ATTENTION!!!
    // Only works properly for real q0,
    // if q0 has an imaginary part calculates the L. function
    // assuming Gamma= 0.
    // Therefore this subroutine provides two different functions
    // depending on whether q0 is real or not!!!!!!!!!!!
    // For the current code q0 is always real
    // If q0 can have an imaginary part then only the real part is used
    // until z and zp are calculated

    // Ref.3, Eq. A4
    if(s > aux3)
    {
        double srot = TMath::Sqrt(s);
        double qcm = TMath::Sqrt(s*s + mpi4 + M4 - 2*(s*mpi2 + s*M2 + mpi2*M2))/2/srot;
        double qcm2 = qcm*qcm;
        gamma = fs2_f2*qcm*qcm2*(M + TMath::Sqrt(M2 + qcm2))/mpi2/12/kPi/srot;
    }
    
    if(sp > aux3)
    {
        double srotp = TMath::Sqrt(sp);
        double qcmp  = TMath::Sqrt(sp*sp + mpi4 + M4 - 2*(sp*mpi2 + sp*M2 + mpi2*M2))/2/srotp;
        double qcmp2 = qcmp*qcmp;
        gammap = fs2_f2*qcmp*qcmp2*(M + TMath::Sqrt(M2 + qcmp2))/mpi2/12/kPi/srotp;
    }
    
    std::complex<double> z ( q0 - dq2/2./MD - wR, gamma/2. );
    std::complex<double> zp(-q0 - dq2/2./MD - wR, gammap/2.);
    z  *= MD/dq/kF;
    zp *= MD/dq/kF;
    
    std::complex<double> pzeta(0, 0);
    std::complex<double> z2(z*z);
    double abs_z = abs(z);
    if(abs_z > 50)
    {
        pzeta = 2.*(1. + 1./5./z2)/3./z;
    }
    else if(abs_z < 1e-2)
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
    else if(abs_zp < 1e-2)
    {
        pzetap = 2.*zp*(1. - zp2/3.) - 1i*kPi*(1. - zp2)/2.;
    }
    else
    {
        pzetap = zp + (1. - zp2)*log((zp + 1.)/(zp - 1.))/2.;
    }
    
    return 2.*rho*MD*(pzeta + pzetap)*fs2_f2/dq/kF/3.;

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
    
    ROOT::Math::IBaseFunctionOneDim * func = new utils::gsl::wrap::NievesQELvcrIntegrand(Rcurr,A,Z);
    const NievesQELCCXSec * integrator = dynamic_cast<const NievesQELCCXSec*>(fXSecIntegrator);
    ROOT::Math::IntegrationOneDim::Type ig_type = utils::gsl::Integration1DimTypeFromString( integrator->Get1DimIntgType() );
    double reltol = integrator->Get1DimRelTol();
    unsigned int nmaxeval = integrator->Get1DimMaxEval();
    ROOT::Math::Integrator ig(*func, ig_type, 0, reltol, nmaxeval);
    double result = ig.Integral(0, Rmax);
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
      Rmax = 3*fR0*std::pow(A, 1./3);
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
  
  int A = target.A();
  int Z = target.Z();
  int N = target.N();
  
  const double k[4] = {neutrinoMom.E(),neutrinoMom.Px(),neutrinoMom.Py(),neutrinoMom.Pz()};
  const double kPrime[4] = {leptonMom.E(),leptonMom.Px(),
                            leptonMom.Py(),leptonMom.Pz()};

  double q2 = qTildeP4.Mag2();

  const double q[4] = {qTildeP4.E(),qTildeP4.Px(),qTildeP4.Py(),qTildeP4.Pz()};

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
  double M2      = M*M;
  double FA2     = FA*FA;
  double F1V2    = F1V*F1V;
  double xiF2V2  = xiF2V*xiF2V;
  double q02     = q[0]*q[0];
  double dq2     = dq*dq;

  double CN(1),CT(1),CL(1);
  if (fRPA)
  {
    CNCTCLimUcalc(qTildeP4, M, r, tgtIsNucleus,
    A, Z, N, CN, CT, CL, assumeFreeNucleon);
  }

  double tulin[4] = {0.,0.,0.,0.};
  double rulin[4][4] = { {0.,0.,0.,0.},
                         {0.,0.,0.,0.},
                         {0.,0.,0.,0.},
                         {0.,0.,0.,0.} };


  // Tulin is the initial nucleon momentum
  tulin[0] = inNucleonMomOnShell.E();
  tulin[1] = inNucleonMomOnShell.Px();
  tulin[2] = inNucleonMomOnShell.Py();
  tulin[3] = inNucleonMomOnShell.Pz();
  
  for(int i=0; i<4; i++)
      for(int j=0; j<4; j++)
        rulin[i][j] = tulin[i]*tulin[j];


  //Additional constants and variables
  double imaginaryPart = 0;

  std::complex<double> sum(0.0,0.0);

  double kPrimek = k[0]*kPrime[0]-k[1]*kPrime[1]-k[2]*kPrime[2]-k[3]*kPrime[3];

  std::complex<double> Lmunu(0.0,0.0),Lnumu(0.0,0.0);
  std::complex<double> Amunu(0.0,0.0),Anumu(0.0,0.0);

  // Calculate LmunuAnumu by iterating over mu and nu
  // In each iteration, add LmunuAnumu to sum if mu=nu, and add
  // LmunuAnumu + LnumuAmunu if mu != nu, since we start nu at mu
  for(int mu=0;mu<4;mu++){
    for(int nu=mu;nu<4;nu++){
      imaginaryPart = 0;
      if(mu == nu){
        //if mu==nu then levi-civita = 0, so imaginary part = 0
        Lmunu = g(mu, mu)*kPrime[mu]*g(nu, nu)*k[nu]+g(nu, nu)*kPrime[nu]*g(mu, mu)*k[mu]-g(mu, nu)*kPrimek;
      }else{
        //if mu!=nu, then g[mu][nu] = 0
        //This same leviCivitaIndex array can be used in the else portion when
        //calculating Anumu
        for(int a=0;a<4;a++){
          for(int b=0;b<4;b++){
            imaginaryPart -= e(mu, nu,a, b)*kPrime[a]*k[b];
          }
        }
        //real(Lmunu) is symmetric, and imag(Lmunu) is antisymmetric
        Lmunu = g(mu, mu)*kPrime[mu]*g(nu, nu)*k[nu]+g(nu, nu)*kPrime[nu]*g(mu, mu)*k[mu] + 1i*imaginaryPart;
        Lnumu = g(nu, nu)*kPrime[nu]*g(mu, mu)*k[mu]+g(mu, mu)*kPrime[mu]*g(nu, nu)*k[nu] - 1i*imaginaryPart;
      } // End Lmunu calculation
      double aux1 = 2*CL*Fp*(Fp*q2 + 4*FA*M);
      if(mu ==0 && nu == 0){
        Amunu = 16.0*F1V2*(2.0*rulin[0][0]*CN+2.0*q[0]*tulin[0]+q2/2.0)+
          2.0*q2*xiF2V2*
          (4.0-4.0*rulin[0][0]/M2-4.0*q[0]*tulin[0]/M2-q02*(4.0/q2+1.0/M2)) +
          4.0*FA2*(2.0*rulin[0][0]+2.0*q[0]*tulin[0]+(q2/2.0-2.0*M2))-
          aux1*q02-16.0*F1V*xiF2V*(-q2+q02)*CN;
        sum += Lmunu*Amunu;
      }else if(mu == 0 && nu == 3){
        Amunu = 16.0*F1V2*((2.0*rulin[0][3]+tulin[0]*dq)*CN+tulin[3]*q[0])+
          -4.0*q2*xiF2V2*(2.0*rulin[0][3]/M2+(dq*tulin[0]+q[0]*tulin[3])/M2+dq*q[0]*(2.0/q2+0.5/M2))+
          4.0*FA2*((2.0*rulin[0][3]+dq*tulin[0])*CL+q[0]*tulin[3])-dq*q[0]*(aux1+16.0*F1V*xiF2V);
        Anumu = Amunu;
        sum += Lmunu*Anumu + Lnumu*Amunu;
      }else if(mu == 3 && nu == 3){
        Amunu = 16.0*F1V2*(2.0*rulin[3][3]+2.0*dq*tulin[3]-q2/2.0)+
          2.0*q2*xiF2V2*(-4.0-4.0*rulin[3][3]/M2-4.0*dq*tulin[3]/M2-dq2*(4.0/q2+1.0/M2))+
          4.0*FA2*(2.0*rulin[3][3]+2.0*dq*tulin[3]-(q2/2.0-2.0*CL*M2))-
          aux1*dq2-16.0*F1V*xiF2V*q02;
        sum += Lmunu*Amunu;
      }else if(mu ==1 && nu == 1){
        Amunu = 16.0*F1V2*(2.0*rulin[1][1]-q2/2.0)+
          2.0*q2*xiF2V2*(-4.0*CT-4.0*rulin[1][1]/M2) +
          4.0*FA2*(2.0*rulin[1][1]-(q2/2.0-2.0*CT*M2))-
          16.0*F1V*xiF2V*CT*q2;
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
        sum += Lmunu*Anumu+Lnumu*Amunu;
      }
      // All other terms will be 0 because the initial nucleus is at rest and
      // qTilde is in the z direction

    } // End loop over nu
  } // End loop over mu

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
              const NievesQELCCPXSec* alg_, const Interaction* interaction_, int mod_):
ROOT::Math::IBaseFunctionOneDim(),
alg(alg_),
interaction(interaction_), 
mod(mod_)
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
   return alg->IntegratedOverMomentum(interaction, rin, mod); 
}
//____________________________________________________________________________
ROOT::Math::IBaseFunctionOneDim *
  utils::gsl::wrap::NievesQELSmithMonizIntegrand::Clone(void) const
{
  return new utils::gsl::wrap::NievesQELSmithMonizIntegrand(alg, interaction, mod);
}
//____________________________________________________________________________
const TVector3 & NievesQELCCPXSec::FinalLeptonPolarization (const Interaction* interaction) const
{
  if (!fIsPreciseLeptonPolarization) return XSecAlgorithmI::FinalLeptonPolarization(interaction);
  
  // Get kinematics and init-state parameters
  const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();
  const Target & target = init_state.Tgt();
  
  double Rmax = MaximalRadius(&target);
  if (Rmax <= 0) return XSecAlgorithmI::FinalLeptonPolarization(interaction);

  const NievesQELCCXSec * integrator = dynamic_cast<const NievesQELCCXSec*>(fXSecIntegrator);
  ROOT::Math::IntegrationOneDim::Type ig_type = utils::gsl::Integration1DimTypeFromString( integrator->Get1DimIntgType() );
  double reltol = integrator->Get1DimRelTol();
  unsigned int nmaxeval = integrator->Get1DimMaxEval();
  
  ROOT::Math::IBaseFunctionOneDim * func = new utils::gsl::wrap::NievesQELSmithMonizIntegrand(this, interaction, 1);
  ROOT::Math::Integrator ig(*func, ig_type, 0, reltol, nmaxeval);
  double R = ig.Integral(0, Rmax);
  delete func;
  if (R == 0)
  {
      fFinalLeptonPolarization = TVector3(0, 0, 0);
      return fFinalLeptonPolarization;
  }
    
  func = new utils::gsl::wrap::NievesQELSmithMonizIntegrand(this, interaction, 2);
  ig.SetFunction(*func);
  double PLR = ig.Integral(0, Rmax);
  delete func;
  
  func = new utils::gsl::wrap::NievesQELSmithMonizIntegrand(this, interaction, 3);
  ig.SetFunction(*func);
  double PPR = ig.Integral(0, Rmax);
  delete func;
  
  double PL = PLR/R;
  double PP = PPR/R;
    
  TLorentzVector * tempNeutrino = init_state.GetProbeP4(kRfLab);
  TLorentzVector   neutrinoMom  = *tempNeutrino;
  delete tempNeutrino;
  const TLorentzVector leptonMom = kinematics.FSLeptonP4();
  
  TVector3 neutrinoMom3 = neutrinoMom.Vect();                                          
  TVector3 leptonMom3   = leptonMom.Vect();
  TVector3 Pz = leptonMom3.Unit();
  TVector3 Px = neutrinoMom3.Cross(leptonMom3).Unit();
  TVector3 Py = Pz.Cross(Px);
  fFinalLeptonPolarization = PP*Py + PL*Pz;
  
  if (fFinalLeptonPolarization.Mag2()>1) return XSecAlgorithmI::FinalLeptonPolarization(interaction);
    
  return fFinalLeptonPolarization;

}
//___________________________________________________________________________________
double NievesQELCCPXSec::IntegratedOverMomentum (const Interaction* interaction, double r, int mod) const
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
  double PlLocal = leptonMom.P();

  bool is_neutrino = pdg::IsNeutrino(init_state.ProbePdg());
  int sign = is_neutrino ? -1 : 1;
  
  double ElLocal = leptonMom.E();
  if ( fCoulomb ) 
  {
    // Coulomb potential
    double Vc = vcr(&target, r);

    // Outgoing lepton energy and momentum including Coulomb potential
    double El      = ElLocal;
    ElLocal = El + sign*Vc;

    if ( ElLocal - ml <= 0 ) return 0;

    // The Coulomb correction factor blows up as pl -> 0. To guard against
    // unphysically huge corrections here, require that the lepton kinetic energy
    // (at infinity) is larger than the magnitude of the Coulomb potential
    // (should be around a few MeV)
    double KEl = El - ml;
    if ( KEl <= TMath::Abs(Vc) ) return 0;

    // Local value of the lepton 3-momentum magnitude for the Coulomb correction
    PlLocal = TMath::Sqrt(ElLocal*ElLocal - ml2);

  }

  double q0Tilde = neutrinoMom.E() - leptonMom.E();
  
  int nucl_pdg_ini = target.HitNucPdg();
  
  double kFi     = fPauliBlocker->GetFermiMomentum(target, nucl_pdg_ini, r);
  double kFf     = fPauliBlocker->GetFermiMomentum(target, interaction->RecoilNucleonPdg(), r);
  double EFi     = TMath::Hypot(M, kFi);
  double EFf     = TMath::Hypot(M, kFf);
  
  bool tgtIsNucleus = target.IsNucleus();
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
  if ( q0Tilde <= 0 ) return 0;

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
  TVector3 leptonMomCoulomb3 = !fCoulomb ? leptonMom3: PlLocal*leptonMom3*(1/leptonMom3.Mag());
  TVector3 q3VecTilde = neutrinoMom3 - leptonMomCoulomb3;
  // Calculate qTilde
  TLorentzVector qTildeP4(0., 0., q3VecTilde.Mag(), q0Tilde);
  double Q2tilde = -qTildeP4.Mag2();
  
  // Check that Q2tilde > 0 (accounting for rounding errors)
  if (Q2tilde < 0) return 0;

  // Store Q2tilde in the kinematic variable representing Q2.
  // This will ensure that the form factors are calculated correctly
  // using the de Forest prescription (Q2tilde instead of Q2).
  interaction->KinePtr()->SetQ2(Q2tilde);
  // Calculate form factors
  fFormFactors.Calculate( interaction );

  // Get the QEL form factors (were calculated before this method was called)
  double F1V   = 0.5*fFormFactors.F1V();
  double xiF2V = 0.5*fFormFactors.xiF2V();
  double FA    = fFormFactors.FA();
  // According to Nieves' paper, Fp = 2.0*M*FA/(kPionMass2-q2), but Llewelyn-
  // Smith uses Fp = 2.0*M^2*FA/(kPionMass2-q2), so I divide by M
  // This gives units of GeV^-1
  double Fp = fFormFactors.Fp()/M;

  double CN(1), CT(1), CL(1);
  if (fRPA) 
  {
    CNCTCLimUcalc(qTildeP4, M, r, tgtIsNucleus, A, Z, N, 
    CN, CT, CL, interaction->TestBit(kIAssumeFreeNucleon) );
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

  double factor  = PlLocal*ElLocal*r*r/dq;

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
  
  double Ep      = M*a7;
  double Ep2     = M2*a4;
  double px2     = 0.5*M2*(a2 - a3); //py2=px2
  double pz      = M*a6;
  double Eppz    = M2*a5;
  double pz2     = M2*a3;
  
  double aux  = 2*CL*Fp*(Fp*q2 + 4*FA*M);
  double W00   = 32*F1V2*(Ep2*CN + Ep*q0 + a1*q2/4)+
                 8*q2*xiF2V2*(a1*(1 - q02*(1/q2 + 1/M2/4)) - Ep2/M2 - Ep*q0/M2) +
                 8*FA2*(Ep2 + Ep*q0 + a1*(q2/4 - M2)) - a1*(aux*q02 + 16*F1V*xiF2V*(q02 - q2)*CN);
  double Wxx   = 32*F1V2*(px2 - a1*q2/4) - 8*q2*xiF2V2*(a1*CT + px2/M2) +
                 8*FA2*(px2 + a1*(CT*M2 - q2/4)) - 16*a1*F1V*xiF2V*CT*q2;
  double Wzz   = 32*F1V2*(pz2 + pz*dq - a1*q2/4)-
                 8*q2*xiF2V2*(a1 + pz2/M2 + pz*dq/M2 + a1*dq2*(1/q2 + 1/M2/4))+
                 8*FA2*(pz2 + pz*dq + a1*(CL*M2 - q2/4)) - a1*(aux*dq2 + 16*F1V*xiF2V*q02);
  double ReW0z = 16*F1V2*((2*Eppz + Ep*dq)*CN + pz*q0)
                -4*q2*xiF2V2*(2*Eppz/M2 + (Ep*dq + pz*q0)/M2 + a1*dq*q0*(2/q2 + 1/M2/2))+
                 4*FA2*((2*Eppz + Ep*dq)*CL + pz*q0) - a1*dq*q0*(aux + 16*F1V*xiF2V);
  double ImWxy = -16*FA*(xiF2V+F1V)*(pz*q0 - Ep*dq*CT);
  
  
  TLorentzVector q4 = neutrinoMom - leptonMom;
  double v  = q4.E();
  double v2 = v*v;
  double qv2 = q4.Vect().Mag2();
  double qv  = TMath::Sqrt(q2);
  
  double W1 = Wxx/2/M;
  double W2 = (W00 + Wxx + v2/qv2*(Wzz - Wxx) - 2*v/qv*ReW0z)/2/M;
  double W3 = -ImWxy/qv;
  double W4 = M/2/qv2*(Wzz - Wxx);
  double W5 = (ReW0z - v/qv*(Wzz - Wxx))/qv;

  double Ev = neutrinoMom.E();
  double El = leptonMom.E();
  double Pl = leptonMom.P();
  double cost = TMath::Cos( neutrinoMom.Angle(leptonMom.Vect()) );
  double sint = TMath::Sqrt(1 - cost*cost);
  
  double auxm  = (El - Pl*cost)/2/M;
  double auxp  = (El + Pl*cost)/2/M;
  double aux1m = (Pl - El*cost)/2/M;
  double aux1p = (Pl + El*cost)/2/M;
  double aux1  = ml2/2/M2;
  double aux2  = (Ev + El)/M;
  
  if (mod == 2) return sign*(2*aux1m*(W1 - aux1*W4) + aux1p*W2 - sign*(aux2*aux1m + aux1*cost)*W3 - aux1*cost*W5); //PL*R
  if (mod == 3) return sign*ml*sint*(2*W1 - W2 -sign*Ev*W3/M - ml2*W4/M2 + El*W5/M)/2/M;  //PP*R
  double R = 2*auxm*(W1 + aux1*W4) + auxp*W2 - sign*(aux2*auxm - aux1)*W3 - aux1*W5;
  if (mod == 1) return R;
  
  double extrafactor  = kGF2*fCos8c2*TMath::Sq(kMw2/(kMw2 - q4.Mag2() ) )/El;
  
  // Calculate xsec
  double xsec = extrafactor*R;
  xsec *= fXSecCCScale ;

  // Number of scattering centers in the target
  int nucpdgc = target.HitNucPdg();
  int NNucl = (pdg::IsProton(nucpdgc)) ? target.Z() : target.N();

  xsec *= NNucl; // nuclear xsec

  return xsec;
  
}
//___________________________________________________________________________________
