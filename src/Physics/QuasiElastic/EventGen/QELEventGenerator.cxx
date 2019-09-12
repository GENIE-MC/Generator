//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Mar 03, 2009 - CA
   Moved into the new QEL package from its previous location (EVGModules)
 @ Mar 05, 2010 - CA
   Added a temprorary SpectralFuncExperimentalCode()
 @ Feb 06, 2013 - CA
   When the value of the differential cross-section for the selected kinematics
   is set to the event, set the corresponding KinePhaseSpace_t value too.
 @ Feb 14, 2013 - CA
   Temporarily disable the kinematical transformation that takes out the
   dipole form from the dsigma/dQ2 p.d.f.
 @ 2015 - AF
   New QELEventgenerator class replaces previous methods in QEL.
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/RunningThreadInfo.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Physics/QuasiElastic/EventGen/QELEventGenerator.h"

#include "Physics/NuclearState/NuclearModelI.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/PrintUtils.h"

using namespace genie;
using namespace genie::controls;
using namespace genie::constants;
using namespace genie::utils;

//___________________________________________________________________________
QELEventGenerator::QELEventGenerator() :
    KineGeneratorWithCache("genie::QELEventGenerator")
{

}
//___________________________________________________________________________
QELEventGenerator::QELEventGenerator(string config) :
    KineGeneratorWithCache("genie::QELEventGenerator", config)
{

}
//___________________________________________________________________________
QELEventGenerator::~QELEventGenerator()
{

}
//___________________________________________________________________________
void QELEventGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
    LOG("QELEvent", pDEBUG) << "Generating QE event kinematics...";

    // Get the random number generators
    RandomGen * rnd = RandomGen::Instance();

    // Access cross section algorithm for running thread
    RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
    const EventGeneratorI * evg = rtinfo->RunningThread();
    fXSecModel = evg->CrossSectionAlg();

    // Get the interaction and check we are working with a nuclear target
    Interaction * interaction = evrec->Summary();
    // Skip if not a nuclear target
    if(interaction->InitState().Tgt().IsNucleus()) {
        // Skip if no hit nucleon is set
        if(! evrec->HitNucleon()) {
            LOG("QELEvent", pFATAL) << "No hit nucleon was set";
            gAbortingInErr = true;
            exit(1);
        }
    } // is nuclear target

    // set the 'trust' bits
    interaction->SetBit(kISkipProcessChk);
    interaction->SetBit(kISkipKinematicChk);

    // Access the hit nucleon and target nucleus entries at the GHEP record
    GHepParticle * nucleon = evrec->HitNucleon();
    GHepParticle * nucleus = evrec->TargetNucleus();
    bool have_nucleus = (nucleus != 0);

    // Store the hit nucleon radius before computing the maximum differential
    // cross section (important when using the local Fermi gas model)
    Target* tgt = interaction->InitState().TgtPtr();
    double hitNucPos = nucleon->X4()->Vect().Mag();
    tgt->SetHitNucPosition( hitNucPos );

    //-- For the subsequent kinematic selection with the rejection method:
    //   Calculate the max differential cross section or retrieve it from the
    //   cache. Throw an exception and quit the evg thread if a non-positive
    //   value is found.
    //   If the kinematics are generated uniformly over the allowed phase
    //   space the max xsec is irrelevant
    double xsec_max = (fGenerateUniformly) ? -1 : this->MaxXSec(evrec);

    // For a composite nuclear target, check to make sure that the
    // final nucleus has a recognized PDG code
    if ( have_nucleus ) {
      // compute A,Z for final state nucleus & get its PDG code
      int nucleon_pdgc = nucleon->Pdg();
      bool is_p  = pdg::IsProton(nucleon_pdgc);
      int Z = (is_p) ? nucleus->Z()-1 : nucleus->Z();
      int A = nucleus->A() - 1;
      TParticlePDG * fnucleus = 0;
      int ipdgc = pdg::IonPdgCode(A, Z);
      fnucleus = PDGLibrary::Instance()->Find(ipdgc);
      if (!fnucleus) {
        LOG("QELEvent", pFATAL)
            << "No particle with [A = " << A << ", Z = " << Z
            << ", pdgc = " << ipdgc << "] in PDGLibrary!";
        exit(1);
      }
    }

    // In the accept/reject loop, each iteration samples a new value of
    //   - the hit nucleon 3-momentum,
    //   - its binding energy (only actually used if fHitNucleonBindingMode == kUseNuclearModel)
    //   - the final lepton scattering angles in the neutrino-and-hit-nucleon COM frame
    //     (measured with respect to the velocity of the COM frame as seen in the lab frame)
    unsigned int iter = 0;
    bool accept = false;
    while (1) {

        iter++;
        LOG("QELEvent", pINFO) << "Attempt #: " << iter;
        if(iter > kRjMaxIterations) {
            LOG("QELEvent", pWARN)
                << "Couldn't select a valid (pNi, Eb, cos_theta_0, phi_0) tuple after "
                << iter << " iterations";
            evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
            genie::exceptions::EVGThreadException exception;
            exception.SetReason("Couldn't select kinematics");
            exception.SwitchOnFastForward();
            throw exception;
        }

        // If the target is a composite nucleus, then sample an initial nucleon
        // 3-momentum and removal energy from the nuclear model.
        if ( tgt->IsNucleus() ) {
          fNuclModel->GenerateNucleon(*tgt, hitNucPos);
        }
        else {
          // Otherwise, just set the nucleon to be at rest in the lab frame and
          // unbound. Use the nuclear model to make these assignments. The call
          // to BindHitNucleon() will apply them correctly below.
          fNuclModel->SetMomentum3( TVector3(0., 0., 0.) );
          fNuclModel->SetRemovalEnergy( 0. );
        }

        // Put the hit nucleon off-shell (if needed) so that we can get the correct
        // value of cos_theta0_max
        genie::utils::BindHitNucleon(*interaction, *fNuclModel,
          fEb, fHitNucleonBindingMode);

        double cos_theta0_max = std::min(1., CosTheta0Max(*interaction));

        // If the allowed range of cos(theta_0) is vanishing, skip doing the
        // full differential cross section calculation (it will be zero)
        if ( cos_theta0_max <= -1. ) continue;

        // Pick a direction
        // NOTE: In the kPSQELEvGen phase space used by this generator,
        // these angles are specified with respect to the velocity of the
        // probe + hit nucleon COM frame as measured in the lab frame. That is,
        // costheta = 1 means that the outgoing lepton's COM frame 3-momentum
        // points parallel to the velocity of the COM frame.
        double costheta = rnd->RndKine().Uniform(-1., cos_theta0_max); // cosine theta
        double phi = rnd->RndKine().Uniform( 2.*kPi ); // phi: [0, 2pi]

        // Set the "bind_nucleon" flag to false in this call to ComputeFullQELPXSec
        // since we've already done that above
        LOG("QELEvent", pDEBUG) << "cth0 = " << costheta << ", phi0 = " << phi;
        double xsec = genie::utils::ComputeFullQELPXSec(interaction, fNuclModel,
          fXSecModel, costheta, phi, fEb, fHitNucleonBindingMode, fMinAngleEM, false);

        // select/reject event
        this->AssertXSecLimits(interaction, xsec, xsec_max);

        double t = xsec_max * rnd->RndKine().Rndm();

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
        LOG("QELEvent", pDEBUG)
            << "xsec= " << xsec << ", Rnd= " << t;
#endif
        accept = (t < xsec);

        // If the generated kinematics are accepted, finish-up module's job
        if(accept) {
            double gQ2 = interaction->KinePtr()->Q2(false);
            LOG("QELEvent", pINFO) << "*Selected* Q^2 = " << gQ2 << " GeV^2";

            // reset bits
            interaction->ResetBit(kISkipProcessChk);
            interaction->ResetBit(kISkipKinematicChk);
            interaction->ResetBit(kIAssumeFreeNucleon);

            // get neutrino energy at struck nucleon rest frame and the
            // struck nucleon mass (can be off the mass shell)
            const InitialState & init_state = interaction->InitState();
            double E  = init_state.ProbeE(kRfHitNucRest);
            double M = init_state.Tgt().HitNucP4().M();
            LOG("QELKinematics", pNOTICE) << "E = " << E << ", M = "<< M;

            // The hadronic inv. mass is equal to the recoil nucleon on-shell mass.
            // For QEL/Charm events it is set to be equal to the on-shell mass of
            // the generated charm baryon (Lamda_c+, Sigma_c+ or Sigma_c++)
            // Similarly for strange baryons
            //
            const XclsTag & xcls = interaction->ExclTag();
            int rpdgc = 0;
            if (xcls.IsCharmEvent()) {
                rpdgc = xcls.CharmHadronPdg();
            } else if (xcls.IsStrangeEvent()) {
                rpdgc = xcls.StrangeHadronPdg();
            } else {
                rpdgc = interaction->RecoilNucleonPdg();
            }
            assert(rpdgc);
            double gW = PDGLibrary::Instance()->Find(rpdgc)->Mass();
            LOG("QELEvent", pNOTICE) << "Selected: W = "<< gW;

            // (W,Q2) -> (x,y)
            double gx=0, gy=0;
            kinematics::WQ2toXY(E,M,gW,gQ2,gx,gy);

            // lock selected kinematics & clear running values
            interaction->KinePtr()->SetQ2(gQ2, true);
            interaction->KinePtr()->SetW (gW,  true);
            interaction->KinePtr()->Setx (gx,  true);
            interaction->KinePtr()->Sety (gy,  true);
            interaction->KinePtr()->ClearRunningValues();

            // set the cross section for the selected kinematics
            evrec->SetDiffXSec(xsec, kPSQELEvGen);

            TLorentzVector lepton(interaction->KinePtr()->FSLeptonP4());
            TLorentzVector outNucleon(interaction->KinePtr()->HadSystP4());
            TLorentzVector x4l(*(evrec->Probe())->X4());

            evrec->AddParticle(interaction->FSPrimLeptonPdg(), kIStStableFinalState,
              evrec->ProbePosition(), -1, -1, -1, interaction->KinePtr()->FSLeptonP4(), x4l);

            GHepStatus_t ist = (tgt->IsNucleus()) ? kIStHadronInTheNucleus : kIStStableFinalState;
            evrec->AddParticle(interaction->RecoilNucleonPdg(), ist, evrec->HitNucleonPosition(),
              -1, -1, -1, interaction->KinePtr()->HadSystP4(), x4l);

            // Store struck nucleon momentum and binding energy
            TLorentzVector p4ptr = interaction->InitStatePtr()->TgtPtr()->HitNucP4();
            LOG("QELEvent",pNOTICE) << "pn: " << p4ptr.X() << ", "
              << p4ptr.Y() << ", " << p4ptr.Z() << ", " << p4ptr.E();
            nucleon->SetMomentum(p4ptr);
            nucleon->SetRemovalEnergy(fEb);

            // add a recoiled nucleus remnant
            this->AddTargetNucleusRemnant(evrec);

            break; // done
        } else { // accept throw
            LOG("QELEvent", pDEBUG) << "Reject current throw...";
        }

    } // iterations - while(1) loop
    LOG("QELEvent", pINFO) << "Done generating QE event kinematics!";
}
//___________________________________________________________________________
void QELEventGenerator::AddTargetNucleusRemnant(GHepRecord * evrec) const
{
    // add the remnant nuclear target at the GHEP record

    LOG("QELEvent", pINFO) << "Adding final state nucleus";

    double Px = 0;
    double Py = 0;
    double Pz = 0;
    double E  = 0;

    GHepParticle * nucleus = evrec->TargetNucleus();
    bool have_nucleus = nucleus != 0;
    if (!have_nucleus) return;

    int A = nucleus->A();
    int Z = nucleus->Z();

    int fd = nucleus->FirstDaughter();
    int ld = nucleus->LastDaughter();

    for(int id = fd; id <= ld; id++) {

        // compute A,Z for final state nucleus & get its PDG code and its mass
        GHepParticle * particle = evrec->Particle(id);
        assert(particle);
        int  pdgc = particle->Pdg();
        bool is_p  = pdg::IsProton (pdgc);
        bool is_n  = pdg::IsNeutron(pdgc);

        if (is_p) Z--;
        if (is_p || is_n) A--;

        Px += particle->Px();
        Py += particle->Py();
        Pz += particle->Pz();
        E  += particle->E();

    }//daughters

    TParticlePDG * remn = 0;
    int ipdgc = pdg::IonPdgCode(A, Z);
    remn = PDGLibrary::Instance()->Find(ipdgc);
    if(!remn) {
        LOG("HadronicVtx", pFATAL)
            << "No particle with [A = " << A << ", Z = " << Z
            << ", pdgc = " << ipdgc << "] in PDGLibrary!";
        assert(remn);
    }

    double Mi = nucleus->Mass();
    Px *= -1;
    Py *= -1;
    Pz *= -1;
    E = Mi-E;

    // Add the nucleus to the event record
    LOG("QELEvent", pINFO)
        << "Adding nucleus [A = " << A << ", Z = " << Z
        << ", pdgc = " << ipdgc << "]";

    int imom = evrec->TargetNucleusPosition();
    evrec->AddParticle(
            ipdgc,kIStStableFinalState, imom,-1,-1,-1, Px,Py,Pz,E, 0,0,0,0);

    LOG("QELEvent", pINFO) << "Done";
    LOG("QELEvent", pINFO) << *evrec;
}
//___________________________________________________________________________
void QELEventGenerator::Configure(const Registry & config)
{
    Algorithm::Configure(config);
    this->LoadConfig();
}
//____________________________________________________________________________
void QELEventGenerator::Configure(string config)
{
    Algorithm::Configure(config);
    this->LoadConfig();
}
//____________________________________________________________________________
void QELEventGenerator::LoadConfig(void)
{
    // Load sub-algorithms and config data to reduce the number of registry
    // lookups
        fNuclModel = 0;

    RgKey nuclkey = "NuclearModel";

    fNuclModel = dynamic_cast<const NuclearModelI *> (this->SubAlg(nuclkey));
    assert(fNuclModel);

    // Safety factor for the maximum differential cross section
    GetParamDef( "MaxXSec-SafetyFactor", fSafetyFactor, 1.6  ) ;

    // Minimum energy for which max xsec would be cached, forcing explicit
    // calculation for lower eneries
    GetParamDef( "Cache-MinEnergy", fEMin, 1.00 ) ;
    
    // Maximum allowed fractional cross section deviation from maxim cross
    // section used in rejection method
    GetParamDef( "MaxXSec-DiffTolerance", fMaxXSecDiffTolerance, 999999. ) ;
    assert(fMaxXSecDiffTolerance>=0);

    // Generate kinematics uniformly over allowed phase space and compute
    // an event weight?
    GetParamDef( "UniformOverPhaseSpace", fGenerateUniformly, false ) ;

    GetParamDef( "SF-MinAngleEMscattering", fMinAngleEM, 0. ) ;

    // Decide how to handle the binding energy of the initial state struck
    // nucleon
    std::string binding_mode;
    GetParamDef( "HitNucleonBindingMode", binding_mode, std::string("UseNuclearModel") );

    fHitNucleonBindingMode = genie::utils::StringToQELBindingMode( binding_mode );

    GetParamDef( "MaxXSecNucleonThrows", fMaxXSecNucleonThrows, 800 );
}
//____________________________________________________________________________
double QELEventGenerator::ComputeMaxXSec(const Interaction * in) const
{
    // Computes the maximum differential cross section in the requested phase
    // space. This method overloads KineGeneratorWithCache::ComputeMaxXSec
    // method and the value is cached at a circular cache branch for retrieval
    // during subsequent event generation.
    // The computed max differential cross section does not need to be the exact
    // maximum. The number used in the rejection method will be scaled up by a
    // safety factor. But it needs to be fast.
    LOG("QELEvent", pINFO) << "Computing maximum cross section to throw against";

    double xsec_max = -1;
    double dummy_Eb = 0.;

    // Clone the input interaction so that we can modify it a bit
    Interaction * interaction = new Interaction( *in );
    interaction->SetBit( kISkipProcessChk );
    interaction->SetBit( kISkipKinematicChk );

    // We'll select the max momentum and zero binding energy.
    // That should give us the nucleon with the highest xsec
    const Target& tgt = interaction->InitState().Tgt();
    // Pick a really slow nucleon to start, but not one at rest,
    // since Prob() for the Fermi gas family of models is zero
    // for a vanishing nucleon momentum
    double max_momentum = 0.010;
    double search_step = 0.1;
    const double STEP_STOP = 1e-6;
    while ( search_step > STEP_STOP ) {
      double pNi_next = max_momentum + search_step;

      // Set the nucleon we're using to be upstream at max energy and unbound
      fNuclModel->SetMomentum3( TVector3(0., 0., -pNi_next) );
      fNuclModel->SetRemovalEnergy( 0. );

      // Set the nucleon total energy to be on-shell with a quick call to
      // BindHitNucleon()
      genie::utils::BindHitNucleon(*interaction, *fNuclModel, dummy_Eb, kOnShell);

      // TODO: document this, won't work for spectral functions
      double dummy_w = -1.;
      double prob = fNuclModel->Prob(pNi_next, dummy_w, tgt,
        tgt.HitNucPosition());

      double costh0_max = genie::utils::CosTheta0Max( *interaction );

      if ( prob > 0. && costh0_max > -1. ) max_momentum = pNi_next;
      else search_step /= 2.;
    }

    {
        // Set the nucleon we're using to be upstream at max energy and unbound
        fNuclModel->SetMomentum3( TVector3(0., 0., -max_momentum) );
        fNuclModel->SetRemovalEnergy( 0. );

        // Set the nucleon total energy to be on-shell with a quick call to
        // BindHitNucleon()
        genie::utils::BindHitNucleon(*interaction, *fNuclModel, dummy_Eb, kOnShell);

        // Just a scoping block for now
        // OK, we're going to scan the COM frame angles to get the point of max xsec
        // We'll bin in solid angle, and find the maximum point
        // Then we'll bin/scan again inside that point
        // Rinse and repeat until the xsec stabilises to within some fraction of our safety factor
        const double acceptable_fraction_of_safety_factor = 0.2;
        const int max_n_layers = 100;
        const int N_theta = 10;
        const int N_phi = 10;
        double phi_at_xsec_max = -1;
        double costh_at_xsec_max = 0;
        double this_nuc_xsec_max = -1;

        double costh_range_min = -1.;
        double costh_range_max = genie::utils::CosTheta0Max( *interaction );
        double phi_range_min = 0.;
        double phi_range_max = 2*TMath::Pi();
        for (int ilayer = 0 ; ilayer < max_n_layers ; ilayer++) {
          double last_layer_xsec_max = this_nuc_xsec_max;
          double costh_increment = (costh_range_max-costh_range_min) / N_theta;
          double phi_increment   = (phi_range_max-phi_range_min) / N_phi;
          // Now scan through centre-of-mass angles coarsely
          for (int itheta = 0; itheta < N_theta; itheta++){
              double costh = costh_range_min + itheta * costh_increment;
              for (int iphi = 0; iphi < N_phi; iphi++) { // Scan around phi
                  double phi = phi_range_min + iphi * phi_increment;
                  // We're after an upper limit on the cross section, so just
                  // put the nucleon on-shell and call it good. The last
                  // argument is false because we've already called
                  // BindHitNucleon() above
                  double xs = genie::utils::ComputeFullQELPXSec(interaction,
                    fNuclModel, fXSecModel, costh, phi, dummy_Eb, kOnShell, fMinAngleEM, false);

                  if (xs > this_nuc_xsec_max){
                      phi_at_xsec_max = phi;
                      costh_at_xsec_max = costh;
                      this_nuc_xsec_max = xs;
                  }
                  //
              } // Done with phi scan
          }// Done with centre-of-mass angles coarsely

          // Calculate the range for the next layer
          costh_range_min = costh_at_xsec_max - costh_increment;
          costh_range_max = costh_at_xsec_max + costh_increment;
          phi_range_min = phi_at_xsec_max - phi_increment;
          phi_range_max = phi_at_xsec_max + phi_increment;

          double improvement_factor = this_nuc_xsec_max/last_layer_xsec_max;
          if (ilayer && (improvement_factor-1) < acceptable_fraction_of_safety_factor * (fSafetyFactor-1)) {
            break;
          }
        }
        if (this_nuc_xsec_max > xsec_max) {
            xsec_max = this_nuc_xsec_max;
            LOG("QELEvent", pINFO) << "best estimate for xsec_max = " << xsec_max;
        }

        delete interaction;
    }
    // Apply safety factor, since value retrieved from the cache might
    // correspond to a slightly different energy
    xsec_max *= fSafetyFactor;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    SLOG("QELEvent", pDEBUG) << interaction->AsString();
    SLOG("QELEvent", pDEBUG) << "Max xsec in phase space = " << max_xsec;
    SLOG("QELEvent", pDEBUG) << "Computed using alg = " << *fXSecModel;
#endif

    LOG("QELEvent", pINFO) << "Computed maximum cross section to throw against - value is " << xsec_max;
    return xsec_max;
}
