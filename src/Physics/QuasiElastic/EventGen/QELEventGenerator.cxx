//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
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
    }//is nucl target

    // set the 'trust' bits
    interaction->SetBit(kISkipProcessChk);
    interaction->SetBit(kISkipKinematicChk);
    // Note: The kinematic generator would be using the free nucleon cross
    // section (even for nuclear targets) so as not to double-count nuclear
    // suppression. This assumes that a) the nuclear suppression was turned
    // on when computing the cross sections for selecting the current event
    // and that b) if the event turns out to be unphysical (Pauli-blocked)
    // the next attempted event will be forced to QEL again.
    // (discussion with Hugh - GENIE/NeuGEN integration workshop - 07APR2006
    interaction->SetBit(kIAssumeFreeNucleon);


    //-- For the subsequent kinematic selection with the rejection method:
    //   Calculate the max differential cross section or retrieve it from the
    //   cache. Throw an exception and quit the evg thread if a non-positive
    //   value is found.
    //   If the kinematics are generated uniformly over the allowed phase
    //   space the max xsec is irrelevant
    double xsec_max = (fGenerateUniformly) ? -1 : this->MaxXSec(evrec);

    //
    // Try to generate (simultaneously):
    //    - Fermi momentum (pF),
    //    - binding energy (w) and
    //    - momentum transfer (Q2)
    //

    unsigned int iter = 0;
    bool accept = false;
    // Access the hit nucleon and target nucleus entries at the GHEP record
    GHepParticle * nucleon = evrec->HitNucleon();
    GHepParticle * nucleus = evrec->TargetNucleus();
    bool have_nucleus = nucleus != 0;

    Target * tgt = interaction->InitState().TgtPtr();
    //TLorentzVector * p4 = tgt->HitNucP4Ptr();
    // Store the hit nucleon radius first
    double hitNucPos = nucleon->X4()->Vect().Mag();
    tgt->SetHitNucPosition(hitNucPos);
    while(1) {
        iter++;
        LOG("QELEvent", pINFO) << "Attempt #: " << iter;
        if(iter > kRjMaxIterations) {
            LOG("QELEvent", pWARN)
                << "Couldn't select a valid (pF, w, Q^2) tuple after "
                << iter << " iterations";
            evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
            genie::exceptions::EVGThreadException exception;
            exception.SetReason("Couldn't select kinematics");
            exception.SwitchOnFastForward();
            throw exception;
        }

        // First, throw Fermi momentum & removal energy from the nuclear model pdfs
        fNuclModel->GenerateNucleon(*tgt, hitNucPos);

        // The nucleon is now accessed in the CalculateXSec method
        //        TVector3 p3 = fNuclModel->Momentum3();
        //        double w    = fNuclModel->RemovalEnergy();
        //
        //        double pF  = p3.Mag();  // (fermi momentum)
        //        double pF2 = p3.Mag2(); // (fermi momentum)^2

        //        LOG("QELEvent", pINFO)
        //            << "Generated nucleon momentum: ("
        //            << p3.Px() << ", " << p3.Py() << ", " << p3.Pz() << "), "
        //            << "|p| = " << pF;
        //        LOG("QELEvent", pINFO)
        //            << "Generated nucleon removal energy: w = " << w;

        if (have_nucleus) {
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

        //
        // Calculate the on-shell and off-shell masses for the struck nucleon
        //
        //         double Mf  = fnucleus -> Mass(); // remnant nucleus mass
        //         double Mi  = nucleus  -> Mass(); // initial nucleus mass
        //       TDatabasePDG *tb = TDatabasePDG::Instance();
        //         double Mn = tb->GetParticle(interaction->InitState().TgtPtr()->HitNucPdg())->Mass();// incoming nucleon mass
        //         double Mp = tb->GetParticle(interaction->RecoilNucleonPdg())->Mass(); // outgoing nucleon mass

        //         double Mn  = nucleon->Mass();
        //         double Mp(0); // outgoing nucleon mass
        //         TDatabasePDG *tb = TDatabasePDG::Instance();
        //         if (nucleon->Pdg() == 2212){
        //           Mp = tb->GetParticle(2112)->Mass();
        //         }
        //         else if (nucleon->Pdg() == 2112){
        //           Mp = tb->GetParticle(2212)->Mass();
        //         }
        //         else{LOG("QELEvent",pDEBUG) << "ERROR - incoming particle not a proton or neutron" << std::endl;}
        //
        //         //double EN_offshell = Mi - TMath::Sqrt(pF2 + Mf*Mf);
        //         double EN_offshell = TMath::Sqrt(Mn*Mn + pF2) - w; // old
        //         double EN_offshell = Mn - w ;
        //         double EN_onshell  = TMath::Sqrt(pF2+Mn*Mn);
        //
        //         double Mn_offshell = TMath::Sqrt(EN_offshell*EN_offshell - pF2);
        //
        //         // Calculate the binding energy
        //         //double Eb = w;
        //         double Eb = EN_onshell - EN_offshell;
        //
        //         LOG("QELEvent", pINFO)
        //          << "Enuc (on-shell) = " << EN_onshell
        //          << " GeV, Enuc (off-shell) = " << EN_offshell
        //          << " GeV, Ebind = " << Eb << " GeV";
        //
        //         // Update the struck nucleon 4-momentum at the interaction summary
        //         // and at the GHEP record
        //         p4->SetPx( p3.Px()    );
        //         p4->SetPy( p3.Py()    );
        //         p4->SetPz( p3.Pz()    );
        //       //p4->SetE ( EN_onshell );
        //         p4->SetE ( EN_offshell);
        //         evrec->Summary()->InitStatePtr()->TgtPtr()->SetHitNucP4(*p4);
        //         nucleon->SetMomentum(*p4); // update GHEP values
        //
        //         // Update the binding energy value at the GHEP record
        //         // nucleon->SetRemovalEnergy(Eb);
        //
        //         //LOG("QELEvent",pDEBUG) << "offshell energy = "<<EN_offshell << ", fRemovalEnergy = "<<w << std::endl;
        //         //LOG("QELEvent",pDEBUG) << "binding energy = "<< Eb << std::endl;
        //
        //        // Sometimes, for interactions near threshold, Fermi momentum might bring
        //        // the neutrino energy in the nucleon rest frame below threshold (for the
        //        // selected interaction). In this case mark, select a new Fermi momentum.
        //        const KPhaseSpace & kps = interaction->PhaseSpace();
        //        if(!kps.IsAboveThreshold()) {
        //            LOG("QELEvent", pNOTICE)
        //                  << "Event below threshold after generating Fermi momentum";
        //            double Ethr = kps.Threshold();
        //            double Ev   = interaction->InitState().ProbeE(kRfHitNucRest);
        //            LOG("QELEvent", pNOTICE)
        //                << "Ev (@ nucleon rest frame) = " << Ev << ", Ethr = " << Ethr;
        //            continue;
        //        }
        //
        //        // Get centre-of-mass energy
        //     const InitialState & init_state = interaction->InitState();
        //        double s = init_state.CMEnergy();
        //        s *= s; //* init_state.CMEnergy(); // centre of mass energy squared

        //        if (TMath::Sqrt(s) < interaction->FSPrimLepton()->Mass() + Mp){ // throw a new event if this one is below threshold
        //          LOG("QELEvent", pINFO) << "Event below threshold, reject throw and try again!";
        //          continue;
        //        }

        // Pick a direction
        double costheta = (rnd->RndKine().Rndm() * 2) - 1; // cosine theta: [-1, 1]
        double phi = 2 * TMath::Pi() * rnd->RndKine().Rndm(); // phi: [0, 2pi]

        //        // Generate the outgoing particles
        //        // with the correct momenta
        //        double lepMass = interaction->FSPrimLepton()->Mass();
        //        double outLeptonEnergy = ( s - Mp*Mp + lepMass*lepMass ) / (2 * TMath::Sqrt(s));
        //        double outMomentum = TMath::Sqrt(outLeptonEnergy*outLeptonEnergy - lepMass*lepMass);
        //
        //        TLorentzVector lepton(outMomentum, 0, 0, TMath::Sqrt(outMomentum*outMomentum + lepMass*lepMass));
        //
        //        lepton.SetTheta(TMath::ACos(costheta));
        //        lepton.SetPhi(phi);
        //
        ////        lepton.SetTheta(0.2);
        ////        lepton.SetPhi(0.7);
        //
        //        TLorentzVector outNucleon(-1*lepton.Px(),-1*lepton.Py(),-1*lepton.Pz(), TMath::Sqrt(outMomentum*outMomentum + Mp*Mp));
        //
        //        // Boost particles
        //        TVector3 beta = this->COMframe2Lab(init_state);
        //
        //        TLorentzVector leptonCOM = TLorentzVector(lepton);
        //
        //        lepton.Boost(beta);
        //        outNucleon.Boost(beta);
        //
        //        // Check if Q2 above Minimum Q2 // important for eA scattering
        //        TLorentzVector qP4 = *(init_state.GetProbeP4()) - lepton;
        //        double Q2 = -1 * qP4.Mag2();
        //        if ( Q2 < kMinQ2Limit){
        //          continue;
        //        }
        //
        //        interaction->KinePtr()->SetFSLeptonP4(lepton);
        //        interaction->KinePtr()->SetHadSystP4(outNucleon);
        //

        double xsec = this->ComputeXSec(interaction, costheta, phi);

        // select/reject event
        this->AssertXSecLimits(interaction, xsec, xsec_max);

        double t = xsec_max * rnd->RndKine().Rndm();
        //        LOG("QELEvent", pNOTICE) << "dsigma/dQ2 (random) = " << t/(1E-38*units::cm2) << " 1E-38 cm^2/GeV^2";

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

            evrec->AddParticle(interaction->FSPrimLeptonPdg(), kIStStableFinalState, evrec->ProbePosition(),-1,-1,-1, interaction->KinePtr()->FSLeptonP4(), x4l);

            GHepStatus_t ist = (tgt->IsNucleus()) ?
                kIStHadronInTheNucleus : kIStStableFinalState;
            evrec->AddParticle(interaction->RecoilNucleonPdg(), ist, evrec->HitNucleonPosition(),-1,-1,-1, interaction->KinePtr()->HadSystP4(), x4l);

            // Store struck nucleon momentum and binding energy
            TLorentzVector p4ptr = interaction->InitStatePtr()->TgtPtr()->HitNucP4();
            LOG("QELEvent",pNOTICE) << "pn: " << p4ptr.X() << ", " <<p4ptr.Y() << ", " <<p4ptr.Z() << ", " <<p4ptr.E();
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
    //  RgAlg nuclalg = fConfig->GetAlgDef(nuclkey, gc->GetAlg(nuclkey));
    //  LOG("FermiMover", pINFO) << "Loading nuclear model: " << nuclalg;

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

    //  fQ2min   = 99999999;
    //  fQ2max   = -1;
    GetParamDef( "SF-MinAngleEMscattering", fMinAngleEM, 0. ) ;

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

    const int nnucthrows = 800;
    double min_energy   = 9E9;
    double max_momentum = -9E9;
    // Loop over tthrown nucleons
    // We'll select the max momentum and the minimum binding energy
    // Which should give us the nucleon with the highest xsec
    for(int inuc = 0; inuc <nnucthrows; inuc++) {

        Interaction * interaction = new Interaction(*in);
        interaction->SetBit(kISkipProcessChk);
        interaction->SetBit(kISkipKinematicChk);
        interaction->SetBit(kIAssumeFreeNucleon);

        // Access the target from the interaction summary
        Target * tgt = interaction->InitState().TgtPtr();
        //        TLorentzVector * p4 = tgt->HitNucP4Ptr();

        // First, throw Fermi momentum & removal energy from the nuclear model pdfs
        // Use r=0. as the radius, since this method should give the max xsec
        // for all possible kinematics
        fNuclModel->GenerateNucleon(*tgt,0.0);

        min_energy   = std::min(min_energy  ,fNuclModel->RemovalEnergy());
        max_momentum = std::max(max_momentum,fNuclModel->Momentum());
        delete interaction;
    } // nucl throws

    { // Just a scoping block for now
        Interaction * interaction = new Interaction(*in);
        interaction->SetBit(kISkipProcessChk);
        interaction->SetBit(kISkipKinematicChk);
        interaction->SetBit(kIAssumeFreeNucleon);

        // Access the target from the interaction summary
        Target * tgt = interaction->InitState().TgtPtr();
        //        TLorentzVector * p4 = tgt->HitNucP4Ptr();

        // Set the nucleon we're using to be upstream at max enegry
        fNuclModel->GenerateNucleon(*tgt,0.0);
        fNuclModel->SetMomentum3(TVector3(0.,0.,-max_momentum));
        fNuclModel->SetRemovalEnergy(min_energy);

        //        TVector3 p3 = fNuclModel->Momentum3();
        //        double w = fNuclModel->RemovalEnergy();
        //
        //        TDatabasePDG *tb = TDatabasePDG::Instance();
        //        double Mn = tb->GetParticle(interaction->InitState().TgtPtr()->HitNucPdg())->Mass();// outgoing nucleon mass
        //        double Mp = tb->GetParticle(interaction->RecoilNucleonPdg())->Mass(); // incoming nucleon mass
        //        double EN_offshell = Mn - w;
        ////      double EN_offshell = TMath::Sqrt(Mn*Mn + p3.Dot(p3)) - w;
        //
        //        p4->SetPx( p3.Px()    );
        //        p4->SetPy( p3.Py()    );
        //        p4->SetPz( p3.Pz()    );
        //        p4->SetE ( EN_offshell );

        // OK, we're going to scan the centre-of-mass angles to get the point of max xsec
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
        double costh_range_max = 1.;
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
                  double xs = this->ComputeXSec(interaction, costh, phi);
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
        if (this_nuc_xsec_max > xsec_max){
            xsec_max = this_nuc_xsec_max;  // this nucleon has the highest xsec!
            LOG("QELEvent", pINFO) << "best estimate for xsec_max = "<<xsec_max;
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
//____________________________________________________________________________
double QELEventGenerator::ComputeXSec( Interaction * interaction, double costheta, double phi) const
{
    Target * tgt = interaction->InitState().TgtPtr();
    TLorentzVector * p4 = tgt->HitNucP4Ptr();

    TVector3 p3 = fNuclModel->Momentum3();
    //  double w = fNuclModel->RemovalEnergy();

    double xsec = 0;
    double pF2 = p3.Mag2(); // (fermi momentum)^2
    double lepMass = interaction->FSPrimLepton()->Mass();

    TDatabasePDG *tb = TDatabasePDG::Instance();
    double Mn = tb->GetParticle(interaction->InitState().TgtPtr()->HitNucPdg())->Mass();// outgoing nucleon mass
    double Mp = tb->GetParticle(interaction->RecoilNucleonPdg())->Mass(); // incoming nucleon mass

    double EN_offshell(0);
    //FermiMoverInteractionType_t interaction_type = fNuclModel->GetFermiMoverInteractionType(); // check the nuclear model essentially
    //if  (interaction_type == kFermiMoveBenharSF){
    //EN_offshell = Mn - w;
    //LOG("QELEvent",pDEBUG) << "Using FermiMoveBenharSF and defined energy" << std::endl;
    //}
    //else if (interaction_type == kFermiMoveDefault){
    //GHepParticle * nucleon = evrec->HitNucleon();
    //GHepParticle * nucleus = evrec->TargetNucleus();
    //assert(nucleon);
    //assert(nucleus);
    //int nucleon_pdgc = nucleon->Pdg();
    //bool is_p  = pdg::IsProton(nucleon_pdgc);
    //int A = nucleus->A() - 1;
    //int Z = (is_p) ? nucleus->Z()-1 : nucleus->Z();
    //TParticlePDG * fnucleus = 0;
    //int ipdgc = pdg::IonPdgCode(A, Z);
    //fnucleus = PDGLibrary::Instance()->Find(ipdgc);

    //double Mf  = fnucleus -> Mass(); // remnant nucleus mass
    //double Mi  = nucleus  -> Mass(); // initial nucleus mass

    //// Get initial and final nucleus masses
    //int nucleus_init = tgt->Pdg();
    //TParticlePDG * p = PDGLibrary::Instance()->Find(nucleus_init);
    //double Mi = p->Mass();
    //
    //int A = tgt->A() - 1;
    //double Mf;
    //if(A>0){
    //  bool is_p = pdg::IsProton(tgt->HitNucPdg());
    //  int Z = (is_p) ? tgt->Z()-1 : tgt->Z();
    //  int nucleus_remnant = pdg::IonPdgCode(A, Z);
    //  p = PDGLibrary::Instance()->Find(nucleus_remnant);
    //  Mf = p->Mass();
    //}else{
    //  Mf = 0.;
    //}
    //
    //EN_offshell = Mi - TMath::Sqrt(pF2 + Mf*Mf);
    //LOG("QELEvent",pDEBUG) << "Using FermiMoveDefault and defined energy" << std::endl;
    //}

    //TESTING: Always leave the nucleon on shell
    EN_offshell = TMath::Sqrt(pF2+Mn*Mn);
    //EN_offshell = TMath::Sqrt(Mn*Mn + pF2) - w;
    double EN_onshell  = TMath::Sqrt(pF2+Mn*Mn);
    p4->SetPx( p3.Px()    );
    p4->SetPy( p3.Py()    );
    p4->SetPz( p3.Pz()    );
    p4->SetE ( EN_offshell );

    fEb = EN_onshell - EN_offshell;

    double s = interaction->InitState().CMEnergy(); // actually sqrt(s)
    s *= s; // now s actually = s


    if (TMath::Sqrt(s) < interaction->FSPrimLepton()->Mass() + Mp){ // throw a new event if this one is below threshold
        LOG("QELEvent", pINFO) << "Event below threshold, reject throw and try again!";
        return 0.;
    }

    double outLeptonEnergy = ( s - Mp*Mp + lepMass*lepMass ) / (2 * TMath::Sqrt(s));

    if(outLeptonEnergy*outLeptonEnergy-lepMass*lepMass < 0.) return 0.;
    double outMomentum = TMath::Sqrt(outLeptonEnergy*outLeptonEnergy - lepMass*lepMass);
    //LOG("QELEvent",pDEBUG) << "calculated root s and outLeptonEnergy" << std::endl;

    TLorentzVector lepton(outMomentum, 0, 0, outLeptonEnergy);

    lepton.SetTheta(TMath::ACos(costheta));
    lepton.SetPhi(phi);
    TLorentzVector outNucleon(-1*lepton.Px(),-1*lepton.Py(),-1*lepton.Pz(),
            TMath::Sqrt(outMomentum*outMomentum + Mp*Mp));

    /*LOG("QELEvent",pDEBUG) << "costheta = " << costheta << ", phi = " << phi << std::endl;
      LOG("QELEvent",pDEBUG) << "lepton = " << utils::print::P4AsString(&lepton) << std::endl;
      LOG("QELEvent",pDEBUG) << "outNucleon = " << utils::print::P4AsString(&outNucleon) << std::endl;
      LOG("QELEvent",pDEBUG) << "inNucleon = " << utils::print::P4AsString(p4) << std::endl;
      LOG("QELEvent",pDEBUG) << "s = " << s << ", outLepE = " << outLeptonEnergy << "outMom = " << outMomentum << std::endl;
      LOG("QELEvent",pDEBUG) << "outMom^2 = " << outLeptonEnergy*outLeptonEnergy - lepMass*lepMass << std::endl;*/

    /*LOG("QELEvent",pDEBUG) << "Lepton COM EvGen:\n";
      lepton.Print();
      LOG("QELEvent",pDEBUG) << "outNucleon COM EvGen:\n";
      outNucleon.Print();*/

    // Boost particles
    TVector3 beta = this->COMframe2Lab(interaction->InitState());
    //LOG("QELEvent", pINFO) << "converted to lab frame";

    TLorentzVector leptonCOM = TLorentzVector(lepton);

    lepton.Boost(beta);
    outNucleon.Boost(beta);
    //LOG("QELEvent", pINFO) << "lepton boosted = " << utils::print::P4AsString(&lepton),
    //LOG("QELEvent", pINFO) << "outNucleon boosted = " << utils::print::P4AsString(&outNucleon);
    /*
       LOG("QELEvent", pINFO) << "neutrino LAB EvGen: ";
       interaction->InitState().GetProbeP4(kRfLab)->Print();
       LOG("QELEvent", pINFO) << "inNucleon LAB EvGen:\n";
       interaction->InitState().Tgt().HitNucP4().Print();
       LOG("QELEvent", pINFO) << "Lepton LAB EvGen:\n";
       lepton.Print();
       LOG("QELEvent", pINFO) << "outNucleon LAB EvGen:\n";
       outNucleon.Print();
       LOG("QELEvent", pINFO) << "beta EvGen:";
       beta.Print();*/


    // Check if event is at a low angle - if so return 0 and stop wasting time
    //double angle = fConfig->GetDoubleDef("MinAngle",  gc->GetDouble("SF-MinAngleEMscattering"));
    //LOG("QELEvent", pINFO) << "min angle = " << fMinAngleEM;
    if (180 * lepton.Theta() / 3.1415 < fMinAngleEM && interaction->ProcInfo().IsEM()){
        return 0;
    }
    // Check if Q2 above Minimum Q2 // important for eA scattering
    TLorentzVector * nuP4 = interaction->InitState().GetProbeP4();
    TLorentzVector qP4 = *nuP4 - lepton;
    delete nuP4;
    double Q2 = -1 * qP4.Mag2();

    interaction->KinePtr()->SetFSLeptonP4(lepton);
    interaction->KinePtr()->SetHadSystP4(outNucleon);
    interaction->KinePtr()->SetQ2(Q2 /*, true */);  // (val, selected) -> probably don't want selected here

    // Compute the QE cross section for the current kinematics ("~" variables)
    interaction->InitStatePtr()->TgtPtr()->HitNucP4Ptr()->SetE(EN_onshell);
    xsec = fXSecModel->XSec(interaction, kPSTnctnBnctl); //

    interaction->InitStatePtr()->TgtPtr()->HitNucP4Ptr()->SetE(EN_offshell);

    // Multiply xsec by Jacobian
    xsec *= 4*kPi*leptonCOM.P()*leptonCOM.P();

    double jac = this->COMJacobian(lepton, leptonCOM, outNucleon, beta);
    xsec *= jac;

    //// BEGIN DEBUG
    //double debug_xsec = fXSecModel->XSec(interaction, kPSQELEvGen);
    //std::cout << "\nDEBUG: xsec = " << xsec << ", debug_xsec = " << debug_xsec
    //  << " xsec / debug_xsec = " << xsec / debug_xsec << '\n';
    //// END DEBUG

    return xsec;
}
//___________________________________________________________________________
TVector3 QELEventGenerator::COMframe2Lab(InitialState initialState) const
{

    TLorentzVector * k4 = initialState.GetProbeP4(kRfLab);
    TLorentzVector * p4 = initialState.TgtPtr()->HitNucP4Ptr();
    TLorentzVector totMom = *k4 + *p4;

    TVector3 beta = totMom.BoostVector();

    delete k4;

    return beta;
}

//___________________________________________________________________________
double QELEventGenerator::COMJacobian(TLorentzVector lepton,
                                      TLorentzVector /* leptonCOM */,
                                      TLorentzVector outNucleon,
                                      TVector3 beta) const
{

    double gamma = 1. / TMath::Sqrt(1. - beta.Dot(beta)); //  gamma factor

    // angle between muon in com frame and com frame velocity
    double theta0 = lepton.Angle(beta);

    // difference in velocity between lepton and outgoing nucleon
    TLorentzVector leptonVel = TLorentzVector(lepton);
    TLorentzVector nucleonVel = TLorentzVector(outNucleon);
    leptonVel *= 1. / lepton.E();
    nucleonVel *= 1. / outNucleon.E();

    double velDiff = (leptonVel - nucleonVel).P(); // centre-of-mass velocity difference between lepton and nucleon

    double jacobian = TMath::Sqrt(1 + (1 - TMath::Cos(theta0)*TMath::Cos(theta0))*(gamma*gamma - 1)) / velDiff;

    return jacobian;
}
