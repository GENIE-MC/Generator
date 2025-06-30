//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Alfonso Garcia <aagarciasoto \at km3net.de>
 IFIC (Valencia)
*/
//____________________________________________________________________________

#include "Physics/Hadronization/LeptoHadPythia6.h"

#ifdef __GENIE_PYTHIA6_ENABLED__
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,6)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif
// the actual PYTHIA call
extern "C" {
  double pyangl_( double *,  double * );
  void   pykfdi_( int *,  int *, int *, int * );
  void   pyzdis_( int *,  int *, double *, double * );
  void   pyrobo_( int *,  int *, double *, double *, double *, double *, double * );
  void   pydecy_( int * );
  void   py2ent_( int *,  int *, int *, double * );
}
#endif

//____________________________________________________________________________
LeptoHadPythia6::LeptoHadPythia6() :
EventRecordVisitorI("genie::LeptoHadPythia6")
{
  this->Initialize();
}
//____________________________________________________________________________
LeptoHadPythia6::LeptoHadPythia6(string config) :
EventRecordVisitorI("genie::LeptoHadPythia6", config)
{
  this->Initialize();
}
//____________________________________________________________________________
LeptoHadPythia6::~LeptoHadPythia6()
{

}
//____________________________________________________________________________
void LeptoHadPythia6::ProcessEventRecord(GHepRecord * event) const
{

  if(!this->Hadronize(event)) {
    LOG("LeptoHad", pWARN) << "Hadronization failed!";
    event->EventFlags()->SetBitNumber(kHadroSysGenErr, true);
    genie::exceptions::EVGThreadException exception;
    exception.SetReason("Could not simulate the hadronic system");
    exception.SwitchOnFastForward();
    throw exception;
    return;
  }


}
//___________________________________________________________________________
bool LeptoHadPythia6::Hadronize(GHepRecord * 
#ifdef __GENIE_PYTHIA6_ENABLED__
  event // avoid unused variable warning if PYTHIA6 is not enabled
#endif
) const
{

#ifdef __GENIE_PYTHIA6_ENABLED__
  // Compute kinematics of hadronic system with energy/momentum conservation
  LongLorentzVector p4v( * event->Probe()->P4()                   );
  LongLorentzVector p4N( * event->HitNucleon()->P4()              );
  LongLorentzVector p4l( * event->FinalStatePrimaryLepton()->P4() );
  LongLorentzVector p4Hadlong( p4v.Px()+p4N.Px()-p4l.Px(), p4v.Py()+p4N.Py()-p4l.Py(), p4v.Pz()+p4N.Pz()-p4l.Pz(), p4v.E()+p4N.E()-p4l.E() );

  LOG("LeptoHad", pDEBUG) << "v [LAB']: " << p4v.E() << " // " << p4v.M2() << " // [ " << p4v.Dx() << " , " << p4v.Dy() << " , " << p4v.Dz() << " ]";
  LOG("LeptoHad", pDEBUG) << "N [LAB']: " << p4N.E() << " // " << p4N.M2() << " // [ " << p4N.Dx() << " , " << p4N.Dy() << " , " << p4N.Dz() << " ]";
  LOG("LeptoHad", pDEBUG) << "l [LAB']: " << p4l.E() << " // " << p4l.M2() << " // [ " << p4l.Dx() << " , " << p4l.Dy() << " , " << p4l.Dz() << " ]";
  LOG("LeptoHad", pDEBUG) << "H [LAB']: " << p4Hadlong.E() << " // " << p4Hadlong.M2() << " // [ " << p4Hadlong.Dx() << " , " << p4Hadlong.Dy() << " , " << p4Hadlong.Dz() << " ]";

  // Translate from long double to double
  const TLorentzVector & vtx = *( event->Probe()->X4());
  TLorentzVector p4Had( (double)p4Hadlong.Px(), (double)p4Hadlong.Py(), (double)p4Hadlong.Pz(), (double)p4Hadlong.E() );
  event->AddParticle(kPdgHadronicSyst, kIStDISPreFragmHadronicState, event->HitNucleonPosition(),-1,-1,-1, p4Had, vtx);

  const Interaction * interaction = event->Summary();
  interaction->KinePtr()->SetHadSystP4(p4Had);
  interaction->KinePtr()->SetW(p4Hadlong.M());

  double W = interaction->Kine().W();
  if(W < fWmin) {
    LOG("LeptoHad", pWARN) << "Low invariant mass, W = " << W << " GeV!!";
    return false;
  }

  const XclsTag &      xclstag    = interaction->ExclTag();
  const Target &       target     = interaction->InitState().Tgt();

  assert(target.HitQrkIsSet());

  bool isp        = pdg::IsProton(target.HitNucPdg());
  int  hit_quark  = target.HitQrkPdg();
  int  frag_quark = xclstag.FinalQuarkPdg();

  LOG("LeptoHad", pDEBUG) << "Hit nucleon pdgc = " << target.HitNucPdg() << ", W = " << W;
  LOG("LeptoHad", pDEBUG) << "Selected hit quark pdgc = " << hit_quark << " // Fragmentation quark = " << frag_quark;

  RandomGen * rnd = RandomGen::Instance();

  //
  // Generate the hadron combination to input PYTHIA
  //

  //If the hit quark is a d we have these options:
  /* uud(->q)     => uu + q */
  /* uud d(->q)db => uu + q (d valence and db sea annihilates)*/
  /* udd(->q)     => ud + q */
  /* udd d(->q)db => ud + q (d valence and db sea annihilates)*/
  if ( pdg::IsDQuark(hit_quark) ) {
    // choose diquark system depending on proton or neutron
    int diquark = 0;
    if (isp) diquark = kPdgUUDiquarkS1;
    else     diquark = rnd->RndHadro().Rndm()>0.75 ? kPdgUDDiquarkS1 : kPdgUDDiquarkS0;
    // Check that the trasnferred energy is higher than the mass of the produced quarks
    double m_frag    = PDGLibrary::Instance()->Find(frag_quark)->Mass();
    double m_diquark = PDGLibrary::Instance()->Find(diquark)->Mass();
    if( W <= m_frag + m_diquark + fMinESinglet ) {
      LOG("LeptoHad", pWARN) << "Low invariant mass, W = " << W << " GeV! Returning a null list";
      LOG("LeptoHad", pWARN) << "frag_quark = " << frag_quark << "    -> m = " << m_frag;
      LOG("LeptoHad", pWARN) << "diquark    = " << diquark    << " -> m = "    << m_diquark;
      return 0;
    }
    // Input the two particles to PYTHIA back to back in the CM frame
    double e_frag    = (W*W - m_diquark*m_diquark + m_frag*m_frag)/2./W;
    double e_diquark = (W*W + m_diquark*m_diquark - m_frag*m_frag)/2./W;
    fPythia->Py1ent( -1, frag_quark, e_frag, 0., 0. ); //k(1,2) = 2
    // If a top quark is produced we decay it because it does not hadronize
    if ( pdg::IsTQuark(frag_quark) ) {
      int ip = 1;
      pydecy_(&ip);
    }
    fPythia->Py1ent( fPythia->GetN()+1, diquark,  e_diquark, fPythia->GetPARU(1), 0. ); //k(2,2) = 1
  }

  //If the hit quark is a u we have these options:
  /* u(->q)ud     => ud + q */
  /* uud u(->q)ub => ud + q (u valence and ub sea annihilates)*/
  /* u(->q)dd     => dd + q */
  /* udd u(->q)ub => dd + q (u valence and ub sea annihilates)*/
  else if ( pdg::IsUQuark(hit_quark) ) {
    // choose diquark system depending on proton or neutron
    int diquark = 0;
    if (isp) diquark = rnd->RndHadro().Rndm()>0.75 ? kPdgUDDiquarkS1 : kPdgUDDiquarkS0;
    else     diquark = kPdgDDDiquarkS1;
    // Check that the trasnferred energy is higher than the mass of the produced quarks.
    double m_frag    = PDGLibrary::Instance()->Find(frag_quark)->Mass();
    double m_diquark = PDGLibrary::Instance()->Find(diquark)->Mass();
    if( W <= m_frag + m_diquark + fMinESinglet ) {
      LOG("LeptoHad", pWARN) << "Low invariant mass, W = " << W << " GeV! Returning a null list";
      LOG("LeptoHad", pWARN) << "frag_quark = " << frag_quark << "    -> m = " << m_frag;
      LOG("LeptoHad", pWARN) << "diquark    = " << diquark    << " -> m = "    << m_diquark;
      return 0;
    }
    // Input the two particles to PYTHIA back to back in the CM frame
    double e_frag    = (W*W - m_diquark*m_diquark + m_frag*m_frag)/2./W;
    double e_diquark = (W*W + m_diquark*m_diquark - m_frag*m_frag)/2./W;
    fPythia->Py1ent( -1, frag_quark, e_frag, 0., 0. ); //k(1,2) = 2
    fPythia->Py1ent( fPythia->GetN()+1, diquark, e_diquark, fPythia->GetPARU(1), 0. ); //k(2,2) = 1
  }


  // If the hit quark is not u or d then is more complicated.
  // We are using the same procedure use in LEPTO (see lqev.F)
  // Our initial systemt will look like this          ->  qqq + hit_q(->frag_q) + rema_q
  // And we have to input PYTHIA something like this  ->  frag_q + rema  + hadron
  // These are the posible combinations               ->  frag_q[q] + meson [qqb]  + diquark [qq]
  //                                                  ->  frag_q[qb] + baryon [qqq] + quark [q]
  else {

    // Remnant of the hit quark (which is from the sea) will be of opposite charge
    int rema_hit_quark = -hit_quark;

    // Check that the trasnfered energy is higher than the mass of the produce quarks plus remnant quark and nucleon
    double m_frag     = PDGLibrary::Instance()->Find(frag_quark)->Mass();
    double m_rema_hit = PDGLibrary::Instance()->Find(rema_hit_quark)->Mass();
    if (W <= m_frag + m_rema_hit + 0.9 + fMinESinglet ) {
      LOG("LeptoHad", pWARN) << "Low invariant mass, W = " << W << " GeV! Returning a null list";
      LOG("LeptoHad", pWARN) << " frag_quark     = " << frag_quark     << " -> m = " << m_frag;
      LOG("LeptoHad", pWARN) << " rema_hit_quark = " << rema_hit_quark << " -> m = " << m_rema_hit;
      return 0;
    }

    //PDG of the two hadronic particles for the final state
    int hadron = 0;
    int rema   = 0;

    int ntwoq = isp ? 2 : 1; //proton two ups & neutron one up
    int counter = 0;

    // Here we select the id and kinematics of the hadron and rema particles
    // Some combinations can be kinematically forbiden so we repeat this process
    // up to 100 times before the event is discarded.
    while( counter<fMaxIterHad ) {

      // Loop to create a combination of hadron + rema. Two options are possible:
      // 1) diquark [qq] + meson [qqb]
      // 2) quark [q] + baryon [qqq]
      while(hadron==0) {
        //choose a valence quark and the remaining will be a diquark system
        int valquark = int(1.+ntwoq/3.+rnd->RndHadro().Rndm());
        int diquark  = 0;
        if ( valquark==ntwoq ) diquark = rnd->RndHadro().Rndm()>0.75 ? kPdgUDDiquarkS1 : kPdgUDDiquarkS0;
        else                   diquark = 1000*ntwoq+100*ntwoq+3;

        // Choose flavours using PYTHIA tool
        int idum;
        if ( rema_hit_quark>0 ) { //create a baryon (qqq)
          pykfdi_(&diquark,&rema_hit_quark,&idum,&hadron);
          rema = valquark;
        }
        else {                    //create a meson (qqbar)
          pykfdi_(&valquark,&rema_hit_quark,&idum,&hadron);
          rema = diquark;
        }
      }

      double m_hadron = PDGLibrary::Instance()->Find(hadron)->Mass();
      double m_rema   = PDGLibrary::Instance()->Find(rema)->Mass();

      // Give balancing pT to hadron and rema particles
      double pT  = fRemnantPT * TMath::Sqrt( -1*TMath::Log( rnd->RndHadro().Rndm() ) );
      double pT2 = TMath::Power(pT,2);
      double pr  = TMath::Power(m_hadron,2)+pT2;
      //to generate the longitudinal scaling variable z in jet fragmentation using PYTHIA function
      // Split energy-momentum of remnant using PYTHIA function
      // z=E-pz fraction for rema forming jet-system with frag_q
      // 1-z=E-pz fraction for hadron
      double z;
      int kfl1 = 1;
      int kfl3 = 0;
      pyzdis_(&kfl1,&kfl3,&pr,&z);

      // Energy of trasnfered to the hadron
      double tm_hadron = pr / z / W;
      double E_hadron   = 0.5 * ( z*W + tm_hadron );  //E_hadron - pz = zW
      double E_pz       = W - tm_hadron;
      double WT         = (1-z) * W * E_pz - pT2;

      // Check if energy in jet system is enough for fragmentation.
      if ( WT > TMath::Power(m_frag+m_rema+fMinESinglet,2) ) {

        // Energy of transfered to the fragmented quark and rema system
        // Applying energy conservation
        WT = TMath::Sqrt( WT + pT2 );
        double tm_rema   = TMath::Power(m_rema,2) + pT2;
        double E_frag    = 0.5 * ( WT + ( TMath::Power(m_frag,2) - tm_rema)/WT ); //E_frag + E_rema = WT
        double E_rema    = 0.5 * ( WT + (-TMath::Power(m_frag,2) + tm_rema)/WT );
        double x_rema    = -1 * TMath::Sqrt( TMath::Power(E_rema,2) - tm_rema );
        double theta_rema;
        theta_rema = pyangl_(&x_rema,&pT);

        // Select a phi angle between between particles randomly
        double phi = 2*kPi*rnd->RndHadro().Rndm();

        double dbez = (E_pz-(1-z)*W)/(E_pz+(1-z)*W);
        double pz_hadron  = -0.5 * ( z*W - tm_hadron );

        // Input the three particles to PYTHIA in the CM frame
        // If a top quark is produced we decay it because it does not hadronize

        fPythia->Py1ent( -1, frag_quark, E_frag, 0.,         0. );           //k(1,2) = 2
        if (TMath::Abs(frag_quark) > 5 ) {
          int ip = 1;
          pydecy_(&ip);
        }
        fPythia->Py1ent( fPythia->GetN()+1, rema, E_rema, theta_rema, phi ); //k(2,2) = 1

        int imin     = 0;
        int imax     = 0;
        double the  = 0.; double ph   = 0.;
        double dbex = 0.; double dbey = 0.; 
        pyrobo_( &imin , &imax, &the, &ph, &dbex, &dbey , &dbez );
        double theta_hadron = pyangl_(&pz_hadron,&pT);

        fPythia->SetMSTU( 10, 1 ); //keep the mass value stored in P(I,5), whatever it is.
        fPythia->SetP( fPythia->GetN()+1, 5, m_hadron );
        fPythia->Py1ent( fPythia->GetN()+1, hadron, E_hadron, theta_hadron, phi + kPi );
        fPythia->SetMSTU( 10, 2 ); //find masses according to mass tables as usual.

        // Target remnants required to go backwards in hadronic cms
        if ( fPythia->GetP(fPythia->GetN()-1,3)<0 && fPythia->GetP(fPythia->GetN(),3)<0 ) break; //quit the while from line 368

        LOG("LeptoHad", pINFO) << "Not backward hadron or rema";
        LOG("LeptoHad", pINFO) << "hadron     = " << hadron     << " -> Pz = " << fPythia->GetP(fPythia->GetN(),3) ;
        LOG("LeptoHad", pINFO) << "rema = " << rema << " -> Pz = " << fPythia->GetP(fPythia->GetN()-1,3) ;
        
      }
      else {
        LOG("LeptoHad", pINFO) << "Low WT value ... ";
        LOG("LeptoHad", pINFO) << "WT = " << TMath::Sqrt(WT) << " // m_frag = " << m_frag << " // m_rema = " << m_rema;
      }

      LOG("LeptoHad", pINFO) << "Hadronization paricles not suitable. Trying again... " << counter;
      counter++;
      if (counter==100) {
        LOG("LeptoHad", pWARN) << "Hadronization particles failed after " << counter << " iterations! Returning a null list";
        return 0;
      }

    }

  }

  // Introduce a primordial kT system
  double pT  = fPrimordialKT * TMath::Sqrt( -1*TMath::Log( rnd->RndHadro().Rndm() ) );
  double phi   = -2*kPi*rnd->RndHadro().Rndm();
  double theta = 0.;

  int imin     = 0;
  int imax     = 0;
  double dbex = 0.; double dbey = 0.; double dbez = 0;
  pyrobo_( &imin , &imax, &theta, &phi, &dbex, &dbey , &dbez );
  phi   = -1 * phi;
  theta = TMath::ATan(2.*pT/W);
  pyrobo_( &imin , &imax, &theta, &phi, &dbex, &dbey , &dbez );

  // Run PYTHIA with the input particles
  fPythia->Pyexec();
  // Use for debugging purposes
  //fPythia->Pylist(3);
  fPythia->GetPrimaries();
  TClonesArray * pythia_particles = (TClonesArray *) fPythia->ImportParticles("All");
  // copy PYTHIA container to a new TClonesArray so as to transfer ownership
  // of the container and of its elements to the calling method
  int np = pythia_particles->GetEntries();
  assert(np>0);
  TClonesArray * particle_list = new TClonesArray("genie::GHepParticle", np);
  particle_list->SetOwner(true);

  // Boost velocity HCM -> LAB
  long double beta = p4Hadlong.P()/p4Hadlong.E();

  //fix numbering for events with top
  bool isTop = false;

  //-- Translate the fragmentation products from TMCParticles to
  //   GHepParticles and copy them to the event record.
  int mom = event->FinalStateHadronicSystemPosition();
  assert(mom!=-1);

  TMCParticle * p = 0;
  TIter particle_iter(pythia_particles);
  while( (p = (TMCParticle *) particle_iter.Next()) ) {

    int pdgc = p->GetKF();
    int ks   = p->GetKS();

    // Final state particles can not be quarks or diquarks but colorless
    if(ks == 1) {
      if( pdg::IsQuark(pdgc) || pdg::IsDiQuark(pdgc) ) {
        LOG("LeptoHad", pERROR) << "Hadronization failed! Bare quark/di-quarks appear in final state!";
        return false;
      }
    }

    // When top quark is produced, it is immidiately decay before hadronization. Then the decayed
    // products are hadronized with the hadron remnants. Therefore, we remove the top quark from
    // the list of particles so that the mother/daugher assigments is at the same level for decayed
    // products and hadron remnants.
    if ( pdg::IsTQuark( TMath::Abs(pdgc) ) ) { isTop=true; continue; }

    // fix numbering scheme used for mother/daughter assignments
    if ( isTop ) {
      (p->GetParent()==0) ? p->SetParent(p->GetParent() - 1) : p->SetParent(p->GetParent() - 2);
      p->SetFirstChild (p->GetFirstChild() - 2);
      p->SetLastChild  (p->GetLastChild()  - 2);
    }
    else  {
      p->SetParent(p->GetParent() - 1);
      p->SetFirstChild (p->GetFirstChild() - 1);
      p->SetLastChild  (p->GetLastChild()  - 1);
    }

    LongLorentzVector p4long( p->GetPx(), p->GetPy(), p->GetPz(), p->GetEnergy()  );
    p4long.BoostZ(beta);
    p4long.Rotate(p4Hadlong);

    // Translate from long double to double
    TLorentzVector p4( (double)p4long.Px(), (double)p4long.Py(), (double)p4long.Pz(), (double)p4long.E() );

    // Somtimes PYTHIA output particles with E smaller than its mass. This is wrong,
    // so we assume that the are at rest.
    double massPDG = PDGLibrary::Instance()->Find(pdgc)->Mass();
    if ( (ks==1 || ks==4) && p4.E()<massPDG ) {
      LOG("LeptoHad", pINFO) << "Putting at rest one stable particle generated by PYTHIA because E < m";
      LOG("LeptoHad", pINFO) << "PDG = " << pdgc << " // State = " << ks;
      LOG("LeptoHad", pINFO) << "E = " << p4.E() << " // |p| = " << p4.P();
      LOG("LeptoHad", pINFO) << "p = [ " << p4.Px() << " , "  << p4.Py() << " , "  << p4.Pz() << " ]";
      LOG("LeptoHad", pINFO) << "m    = " << p4.M() << " // mpdg = " << massPDG;
      p4.SetXYZT(0,0,0,massPDG);
    }

    // copy final state particles to the event record
    GHepStatus_t ist = (ks==1 || ks==4) ? kIStStableFinalState : kIStDISPreFragmHadronicState;

    int im  = mom + 1 + p->GetParent();
    int ifc = (p->GetFirstChild() <= -1) ? -1 : mom + 1 + p->GetFirstChild();
    int ilc = (p->GetLastChild()  <= -1) ? -1 : mom + 1 + p->GetLastChild();

    double vx = vtx.X() + p->GetVx()*1e12; //pythia gives position in [mm] while genie uses [fm]
    double vy = vtx.Y() + p->GetVy()*1e12;
    double vz = vtx.Z() + p->GetVz()*1e12;
    double vt = vtx.T() + p->GetTime()*(units::millimeter/units::second);
    TLorentzVector pos( vx, vy, vz, vt );

    event->AddParticle( pdgc, ist, im,-1, ifc, ilc, p4, pos );

  }

  return true;
#else
  return false;
#endif

}
//____________________________________________________________________________
void LeptoHadPythia6::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void LeptoHadPythia6::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void LeptoHadPythia6::LoadConfig(void)
{

#ifdef __GENIE_PYTHIA6_ENABLED__
  GetParam("MaxIter-Had", fMaxIterHad ) ;

  // Width of Gaussian distribution for transverse momentums
  // Define in LEPTO with PARL(3) and PARL(14)
  GetParam("Primordial-kT", fPrimordialKT ) ;
  GetParam("Remnant-pT",    fRemnantPT ) ;

  // It is, with quark masses added, used to define the minimum allowable energy of a colour-singlet parton system.
  GetParam( "Energy-Singlet", fMinESinglet ) ;

  GetParam( "Xsec-Wmin", fWmin ) ;

  GetParam( "SSBarSuppression",       fSSBarSuppression       );
  GetParam( "GaussianPt2",            fGaussianPt2            );
  GetParam( "NonGaussianPt2Tail",     fNonGaussianPt2Tail     );
  GetParam( "RemainingEnergyCutoff",  fRemainingECutoff       );
  GetParam( "DiQuarkSuppression",     fDiQuarkSuppression     );
  GetParam( "LightVMesonSuppression", fLightVMesonSuppression );
  GetParam( "SVMesonSuppression",     fSVMesonSuppression     );
  GetParam( "Lunda",                  fLunda                  );
  GetParam( "Lundb",                  fLundb                  );
  GetParam( "LundaDiq",               fLundaDiq               );

  int warnings;       GetParam( "Warnings",      warnings ) ;
  int errors;         GetParam( "Errors",        errors ) ;
  int qrk_mass;       GetParam( "QuarkMass",     qrk_mass ) ;

  // PYTHIA6 parameters only valid for HEDIS
  fPythia->SetPARP(2,  fWmin);     // (D = 10. GeV) lowest c.m. energy for the event as a whole that the program will accept to simulate. (bellow 2GeV pythia crashes)
  fPythia->SetMSTU(26, warnings);  // (Default=10) maximum number of warnings that are printed
  fPythia->SetMSTU(22, errors);    // (Default=10) maximum number of errors that are printed
  fPythia->SetMSTJ(93, qrk_mass);  // light (d, u, s, c, b) quark masses are taken from PARF(101) - PARF(105) rather than PMAS(1,1) - PMAS(5,1). Diquark masses are given as sum of quark masses, without spin splitting term.
  fPythia->SetPMAS(24,1,kMw);      // mass of the W boson (pythia=80.450 // genie=80.385)
  fPythia->SetPMAS(24,2,0.);       // set to 0 the width of the W boson to avoid problems with energy conservation
  fPythia->SetPMAS(6,2,0.);        // set to 0 the width of the top to avoid problems with energy conservation
  fPythia->SetMDME(192,1,0);   // W->dbar+t decay off 
  fPythia->SetMDME(196,1,0);   // W->cbar+t decay off 
  fPythia->SetMDME(200,1,0);   // W->cbar+t decay off 
  
  // PYTHIA tuned parameters
  fPythia->SetPARJ(2,  fSSBarSuppression       );
  fPythia->SetPARJ(21, fGaussianPt2            );
  fPythia->SetPARJ(23, fNonGaussianPt2Tail     );
  fPythia->SetPARJ(33, fRemainingECutoff       );
  fPythia->SetPARJ(1,  fDiQuarkSuppression     );
  fPythia->SetPARJ(11, fLightVMesonSuppression );
  fPythia->SetPARJ(12, fSVMesonSuppression     );
  fPythia->SetPARJ(41, fLunda                  );
  fPythia->SetPARJ(42, fLundb                  );
  fPythia->SetPARJ(45, fLundaDiq               );
#endif

}
//____________________________________________________________________________
void LeptoHadPythia6::Initialize(void) const
{
#ifdef __GENIE_PYTHIA6_ENABLED__
  fPythia = TPythia6::Instance();

  // sync GENIE/PYTHIA6 seed number
  RandomGen::Instance();
#endif

}
