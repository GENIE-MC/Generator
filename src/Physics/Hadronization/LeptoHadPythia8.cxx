//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Alfonso Garcia <aagarciasoto \at km3net.de>
 IFIC (Valencia)
*/
//____________________________________________________________________________

#include "Physics/Hadronization/LeptoHadPythia8.h"

int LeptoHadPythia8::getMeson(int q, int qb, double rnd) const {
    // options: q -> d=1 / u=1
    // options: qb -> sb=-3 / cb=-4 / bb=-5
    int aqb = abs(qb); 
    int MultipletCode = ( mesonRateSum[aqb-3]*rnd>1 ) ? 1: 3;
    int idMeson = 100*aqb + 10*q + MultipletCode;
    if (qb==-4) return -1*idMeson;
    else        return    idMeson;
}

int LeptoHadPythia8::getBaryon(int qq, int q, double rnd) const {
    // options: qq -> dd1=1103 / ud0=2101 / ud1=2103 / uu1=2203 
    // options: q -> d=1 / u=2 / s=3 / c=4 / b=5
    int id1 = qq / 1000;
    int id2 = (qq / 100) % 10;
    int id3 = qq % 10;
    int spin = id3 - 1;
    if (spin==2 && id1!=id2) spin = 4;
    if (q!=id1 && q!=id2) spin++;
    int o1 = TMath::Max( q, TMath::Max( id1, id2) );
    int o3 = TMath::Min( q, TMath::Min( id1, id2) );
    int o2 = q + id1 + id2 - o1 - o3;
    int spinBar = ( CGSum[spin]*rnd<CGOct[spin] ) ? 2 : 4;
    if ( spinBar==2 && o1>o2 && o2>o3 && id3==1 ) return 1000*o1 + 100*o3 + 10*o2 + spinBar; //lambdas
    else                                          return 1000*o1 + 100*o2 + 10*o3 + spinBar;
}

double LeptoHadPythia8::getRandomZ( double a, double b) const {
  // fragmentation function f(x) = ((1-x)^a*exp(-b/x))/x
  // not optimal if a->0 or a->1
  double zpeak = (b+1.-sqrt(pow(b-1.,2)+4.*a*b))/(1.-a)/2.;
  if ( zpeak>0.9999 && b>100. ) zpeak = TMath::Min(zpeak, 1.-a/b);

  bool closeto0 = ( zpeak<0.1);
  bool closeto1 = ( zpeak>0.85 && b>1. );

  double flw = 1.;
  double fup = 1.;
  double frn = 2.;
  double wdt = 0.5;

  if ( closeto0 ) {
    wdt = 2.75*zpeak;
    flw = wdt;
    fup = -wdt*log(wdt);
    frn = flw+fup;
  } 
  else if ( closeto1 ) {
    double rcb = sqrt(4.+pow(1./b,2));
    wdt = rcb - 1./zpeak - (1./b)*log(zpeak*(rcb+1./b)/2. );
    wdt += (a/b)*log(1.-zpeak);
    wdt = TMath::Min(zpeak,TMath::Max(0.,wdt));
    flw = 1./b;
    fup = 1.-wdt;
    frn = flw+fup;
  }

  RandomGen * rnd = RandomGen::Instance();

  double z,prel,val;  
  while(1){
    z = rnd->RndHadro().Rndm();
    prel = 1.;
    if (closeto0) {
      if ( frn*rnd->RndHadro().Rndm()<flw ) z = wdt*z;
      else { z = pow(wdt,z); prel = wdt/z; }
    } else if (closeto1) {
      if ( frn*rnd->RndHadro().Rndm()<flw) { z = wdt+log(z)/b; prel = exp(b*(z-wdt)); } 
      else z = wdt+(1.-wdt)*z;
    }
    if ( z>0 && z<1 ) {
      double pw = b*(1./zpeak-1./z)+ log(zpeak/z) + a*log((1.-z)/(1.-zpeak));
      val = exp(TMath::Max(-50.,TMath::Min(50.,pw)));
    } 
    else val = 0.;
    if ( val >= rnd->RndHadro().Rndm() * prel) break;
  }

  return z;
}


//____________________________________________________________________________
LeptoHadPythia8::LeptoHadPythia8() :
EventRecordVisitorI("genie::LeptoHadPythia8")
{
  this->Initialize();
}
//____________________________________________________________________________
LeptoHadPythia8::LeptoHadPythia8(string config) :
EventRecordVisitorI("genie::LeptoHadPythia8", config)
{
  this->Initialize();
}
//____________________________________________________________________________
LeptoHadPythia8::~LeptoHadPythia8()
{

}
//____________________________________________________________________________
void LeptoHadPythia8::ProcessEventRecord(GHepRecord * event) const
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
bool LeptoHadPythia8::Hadronize(GHepRecord * 
#ifdef __GENIE_PYTHIA8_ENABLED__
  event // avoid unused variable warning if PYTHIA8 is not enabled
#endif
) const
{

#ifdef __GENIE_PYTHIA8_ENABLED__
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
  fPythia->event.reset();

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
    double pz_cm = Pythia8::sqrtpos( e_frag*e_frag - m_frag*m_frag );
    fPythia->event.append(frag_quark, 23, 101, 0, 0., 0.,  pz_cm, e_frag, m_frag);
    fPythia->event.append(diquark,    23, 0, 101, 0., 0., -pz_cm, e_diquark, m_diquark);
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
    double pz_cm = Pythia8::sqrtpos( e_frag*e_frag - m_frag*m_frag );
    fPythia->event.append(frag_quark, 23, 101, 0, 0., 0.,  pz_cm, e_frag, m_frag);
    fPythia->event.append(diquark,    23, 0, 101, 0., 0., -pz_cm, e_diquark, m_diquark);
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
        if ( rema_hit_quark>0 ) { //create a baryon (qqq)
          hadron = getBaryon(diquark,rema_hit_quark,rnd->RndHadro().Rndm());
          rema = valquark;
        }
        else {                    //create a meson (qqbar)
          hadron = getMeson(valquark,rema_hit_quark,rnd->RndHadro().Rndm());
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
      double z = getRandomZ(Afrag,Bfrag*pr);

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
        theta_rema = TMath::ATan2(pT,x_rema);

        // Select a phi angle between between particles randomly
        double phi = 2*kPi*rnd->RndHadro().Rndm();

        double dbez = (E_pz-(1-z)*W)/(E_pz+(1-z)*W);
        double pz_hadron  = -0.5 * ( z*W - tm_hadron );

        // Input the three particles to PYTHIA in the CM frame
        // If a top quark is produced we decay it because it does not hadronize
        fPythia->event.append(frag_quark, 23, 101, 0, 0., 0., sqrt(E_frag*E_frag-m_frag*m_frag), E_frag, m_frag);

        double p_rema  = sqrt(E_rema*E_rema-m_rema*m_rema);
        fPythia->event.append(rema, 23, 0, 101, p_rema*sin(theta_rema)*sin(phi), p_rema*sin(theta_rema)*cos(phi), p_rema*cos(theta_rema), E_rema, m_rema);

        fPythia->event.bst(0,0,dbez);

        double theta_hadron = TMath::ATan2(pT,pz_hadron);

        double p_hadron = sqrt(E_hadron*E_hadron-m_hadron*m_hadron);
        fPythia->event.append(hadron, 23, 0, 0, p_hadron*sin(theta_hadron)*sin(phi+kPi), p_hadron*sin(theta_hadron)*cos(phi+kPi), p_hadron*cos(theta_hadron), E_hadron, m_hadron);

        // Target remnants required to go backwards in hadronic cms
        int nsize = fPythia->event.size();
        if ( fPythia->event[nsize-1].pz()<0 && fPythia->event[nsize-2].pz()<0 ) break; //quit the while( counter<fMaxIterHad )

        // break;

        LOG("LeptoHad", pINFO) << "Not backward hadron or rema";
        LOG("LeptoHad", pINFO) << "hadron     = " << hadron     << " -> Pz = " << fPythia->event[nsize-1].pz() ;
        LOG("LeptoHad", pINFO) << "rema = " << rema << " -> Pz = " << fPythia->event[nsize-2].pz() ;
        fPythia->event.reset();
        
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

  fPythia->event.rot(theta,phi);
  phi   = -1 * phi;
  theta = TMath::ATan(2.*pT/W);
  fPythia->event.rot(theta,phi);

  // fPythia->event.list();
  // fPythia->stat();

  // Run PYTHIA with the input particles
  fPythia->next();
  // fPythia->event.list();
  // fPythia->stat();
  Pythia8::Event &fEvent = fPythia->event;
  int np = fEvent.size();
  assert(np>0);

  // Boost velocity HCM -> LAB
  long double beta = p4Hadlong.P()/p4Hadlong.E();

  //fix numbering for events with top
  bool isTop = false;

  //-- Translate the fragmentation products from TMCParticles to
  //   GHepParticles and copy them to the event record.
  int mom = event->FinalStateHadronicSystemPosition();
  assert(mom!=-1);

  for (int i = 1; i < np; i++) { // ignore first entry -> (system) pseudoparticle

    int pdgc = fEvent[i].id();
    if (!PDGLibrary::Instance()->Find(pdgc)) continue; // some intermediate particles not part of genie tables

    int ks   = fEvent[i].status();

    // Final state particles can not be quarks or diquarks but colorless
    if(ks > 0) {
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
      (fEvent[i].mother1()==0) ? fEvent[i].mothers(fEvent[i].mother1()-1,fEvent[i].mother2()) : fEvent[i].mothers(fEvent[i].mother1()-2,fEvent[i].mother2());
      fEvent[i].daughters(fEvent[i].daughter1()-2,fEvent[i].daughter2()-2);
    }
    else  {
      fEvent[i].mothers(fEvent[i].mother1()-1,fEvent[i].mother2());
      fEvent[i].daughters(fEvent[i].daughter1()-1,fEvent[i].daughter2()-1);
    }

    LongLorentzVector p4long( fEvent[i].px(), fEvent[i].py(), fEvent[i].pz(), fEvent[i].e()  );
    p4long.BoostZ(beta);
    p4long.Rotate(p4Hadlong);

    // Translate from long double to double
    TLorentzVector p4( (double)p4long.Px(), (double)p4long.Py(), (double)p4long.Pz(), (double)p4long.E() );

    // Somtimes PYTHIA output particles with E smaller than its mass. This is wrong,
    // so we assume that the are at rest.
    double massPDG = PDGLibrary::Instance()->Find(pdgc)->Mass();
    if ( ks>0 && p4.E()<massPDG ) {
      LOG("LeptoHad", pINFO) << "Putting at rest one stable particle generated by PYTHIA because E < m";
      LOG("LeptoHad", pINFO) << "PDG = " << pdgc << " // State = " << ks;
      LOG("LeptoHad", pINFO) << "E = " << p4.E() << " // |p| = " << p4.P();
      LOG("LeptoHad", pINFO) << "p = [ " << p4.Px() << " , "  << p4.Py() << " , "  << p4.Pz() << " ]";
      LOG("LeptoHad", pINFO) << "m    = " << p4.M() << " // mpdg = " << massPDG;
      p4.SetXYZT(0,0,0,massPDG);
    }

    // copy final state particles to the event record
    GHepStatus_t ist = (ks>0) ? kIStStableFinalState : kIStDISPreFragmHadronicState;

    int im  = mom + 1 + fEvent[i].mother1();
    int ifc = (fEvent[i].daughter1() <= -1) ? -1 : mom + 1 + fEvent[i].daughter1();
    int ilc = (fEvent[i].daughter2()  <= -1) ? -1 : mom + 1 + fEvent[i].daughter2();

    double vx = vtx.X() + fEvent[i].xProd()*1e12; //pythia gives position in [mm] while genie uses [fm]
    double vy = vtx.Y() + fEvent[i].yProd()*1e12;
    double vz = vtx.Z() + fEvent[i].zProd()*1e12;
    double vt = vtx.T() + fEvent[i].tProd()*(units::millimeter/units::second);
    TLorentzVector pos( vx, vy, vz, vt );

    event->AddParticle( pdgc, ist, im,-1, ifc, ilc, p4, pos );

  }

  return true;
#else
  return false;
#endif

}
//____________________________________________________________________________
void LeptoHadPythia8::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void LeptoHadPythia8::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void LeptoHadPythia8::LoadConfig(void)
{

#ifdef __GENIE_PYTHIA8_ENABLED__
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

  fPythia->settings.parm("StringFlav:probStoUD",         fSSBarSuppression);
  fPythia->settings.parm("Diffraction:primKTwidth",      fGaussianPt2);
  fPythia->settings.parm("StringPT:enhancedFraction",    fNonGaussianPt2Tail);
  fPythia->settings.parm("StringFragmentation:stopMass", fRemainingECutoff);
  fPythia->settings.parm("StringFlav:probQQtoQ",         fDiQuarkSuppression);
  fPythia->settings.parm("StringFlav:mesonUDvector",     fLightVMesonSuppression);
  fPythia->settings.parm("StringFlav:mesonSvector",      fSVMesonSuppression);
  fPythia->settings.parm("StringZ:aLund",                fLunda);
  fPythia->settings.parm("StringZ:bLund",                fLundb);
  fPythia->settings.parm("StringZ:aExtraDiquark",        fLundaDiq);

  // Same default mass of the W boson in pythia8 and genie, so no need to change in pythia8
  // No problem with energy conservation W and top decays, so no need to set the width to 0

  Afrag = fPythia->settings.parm("StringZ:aLund");
  Bfrag = fPythia->settings.parm("StringZ:bLund");

  bool isAvalid = true;
  if (Afrag<0.02 || abs(Afrag-1)==0.01 ) isAvalid = false;
  if (!isAvalid) {
    LOG("LeptoHad", pFATAL) << "Invalid A factor for fragmenation function" ;
    LOG("LeptoHad", pFATAL) << "A must be (>0.02 & <0.99) | >1.01" ;
    exit(1);
  }

  mesonRateSum[0] = 1. + fPythia->settings.parm("StringFlav:mesonSvector"); //0.55 
  mesonRateSum[1] = 1. + fPythia->settings.parm("StringFlav:mesonCvector"); //0.88 
  mesonRateSum[2] = 1. + fPythia->settings.parm("StringFlav:mesonBvector"); //2.20 

  double decupletSup = fPythia->settings.parm("StringFlav:decupletSup");;
  for (int i = 0; i < 6; ++i) CGSum[i] = CGOct[i] + decupletSup*CGDec[i];

  LOG("LeptoHad", pINFO) << "Initialising PYTHIA..." ;
  fPythia->init(); 
#endif

}
//____________________________________________________________________________
void LeptoHadPythia8::Initialize(void) const
{

#ifdef __GENIE_PYTHIA8_ENABLED__
  fPythia = Pythia8Singleton::Instance()->Pythia8();

  fPythia->readString("Print:quiet = on");

  // sync GENIE and PYTHIA8 seeds
  RandomGen * rnd = RandomGen::Instance();
  long int seed = rnd->GetSeed();
  fPythia->readString("Random:setSeed = on");
  fPythia->settings.mode("Random:seed", seed);
  LOG("LeptoHad", pINFO) << "PYTHIA8  seed = " << fPythia->settings.mode("Random:seed");

  //needed to only do hadronization
  fPythia->readString("ProcessLevel:all = off");
#endif

}

