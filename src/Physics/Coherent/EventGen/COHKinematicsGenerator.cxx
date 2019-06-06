//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 09, 2009 - CA
   Moved into the new Coherent package from its previous location  (EVGModules 
   package)
 @ Mar 03, 2009 - CA
   Renamed COHPiKinematicsGenerator -> COHKinematicsGenerator in
   anticipation of reusing the code for simulating coherent production of
   vector mesons.
 @ May 06, 2009 - CA
   Fix a problem with the search for the max cross section over the allowed
   phase space which prevented kinematics to be generated for events near the 
   energy threshold.
 @ Feb 06, 2013 - CA
   When the value of the differential cross-section for the selected kinematics
   is set to the event, set the corresponding KinePhaseSpace_t value too.

*/
//____________________________________________________________________________

#include <cstdlib>

#include <TROOT.h>
#include <TMath.h>
#include <TF2.h>
#include "Math/Minimizer.h"
#include "Math/Factory.h"

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Physics/Coherent/EventGen/COHKinematicsGenerator.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/RunningThreadInfo.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;
using namespace genie::utils;

//___________________________________________________________________________
COHKinematicsGenerator::COHKinematicsGenerator() :
  KineGeneratorWithCache("genie::COHKinematicsGenerator")
{
  fEnvelope = 0;
}
//___________________________________________________________________________
COHKinematicsGenerator::COHKinematicsGenerator(string config) :
  KineGeneratorWithCache("genie::COHKinematicsGenerator", config)
{
  fEnvelope = 0;
}
//___________________________________________________________________________
COHKinematicsGenerator::~COHKinematicsGenerator()
{
  if(fEnvelope) delete fEnvelope;
}
//___________________________________________________________________________
void COHKinematicsGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
  if(fGenerateUniformly) {
    LOG("COHKinematics", pNOTICE)
      << "Generating kinematics uniformly over the allowed phase space";
  }

  //-- Access cross section algorithm for running thread
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  fXSecModel = evg->CrossSectionAlg();
  if (fXSecModel->Id().Name() == "genie::ReinSehgalCOHPiPXSec") {
    CalculateKin_ReinSehgal(evrec);
  } else if (fXSecModel->Id().Name() == "genie::BergerSehgalCOHPiPXSec2015") {
    CalculateKin_BergerSehgal(evrec);
  } else if (fXSecModel->Id().Name() == "genie::BergerSehgalFMCOHPiPXSec2015") {
    CalculateKin_BergerSehgalFM(evrec);
  } else if ((fXSecModel->Id().Name() == "genie::AlvarezRusoCOHPiPXSec")) {
    CalculateKin_AlvarezRuso(evrec);
  }
  else {
    LOG("COHKinematicsGenerator",pFATAL) <<
      "ProcessEventRecord >> Cannot calculate kinematics for " <<
      fXSecModel->Id().Name();
  }
}
//___________________________________________________________________________
void COHKinematicsGenerator::CalculateKin_BergerSehgal(GHepRecord * evrec) const
{
  // Get the Primary Interacton object
  Interaction * interaction = evrec->Summary();
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);    // TODO: Why turn this off?

  // Initialise a random number generator 
  RandomGen * rnd = RandomGen::Instance();

  //-- For the subsequent kinematic selection with the rejection method:
  //   Calculate the max differential cross section or retrieve it from the
  //   cache. Throw an exception and quit the evg thread if a non-positive
  //   value is found.
  //
  //   TODO: We are not offering the "fGenerateUniformly" option here.
  double xsec_max = this->MaxXSec(evrec);

  //-- Get the kinematical limits for the generated x,y
  const KPhaseSpace & kps = interaction->PhaseSpace();
  Range1D_t y  = kps.YLim();
  Range1D_t Q2(fQ2Min,fQ2Max);
  assert(y.min>0. && y.max>0. && y.min<1. && y.max<1. && y.min<y.max);

  const double ymin  = y.min + kASmallNum;
  const double ymax  = y.max - kASmallNum;
  const double dy    = ymax - ymin;
  const double Q2min = Q2.min + kASmallNum;
  const double Q2max = Q2.max - kASmallNum; 
  const double dQ2   = Q2max - Q2min;

  unsigned int iter = 0;
  bool accept=false;
  double xsec=-1, gy=-1, gQ2=-1;

  while(1) {
    iter++;
    if(iter > kRjMaxIterations) this->throwOnTooManyIterations(iter,evrec);

    //-- Select unweighted kinematics using importance sampling method. 
    // TODO: The importance sampling envelope is not used. Currently, 
    // we just employ a standard rejection-method approach.

    gy  = ymin  + dy  * rnd->RndKine().Rndm(); 
    gQ2 = Q2min + dQ2 * rnd->RndKine().Rndm(); 

    LOG("COHKinematics", pINFO) << 
      "Trying: Q^2 = " << gQ2 << ", y = " << gy; /* << ", t = " << gt; */

    interaction->KinePtr()->Sety(gy);
    interaction->KinePtr()->SetQ2(gQ2);
    kinematics::UpdateXFromQ2Y(interaction);

    // computing cross section for the current kinematics
    xsec = fXSecModel->XSec(interaction, kPSQ2yfE);

    //-- decide whether to accept the current kinematics
    accept = (xsec_max * rnd->RndKine().Rndm() < xsec);

    //-- If the generated kinematics are accepted, finish-up module's job
    if(accept) {
      LOG("COHKinematics", pNOTICE)
        << "Selected: Q^2 = " << gQ2 << ", y = " << gy; /* << ", t = " << gt; */

      // the Berger-Sehgal COH cross section should be a triple differential cross section
      // d^2xsec/dQ2dydt where t is the the square of the 4p transfer to the
      // nucleus. The cross section used for kinematical selection should have
      // the t-dependence integrated out. The t-dependence is of the form
      // ~exp(-bt). Now that the x,y kinematical variables have been selected
      // we can generate a t using the t-dependence as a PDF.
      const InitialState & init_state = interaction->InitState();
      double Ev    = init_state.ProbeE(kRfLab);
      double Epi   = gy*Ev; // pion energy
      double Epi2  = TMath::Power(Epi,2);
      double pme2  = kPionMass2/Epi2;   
      double gx    = interaction->KinePtr()->x(); // utils::kinematics::Q2YtoX(Ev,kNucleonMass,gQ2,gy); // gQ2 / ( 2. * gy * kNucleonMass * Ev); 
      double xME   = kNucleonMass*gx/Epi;
      double tA    = 1. + xME - 0.5*pme2;
      /* Range1D_t tl = kps.TLim();   // this becomes the bounds for t */
      double tB    = TMath::Sqrt(1.+ 2*xME) * TMath::Sqrt(1.-pme2);
      double tmin  = 2*Epi2 * (tA-tB);
      double tmax  = 2*Epi2 * (tA+tB);
      double A     = (double) init_state.Tgt().A(); 
      double A13   = TMath::Power(A,1./3.);
      double R     = fRo * A13 * units::fermi; // nuclear radius
      double R2    = TMath::Power(R,2.);
      double b     = 0.33333 * R2;
      double tsum  = (TMath::Exp(-b*tmin) - TMath::Exp(-b*tmax))/b; 
      double rt    = tsum * rnd->RndKine().Rndm();
      double gt    = -1.*TMath::Log(-1.*b*rt + TMath::Exp(-1.*b*tmin))/b;

      // TODO: If we re-install the fGenerateUniformly option, we 
      // would compute the event weight here.

      // reset bits
      interaction->ResetBit(kISkipProcessChk);
      interaction->ResetBit(kISkipKinematicChk);

      // lock selected kinematics & clear running values
      interaction->KinePtr()->Setx(gx, true);
      interaction->KinePtr()->Sety(gy, true);
      interaction->KinePtr()->Sett(gt, true);
      interaction->KinePtr()->SetW(this->pionMass(interaction), true);
      interaction->KinePtr()->SetQ2(gQ2, true);
      interaction->KinePtr()->ClearRunningValues();

      // set the cross section for the selected kinematics
      evrec->SetDiffXSec(xsec * TMath::Exp(-b * gt) / tsum, kPSQ2yfE);

      return;
    }
  }// iterations
}
//___________________________________________________________________________
void COHKinematicsGenerator::CalculateKin_BergerSehgalFM(GHepRecord * evrec) const
{
  // Get the Primary Interacton object
  Interaction * interaction = evrec->Summary();
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

  // Initialise a random number generator 
  RandomGen * rnd = RandomGen::Instance();

  //-- For the subsequent kinematic selection with the rejection method:
  //   Calculate the max differential cross section or retrieve it from the
  //   cache. Throw an exception and quit the evg thread if a non-positive
  //   value is found.
  //
  //   TODO: We are not offering the "fGenerateUniformly" option here.
  double xsec_max = this->MaxXSec(evrec);

  //-- Get the kinematical limits for the generated x,y
  const KPhaseSpace & kps = interaction->PhaseSpace();
  Range1D_t y  = kps.YLim();
  Range1D_t Q2(fQ2Min,fQ2Max);
  assert(y.min>0. && y.max>0. && y.min<1. && y.max<1. && y.min<y.max);

  const double ymin  = y.min + kASmallNum;
  const double ymax  = y.max - kASmallNum;
  const double dy    = ymax - ymin;
  const double Q2min = Q2.min + kASmallNum;
  const double Q2max = Q2.max - kASmallNum; 
  const double dQ2   = Q2max - Q2min;
  const double tmin  = kASmallNum;
  const double tmax  = fTMax - kASmallNum; // TODO: Choose realistic t bounds
  const double dt    = tmax - tmin;

  //-- Try to select a valid (Q^2,y,t) triple.

  unsigned int iter = 0;
  bool accept=false;
  double xsec=-1, gy=-1, gt=-1, gQ2=-1;

  while(1) {
    iter++;
    if(iter > kRjMaxIterations) this->throwOnTooManyIterations(iter,evrec);

    //-- Select unweighted kinematics using importance sampling method. 
    // TODO: The importance sampling envelope is not used. Currently, 
    // we just employ a standard rejection-method approach.

    gy  = ymin  + dy  * rnd->RndKine().Rndm(); 
    gt  = tmin  + dt  * rnd->RndKine().Rndm(); 
    gQ2 = Q2min + dQ2 * rnd->RndKine().Rndm(); 

    LOG("COHKinematics", pINFO) << 
      "Trying: Q^2 = " << gQ2 << ", y = " << gy << ", t = " << gt;

    interaction->KinePtr()->Sety(gy);
    interaction->KinePtr()->Sett(gt);  
    interaction->KinePtr()->SetQ2(gQ2);

    // computing cross section for the current kinematics
    xsec = fXSecModel->XSec(interaction, kPSxyfE);

    //-- decide whether to accept the current kinematics
    accept = (xsec_max * rnd->RndKine().Rndm() < xsec);

    //-- If the generated kinematics are accepted, finish-up module's job
    if(accept) {
      LOG("COHKinematics", pNOTICE)
        << "Selected: Q^2 = " << gQ2 << ", y = " << gy << ", t = " << gt; 

      // TODO: If we re-install the fGenerateUniformly option, we 
      // would compute the event weight here.

      // reset bits
      interaction->ResetBit(kISkipProcessChk);
      interaction->ResetBit(kISkipKinematicChk);

      // lock selected kinematics & clear running values
      interaction->KinePtr()->SetQ2(gQ2, true);
      interaction->KinePtr()->Sety(gy, true);
      interaction->KinePtr()->Sett(gt, true);
      interaction->KinePtr()->SetW(this->pionMass(interaction), true); 
      interaction->KinePtr()->ClearRunningValues();

      // set the cross section for the selected kinematics
      evrec->SetDiffXSec(xsec, kPSxytfE);

      return;
    }
  }// iterations
}
//___________________________________________________________________________
void COHKinematicsGenerator::CalculateKin_ReinSehgal(GHepRecord * evrec) const
{
  // Get the Primary Interacton object
  Interaction * interaction = evrec->Summary();
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

  //-- Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  //-- For the subsequent kinematic selection with the rejection method:
  //   Calculate the max differential cross section or retrieve it from the
  //   cache. Throw an exception and quit the evg thread if a non-positive
  //   value is found.
  //   If the kinematics are generated uniformly over the allowed phase
  //   space the max xsec is irrelevant
  double xsec_max = (fGenerateUniformly) ? -1 : this->MaxXSec(evrec);

  //-- Get the kinematical limits for the generated x,y
  const KPhaseSpace & kps = interaction->PhaseSpace();
  Range1D_t y = kps.YLim();
  assert(y.min>0. && y.max>0. && y.min<1. && y.max<1. && y.min<y.max);

  const double xmin = kASmallNum;
  const double xmax = 1.- kASmallNum;
  const double ymin = y.min + kASmallNum;
  const double ymax = y.max - kASmallNum;
  const double dx   = xmax - xmin;
  const double dy   = ymax - ymin;

  //------ Try to select a valid x,y pair

  unsigned int iter = 0;
  bool accept=false;
  double xsec=-1, gx=-1, gy=-1;

  while(1) {
    iter++;
    if(iter > kRjMaxIterations) this->throwOnTooManyIterations(iter,evrec);

    if(fGenerateUniformly) {
      //-- Generate a x,y pair uniformly in the kinematically allowed range.
      gx = xmin + dx * rnd->RndKine().Rndm();
      gy = ymin + dy * rnd->RndKine().Rndm();

    } else {
      //-- Select unweighted kinematics using importance sampling method. 

      if(iter==1) {
        LOG("COHKinematics", pNOTICE) << "Initializing the sampling envelope";
        double Ev = interaction->InitState().ProbeE(kRfLab);
        fEnvelope->SetRange(xmin,ymin,xmax,ymax);
        fEnvelope->SetParameter(0, xsec_max);  
        fEnvelope->SetParameter(1, Ev);        
      }

      // Generate W,QD2 using the 2-D envelope as PDF
      fEnvelope->GetRandom2(gx,gy);
    }

    LOG("COHKinematics", pINFO) << "Trying: x = " << gx << ", y = " << gy;

    interaction->KinePtr()->Setx(gx);
    interaction->KinePtr()->Sety(gy);

    // computing cross section for the current kinematics
    xsec = fXSecModel->XSec(interaction, kPSxyfE);

    //-- decide whether to accept the current kinematics
    if(!fGenerateUniformly) {
      double max = fEnvelope->Eval(gx, gy);
      double t   = max * rnd->RndKine().Rndm();

      this->AssertXSecLimits(interaction, xsec, max);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
      LOG("COHKinematics", pDEBUG) 
        << "xsec= " << xsec << ", J= 1, Rnd= " << t;
#endif
      accept = (t<xsec);
    }
    else { 
      accept = (xsec>0);
    }

    //-- If the generated kinematics are accepted, finish-up module's job
    if(accept) {
      LOG("COHKinematics", pNOTICE) << "Selected: x = "<< gx << ", y = "<< gy;

      // the Rein-Sehgal COH cross section should be a triple differential cross section
      // d^2xsec/dxdydt where t is the the square of the 4p transfer to the
      // nucleus. The cross section used for kinematical selection should have
      // the t-dependence integrated out. The t-dependence is of the form
      // ~exp(-bt). Now that the x,y kinematical variables have been selected
      // we can generate a t using the t-dependence as a PDF.
      const InitialState & init_state = interaction->InitState();
      double Ev    = init_state.ProbeE(kRfLab);
      double Epi   = gy*Ev; // pion energy
      double Epi2  = TMath::Power(Epi,2);
      double pme2  = kPionMass2/Epi2;   
      double xME   = kNucleonMass*gx/Epi;
      double tA    = 1. + xME - 0.5*pme2;
      double tB    = TMath::Sqrt(1.+ 2*xME) * TMath::Sqrt(1.-pme2);
      double tmin  = 2*Epi2 * (tA-tB);
      double tmax  = 2*Epi2 * (tA+tB);
      double A     = (double) init_state.Tgt().A(); 
      double A13   = TMath::Power(A,1./3.);
      double R     = fRo * A13 * units::fermi; // nuclear radius
      double R2    = TMath::Power(R,2.);
      double b     = 0.33333 * R2;
      double tsum  = (TMath::Exp(-b*tmin) - TMath::Exp(-b*tmax))/b; 
      double rt    = tsum * rnd->RndKine().Rndm();
      double gt    = -1.*TMath::Log(-1.*b*rt + TMath::Exp(-1.*b*tmin))/b;

      LOG("COHKinematics", pNOTICE)
        << "Selected: t = "<< gt << ", from ["<< tmin << ", "<< tmax << "]";

      // for uniform kinematics, compute an event weight as
      // wght = (phase space volume)*(differential xsec)/(event total xsec)
      if(fGenerateUniformly) {
        double vol     = y.max-y.min; // dx=1, dt: irrelevant
        double totxsec = evrec->XSec();
        double wght    = (vol/totxsec)*xsec;
        LOG("COHKinematics", pNOTICE)  << "Kinematics wght = "<< wght;

        // apply computed weight to the current event weight
        wght *= evrec->Weight();
        LOG("COHKinematics", pNOTICE) << "Current event wght = " << wght;
        evrec->SetWeight(wght);
      }

      // reset bits
      interaction->ResetBit(kISkipProcessChk);
      interaction->ResetBit(kISkipKinematicChk);

      // lock selected kinematics & clear running values
      interaction->KinePtr()->Setx(gx, true);
      interaction->KinePtr()->Sety(gy, true);
      interaction->KinePtr()->Sett(gt, true);
      interaction->KinePtr()->SetW(kPionMass, true);
      interaction->KinePtr()->SetQ2(2*kNucleonMass*gx*gy*Ev, true);
      interaction->KinePtr()->ClearRunningValues();

      // set the cross section for the selected kinematics
      evrec->SetDiffXSec(xsec * TMath::Exp(-b * gt) / tsum, kPSxytfE);

      return;
    }
  }// iterations
}
//___________________________________________________________________________
void COHKinematicsGenerator::CalculateKin_AlvarezRuso(GHepRecord * evrec) const
{

  LOG("COHKinematics", pNOTICE) << "Using AlvarezRuso Model";
  // Get the Primary Interacton object
  Interaction * interaction = evrec->Summary();
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

  // Initialise a random number generator 
  RandomGen * rnd = RandomGen::Instance();

  //-- For the subsequent kinematic selection with the rejection method:
  //   Calculate the max differential cross section or retrieve it from the
  //   cache. Throw an exception and quit the evg thread if a non-positive
  //   value is found.
  //   If the kinematics are generated uniformly over the allowed phase
  //   space the max xsec is irrelevant
  double xsec_max = (fGenerateUniformly) ? -1 : this->MaxXSec(evrec);

  //Set up limits of integration variables
  // Primary lepton energy
  const double E_l_min = interaction->FSPrimLepton()->Mass();
  const double E_l_max = interaction->InitStatePtr()->GetProbeP4(kRfLab)->E() - kPionMass;
  // Primary lepton angle with respect to the beam axis
  const double ctheta_l_min = 0.4;
  const double ctheta_l_max = 1.0 - kASmallNum;
  // Pion angle with respect to the beam axis
  const double ctheta_pi_min = 0.4;
  const double ctheta_pi_max = 1.0 - kASmallNum;
  // Pion angle transverse to the beam axis
  const double phi_min = kASmallNum;
  const double phi_max = (2.0 * kPi) - kASmallNum;
  // 
  const double d_E_l = E_l_max - E_l_min;
  const double d_ctheta_l  = ctheta_l_max  - ctheta_l_min;
  const double d_ctheta_pi = ctheta_pi_max - ctheta_pi_min;
  const double d_phi = phi_max - phi_min;

  //------ Try to select a valid set of kinematics
  unsigned int iter = 0;
  bool accept=false;
  double xsec=-1, g_E_l=-1, g_theta_l=-1, g_phi_l=-1, g_theta_pi=-1, g_phi_pi=-1;
  double g_ctheta_l, g_ctheta_pi;

  while(1) {
    iter++;
    if(iter > kRjMaxIterations) this->throwOnTooManyIterations(iter,evrec);

    //Select kinematic point
    g_E_l = E_l_min + d_E_l * rnd->RndKine().Rndm();
    g_ctheta_l  = ctheta_l_min  + d_ctheta_l  * rnd->RndKine().Rndm();
    g_ctheta_pi = ctheta_pi_min + d_ctheta_pi * rnd->RndKine().Rndm();
    g_phi_l = phi_min + d_phi * rnd->RndKine().Rndm();
    // random phi is relative to phi_l
    g_phi_pi = g_phi_l + (phi_min + d_phi * rnd->RndKine().Rndm()); 
    g_theta_l = TMath::ACos(g_ctheta_l);
    g_theta_pi = TMath::ACos(g_ctheta_pi);

    LOG("COHKinematics", pINFO) << "Trying: Lep(" <<g_E_l << ", " << 
      g_theta_l << ", " << g_phi_l << ") Pi(" << 
      g_theta_pi << ",     " << g_phi_pi << ")";

    this->SetKinematics(g_E_l, g_theta_l, g_phi_l, g_theta_pi, g_phi_pi, 
                        interaction, interaction->KinePtr());

    // computing cross section for the current kinematics
    xsec = fXSecModel->XSec(interaction,kPSElOlOpifE) / (1E-38 * units::cm2);

    if (!fGenerateUniformly) {
      //-- decide whether to accept the current kinematics
      double t   = xsec_max * rnd->RndKine().Rndm();

      LOG("COHKinematics", pINFO) << "Got: xsec = " << xsec << ", t = " << 
        t << " (max_xsec = " << xsec_max << ")";

      this->AssertXSecLimits(interaction, xsec, xsec_max);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
      LOG("COHKinematics", pDEBUG)
        << "xsec= " << xsec << ", J= 1, Rnd= " << t;
#endif
      accept = (t<xsec);
    }
    else {
      accept = (xsec>0);
    }

    //-- If the generated kinematics are accepted, finish-up module's job
    if(accept) {
      LOG("COHKinematics", pNOTICE) << "Selected: Lepton(" << 
        g_E_l << ", " << g_theta_l << ", " << 
        g_phi_l << ") Pion(" << g_theta_pi << ", " << g_phi_pi << ")";

      double E_l = g_E_l;
      double theta_l = g_theta_l;
      double theta_pi = g_theta_pi;
      double phi_l = g_phi_l;
      double phi_pi = g_phi_pi;
      const TLorentzVector P4_nu = *(interaction->InitStatePtr()->GetProbeP4(kRfLab));
      double E_nu       = P4_nu.E();
      double E_pi= E_nu-E_l;
      double m_l = interaction->FSPrimLepton()->Mass();
      double m_pi = this->pionMass(interaction);

      double p_l = TMath::Sqrt(E_l*E_l - m_l*m_l);
      TVector3 lepton_3vector = TVector3(0,0,0);
      lepton_3vector.SetMagThetaPhi(p_l,theta_l,phi_l);
      TLorentzVector P4_lep    = TLorentzVector(lepton_3vector , E_l );

      double p_pi = TMath::Sqrt(E_pi*E_pi - m_pi*m_pi);
      TVector3 pion_3vector = TVector3(0,0,0);
      pion_3vector.SetMagThetaPhi(p_pi,theta_pi,phi_pi);
      TLorentzVector P4_pion   = TLorentzVector(pion_3vector   , E_pi);

      TLorentzVector q = P4_nu - P4_lep;
      double Q2 = -q.Mag2();
      double x = Q2/(2*E_pi*constants::kNucleonMass);
      double y = E_pi/E_nu;

      double t = TMath::Abs( (q - P4_pion).Mag2() );

      // for uniform kinematics, compute an event weight as
      // wght = (phase space volume)*(differential xsec)/(event total xsec)
      if(fGenerateUniformly) {
        // Phase space volume needs checking
        double vol     = d_E_l*d_ctheta_l*d_phi*d_ctheta_pi*d_phi;
        double totxsec = evrec->XSec();
        double wght    = (vol/totxsec)*xsec;
        LOG("COHKinematics", pNOTICE)  << "Kinematics wght = "<< wght;

        // apply computed weight to the current event weight
        wght *= evrec->Weight();
        LOG("COHKinematics", pNOTICE) << "Current event wght = " << wght;
        evrec->SetWeight(wght);
      }

      // reset bits
      interaction->ResetBit(kISkipProcessChk);
      interaction->ResetBit(kISkipKinematicChk);
      // lock selected kinematics & clear running values
      interaction->KinePtr()->Setx(x, true);
      interaction->KinePtr()->Sety(y, true);
      interaction->KinePtr()->Sett(t, true);
      interaction->KinePtr()->SetW(kPionMass, true);
      interaction->KinePtr()->SetQ2(2*kNucleonMass*x*y*E_nu, true);
      interaction->KinePtr()->ClearRunningValues();
      // set the cross section for the selected kinematics
      evrec->SetDiffXSec(xsec,kPSElOlOpifE);
      return;
    }
  }//while
}
//___________________________________________________________________________
void COHKinematicsGenerator::SetKinematics(const double E_l,
                                           const double theta_l,
                                           const double phi_l,
                                           const double theta_pi,
                                           const double phi_pi,
                                           const Interaction* interaction,
                                           Kinematics* kinematics) const
{
  const TLorentzVector P4_nu = *(interaction->InitStatePtr()->GetProbeP4(kRfLab));
  double E_nu       = P4_nu.E();
  double E_pi= E_nu-E_l;
  double m_l = interaction->FSPrimLepton()->Mass();
  double m_pi;
  if ( interaction->ProcInfo().IsWeakCC() ) {
    m_pi = constants::kPionMass;
  } else {
    m_pi = constants::kPi0Mass;
  }
  double p_l=0.0;
  if (E_l > m_l) {
    p_l = TMath::Sqrt(E_l*E_l - m_l*m_l);
  }
  TVector3 lepton_3vector = TVector3(0,0,0);
  lepton_3vector.SetMagThetaPhi(p_l,theta_l,phi_l);
  TLorentzVector P4_lep    = TLorentzVector(lepton_3vector , E_l );

  double p_pi=0.0;
  if (E_pi > m_pi) {
    p_pi = TMath::Sqrt(E_pi*E_pi - m_pi*m_pi);
  }
  TVector3 pion_3vector = TVector3(0,0,0);
  pion_3vector.SetMagThetaPhi(p_pi,theta_pi,phi_pi);
  TLorentzVector P4_pion   = TLorentzVector(pion_3vector   , E_pi);

  double Q2 = -(P4_nu-P4_lep).Mag2();
  double x = Q2/(2*E_pi*constants::kNucleonMass);
  double y = E_pi/E_nu;

  kinematics->Setx(x);
  kinematics->Sety(y);
  kinematics::UpdateWQ2FromXY(interaction);

  kinematics->SetFSLeptonP4(P4_lep );
  kinematics->SetHadSystP4 (P4_pion); // use Hadronic System variable to store pion momentum
}
//___________________________________________________________________________
bool COHKinematicsGenerator::CheckKinematics(const double E_l,
                                             const double /*  theta_l    */ ,  
                                             const double /*  phi_l      */ ,
                                             const double /*  theta_pi   */ ,
                                             const double /*  phi_pi     */ ,
                                             const Interaction* interaction) const
{
  const TLorentzVector P4_nu = *(interaction->InitStatePtr()->GetProbeP4(kRfLab));
  double E_nu       = P4_nu.E();
  double E_pi= E_nu-E_l;
  double m_l = interaction->FSPrimLepton()->Mass();
  double m_pi;
  if ( interaction->ProcInfo().IsWeakCC() ) {
    m_pi = constants::kPionMass;
  }
  else {
    m_pi = constants::kPi0Mass;
  }
  if (E_l <= m_l) {
    return false;
  }
  if (E_pi <= m_pi) {
    return false;
  }
  return true;
}
//___________________________________________________________________________
double COHKinematicsGenerator::ComputeMaxXSec(const Interaction * in) const
{
  // Computes the maximum differential cross section in the requested phase
  // space. This method overloads KineGeneratorWithCache::ComputeMaxXSec
  // method and the value is cached at a circular cache branch for retrieval
  // during subsequent event generation.

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  SLOG("COHKinematics", pDEBUG)
    << "Scanning the allowed phase space {K} for the max(dxsec/d{K})";
#endif
  double max_xsec = 0.;
  if (fXSecModel->Id().Name() == "genie::ReinSehgalCOHPiPXSec") {
    max_xsec = MaxXSec_ReinSehgal(in);
  } else if ((fXSecModel->Id().Name() == "genie::BergerSehgalCOHPiPXSec2015")) {
    max_xsec = MaxXSec_BergerSehgal(in);
  } else if ((fXSecModel->Id().Name() == "genie::BergerSehgalFMCOHPiPXSec2015")) {
    max_xsec = MaxXSec_BergerSehgalFM(in);
  } else if ((fXSecModel->Id().Name() == "genie::AlvarezRusoCOHPiPXSec")) {
    max_xsec = MaxXSec_AlvarezRuso(in);
  }
  else {
    LOG("COHKinematicsGenerator",pFATAL) <<
      "ComputeMaxXSec >> Cannot calculate max cross-section for " <<
      fXSecModel->Id().Name();
  }

  // Apply safety factor, since value retrieved from the cache might
  // correspond to a slightly different energy.
  max_xsec *= fSafetyFactor;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  SLOG("COHKinematics", pDEBUG) << in->AsString();
  SLOG("COHKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  SLOG("COHKinematics", pDEBUG) << "Computed using alg = " << fXSecModel->Id();
#endif

  return max_xsec;
}
//___________________________________________________________________________
double COHKinematicsGenerator::MaxXSec_BergerSehgal(const Interaction * in) const
{
  double max_xsec = 0;
  const int NQ2   = 50;
  const int Ny    = 50;

  const KPhaseSpace & kps = in->PhaseSpace();
  Range1D_t Q2r = kps.Q2Lim();
  Q2r.max = fQ2Max;

  const double logQ2min = TMath::Log10(Q2r.min + kASmallNum);
  const double logQ2max = TMath::Log10(Q2r.max); 
  const double dlogQ2   = (logQ2max - logQ2min) /(NQ2-1);

  for(int i=0; i<NQ2; i++) {
    double Q2 = TMath::Power(10, logQ2min + i * dlogQ2);
    in->KinePtr()->SetQ2(Q2);

    Range1D_t yr = kps.YLim();
    if ((yr.max < 0) || (yr.max < yr.min) || 
        (yr.max > 1) || (yr.min < 0)) { // forbidden kinematics
      continue;
    }
    const double logymin  = TMath::Log10(yr.min);
    const double logymax  = TMath::Log10(yr.max);
    const double dlogy    = (logymax - logymin) /(Ny-1);

    for(int j=0; j<Ny; j++) {
      double gy = TMath::Power(10, logymin + j * dlogy);
      in->KinePtr()->Sety(gy);

      /* Range1D_t tl = kps.TLim();   // TESTING! - this becomes a loop over t */
      kinematics::UpdateXFromQ2Y(in);

      // Note: We're not stepping through log Q^2, log y - we "unpacked"
      double xsec = fXSecModel->XSec(in, kPSQ2yfE);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
      LOG("COHKinematics", pDEBUG)  
        << "xsec(Q2= " << Q2 << ", y= " << gy << ", t = " << gt << ") = " << xsec;
#endif
      max_xsec = TMath::Max(max_xsec, xsec);

    } // y
  } // Q2
  return max_xsec;
}
//___________________________________________________________________________
double COHKinematicsGenerator::MaxXSec_BergerSehgalFM(const Interaction * in) const
{
  double max_xsec = 0;
  // How many sampling bins in each variable for max xsec calculation?
  const int NQ2   = 50;
  const int Ny    = 50;
  const int Nt    = 50;

  const KPhaseSpace & kps = in->PhaseSpace();
  Range1D_t Q2r = kps.Q2Lim();
  Q2r.max = fQ2Max;

  const double logQ2min = TMath::Log10(Q2r.min + kASmallNum);
  const double logQ2max = TMath::Log10(Q2r.max); 
  const double logtmin  = TMath::Log10(kASmallNum);
  const double logtmax  = TMath::Log10(fTMax - kASmallNum); 
  const double dlogQ2   = (logQ2max - logQ2min) /(NQ2-1);
  const double dlogt    = (logtmax - logtmin) /(Nt-1);

  for(int i=0; i<NQ2; i++) {
    double Q2 = TMath::Power(10, logQ2min + i * dlogQ2);
    in->KinePtr()->SetQ2(Q2);

    Range1D_t yr = kps.YLim();
    if ((yr.max < 0) || (yr.max < yr.min) || 
        (yr.max > 1) || (yr.min < 0)) { // forbidden kinematics
      continue;
    }
    const double logymin  = TMath::Log10(yr.min);
    const double logymax  = TMath::Log10(yr.max);
    const double dlogy    = (logymax - logymin) /(Ny-1);

    for(int j=0; j<Ny; j++) {
      double gy = TMath::Power(10, logymin + j * dlogy);
      
      for(int k=0; k<Nt; k++) {
        double gt = TMath::Power(10, logtmin + k * dlogt);

        in->KinePtr()->Sety(gy);
        in->KinePtr()->Sett(gt);

        double xsec = fXSecModel->XSec(in, kPSxyfE);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
        LOG("COHKinematics", pDEBUG)  
          << "xsec(Q2= " << Q2 << ", y= " << gy << ", t = " << gt << ") = " << xsec;
#endif
        max_xsec = TMath::Max(max_xsec, xsec);

      } // t
    } // y
  } // Q2
  return max_xsec;
}
//___________________________________________________________________________
double COHKinematicsGenerator::MaxXSec_ReinSehgal(const Interaction * in) const
{
  double max_xsec = 0;
  double Ev = in->InitState().ProbeE(kRfLab);

  const int Nx = 50;
  const int Ny = 50;

  const KPhaseSpace & kps = in->PhaseSpace();
  Range1D_t y = kps.YLim();

  const double logxmin = TMath::Log10(1E-5);
  const double logxmax = TMath::Log10(1.0);
  const double logymin = TMath::Log10(y.min);
  const double logymax = TMath::Log10(y.max);

  const double dlogx   = (logxmax - logxmin) /(Nx-1);
  const double dlogy   = (logymax - logymin) /(Ny-1);

  for(int i=0; i<Nx; i++) {
    double gx = TMath::Power(10, logxmin + i * dlogx);
    for(int j=0; j<Ny; j++) {
      double gy = TMath::Power(10, logymin + j * dlogy);

      double Q2 = 2*kNucleonMass*gx*gy*Ev;
      if(Ev>1.0 && Q2>0.01) continue;

      in->KinePtr()->Setx(gx);
      in->KinePtr()->Sety(gy);

      double xsec = fXSecModel->XSec(in, kPSxyfE);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
      LOG("COHKinematics", pDEBUG)  
        << "xsec(x= " << gx << ", y= " << gy << ") = " << xsec;
#endif
      max_xsec = TMath::Max(max_xsec, xsec);

    }//y
  }//x
  return max_xsec;
}
//___________________________________________________________________________
double COHKinematicsGenerator::MaxXSec_AlvarezRuso(const Interaction * in) const
{
  // Computes the maximum differential cross section in the requested phase
  // space. This method overloads KineGeneratorWithCache::ComputeMaxXSec
  // method and the value is cached at a circular cache branch for retrieval
  // during subsequent event generation.

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  SLOG("COHKinematics", pDEBUG)
    << "Scanning the allowed phase space {K} for the max(dxsec/d{K})";
#endif
  double max_xsec = 0.;
  double Ev = in->InitState().ProbeE(kRfLab);

  const KPhaseSpace & kps = in->PhaseSpace();
  Range1D_t y = kps.YLim();

  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit2");
  gsl::d4Xsec_dEldThetaldOmegapi f(fXSecModel,in);
  f.SetFactor(-1.); // Make it return negative of cross-section so we can minimize

  min->SetFunction( f );
  min->SetMaxFunctionCalls(10000);
  min->SetTolerance(0.05);

  const double min_el = in->FSPrimLepton()->Mass();
  const double max_el = Ev - kPionMass;
  const unsigned int n_el = 100;
  const double d_el = (max_el - min_el) / double(n_el - 1);

  const double min_thetal = kASmallNum;
  const double max_thetal = kPi / 4.0;
  const unsigned int n_thetal = 10;
  const double d_thetal = (max_thetal - min_thetal) / double(n_thetal - 1);

  const double min_thetapi = kASmallNum;
  const double max_thetapi = kPi / 2.0;
  const unsigned int n_thetapi = 10;
  const double d_thetapi = (max_thetapi - min_thetapi) / double(n_thetapi - 1);

  //~ const double min_phipi = kPi;
  //~ const double max_phipi = 0.5 * kPi;
  const double min_phipi = kASmallNum;
  const double max_phipi = 2*kPi-kASmallNum;
  const unsigned int n_phipi = 10;
  const double d_phipi = (max_phipi - min_phipi) / double(n_phipi - 1);

  min->SetLimitedVariable ( 0 ,"E_lep"    , max_el     -kASmallNum , d_el      , min_el     , max_el      );
  min->SetLimitedVariable ( 1 ,"theta_l"  , min_thetal +kASmallNum , d_thetal  , min_thetal , max_thetal  );
  min->SetLimitedVariable ( 2 ,"theta_pi" , min_thetapi+kASmallNum , d_thetapi , min_thetapi, max_thetapi );
  min->SetLimitedVariable ( 3 ,"phi_pi"   , min_phipi  +kASmallNum , d_phipi   , min_phipi  , max_phipi   );

  min->Minimize();
  max_xsec = -min->MinValue(); //back to positive xsec

  // Apply safety factor, since value retrieved from the cache might
  // correspond to a slightly different energy.
  max_xsec *= fSafetyFactor;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  SLOG("COHKinematics", pDEBUG) << in->AsString();
  SLOG("COHKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  SLOG("COHKinematics", pDEBUG) << "Computed using alg = " << fXSecModel->Id();
#endif

  delete min;

  return max_xsec;
}
//___________________________________________________________________________
double COHKinematicsGenerator::Energy(const Interaction * interaction) const
{
  // Override the base class Energy() method to cache the max xsec for the
  // neutrino energy in the LAB rather than in the hit nucleon rest frame.

  const InitialState & init_state = interaction->InitState();
  double E = init_state.ProbeE(kRfLab);
  return E;
}
//___________________________________________________________________________
double COHKinematicsGenerator::pionMass(const Interaction* in) const
{
  double m_pi = 0.0;
  if ( in->ProcInfo().IsWeakCC() ) {
    m_pi = constants::kPionMass;
  } else {
    m_pi = constants::kPi0Mass;
  }
  return m_pi;
}
//___________________________________________________________________________
void COHKinematicsGenerator::throwOnTooManyIterations(unsigned int iters,
                                                      GHepRecord* evrec) const
{
  LOG("COHKinematics", pWARN)
    << "*** Could not select valid kinematics after "
    << iters << " iterations";
  evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
  genie::exceptions::EVGThreadException exception;
  exception.SetReason("Couldn't select kinematics");
  exception.SwitchOnFastForward();
  throw exception;
}
//___________________________________________________________________________
void COHKinematicsGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHKinematicsGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHKinematicsGenerator::LoadConfig(void)
{
  //-- COH model parameter Ro
  GetParam( "COH-Ro", fRo );
  //-- COH model parameter t_max for t = (q - p_pi)^2
  GetParam( "COH-t-max", fTMax ) ;
  //-- COH model bounds of integration for Q^2
  GetParam( "COH-Q2-min", fQ2Min ) ;
  GetParam( "COH-Q2-max", fQ2Max ) ;

  //-- max xsec safety factor (for rejection method) and min cached energy
  GetParamDef( "MaxXSec-SafetyFactor", fSafetyFactor, 1.6 ) ;
  GetParamDef( "Cache-MinEnergy", fEMin,  -1.0 ) ;

  //-- Generate kinematics uniformly over allowed phase space and compute
  //   an event weight?
  GetParamDef( "UniformOverPhaseSpace", fGenerateUniformly, false ) ;

  //-- Maximum allowed fractional cross section deviation from maxim cross
  //   section used in rejection method
  GetParamDef( "MaxXSec-DiffTolerance", fMaxXSecDiffTolerance, 999999. ) ;
    assert(fMaxXSecDiffTolerance>=0);

  //-- Envelope employed when importance sampling is used 
  //   (initialize with dummy range)
  if(fEnvelope) delete fEnvelope;
  fEnvelope = new TF2("CohKinEnvelope",
                      kinematics::COHImportanceSamplingEnvelope,0.,1,0.,1,2);
  // stop ROOT from deleting this object of its own volition
  gROOT->GetListOfFunctions()->Remove(fEnvelope);
}
//____________________________________________________________________________

