//____________________________________________________________________________
/*!

\namespace genie::intranuke

\brief     INTRANUKE utilities

\author    Jim Dobson <j.dobson07 \at imperial.ac.uk>
           Imperial College London

           Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
           University of Liverpool & STFC Rutherford Appleton Lab

	   Aaron Meyer <asm58 \at pitt.edu>
	   Pittsburgh University

\created   Mar 03, 2009

\cpright   Copyright (c) 2003-2019, The GENIE Collaboration
           For the full text of the license visit http://copyright.genie-mc.org
           or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _INTRANUKE_UTILS_2018_H_
#define _INTRANUKE_UTILS_2018_H_

#include <TGenPhaseSpace.h>

#include "Physics/NuclearState/NuclearModelI.h"
#include "Physics/HadronTransport/INukeHadroFates2018.h"
#include "Physics/HadronTransport/INukeMode.h"
#include "Physics/HadronTransport/INukeNucleonCorr.h"

class TLorentzVector;

namespace genie {

class GHepRecord;
class GHepParticle;
class PDGCodeList;

namespace utils {
namespace intranuke2018
{
  //! Hadron survival probability
  double ProbSurvival(
    int pdgc, const TLorentzVector & x4, const TLorentzVector & p4, double A,
    double Z, double mfp_scale_factor=1.0,
    double nRpi=0.5, double nRnuc=1.0, double NR=3, double R0=1.4);

  //! Mean free path (pions, nucleons)
  double MeanFreePath(
    int pdgc, const TLorentzVector & x4, const TLorentzVector & p4, double A,
    double Z, double nRpi=0.5, double nRnuc=1.0, const bool useOset = false, const bool altOset = false, const bool xsecNNCorr = false, string INukeMode = "XX2018");
 
  //! Mean free path (Delta++ **test**)
  double MeanFreePath_Delta(
			    int pdgc, const TLorentzVector & x4, const TLorentzVector & p4, double A );

  //! Distance to exit
  double Dist2Exit(
    const TLorentzVector & x4, const TLorentzVector & p4, 
    double A, double NR=3, double R0=1.4);

  //! Distance to exit
  double Dist2ExitMFP(
   int pdgc, const TLorentzVector & x4, const TLorentzVector & p4, 
    double A, double Z, double NR=3, double R0=1.4);

  //! Step particle
  void StepParticle(
    GHepParticle * p, double step, double nuclear_radius=-1.);


  //! Intranuke utility functions

  bool TwoBodyCollision(
    GHepRecord* ev, int pcode, int tcode, int scode, int s2code, double C3CM, GHepParticle* p,
    GHepParticle* t, int &RemnA, int &RemnZ, TLorentzVector &RemnP4, EINukeMode mode=kIMdHA);

  bool TwoBodyKinematics(
    double M3, double M4, TLorentzVector tP1L, TLorentzVector tP2L, 
    TLorentzVector &tP3L, TLorentzVector &tP4L, double C3CM, TLorentzVector &RemnP4, double bindE=0);

  bool ThreeBodyKinematics(
    GHepRecord* ev, GHepParticle* p, int tcode, GHepParticle* s1, GHepParticle* s2, GHepParticle* s3,
    bool DoFermi=false, double FermiFac=0, double FermiMomentum=0, const NuclearModelI* Nuclmodel=(const NuclearModelI*)0);

  bool PionProduction(
    GHepRecord* ev, GHepParticle* p, GHepParticle* s1, GHepParticle* s2, GHepParticle* s3, int &RemnA, int &RemnZ,
    TLorentzVector &RemnP4,bool DoFermi, double FermiFac, double FermiMomentum, const NuclearModelI* Nuclmodel);

  double CalculateEta(
    double Minc, double ke, double Mtarg, double Mtwopart, double Mpi);

  void Equilibrium(
    GHepRecord* ev, GHepParticle* p, int &RemnA, int &RemnZ, TLorentzVector &RemnP4, bool DoFermi,
    double FermiFac, const NuclearModelI* Nuclmodel, double NucRmvE, EINukeMode mode=kIMdHN);

  void PreEquilibrium(
    GHepRecord* ev, GHepParticle* p, int &RemnA, int &RemnZ, TLorentzVector &RemnP4, bool DoFermi,
    double FermiFac, const NuclearModelI* Nuclmodel, double NucRmvE, EINukeMode mode=kIMdHN);


  //! general phase space decay method
  bool PhaseSpaceDecay (
    GHepRecord* ev, GHepParticle* p, const PDGCodeList & pdgv, TLorentzVector &RemnP4,
    double NucRmvE, EINukeMode mode=kIMdHA);

  // calculate pion-nucleon total cross section based on Oset model
  // use only for pion with kinetic energy up to 350 MeV
  double sigmaTotalOset (const double &pionKineticEnergy, const double &density,
                         const int &pionPDG, const double &protonFraction,
                         const bool &isTableChosen = true
                         );

}      // intranuke namespace
}      // utils     namespace
}      // genie     namespace


#endif // _INTRANUKE_UTILS_2018_H_
