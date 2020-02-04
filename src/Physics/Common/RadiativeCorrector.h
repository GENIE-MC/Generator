//____________________________________________________________________________
/*!

\class    genie::RadiativeCorrector

\brief    

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  June 16, 2007

\cpright  Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _RADIATIVE_CORRECTOR_H_
#define _RADIATIVE_CORRECTOR_H_

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/Utils/Range1.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Interaction/InitialState.h"
#include "Framework/Ntuple/NtpWriter.h"
#include "Framework/Ntuple/NtpMCFormat.h"
#include "Framework/GHEP/GHepRecordHistory.h"

namespace genie {

class GHepParticle;
class Interaction;
class InitialState;
class Target;
class NtpMCFormat;

class RadiativeCorrector : public EventRecordVisitorI {

public :
  RadiativeCorrector();
  RadiativeCorrector(string config);
 ~RadiativeCorrector();
 
  //double r_eprime;
  //double r_etheta;
  //double r_e;

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord (GHepRecord * event_rec) const;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);
  void Configure(const InitialState & is);
  void SetISR(bool isr); 
  void SetModel(std::string model); 
  void SetQ2(double Q2); 
  void SetP4l(TLorentzVector p4l); 
  void SetCutoff(double cutoff); 
  void SetThickness(double thickness); 
  void SetDoInternalRad(bool doInternal);


  //int           rad_kDefOptNevents;       // n-events to generate
  //NtpMCFormat_t rad_kDefOptNtpFormat; // ntuple format
  //Long_t        rad_kDefOptRunNu;       // default run number
  //Long_t        rad_gOptRunNu;        // run number


private:

  void  LoadConfig (void);
  void  BuildInitialState(const InitialState & init_state);
  bool  ToBeDecayed       (GHepParticle * particle) const;

  bool 		 fISR;                 ///< to distinguish between ISR and FSR
  std::string	 fModel;               ///< to distinguish between differnt models, right now simc / vanderhagen
  //double         fRadiatedEnergyLimit; 
  //map<int, double> fThicknesses;       ///< thicnknesses of targets in CLAS in radiation length
  double fThickness;       ///< thicnknesses of targets in CLAS in radiation length

  double         fQ2; 
  TLorentzVector fP4l; 
  InitialState * fInitState;       ///< initial state information for changing probe
  double         fCutoff;
  bool		 fDoInternal;
  mutable GHepRecordHistory             fRecHistory;     ///< event record history 

};

}      // genie namespace
#endif // _RADIATIVE_CORRECTOR_H_
