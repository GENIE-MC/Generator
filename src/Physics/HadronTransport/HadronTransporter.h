//____________________________________________________________________________
/*!

\class    genie::HadronTransporter

\brief    Intranuclear hadronic transport module. 
          It is being used to transfer all hadrons outside the nucleus without
          rescattering -if rescattering is switched off- or to call one of the 
          supported hadron transport MCs -if rescattering is switched on-
         
\author   Costas Andreopoulos <c.andreopoulos \at cern.ch> STFC, Rutherford Lab

\created  September 14, 2006

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          
*/
//____________________________________________________________________________

#ifndef _HADRON_TRANSPORTER_H_
#define _HADRON_TRANSPORTER_H_

#include "Framework/EventGen/EventRecordVisitorI.h"

namespace genie {

class HadronTransporter : public EventRecordVisitorI {

public :
  HadronTransporter();
  HadronTransporter(string config);
 ~HadronTransporter();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  // override the Algorithm::Configure methods to load configuration
  // data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:
  void  LoadConfig                (void);
  void  TransportInTransparentNuc (GHepRecord * ev) const;

  bool  fEnabled;                              ///< hadron transport enabled?
  const EventRecordVisitorI * fHadTranspModel; ///< hadron transport MC to use
  const EventRecordVisitorI * fCascadeReweight; ///< Cascade reweight member

};

}      // genie namespace
#endif // _HADRON_TRANSPORTER_H_
