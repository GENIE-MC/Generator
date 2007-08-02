//____________________________________________________________________________
/*!

\class    genie::GIBUU

\brief    A GENIE interface to the GIBUU hadron transport code.
          Is a concerete implementation of the EventRecordVisitorI interface.
          Note: To use this event generation module you need to obtain GIBUU
          from its official distribution point and enable it during the GENIE
          installation.

\ref      http://tp8.physik.uni-giessen.de:8080/GiBUU/

          GIBUU team: 
          U.Model, O.Bub, T.Gaitanos, K.Gallmeister, D.Kalok, M.Kaskulov,
          A.Larionov, T.Leitner, B.Steinmu¼lle, L. Alvarez-Ruso, P. Mu¼hlich 

          Please cite GIBUU separately if you include this event generation
          module in your event generation threads.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk> STFC, Rutherford Lab
          Tina Leitner <Tina.J.Leitner@theo.physik.uni-giessen.de> Giessen Univ.

\created  December 13, 2005

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GIBUU_H_
#define _GIBUU_H_

#include "Conventions/GBuild.h"
#include "EVGCore/EventRecordVisitorI.h"

#ifdef __GENIE_GIBUU_ENABLED__

// fortran90 subroutine/function name aliases
//
#define SetNucleus              __gibuu_genie__setnucleus 
#define AddHadron               __gibuu_genie__addhadron
#define HadronTransportForGENIE __gibuu_genie__hadrontransportforgenie
#define NumOfFinalStateHadrons  __gibuu_genie__numoffinalstatehadrons
#define FinalStateHadronPDG     __gibuu_genie__finalstatehadronpdg
#define FinalStateHadronPX      __gibuu_genie__finalstatehadronpx
#define FinalStateHadronPY      __gibuu_genie__finalstatehadronpy
#define FinalStateHadronPZ      __gibuu_genie__finalstatehadronpz
#define FinalStateHadronE       __gibuu_genie__finalstatehadrone
#define FinalStateHadronX       __gibuu_genie__finalstatehadronx
#define FinalStateHadronY       __gibuu_genie__finalstatehadrony
#define FinalStateHadronZ       __gibuu_genie__finalstatehadronz
#define FinalStateHadronT       __gibuu_genie__finalstatehadront

// C bindings for fortran90 GIBUU/GENIE module subroutines/functions
//
extern "C" {
  void  SetNucleus              (int*);
  void  AddHadron               (int*, int*,float*,float*,float*,float*,float*,float*,float*,float*);
  void  HadronTransportForGENIE (void);
  int   NumOfFinalStateHadrons  (int*);
  int   FinalStateHadronPDG     (int*,int*);
  float FinalStateHadronPX      (int*,int*);
  float FinalStateHadronPY      (int*,int*);
  float FinalStateHadronPZ      (int*,int*);
  float FinalStateHadronE       (int*,int*);
  float FinalStateHadronX       (int*,int*);
  float FinalStateHadronY       (int*,int*);
  float FinalStateHadronZ       (int*,int*);
  float FinalStateHadronT       (int*,int*);
}
#endif // __GENIE_GIBUU_ENABLED__

namespace genie {

class GIBUU : public EventRecordVisitorI {

public:
  GIBUU();
  GIBUU(string config);
 ~GIBUU();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //-- members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void LoadConfig (void);
};

}      // genie namespace
#endif // _GIBUU_H_
