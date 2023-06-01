//____________________________________________________________________________
/*!

\class    genie::ElectronVelocityMap

\brief    This class is a hook for  Electron velocity distributions and allows associating each
          one of them with specific nuclei.
          Is a concrete implementation of the ElectronVelocity interface.

 \author   Marco Roda <mroda \at liverpool.ac.uk>
          University of Liverpool
  
  \created March 28, 2023

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          
*/
//____________________________________________________________________________

#ifndef _ELECTRON_VELOCITY_MAP_H_
#define _ELECTRON_VELOCITY_MAP_H_

#include "Physics/NuclearState/ElectronVelocity.h"

namespace genie {

class ElectronVelocityMap : public ElectronVelocity {

public :

  ElectronVelocityMap();
  ElectronVelocityMap(string config);
  
  ElectronVelocityMap(const ElectronVelocityMap & ) = delete;
  ElectronVelocityMap(ElectronVelocityMap && ) = delete;
  

  virtual ~ElectronVelocityMap() {;}
  
  //-- overload the ElectronVelocity::Configure() methods 
  // to load data from ModelConfig.xml
  void Configure(string config) override final ;

  void InitializeVelocity(Interaction & interaction) const override; //Give initial velocity

protected:
  void LoadConfig () override;
  const ElectronVelocity & SelectModel(const Target &) const;

private:
  const ElectronVelocity * fDefGlobalVelocity = nullptr;
  map<int, const ElectronVelocity *> fSpecificModels ;
  // the key is the Z of the atom

};

}      // genie namespace
#endif // _ELECTRON_VELOCITY_MAP_H_
