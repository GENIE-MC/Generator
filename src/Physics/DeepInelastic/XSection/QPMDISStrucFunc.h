//____________________________________________________________________________
/*!

\class    genie::QPMDISStrucFunc

\brief    Standard Quark Parton Model (QPM) Deep Inelastic Scatering (DIS)
          Structure Functions (SF)

\author   Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

\created  May 03, 2004

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _QPM_DIS_STRUC_FUNC_H_
#define _QPM_DIS_STRUC_FUNC_H_

#include "Physics/DeepInelastic/XSection/QPMDISStrucFuncBase.h"

namespace genie {

class QPMDISStrucFunc : public QPMDISStrucFuncBase {

public:
  QPMDISStrucFunc();
  QPMDISStrucFunc(string config);
  virtual ~QPMDISStrucFunc();
};

}         // genie namespace
#endif    // _QPM_DIS_STRUC_FUNC_H_
