//____________________________________________________________________________
/*!

\class    genie::QPMDISStrucFunc

\brief    Standard Quark Parton Model (QPM) Deep Inelastic Scatering (DIS)
          Structure Functions (SF)

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 03, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
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

