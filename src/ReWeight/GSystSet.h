//____________________________________________________________________________
/*!

\class    genie::rew::GSystSet

\brief    Set of systematics to be considered 

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Aug 1, 2009

\cpright  Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _G_SET_OF_SYSTEMATICS_H_
#define _G_SET_OF_SYSTEMATICS_H_

#include <string>
#include <map>
#include <vector>

#include "ReWeight/GSyst.h"
#include "ReWeight/GSystType.h"

using std::string;
using std::map;
using std::vector;

namespace genie {
namespace rew   {

class GSystInfo;

class GSystSet {

public:  
  GSystSet();
  GSystSet(const GSystSet & syst_set);
 ~GSystSet();

  void    Include      (GSyst_t syst);
  void    Remove       (GSyst_t syst);
  int     NIncluded    (void) const;
  bool    IsIncluded   (GSyst_t syst) const;
  double  CurValue     (GSyst_t syst) const;
  double  DefValue     (GSyst_t syst) const;
  double  MinValue     (GSyst_t syst) const;
  double  MaxValue     (GSyst_t syst) const;
  void    SetCurValue  (GSyst_t syst, double val);
  void    SetDefValue  (GSyst_t syst, double val);
  void    SetRange     (GSyst_t syst, double min, double max);
  void    PrintSummary (void);
  void    Copy         (const GSystSet & syst_set);

  vector<genie::rew::GSyst_t> AllIncluded (void);

private:
  
  map<GSyst_t, GSystInfo *>  fSystematics;  
};

class GSystInfo {

public:
  GSystInfo() : 
     CurValue(0), DefValue(0), MinValue(0), MaxValue(0) 
  { 

  }
  GSystInfo(double val, double def, double min, double max) : 
     CurValue(val), DefValue(def), MinValue(min), MaxValue(max) 
  {
 
  }
 ~GSystInfo() 
  { 

  }

  double CurValue;
  double DefValue;
  double MinValue;
  double MaxValue;
};

} // rew   namespace
} // genie namespace

#endif 

