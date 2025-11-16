//____________________________________________________________________________
/*!

\namespace  genie::utils::style

\brief      GENIE style!

\author     Costas Andreopoulos <c.andreopoulos \at cern.ch>
            University of Liverpool

\created    July 29, 2010

\cpright    Copyright (c) 2003-2025, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org            
*/
//____________________________________________________________________________

#ifndef _STYLE_UTILS_H_
#define _STYLE_UTILS_H_

class TGraph;
class TH1;

namespace genie {
namespace utils {

namespace style
{

  void SetDefaultStyle(bool black_n_white=false);

  void Format(TGraph* gr,
        int lcol, int lsty, int lwid, int mcol, int msty, double msiz);
  void Format(TH1* hst,
        int lcol, int lsty, int lwid, int mcol, int msty, double msiz);

} // style namespace
} // utils namespace
} // genie namespace

#endif // _SYST_UTILS_H_
