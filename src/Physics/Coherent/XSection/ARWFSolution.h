//____________________________________________________________________________
/*!

\class    genie::alvarezruso::ARWFSolution

\brief    Abstract base class for Alvarez-Ruso wavefunction solution.

\ref      

\author   Steve Dennis
          University of Warwick, Rutherford Appleton Laboratory

\created  05/12/2013

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _AR_WF_SOLUTION_H_
#define _AR_WF_SOLUTION_H_

#include <complex>

namespace genie
{
namespace alvarezruso
{

enum pUnits_t{ kInMeV=0, kInNatural };

class ARWFSolution
{
  public:
    
    ARWFSolution(bool debug = false);
    virtual ~ARWFSolution();
    virtual std::complex<double>  Element(const double radius, const double cosine_rz, const double e_pion) = 0;
    virtual void Solve() = 0;
    bool debug_;
};

} //namespace alvarezruso
} //namespace genie

#endif
