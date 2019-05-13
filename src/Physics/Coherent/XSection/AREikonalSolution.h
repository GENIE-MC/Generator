//____________________________________________________________________________
/*!

\class    genie::alvarezruso::AREikonalSolution

\brief    Eikonal wavefunction solution for Alvarez-Ruso Coherent Pion Production xsec

\ref      

\author   Steve Dennis
          University of Warwick, Rutherford Appleton Laboratory

\created  05/12/2013

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________
#ifndef _AR_EIKONAL_SOLUTION_H_
#define _AR_EIKONAL_SOLUTION_H_

#include <complex>

#include "Physics/Coherent/XSection/AlvarezRusoCOHPiPDXSec.h"
#include "Physics/Coherent/XSection/ARSampledNucleus.h"
#include "Physics/Coherent/XSection/ARConstants.h"
#include "Physics/Coherent/XSection/ARWFSolution.h"

namespace genie
{
namespace alvarezruso
{

class AREikonalSolution: public ARWFSolution
{
  public:
  
    AREikonalSolution(bool debug, AlvarezRusoCOHPiPDXSec* parent);
    AREikonalSolution(bool debug, ARSampledNucleus* nucl);

    virtual ~AREikonalSolution();
    virtual std::complex<double>  Element(const double radius, const double cosine_rz, 
                 const double e_pion);
    void Solve();
  
  private:
  
    AlvarezRusoCOHPiPDXSec* Parent() { return this->parent_; }
    ARSampledNucleus* Nucleus() { return fNucleus; }
    ARConstants* Con() { return this->constants_; }
    
    std::complex<double>  PionSelfEnergy(const double rhop_cent, const double rhon_cent, 
                const double omepi, const double ppim);
    void Deltamed(const double sdel, const double pf, const double rat, double& gamdpb, 
            double& imsig, const double ppim, const double omepi);
    double Cc(const double a, const double b, const double c, const double ome);
    double Gamd(const double s);
    double Qcm(const double s);
    
    AlvarezRusoCOHPiPDXSec* parent_;
    ARSampledNucleus* fNucleus;
    ARConstants* constants_;
    
    bool owns_constants;
    
};


} //namespace alvarezruso
} //namespace genie

#endif
