//____________________________________________________________________________
/*!

\class    genie::alvarezruso::ARWavefunction

\brief    Wave function class for AlvarezRuso Coherent pion production xsec

\ref      

\author   Steve Dennis
          University of Warwick, Rutherford Appleton Laboratory

\created  05/12/2013

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________
#ifndef _AR_WAVEFUNCTION_H_
#define _AR_WAVEFUNCTION_H_

#include <string>
#include <complex>

namespace genie
{
namespace alvarezruso
{

class ARWavefunction
{
  public:
    
    ARWavefunction(unsigned int sampling_in, bool debug = false);
    
    ~ARWavefunction();
    
    std::string print() const;
    
    const std::vector<std::complex<double> >& operator[] (unsigned int i) const;
    
    const std::complex<double> & operator() (unsigned int i, unsigned int j) const;
    
    std::complex<double>  get(unsigned int i, unsigned int j) const;
    
    void set(unsigned int i, unsigned int j, const std::complex<double> & value);
    
    unsigned int sampling() const;
    
    private:
      
      bool debug_;
      unsigned int sampling_;
      std::vector< std::vector<std::complex<double> > > wavefunction_;
}; // class ARWavefunction

} //namespace alvarezruso
} //namespace genie

#endif
