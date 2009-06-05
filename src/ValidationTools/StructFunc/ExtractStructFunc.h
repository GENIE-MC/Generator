//____________________________________________________________________________
/*!

\class   ExtractStructFunc

\brief   Extract structure functions treating the generator as a 'black box'
         using the method described in:

         H.Gallagher, Nucl.Phys.Proc.Suppl.159:229-234,2006

         In a nutshell:
          - Find 3 pairs of {y,E} such as {x,Q2} is constant
          - Calculate the total d2sigma/dxdy For each pair
          - Write down a system of 3 equations:
            (d2sig/dxdy) / (G2*M*E/pi) = (xy^2)       * F1 + 
                                         (1-y-Mxy/2E) * F2 +/- 
                                         y(1-y/2)     * xF3
          - Solve the system to determine F1,F2,xF3 at the input x,Q2
 
\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  May 30, 2009

\cpright  Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _EXTRACT_SF_H_
#define _EXTRACT_SF_H_

namespace genie {

class XSecAlgorithmI;
class Interaction;

namespace vld_structfunc { 
   
  void ExtractStructFunc(
     double x, double Q2, int lepton, int nucleon, double & F1, double & F2, double & xF3);
  
} // vld_structfunc
} // genie

#endif
