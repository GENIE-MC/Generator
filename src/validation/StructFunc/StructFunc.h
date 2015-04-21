//____________________________________________________________________________
/*!

\class    genie::mc_vs_data::StructFunc

\brief    A utility class to extract F1, F2 and xF3.

          Doesn't use the structure function model used for constructing the 
          DIS cross-section. Rather, it treats the generator as a black box
          an extracts the values of structure functions at, specific x, Q2, 
          from calculated cross-section values.
          
          The method is described in:        
          H.Gallagher, Nucl.Phys.Proc.Suppl.159:229-234,2006

          In a nutshell, the methodology is:
          - Find 3 pairs of {y,E} such as {x,Q2} is constant
          - Calculate the total d2sigma/dxdy For each pair
          - Write down a system of 3 equations:
             (d2sig/dxdy) / (G2*M*E/pi) = (xy^2)       * F1 + 
                                          (1-y-Mxy/2E) * F2 +/- 
                                          y(1-y/2)     * xF3
          - Solve the system to determine F1,F2,xF3 at the input x,Q2

          The obvious advantage here is that resonance contributions can be
          added to the DIS d2sig/dxdy.
 
\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

          Jelena Ilic <jelena.ilic \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Apr 25, 2012

\cpright  Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _STRUCT_FUNC_H_
#define _STRUCT_FUNC_H_

namespace genie {

class XSecAlgorithmI;

namespace mc_vs_data {

class StructFunc
{
public:
  StructFunc();
 ~StructFunc() {}

  void SetResonanceXSecModel (const XSecAlgorithmI * xsec) { fRESXSecModel      = xsec; }
  void SetDISXSecModel       (const XSecAlgorithmI * xsec) { fDISXSecModel      = xsec; }
  void SetDISCharmXSecModel  (const XSecAlgorithmI * xsec) { fDISCharmXSecModel = xsec; }

  bool ExtractF1F2xF3(
          double x, double Q2, int lepton_pdg, int nucleon_pdg, 
          double & F1, double & F2, double & xF3) const;

  double d2SigTot_dxdy(
          double x, double y, double E, int lepton_pdg, int nucleon_pdg) const;

private:

  const XSecAlgorithmI * fRESXSecModel;
  const XSecAlgorithmI * fDISXSecModel;
  const XSecAlgorithmI * fDISCharmXSecModel;
};

} // mc_vs_data namespace
} // genie namespace

#endif
