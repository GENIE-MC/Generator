//____________________________________________________________________________
/*!

\class    genie::NuclearData

\brief    

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  Nov 20, 2009

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          
*/
//____________________________________________________________________________

#ifndef _NUCLEAR_DATA_H_
#define _NUCLEAR_DATA_H_

namespace genie {

class Spline;

class NuclearData
{
public:
  static NuclearData * Instance (void);

  double DeuteriumSuppressionFactor(double Q2);

private:
  NuclearData();
  NuclearData(const NuclearData & nd);
 ~NuclearData();

  void Load(void); 

  static NuclearData * fInstance;

  // Loaded data
  Spline * fNuclSupprD2;

  // Sinleton cleaner
  struct Cleaner {
      void DummyMethodAndSilentCompiler() { }
      ~Cleaner() {
         if (NuclearData::fInstance !=0) {
            delete NuclearData::fInstance;
            NuclearData::fInstance = 0;
         }
      }
  };
  friend struct Cleaner;
};

}      // genie namespace
#endif //_NUCLEAR_DATA_H_

