//____________________________________________________________________________
/*!

\class    genie::NaturalIsotopes

\brief    Singleton class to load & serve tables of natural occurring isotopes

\author   Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

	        Jim Dobson <j.dobson07@imperial.ac.uk>
          Imperial College London

\created  May 30, 2008

\cpright   Copyright (c) 2003-2020, The GENIE Collaboration
           For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _NATURAL_ISOTOPES_H_
#define _NATURAL_ISOTOPES_H_

#include <map>
#include <vector>

using std::map;
using std::vector;

namespace genie {

class NaturalIsotopeElementData;
class NaturalIsotopes
{
public:
  static NaturalIsotopes * Instance (void);

  int NElements(int Z) const;
  const NaturalIsotopeElementData * ElementData (int Z, int ielement) const;

private:
  NaturalIsotopes();
  NaturalIsotopes(const NaturalIsotopes &);
  virtual ~NaturalIsotopes();

  bool LoadTable(void);

  static NaturalIsotopes * fInstance;

  map<int, vector<NaturalIsotopeElementData*> > fNaturalIsotopesTable;

  struct Cleaner {
      void DummyMethodAndSilentCompiler() { }
      ~Cleaner() {
         if (NaturalIsotopes::fInstance !=0) {
            delete NaturalIsotopes::fInstance;
            NaturalIsotopes::fInstance = 0;
         }
      }
  };
  friend struct Cleaner;
};

class NaturalIsotopeElementData {
public:
  NaturalIsotopeElementData()                           : fPdgCode(0),    fAbundance(0)         { }
  NaturalIsotopeElementData(int code, double abundance) : fPdgCode(code), fAbundance(abundance) { }
 ~NaturalIsotopeElementData() { }

  int    PdgCode   (void) const { return fPdgCode;   }
  double Abundance (void) const { return fAbundance; }

private:
  int    fPdgCode;
  double fAbundance;
};

}      // genie namespace
#endif // _NATURAL_ISOTOPES_H_
