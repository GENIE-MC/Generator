//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool

 Jim Dobson <j.dobson07@imperial.ac.uk>
 Imperial College London
*/
//____________________________________________________________________________

#include <fstream>
#include <string>

#include <TSystem.h>

#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/NaturalIsotopes.h"

using std::string;
using std::cout;
using std::endl;

using namespace genie;

//____________________________________________________________________________
NaturalIsotopes * NaturalIsotopes::fInstance = 0;
//____________________________________________________________________________
NaturalIsotopes::NaturalIsotopes()
{
  if( ! this->LoadTable() ) {
     LOG("NatIsotop", pERROR) << "NaturalIsotopes initialization failed!";
  }
  fInstance =  0;
}
//____________________________________________________________________________
NaturalIsotopes::~NaturalIsotopes()
{
  cout << "NaturalIsotopes singleton dtor: "
       << "Deleting natural isotope data tables" << endl;

  map<int, vector<NaturalIsotopeElementData*> >::iterator miter;
  vector<NaturalIsotopeElementData*>::iterator viter;

  for(miter = fNaturalIsotopesTable.begin();
                 miter != fNaturalIsotopesTable.end(); ++miter) {
     vector<NaturalIsotopeElementData*> vec = miter->second;
     for(viter = vec.begin(); viter != vec.end(); ++viter) {
	NaturalIsotopeElementData * element_data = *viter;
        if(element_data) {
	    delete element_data;
	    element_data = 0;
        }
     }
     vec.clear();
  }
  fNaturalIsotopesTable.clear();
  fInstance = 0;
}
//____________________________________________________________________________
NaturalIsotopes * NaturalIsotopes::Instance()
{
  if(fInstance == 0) {
    LOG("NatIsotop", pINFO) << "NaturalIsotopes late initialization";

    static NaturalIsotopes::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();

    fInstance = new NaturalIsotopes;
  }
  return fInstance;
}
//____________________________________________________________________________
int NaturalIsotopes::NElements(int Z) const
{
  map<int, vector<NaturalIsotopeElementData*> >::const_iterator miter;

  if( (miter=fNaturalIsotopesTable.find(Z)) == fNaturalIsotopesTable.end()) {
    LOG("NatIsotop", pWARN)
       << "Table has no elements for natural isotope  Z = " << Z;
    return 0;
  }
  vector<NaturalIsotopeElementData*> vec = miter->second;
  return vec.size();
}
//____________________________________________________________________________
const NaturalIsotopeElementData *
         NaturalIsotopes::ElementData(int Z, int ielement) const
{
  map<int, vector<NaturalIsotopeElementData*> >::const_iterator miter;

  if( (miter=fNaturalIsotopesTable.find(Z)) == fNaturalIsotopesTable.end()) {
    LOG("NatIsotop", pWARN)
       << "Table has no elements for natural isotope  Z = " << Z;
    return 0;
  }
  vector<NaturalIsotopeElementData*> vec = miter->second;
  if(ielement >= (int)vec.size() || ielement < 0) {
    LOG("NatIsotop", pWARN)
      << "Natural isotope Z = " << Z << " has " << vec.size() << " elements"
      << " (element = " << ielement << " was requested)";
    return 0;
  }
  return vec[ielement];
}
//____________________________________________________________________________
const NaturalIsotopeElementData *
         NaturalIsotopes::ElementDataPdg(int Z, int pdgcode) const
{
  map<int, vector<NaturalIsotopeElementData*> >::const_iterator miter;

  if( (miter=fNaturalIsotopesTable.find(Z)) == fNaturalIsotopesTable.end()) {
    LOG("NatIsotop", pWARN)  
       << "Table has no elements for natural isotope  Z = " << Z;
    return 0;
  }
  
  vector<NaturalIsotopeElementData*> vec = miter->second;
  for (int i; i<vec.size(); i++) {
    if (vec[i]->PdgCode()==pdgcode) return vec[i];
  }

  LOG("NatIsotop", pWARN) 
    << "Natural isotope Z = " << Z << " has " << vec.size() << " elements"
    << " (pdgcode = " << pdgcode << " was requested)";
  return 0;

}
//____________________________________________________________________________
bool NaturalIsotopes::LoadTable(void)
{
  // get the natural isotopes table filename
  string filename = string(gSystem->Getenv("GENIE")) +
                    string("/data/evgen/catalogues/iso/natural-isotopes.data");

  LOG("NatIsotop", pINFO)
    << "Loading natural occurring isotope table from file: " << filename;

  bool is_accessible = ! (gSystem->AccessPathName( filename.c_str() ));
  if (!is_accessible) {
     LOG("NatIsotop", pWARN) << "Can not read file: " << filename;
     return false;
  }

  // load the natural isotopes .txt file
  string input_buf;
  std::ifstream input(filename.c_str());
  if (input.is_open()){

    //skip first 8 lines (comments)
    for(int i=0; i<8; i++){
      string buffer;
      getline(input, buffer);
    }

    int Z = -1, Z_previous = -1, nelements = 0, pdgcode = 0;
    double atomicmass = 0, abundance = 0;
    string elementname, subelementname;

    while( !input.eof() ) {

      //read in naturally occuring element info
      input >> Z;
      input >> elementname;
      input >> nelements;

      vector<NaturalIsotopeElementData *> vec;
      NaturalIsotopeElementData * data = 0;

      // check not re-reading same element
      if(Z!=Z_previous){
	LOG("NatIsotop", pDEBUG) << "Reading entry for Z = " << Z;
	for(int n=0 ; n < nelements; n++){
	  input >> subelementname;
	  input >> pdgcode;
	  input >> atomicmass;
	  input >> abundance;
  	  LOG("NatIsotop", pDEBUG)
            << " - Element: " << n << ", pdg = " << pdgcode
            << ", A = " << atomicmass << ", abundance = " << abundance;
          data = new NaturalIsotopeElementData(pdgcode, abundance,atomicmass);
  	  vec.push_back(data);
	}
	fNaturalIsotopesTable.insert(
           map<int,vector<NaturalIsotopeElementData*> >::value_type(Z,vec));
      }
      Z_previous = Z;
    } //!eof

  } else {
    return false;
  } //open?

  return true;
}
//____________________________________________________________________________
