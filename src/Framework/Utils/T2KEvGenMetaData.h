//____________________________________________________________________________
/*!

\class    genie::utils::T2KEvGenMetaData

\brief    Utility class to store MC job meta-data

\author   Jim Dobson
          Imperial College London

\created  Mar 04, 2010

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _T2KEVGEN_METADATA_H_
#define _T2KEVGEN_METADATA_H_

#include <iostream>
#include <string>
#include <map>

#include <TObject.h>
#include <TH1D.h>

using std::ostream;
using std::string;
using std::map;

namespace genie {
namespace utils {

class T2KEvGenMetaData;
ostream & operator << (ostream & stream, const T2KEvGenMetaData & md);

class T2KEvGenMetaData: public TObject
{
public:
  using TObject::Print; // suppress clang 'hides overloaded virtual function [-Woverloaded-virtual]' warnings

  T2KEvGenMetaData() :
      jnubeam_version(""), 
      jnubeam_file(""), 
      detector_location(""), 
      geom_file(""), 
      geom_top_volume(""), 
      geom_length_units(1.), 
      geom_density_units(1.),
      using_root_geom(false), 
      using_hist_flux(false) 
  { 
  }

  ~T2KEvGenMetaData() 
  { 
  }

  void Print(ostream & stream) const;

  friend ostream & operator << (ostream & stream, const T2KEvGenMetaData & md);

  string           jnubeam_version;
  string           jnubeam_file;
  string           detector_location;
  string           geom_file;
  string           geom_top_volume;
  double           geom_length_units;
  double           geom_density_units;
  bool             using_root_geom;
  bool             using_hist_flux;
  map<int, double> target_mix;
  map<int, TH1D*>  flux_hists;

  ClassDef(T2KEvGenMetaData,1)

};


} // utils namespace
} // genie namespace

#endif // _T2KEVGEN_METADATA_H_
