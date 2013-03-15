//____________________________________________________________________________
/*!

\class    genie::nuvld::XmlExperimentInfo

\brief

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Aug, 2003

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _EXPERIMENT_INFO_H_
#define _EXPERIMENT_INFO_H_

#include <string>
#include <iostream>

using std::string;
using std::ostream;
using std::cout;
using std::cerr;
using std::endl;

namespace genie {
namespace nuvld {
  
const int c_exp_info_ntags = 7;

const string c_exp_info_tag[c_exp_info_ntags] = {
   string("facility"),
   string("detector"),
   string("target"),
   string("exposure"),
   string("year"),
   string("beam"),
   string("E")
};

class XmlExperimentInfo
{
public:

  XmlExperimentInfo();
  XmlExperimentInfo(const XmlExperimentInfo & info);
  virtual ~XmlExperimentInfo() { }

  void Add(string key, string value);

  virtual const string Name          (void) const { return _name;           }
  virtual const string Comment       (void) const { return _comment;        }
  virtual const string Facility      (void) const { return _facility;       }
  virtual const string Detector      (void) const { return _detector;       }
  virtual const string Beam          (void) const { return _beam;           }
  virtual const string Target        (void) const { return _target;         }
  virtual const string YearBegin     (void) const { return _year_begin;     }
  virtual const string YearEnd       (void) const { return _year_end;       }
  virtual const string Exposure      (void) const { return _exposure;       }
  virtual const string ExposureUnits (void) const { return _exposure_units; }
  virtual const string EnergyMin     (void) const { return _energy_min;     }
  virtual const string EnergyMax     (void) const { return _energy_max;     }
  virtual const string EnergyUnits   (void) const { return _energy_units;   }
  virtual const string EnergyFrame   (void) const { return _energy_frame;   }

  friend ostream & operator <<(ostream & stream, const XmlExperimentInfo & info);

protected:

  string  _name;
  string  _comment;
  string  _facility;
  string  _detector;
  string  _beam;
  string  _target;
  string  _exposure;
  string  _exposure_units;
  string  _year_begin;
  string  _year_end;
  string  _energy_min;
  string  _energy_max;
  string  _energy_units;   
  string  _energy_frame;    

private:

  void ProcessYearField   (string value);
  void ProcessEnergyField (string value);
};

} // nuvld namespace
} // genie namespace

#endif // _EXPERIMENT_INFO_H_
