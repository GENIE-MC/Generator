//____________________________________________________________________________
/*!

\class    genie::nuvld::XmlMeasurementHeader

\brief    NuVld DB measurement header

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Aug, 2003

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _MEASUREMENT_HEADER_H_
#define _MEASUREMENT_HEADER_H_

#include <string>
#include <vector>
#include <iostream>

#include "ValidationTools/NuVld/XmlCitation.h"

using std::string;
using std::ostream;
using std::endl;
using std::cerr;
using std::vector;

namespace genie {
namespace nuvld {

const int c_meas_header_ntags = 9;
const string c_meas_header_tag[c_meas_header_ntags] = {
   string("observable"),
   string("target"),
   string("reaction"),
   string("A"),          
   string("exposure"),
   string("exposure_units"),
   string("data_source"),
   string("npoints"),
   string("err_status"),
};

class XmlMeasurementHeader
{
public:

  XmlMeasurementHeader();
  XmlMeasurementHeader(const XmlMeasurementHeader & header);
  virtual ~XmlMeasurementHeader() { }

  void Add(string key, string value = "");
  void Add(XmlCitation * ref);

  virtual const string Tag           (void) const { return fTag;            }
  virtual const string Observable    (void) const { return fObservable;     }
  virtual const string Target        (void) const { return fTarget;         }
  virtual const string Reaction      (void) const { return fReaction;       }
  virtual const string A             (void) const { return fA;              }
  virtual const string Exposure      (void) const { return fExposure;       }
  virtual const string ExposureUnits (void) const { return fExposureUnits;  }
  virtual const string DataSource    (void) const { return fDataSource;     }
  virtual const string NPoints       (void) const { return fNpoints;        }
  virtual const string ErrorStatus   (void) const { return fErrStatus;      }
  virtual const string Comment       (void) const { return fComment;        }

  const vector<XmlCitation *> & GetRefs(void) const { return * fRefs; }

  friend ostream & operator << (ostream & stream, const XmlMeasurementHeader & header);

protected:

  vector<XmlCitation *> * fRefs;

  string  fTag;
  string  fObservable;
  string  fTarget;
  string  fReaction;
  string  fA;
  string  fExposure;
  string  fExposureUnits;
  string  fDataSource;
  string  fNpoints;
  string  fErrStatus;
  string  fComment;
};

} // nuvld namespace
} // genie namespace

#endif // _MEASUREMENT_HEADER_H_
