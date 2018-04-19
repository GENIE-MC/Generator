//____________________________________________________________________________
/*!

\class    genie::TuneId

\brief    GENIE tune ID

\author   Marco Roda <Marco.Roda \at liverpool.ac.uk>
          University of Liverpool

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  April 19, 2018

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _TUNE_ID_H_
#define _TUNE_ID_H_

#include <string>
#include <iostream>

#include <TObject.h>

using std::string;
using std::ostream;

namespace genie {

class TuneId;
ostream & operator << (ostream & stream, const TuneId & id);
bool      operator == (const TuneId & id1, const TuneId & id2);
bool      operator != (const TuneId & id1, const TuneId & id2);

class TuneId: public TObject {

public:
  using TObject::Copy;    // Suppress clang 'hides overloaded virtual function [-Woverloaded-virtual]' warnings
  using TObject::Compare; // Ditto
  using TObject::Print;   // Ditto

  TuneId();
  TuneId(string id_str);
  TuneId(const TuneId & id);
 ~TuneId();

  // The typical tune name in neutrino mode is : Gdd_MMv_PP_xxx
  string Name            (void) const { return fName;                         } // Gdd_MMv_PP_xxx
  string Prefix          (void) const { return fPrefix;                       } // G
  string Year            (void) const { return fYear;                         } // dd
  string ModelId         (void) const { return fMajorModelId + fMinorModelId; } // MMv
  string MajorModelId    (void) const { return fMajorModelId;                 } // MM
  string MinorModelId    (void) const { return fMinorModelId;                 } // v
  string TunedParamSetId (void) const { return fTunedParamSetId;              } // PP
  string FitDataSetId    (void) const { return fFitDataSetId;                 } // xxx
 
  void   Decode  (string id_str);
  void   Copy    (const TuneId & id);
  bool   Compare (const TuneId & id) const;
  void   Print   (ostream & stream) const;

  friend ostream & operator << (ostream & stream, const TuneId & id);

private:

  void Init (void);

  string fName;
  string fPrefix;
  string fYear;
  string fModelId;
  string fMajorModelId;
  string fMinorModelId;
  string fTunedParamSetId; 
  string fFitDataSetId;
};

}       // genie namespace

#endif  // _TUNE_ID_H_
