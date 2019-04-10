//__________________________________________________________________________
/*!

\class    GSimFiles

\brief    Holds GENIE simulation outputs (cross-section ROOT files, simulated 
          event samples in GHEP, GST or other format) typically used as inputs
          in physics validation / tuning apps.

          The file lists are stored in XML format as shown below.

          Simulation results from multiple version of the code or from multiple
          physics models can be stored.

          <?xml version="1.0" encoding="ISO-8859-1"?>
          <genie_simulation_outputs>
            <model name="a_model_name">
               <xsec_file>             /path/model_1/xsec.root     </xsec_file>
               <evt_file format="gst"> /path/model_1/evtfile0.root </evt_file>
               <evt_file format="gst"> /path/model_1/evtfile1.root </evt_file>
               <evt_file format="gst"> /path/model_1/evtfile2.root </evt_file>
               ...
            </model>
            <model name="another_model_name">
               <xsec_file>             /path/model_2/xsec.root     </xsec_file>
               <evt_file format="gst"> /path/model_2/evtfile0.root </evt_file>
               <evt_file format="gst"> /path/model_2/evtfile1.root </evt_file>
               <evt_file format="gst"> /path/model_2/evtfile2.root </evt_file>
               ...
             </model>
             ...
          </genie_simulation_outputs>

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Laboratory

\created  Oct 12, 2009

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//__________________________________________________________________________

#ifndef _GSIM_FILES_H_
#define _GSIM_FILES_H_

#include <iostream>
#include <string>
#include <vector>

#include <TChain.h>
#include <TTree.h>
#include <TFile.h>

using std::ostream;
using std::string;
using std::vector;

namespace genie {

class GSimFiles;
ostream & operator << (ostream & stream, const GSimFiles & gsimf);  

class GSimFiles
{
public:
  GSimFiles(bool chain=true, const int nmaxmodels=10);
 ~GSimFiles(void);

  int              NModels       (void)             const;
  int              FindModelID   (string tag)       const;
  string           ModelTag      (int imodel)       const;
  TFile *          XSecFile      (int imodel)       const;
  string           XSecFileName  (int imodel)       const;
  TChain *         EvtChain      (int imodel)       const;
  vector<string> & EvtFileNames  (int imodel)       const;
  const string   & PathToXMLFile(void)              const;             
  void             Print         (ostream & stream) const;
  bool             LoadFromFile  (string xmlfile);

  friend ostream & operator << (ostream & stream, const GSimFiles & gsimf);  

private:

  void Init    (const int nmaxmodels);
  void CleanUp (void);

  bool                      fDoChain;
  int                       fNModels;
  vector<string>          * fModelTag;
  vector<TFile*>          * fXSecFile;
  vector<string>          * fXSecFileName;
  vector<TChain*>         * fEvtChain;
  vector<vector<string> > * fEvtFileNames;
  string                    fPath2XMLFile;
};

}      // genie namespace

#endif // _GSIM_FILES_H_

