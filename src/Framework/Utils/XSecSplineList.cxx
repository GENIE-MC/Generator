//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab
*/
//____________________________________________________________________________

#include <fenv.h>  //provides: int feenableexcept(int excepts);
#include <cmath>   //provides: std::isnan()

#include <fstream>
#include <cstdlib>

#include "libxml/parser.h"
#include "libxml/xmlmemory.h"
#include "libxml/xmlreader.h"

#include <TMath.h>
#include <TLorentzVector.h>

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Utils/XSecSplineList.h"
#include "Framework/Utils/XmlParserUtils.h"

using std::ofstream;
using std::endl;

namespace genie {

//____________________________________________________________________________
ostream & operator << (ostream & stream, const XSecSplineList & list)
{
  list.Print(stream);
  return stream;
}
//____________________________________________________________________________
XSecSplineList * XSecSplineList::fInstance = 0;
//____________________________________________________________________________
XSecSplineList::XSecSplineList()
{
  fInstance    =  0;
  fCurrentTune = "";
  fUseLogE     = true;
  fNKnots      = 100;
  fEmin        =   0.01; // GeV
  fEmax        = 100.00; // GeV
}
//____________________________________________________________________________
XSecSplineList::~XSecSplineList()
{
// Clean up.

  map<string,  map<string, Spline *> >::iterator mm_iter = fSplineMap.begin();
  for( ; mm_iter != fSplineMap.end(); ++mm_iter) {
    // loop over splines for given tune
    map<string, Spline *> & spl_map_curr_tune = mm_iter->second;
    map<string, Spline *>::iterator m_iter = spl_map_curr_tune.begin();
    for( ; m_iter != spl_map_curr_tune.end(); ++m_iter) {
      Spline * spline = m_iter->second;
      delete spline;
      spline = 0;
    }
    spl_map_curr_tune.clear();
  }
  fSplineMap.clear();
  fInstance = 0;
}
//____________________________________________________________________________
XSecSplineList * XSecSplineList::Instance()
{
  if(fInstance == 0) {
    static XSecSplineList::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new XSecSplineList;
  }
  return fInstance;
}
//____________________________________________________________________________
bool XSecSplineList::SplineExists(
            const XSecAlgorithmI * alg, const Interaction * interaction) const
{
  string key = this->BuildSplineKey(alg,interaction);
  return this->SplineExists(key);
}
//____________________________________________________________________________
bool XSecSplineList::SplineExists(string key) const
{

  if ( fCurrentTune.size() == 0 ) {
    SLOG("XSecSplLst", pERROR) << "Spline requested while CurrentTune not set" ;
    return false ;
  }

  SLOG("XSecSplLst", pDEBUG)
    << "Checking for spline: " << key << " in tune: " << fCurrentTune;

  map<string,  map<string, Spline *> >::const_iterator //
  mm_iter = fSplineMap.find(fCurrentTune);
  if(mm_iter == fSplineMap.end()) {
    SLOG("XSecSplLst", pWARN)
       << "No splines for tune " << fCurrentTune << " were found!";
    return false;
  }
  const map<string, Spline *> & spl_map_curr_tune = mm_iter->second;
  bool exists = (spl_map_curr_tune.count(key) == 1);
  SLOG("XSecSplLst", pDEBUG)
    << "Spline found?...." << utils::print::BoolAsYNString(exists);
  return exists;
}
//____________________________________________________________________________
const Spline * XSecSplineList::GetSpline(
            const XSecAlgorithmI * alg, const Interaction * interaction) const
{
  string key = this->BuildSplineKey(alg,interaction);
  return this->GetSpline(key);
}
//____________________________________________________________________________
const Spline * XSecSplineList::GetSpline(string key) const
{

  if ( fCurrentTune.size() == 0 ) {
    SLOG("XSecSplLst", pFATAL) << "Spline requested while CurrentTune not set" ;
    exit(0) ;
  }

  SLOG("XSecSplLst", pDEBUG)
    << "Getting spline: " << key << " in tune: " << fCurrentTune;

  map<string,  map<string, Spline *> >::const_iterator //\/
  mm_iter = fSplineMap.find(fCurrentTune);
  if(mm_iter == fSplineMap.end()) {
    SLOG("XSecSplLst", pWARN)
       << "No splines for tune " << fCurrentTune << " were found!";
    return 0;
  }
  const map<string, Spline *> & spl_map_curr_tune = mm_iter->second;
  map<string, Spline *>::const_iterator //\/
  m_iter = spl_map_curr_tune.find(key);
  if(m_iter == spl_map_curr_tune.end()) {
    SLOG("XSecSplLst", pWARN)
      << "Couldn't find spline: " << key << " in tune: " << fCurrentTune;
    return 0;
  }
  return m_iter->second;
}
//____________________________________________________________________________
void XSecSplineList::CreateSpline(const XSecAlgorithmI * alg,
        const Interaction * interaction, int nknots, double e_min, double e_max)
{
// Build a cross section spline for the input interaction using the input
// cross section algorithm and store in the list.
// For building this specific entry of the spline list, the user is allowed
// to override the list-wide nknots,e_min,e_max

  // FE_ALL_EXCEPT FE_INEXACT FE_UNDERFLOW
  // FE_DIVBYZERO  FE_INVALID FE_OVERFLO
  // rwh -- uncomment to catch NaN
  // feenableexcept(FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW);

  double xsec[nknots];
  double E   [nknots];

  SLOG("XSecSplLst", pNOTICE)
     << "Creating cross section spline using the algorithm: " << *alg;

  string key = this->BuildSplineKey(alg,interaction);

  // If any of the nknots,e_min,e_max was not set or its value is not acceptable
  // use the list values
  //
  if (e_min   < 0.) e_min = this->Emin();
  if (e_max   < 0.) e_max = this->Emax();
  if (nknots <= 2) nknots = this->NKnots();
  assert( e_min < e_max );

  // Distribute the knots in the energy range (e_min,e_max) :
  // - Will use 5 knots linearly spaced below the energy thresholds so that the
  //   spline behaves correctly in (e_min,Ethr)
  // - Place 1 knot exactly on the input interaction threshold
  // - Place the remaining n-6 knots spaced either linearly or logarithmically
  //   above the input interaction threshold
  // The above scheme schanges appropriately if Ethr<e_min (i.e. no knots
  // are computed below threshold)
  //
  double Ethr = interaction->PhaseSpace().Threshold();
  SLOG("XSecSplLst", pNOTICE)
    << "Energy threshold for current interaction = " << Ethr << " GeV";

  int nkb = (Ethr>e_min) ? 5 : 0; // number of knots <  threshold
  int nka = nknots-nkb;           // number of knots >= threshold

  // knots < energy threshold
  double dEb =  (Ethr>e_min) ? (Ethr - e_min) / nkb : 0;
  for(int i=0; i<nkb; i++) {
     E[i] = e_min + i*dEb;
  }
  // knots >= energy threshold
  double E0  = TMath::Max(Ethr,e_min);
  double dEa = 0;
  if(this->UseLogE())
    dEa = (TMath::Log10(e_max) - TMath::Log10(E0)) /(nka-1);
  else
    dEa = (e_max-E0) /(nka-1);

  for(int i=0; i<nka; i++) {
     if(this->UseLogE())
       E[i+nkb] = TMath::Power(10., TMath::Log10(E0) + i * dEa);
     else
       E[i+nkb] = E0 + i * dEa;
  }
  // force last point to avoid floating point cumulative slew
  E[nknots-1] = e_max;

  // Compute cross sections for the input interaction at the selected
  // set of energies
  //
  double pr_mass = interaction->InitStatePtr()->Probe()->Mass();
  for (int i = 0; i < nknots; i++) {
    TLorentzVector p4(0,0,E[i],E[i]);
    if (pr_mass > 0.) {
      double pz = TMath::Max(0.,E[i]*E[i] - pr_mass*pr_mass);
      pz = TMath::Sqrt(pz);
      p4.SetPz(pz);
    }
    interaction->InitStatePtr()->SetProbeP4(p4);
    xsec[i] = alg->Integral(interaction);
    SLOG("XSecSplLst", pNOTICE)
                       << "xsec(E = " << E[i] << ") =  "
                       << (1E+38/units::cm2)*xsec[i] << " x 1E-38 cm^2";
    if ( std::isnan(xsec[i]) ) {
      // this sometimes happens near threshold, warn and move on
      SLOG("XSecSplLst", pWARN)
                       << "xsec(E = " << E[i] << ") =  "
                       << (1E+38/units::cm2)*xsec[i] << " x 1E-38 cm^2"
                       << " : converting NaN to 0.0";
      xsec[i] = 0.0;
    }

  }

  // Warn about odd case of decreasing cross section
  //    but allow for small variation due to integration errors
  const double eps_xsec = 1.0e-5;
  const double xsec_scale = (1.0-eps_xsec);
  if ( xsec[nknots-1] < xsec[nknots-2]*xsec_scale ) {
    SLOG("XSecSplLst", pWARN)
      << "Last point oddity: " << key <<  " has "
      << " xsec[nknots-1] " << xsec[nknots-1] << " < "
      << " xsec[nknots-2] " << xsec[nknots-2];
  }

  // Build
  //
  Spline * spline = new Spline(nknots, E, xsec);

  // Save
  //
  map<string,  map<string, Spline *> >::iterator //\/
  mm_iter = fSplineMap.find(fCurrentTune);
  if(mm_iter == fSplineMap.end()) {
    map<string, Spline *> spl_map_curr_tune;
    fSplineMap.insert( map<string, map<string, Spline *> >::value_type(
      fCurrentTune, spl_map_curr_tune) );
    mm_iter = fSplineMap.find(fCurrentTune);
  }
  map<string, Spline *> & spl_map_curr_tune = mm_iter->second;
  spl_map_curr_tune.insert( map<string, Spline *>::value_type(key, spline) );
}
//____________________________________________________________________________
int XSecSplineList::NSplines(void) const
{
  map<string,  map<string, Spline *> >::const_iterator //
  mm_iter = fSplineMap.find(fCurrentTune);
  if(mm_iter == fSplineMap.end()) {
    SLOG("XSecSplLst", pWARN)
       << "No splines for tune " << fCurrentTune << " were found!";
    return 0;
  }
  const map<string, Spline *> & spl_map_curr_tune = mm_iter->second;
  return (int) spl_map_curr_tune.size();
}
//____________________________________________________________________________
bool XSecSplineList::IsEmpty(void) const
{
  int n = this->NSplines();
  return (n == 0);
}
//____________________________________________________________________________
void XSecSplineList::SetLogE(bool on)
{
  fUseLogE = on;
}
//____________________________________________________________________________
void XSecSplineList::SetNKnots(int nk)
{
  fNKnots = nk;
  if(fNKnots<10) fNKnots = 10; // minimum acceptable number of knots
}
//____________________________________________________________________________
void XSecSplineList::SetMinE(double Ev)
{
  if(Ev>0) fEmin = Ev;
}
//____________________________________________________________________________
void XSecSplineList::SetMaxE(double Ev)
{
  if(Ev>0) fEmax = Ev;
}
//____________________________________________________________________________
void XSecSplineList::SaveAsXml(const string & filename, bool save_init) const
{
//! Save XSecSplineList to XML file

  SLOG("XSecSplLst", pNOTICE)
       << "Saving XSecSplineList as XML in file: " << filename;

  ofstream outxml(filename.c_str());
  if(!outxml.is_open()) {
    SLOG("XSecSplLst", pERROR) << "Couldn't create file = " << filename;
    return;
  }
  outxml << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>";
  outxml << endl << endl;
  outxml << "<!-- generated by genie::XSecSplineList::SaveSplineList() -->";
  outxml << endl << endl;

  int uselog = (fUseLogE ? 1 : 0);
  outxml << "<genie_xsec_spline_list "
                << "version=\"3.00\" uselog=\"" << uselog << "\">";
  outxml << endl << endl;

  // loop over tunes
  map<string,  map<string, Spline *> >::const_iterator //\/
  mm_iter = fSplineMap.begin();
  for( ; mm_iter != fSplineMap.end(); ++mm_iter) {

    string tune_name = mm_iter->first;
    outxml << "  <genie_tune name=\"" << tune_name << "\">";
    outxml << endl << endl;

    // loop over splines for given tune
    const map<string, Spline *> & spl_map_curr_tune = mm_iter->second;
    map<string, Spline *>::const_iterator //\/
    m_iter = spl_map_curr_tune.begin();
    for( ; m_iter != spl_map_curr_tune.end(); ++m_iter) {
      string key = m_iter->first;

      // If current spline is from the initial loaded set,
      // look-up input option to decide whether to write out in
      // new output file or not
      bool from_init_set = false;
      map<string, set<string> >::const_iterator //\/
      it = fLoadedSplineSet.find(tune_name);
      if(it != fLoadedSplineSet.end()) {
         const set<string> & init_set_curr_tune = it->second;
         from_init_set = (init_set_curr_tune.count(key) == 1);
      }
      if(from_init_set && !save_init) continue;

      // Add current spline to output file
      Spline * spline = m_iter->second;
      spline->SaveAsXml(outxml,"E","xsec", key);
    }//spline loop

    outxml << "  </genie_tune>" << endl;
  }//tune loop

  outxml << "</genie_xsec_spline_list>" << endl;

  outxml.close();
}
//____________________________________________________________________________
XmlParserStatus_t XSecSplineList::LoadFromXml(const string & filename, bool keep)
{
//! Load XSecSplineList from ROOT file. If keep = true, then the loaded splines
//! are added to the existing list. If false, then the existing list is reset
//! before loading the splines.

  SLOG("XSecSplLst", pNOTICE)
    << "Loading splines from: " << filename;
  SLOG("XSecSplLst", pINFO)
    << "Option to keep pre-existing splines is switched "
    << ( (keep) ? "ON" : "OFF" );

  if(!keep) fSplineMap.clear();

  const int kNodeTypeStartElement = 1;
  const int kNodeTypeEndElement   = 15;
  const int kKnotX                = 0;
  const int kKnotY                = 1;

  xmlTextReaderPtr reader;

  int ret = 0, val_type = -1, iknot = 0, nknots = 0;
  double * E = 0, * xsec = 0;
  string spline_name = "";
  string temp_tune ;

  reader = xmlNewTextReaderFilename(filename.c_str());
  if (reader != NULL) {
        ret = xmlTextReaderRead(reader);
        while (ret == 1) {
            xmlChar * name  = xmlTextReaderName     (reader);
            xmlChar * value = xmlTextReaderValue    (reader);
            int       type  = xmlTextReaderNodeType (reader);
            int       depth = xmlTextReaderDepth    (reader);

            if(depth==0 && type==kNodeTypeStartElement) {
               LOG("XSecSplLst", pDEBUG) << "Root element = " << name;
               if(xmlStrcmp(name, (const xmlChar *) "genie_xsec_spline_list")) {
                   LOG("XSecSplLst", pERROR)
                     << "\nXML doc. has invalid root element! [filename: " << filename << "]";
                   return kXmlInvalidRoot;
               }

               xmlChar * xvrs   = xmlTextReaderGetAttribute(reader,(const xmlChar*)"version");
               xmlChar * xinlog = xmlTextReaderGetAttribute(reader,(const xmlChar*)"uselog");
               string svrs      = utils::str::TrimSpaces((const char *)xvrs);
               string sinlog    = utils::str::TrimSpaces((const char *)xinlog);

               LOG("XSecSplLst", pNOTICE)
                   << "Input x-section spline XML file format version: " << svrs;

               if (atoi(sinlog.c_str()) == 1) this->SetLogE(true);
               else this->SetLogE(false);

               xmlFree(xvrs);
               xmlFree(xinlog);
            }

            if( (!xmlStrcmp(name, (const xmlChar *) "genie_tune")) && type==kNodeTypeStartElement) {
               xmlChar * xtune = xmlTextReaderGetAttribute(reader,(const xmlChar*)"name");
               temp_tune    = utils::str::TrimSpaces((const char *)xtune);
               SLOG("XSecSplLst", pNOTICE) << "Loading x-section splines for GENIE tune: " << temp_tune;
               xmlFree(xtune);
            }

            if( (!xmlStrcmp(name, (const xmlChar *) "spline")) && type==kNodeTypeStartElement) {
               xmlChar * xname = xmlTextReaderGetAttribute(reader,(const xmlChar*)"name");
               xmlChar * xnkn  = xmlTextReaderGetAttribute(reader,(const xmlChar*)"nknots");
               string sname    = utils::str::TrimSpaces((const char *)xname);
               string snkn     = utils::str::TrimSpaces((const char *)xnkn);

               spline_name = sname;
               SLOG("XSecSplLst", pNOTICE) << "Loading spline: " << spline_name;

               nknots = atoi( snkn.c_str() );
               iknot=0;
               E     = new double[nknots];
               xsec  = new double[nknots];

               xmlFree(xname);
               xmlFree(xnkn);
            }

            if( (!xmlStrcmp(name, (const xmlChar *) "E"))    && type==kNodeTypeStartElement) { val_type = kKnotX; }
            if( (!xmlStrcmp(name, (const xmlChar *) "xsec")) && type==kNodeTypeStartElement) { val_type = kKnotY; }

            if( (!xmlStrcmp(name, (const xmlChar *) "#text")) && depth==5) {
                if      (val_type==kKnotX) E   [iknot] = atof((const char *)value);
                else if (val_type==kKnotY) xsec[iknot] = atof((const char *)value);
            }
            if( (!xmlStrcmp(name, (const xmlChar *) "knot")) && type==kNodeTypeEndElement) {
               iknot++;
            }
            if( (!xmlStrcmp(name, (const xmlChar *) "spline")) && type==kNodeTypeEndElement) {
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
               LOG("XSecSplLst", pINFO) << "Done with current spline";
               for(int i=0; i<nknots; i++) {
                  LOG("XSecSplLst", pINFO) << "xsec[E = " << E[i] << "] = " << xsec[i];
               }
#endif
               // done looping over knots - build the spline
               Spline * spline = new Spline(nknots, E, xsec);
               delete [] E;
               delete [] xsec;

               // insert the spline to the map
               map<string,  map<string, Spline *> >::iterator //\/
               mm_iter = fSplineMap.find( temp_tune );
               if(mm_iter == fSplineMap.end()) {
                 map<string, Spline *> spl_map_curr_tune;
                 fSplineMap.insert( map<string, map<string, Spline *> >::value_type(
                    temp_tune, spl_map_curr_tune) );
                 mm_iter = fSplineMap.find( temp_tune );
               }
               map<string, Spline *> & spl_map_curr_tune = mm_iter->second;
               spl_map_curr_tune.insert(
                  map<string, Spline *>::value_type(spline_name, spline) );
               fLoadedSplineSet[temp_tune].insert(spline_name);
            }
            xmlFree(name);
            xmlFree(value);
            ret = xmlTextReaderRead(reader);
        }
        xmlFreeTextReader(reader);
        if (ret != 0) {
          LOG("XSecSplLst", pERROR)
            << "\nXML file could not be parsed! [filename: " << filename << "]";
          return kXmlNotParsed;
        }
  } else {
    LOG("XSecSplLst", pERROR)
          << "\nXML file could not be found! [filename: " << filename << "]";
  }

  return kXmlOK;
}
//____________________________________________________________________________
string XSecSplineList::BuildSplineKey(
            const XSecAlgorithmI * alg, const Interaction * interaction) const
{
  if(!alg) {
    LOG("XSecSplLst", pWARN)
            << "Null XSecAlgorithmI - Returning empty spline key";
    return "";
  }

  if(!interaction) {
    LOG("XSecSplLst", pWARN)
            << "Null Interaction - Returning empty spline key";
    return "";
  }

  string alg_name  = alg->Id().Name();
  string param_set = alg->Id().Config();
  string intkey    = interaction->AsString();

  string key = alg_name + "/" + param_set + "/" + intkey;

  return key;
}
//____________________________________________________________________________
const vector<string> * XSecSplineList::GetSplineKeys(void) const
{
  map<string,  map<string, Spline *> >::const_iterator //\/
  mm_iter = fSplineMap.find(fCurrentTune);
  if(mm_iter == fSplineMap.end()) {
    SLOG("XSecSplLst", pWARN)
       << "No splines for tune " << fCurrentTune << " were found!";
    return 0;
  }
  const map<string, Spline *> & spl_map_curr_tune = mm_iter->second;
  vector<string> * keyv = new vector<string>(spl_map_curr_tune.size());
  unsigned int i=0;
  map<string, Spline *>::const_iterator m_iter = spl_map_curr_tune.begin();
  for( ; m_iter != spl_map_curr_tune.end(); ++m_iter) {
    string key = m_iter->first;
    (*keyv)[i++]=key;
  }
  return keyv;
}
//____________________________________________________________________________
void XSecSplineList::Print(ostream & stream) const
{
  stream << "\n ******************* XSecSplineList *************************";
  stream << "\n [-] Options:";
  stream << "\n  |";
  stream << "\n  |-----o  UseLogE..................." << fUseLogE;
  stream << "\n  |-----o  Spline Emin..............." << fNKnots;
  stream << "\n  |-----o  Spline Emax..............." << fEmin;
  stream << "\n  |-----o  Spline NKnots............." << fEmax;
  stream << "\n  |";

  map<string, map<string, Spline *> >::const_iterator mm_iter;
  for(mm_iter = fSplineMap.begin(); mm_iter != fSplineMap.end(); ++mm_iter) {

    string curr_tune = mm_iter->first;
    stream << "\n [-] Available x-section splines for tune: " << curr_tune ;
    stream << "\n  |";

    const map<string, Spline *> & spl_map_curr_tune = mm_iter->second;
    map<string, Spline *>::const_iterator m_iter = spl_map_curr_tune.begin();
    for( ; m_iter != spl_map_curr_tune.end(); ++m_iter) {
      string key = m_iter->first;
      stream << "\n  |-----o  " << key;
    }
    stream << "\n";
  }
}
//___________________________________________________________________________

} // genie namespace
