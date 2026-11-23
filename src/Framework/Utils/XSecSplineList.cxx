//____________________________________________________________________________
/*
 Copyright (c) 2003-2026, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool
*/
//____________________________________________________________________________

#include <fenv.h>  //provides: int feenableexcept(int excepts);
#include <cmath>   //provides: std::isnan()
#include <string>
#include <vector>

#include <algorithm>
#include <chrono>

#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <sstream>

#include "libxml/parser.h"
#include "libxml/xmlmemory.h"
#include "libxml/xmlreader.h"

#include <TMath.h>
#include <TLorentzVector.h>

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Interaction/KPhaseSpaceCuts.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Utils/XSecSplineList.h"
#include "Framework/Utils/XmlParserUtils.h"

using std::ofstream;
using std::endl;

namespace genie {

  using namespace std::chrono ;

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

  map<string, SplineMap>::const_iterator //
  mm_iter = fSplineMap.find(fCurrentTune);
  if(mm_iter == fSplineMap.end()) {
    SLOG("XSecSplLst", pWARN)
       << "No splines for tune " << fCurrentTune << " were found!";
    return false;
  }
  const SplineMap & spl_map_curr_tune = mm_iter->second;
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
  const SplineEntry * entry = this->FindSpline(key);
  if(!entry) return 0;
  if(fValidateQ2Cuts && !entry->validated) this->ValidateQ2Cut(key, interaction, *entry);
  return entry->spline.get();
}
//____________________________________________________________________________
const Spline * XSecSplineList::GetSpline(string key) const
{
  const SplineEntry * entry = this->FindSpline(key);
  return entry ? entry->spline.get() : 0;
}
//____________________________________________________________________________
const XSecSplineList::SplineEntry * XSecSplineList::FindSpline(const string & key) const
{
  if(fCurrentTune.empty()) {
    SLOG("XSecSplLst", pFATAL) << "Spline requested while CurrentTune not set";
    exit(0);
  }

  auto tune = fSplineMap.find(fCurrentTune);
  if(tune == fSplineMap.end()) {
    SLOG("XSecSplLst", pWARN) << "No splines for tune " << fCurrentTune << " were found!";
    return 0;
  }
  auto spline = tune->second.find(key);
  if(spline == tune->second.end()) {
    SLOG("XSecSplLst", pWARN) << "Couldn't find spline: " << key << " in tune: " << fCurrentTune;
    return 0;
  }
  return &spline->second;
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

  std::vector<double> xsec( nknots, 0. ) ;
  std::vector<double> E( nknots, 0. ) ;

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

  if (Ethr>e_max) {
    SLOG("XSecSplLst", pFATAL) << "Energy threshold higher than maximum.";
    SLOG("XSecSplLst", pFATAL) << "Energy threshold = " << Ethr << " GeV";
    SLOG("XSecSplLst", pFATAL) << "Energy maximum = " << e_max << " GeV";
    return;
  }

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

    steady_clock::time_point start = steady_clock::now();

    xsec[i] = TMath::Max(alg->Integral(interaction), 0.);

    steady_clock::time_point end = steady_clock::now();

    duration<double> time_span = duration_cast<duration<double>>(end - start);

    SLOG("XSecSplLst", pNOTICE)
                       << "xsec(E = " << E[i] << ") =  "
                       << (1E+38/units::cm2)*xsec[i] << " x 1E-38 cm^2, evaluated in " << time_span.count() << " s";
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

  SplineEntry entry;
  entry.spline.reset(new Spline(nknots, E.data(), xsec.data()));
  entry.cuts.known = true;
  KPhaseSpaceCuts * cuts = KPhaseSpaceCuts::Instance();
  if(cuts->HasQ2MinCut(interaction)) {
    entry.cuts.Q2_min = cuts->Q2MinCut(interaction, 0.);
    entry.cuts.Q2_source = "CommonPhaseSpaceCuts.xml";
  }
  fSplineMap[fCurrentTune].emplace(key, std::move(entry));
}
//____________________________________________________________________________
int XSecSplineList::NSplines(void) const
{
  map<string, SplineMap>::const_iterator //
  mm_iter = fSplineMap.find(fCurrentTune);
  if(mm_iter == fSplineMap.end()) {
    SLOG("XSecSplLst", pWARN)
       << "No splines for tune " << fCurrentTune << " were found!";
    return 0;
  }
  const SplineMap & spl_map_curr_tune = mm_iter->second;
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
namespace {
  bool SameQ2Cut(double left, double right)
  {
    if(left == 0. || right == 0.) return left == right;
    double scale = std::max(1., std::max(std::fabs(left), std::fabs(right)));
    return std::fabs(left - right) <= 1E-12 * scale;
  }

  string ReadXMLAttribute(xmlTextReaderPtr reader, const char * name)
  {
    xmlChar * value = xmlTextReaderGetAttribute(reader, BAD_CAST name);
    if(!value) return "";
    string result = utils::str::TrimSpaces(reinterpret_cast<const char *>(value));
    xmlFree(value);
    return result;
  }
}
//____________________________________________________________________________
bool XSecSplineList::GetSplineQ2MinKinematics(const string & key, double & Q2min) const
{
  const SplineEntry * entry = this->FindSpline(key);
  if(!entry || !entry->cuts.known) return false;
  Q2min = entry->cuts.Q2_min;
  return true;
}
//____________________________________________________________________________
void XSecSplineList::ValidateQ2Cut(
  const string & key, const Interaction * interaction, const SplineEntry & entry) const
{
  KPhaseSpaceCuts * cuts = KPhaseSpaceCuts::Instance();
  double Q2min = cuts->HasQ2MinCut(interaction) ? cuts->Q2MinCut(interaction, 0.) : 0.;
  if(!entry.cuts.known && Q2min > 0.) {
    LOG("XSecSplLst", pFATAL)
      << "Unknown Q2 minimum cut for spline [" << key << "] in tune " << fCurrentTune
      << ". Regenerate it with per-spline Q2 metadata to verify the configured cut of "
      << Q2min << " GeV^2.";
    gAbortingInErr = true;
    exit(1);
  }
  if(entry.cuts.known && !SameQ2Cut(Q2min, entry.cuts.Q2_min)) {
    LOG("XSecSplLst", pFATAL)
      << "Q2 minimum cut mismatch for spline [" << key << "] in tune " << fCurrentTune
      << ". Configured cut: " << Q2min << " GeV^2. Spline cut: "
      << entry.cuts.Q2_min << " GeV^2. Zero denotes no configured cut.";
    gAbortingInErr = true;
    exit(1);
  }
  entry.validated = true;
}
//____________________________________________________________________________
bool XSecSplineList::CheckEMCutForMerge(const string & tune, const string & key,
  const KinematicsMetadata & cuts, const string & filename)
{
  // Interaction::AsString() includes the interaction type as ";proc:TYPE,...".
  if(!cuts.known || key.find(";proc:EM,") == string::npos) return true;

  const auto first = fMergedEMCuts.emplace(tune, EMCutOrigin{cuts.Q2_min, filename}).first;
  if(SameQ2Cut(cuts.Q2_min, first->second.Q2_min)) return true;
  LOG("XSecSplLst", pERROR)
    << "Conflicting EM Q2 minimum cuts for tune " << tune
    << ": [" << first->second.filename << "] has " << first->second.Q2_min
    << " GeV^2, while [" << filename << "] has " << cuts.Q2_min
    << " GeV^2 for spline [" << key << "]. Zero denotes no configured cut.";
  return false;
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
  map<string, SplineMap>::const_iterator //\/
  mm_iter = fSplineMap.begin();
  for( ; mm_iter != fSplineMap.end(); ++mm_iter) {

    string tune_name = mm_iter->first;
    outxml << "  <genie_tune name=\"" << tune_name << "\">";
    outxml << endl << endl;

    for(const auto & item : mm_iter->second) {
      const SplineEntry & entry = item.second;
      if(entry.loaded && !save_init) continue;

      std::ostringstream header;
      if(entry.cuts.known) {
        header << "      <kinematics>\n";
        if(entry.cuts.Q2_min > 0.) {
          xmlChar * source = xmlEncodeSpecialChars(nullptr, BAD_CAST entry.cuts.Q2_source.c_str());
          header << "        <q2 min=\"" << std::setprecision(std::numeric_limits<double>::max_digits10)
                 << entry.cuts.Q2_min << "\" unit=\"GeV^2\" source=\""
                 << (source ? reinterpret_cast<const char *>(source) : "") << "\"/>\n";
          xmlFree(source);
        } else {
          header << "        <q2 status=\"unconfigured\"/>\n";
        }
        header << "      </kinematics>\n";
      }

      entry.spline->SaveAsXml(outxml, "E", "xsec", item.first, header.str());
    }
    outxml << "  </genie_tune>" << endl;
  }//tune loop

  outxml << "</genie_xsec_spline_list>" << endl;

  outxml.close();
}
//____________________________________________________________________________
XmlParserStatus_t XSecSplineList::LoadFromXml(const string & filename, bool keep)
{
  return this->LoadFromXml(filename, keep, false);
}
//____________________________________________________________________________
XmlParserStatus_t XSecSplineList::LoadFromXml(
  const string & filename, bool keep, bool require_consistent_em_cuts)
{
// If keep is true, add the loaded splines to the existing list.

  SLOG("XSecSplLst", pNOTICE) << "Loading splines from: " << filename;
  if(!keep) {
    fSplineMap.clear();
    fMergedEMCuts.clear();
  }

  std::unique_ptr<_xmlTextReader, decltype(&xmlFreeTextReader)> reader(
    xmlNewTextReaderFilename(filename.c_str()), xmlFreeTextReader);
  if(!reader) {
    LOG("XSecSplLst", pERROR) << "Could not open spline XML file [" << filename << "]";
    return kXmlNotParsed;
  }

  string tune, key;
  vector<double> E, xsec;
  int nknots = 0;
  vector<double> * values = nullptr;
  bool in_kinematics = false;
  KinematicsMetadata cuts;

  int result = 0;
  while((result = xmlTextReaderRead(reader.get())) == 1) {
    const xmlChar * name = xmlTextReaderConstName(reader.get());
    int type = xmlTextReaderNodeType(reader.get());
    int depth = xmlTextReaderDepth(reader.get());

    if(type == XML_READER_TYPE_ELEMENT) {
      if(depth == 0) {
        if(!xmlStrEqual(name, BAD_CAST "genie_xsec_spline_list")) return kXmlInvalidRoot;
        this->SetLogE(ReadXMLAttribute(reader.get(), "uselog") == "1");
        continue;
      }
      if(depth == 1 && xmlStrEqual(name, BAD_CAST "genie_tune")) {
        tune = ReadXMLAttribute(reader.get(), "name");
        continue;
      }
      if(depth == 2 && xmlStrEqual(name, BAD_CAST "spline")) {
        key = ReadXMLAttribute(reader.get(), "name");
        nknots = atoi(ReadXMLAttribute(reader.get(), "nknots").c_str());
        E.clear();
        xsec.clear();
        values = nullptr;
        in_kinematics = false;
        cuts = KinematicsMetadata{};
        continue;
      }
      if(depth == 3 && xmlStrEqual(name, BAD_CAST "kinematics")) {
        // Tune-level stamps from older files cannot attest to any individual spline.
        in_kinematics = !key.empty() && !xmlTextReaderIsEmptyElement(reader.get());
        continue;
      }
      if(depth == 4 && in_kinematics && xmlStrEqual(name, BAD_CAST "q2")) {
        string minimum = ReadXMLAttribute(reader.get(), "min");
        string unit = ReadXMLAttribute(reader.get(), "unit");
        string status = ReadXMLAttribute(reader.get(), "status");
        char * end = nullptr;
        double Q2min = strtod(minimum.c_str(), &end);
        bool configured = !minimum.empty() && *end == '\0' && std::isfinite(Q2min)
                       && Q2min > 0. && (unit.empty() || unit == "GeV^2") && status.empty();
        bool unconfigured = minimum.empty() && status == "unconfigured";
        if(cuts.known || (!configured && !unconfigured)) {
          LOG("XSecSplLst", pERROR) << "Invalid or duplicate Q2 cut stamp for spline ["
            << key << "] in file [" << filename << "]";
          return kXmlNotParsed;
        }
        cuts.known = true;
        cuts.Q2_min = configured ? Q2min : 0.;
        cuts.Q2_source = ReadXMLAttribute(reader.get(), "source");
        continue;
      }
      if(depth == 4) {
        if(xmlStrEqual(name, BAD_CAST "E")) values = &E;
        if(xmlStrEqual(name, BAD_CAST "xsec")) values = &xsec;
      }
    } else if(type == XML_READER_TYPE_TEXT && depth == 5 && values) {
      values->push_back(atof(reinterpret_cast<const char *>(xmlTextReaderConstValue(reader.get()))));
      values = nullptr;
    } else if(type == XML_READER_TYPE_END_ELEMENT) {
      if(depth == 3 && xmlStrEqual(name, BAD_CAST "kinematics")) in_kinematics = false;
      if(depth != 2 || !xmlStrEqual(name, BAD_CAST "spline")) continue;
      if(tune.empty() || key.empty() || nknots <= 0
         || E.size() != static_cast<size_t>(nknots) || xsec.size() != E.size()) {
        LOG("XSecSplLst", pERROR) << "Invalid knots for spline [" << key
          << "] in file [" << filename << "]";
        return kXmlNotParsed;
      }

      // Compare input stamps before duplicate keys can discard a spline.
      if(require_consistent_em_cuts && !this->CheckEMCutForMerge(tune, key, cuts, filename)) {
        return kXmlNotParsed;
      }
      auto & splines = fSplineMap[tune];
      if(splines.count(key) == 0) {
        SplineEntry entry;
        entry.spline.reset(new Spline(nknots, E.data(), xsec.data()));
        entry.cuts = std::move(cuts);
        entry.loaded = true;
        splines.emplace(key, std::move(entry));
      }
      key.clear();
    }
  }

  if(result != 0) {
    LOG("XSecSplLst", pERROR) << "Could not parse spline XML file [" << filename << "]";
    return kXmlNotParsed;
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
  map<string, SplineMap>::const_iterator //\/
  mm_iter = fSplineMap.find(fCurrentTune);
  if(mm_iter == fSplineMap.end()) {
    SLOG("XSecSplLst", pWARN)
       << "No splines for tune " << fCurrentTune << " were found!";
    return 0;
  }
  const SplineMap & spl_map_curr_tune = mm_iter->second;
  vector<string> * keyv = new vector<string>(spl_map_curr_tune.size());
  unsigned int i=0;
  SplineMap::const_iterator m_iter = spl_map_curr_tune.begin();
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
  stream << "\n  |-----o  Spline Emin..............." << fEmin;
  stream << "\n  |-----o  Spline Emax..............." << fEmax;
  stream << "\n  |-----o  Spline NKnots............." << fNKnots;
  stream << "\n  |";

  map<string, SplineMap>::const_iterator mm_iter;
  for(mm_iter = fSplineMap.begin(); mm_iter != fSplineMap.end(); ++mm_iter) {

    string curr_tune = mm_iter->first;
    stream << "\n [-] Available x-section splines for tune: " << curr_tune ;
    stream << "\n  |";

    const SplineMap & spl_map_curr_tune = mm_iter->second;
    SplineMap::const_iterator m_iter = spl_map_curr_tune.begin();
    for( ; m_iter != spl_map_curr_tune.end(); ++m_iter) {
      string key = m_iter->first;
      stream << "\n  |-----o  " << key;
    }
    stream << "\n";
  }
}
//___________________________________________________________________________
void XSecSplineList::SetInterpolationType(string type)
{
  for(auto & tune : fSplineMap) {
    for(auto & entry : tune.second) {
      entry.second.spline->SetType(type);
    }
  }
}
//___________________________________________________________________________
std::vector< std::string > XSecSplineList::GetLoadedTunes() const {
  std::vector< std::string > tunes;
  for ( auto it = fSplineMap.begin(); it != fSplineMap.end(); ++it )
    tunes.push_back( it->first );
  return tunes;
}
//___________________________________________________________________________
} // genie namespace
