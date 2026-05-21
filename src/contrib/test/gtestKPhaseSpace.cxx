//____________________________________________________________________________
/*!

\program gtestKPhaseSpace

\brief   Program used for testing / debugging the kinematic phase space calc

\author  Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool

\created June 20, 2004

\cpright Copyright (c) 2003-2025, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         
*/
//____________________________________________________________________________

#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sys/wait.h>
#include <unistd.h>

#include "Framework/Conventions/Controls.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Interaction/KPhaseSpaceCuts.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Utils/XSecSplineList.h"

using namespace genie;

void PrintLimits(const Interaction * interaction);
int  RunDefaultChecks(const char * self);
int  RunMissingEMCutCheck(void);
int  RunWeakConfigCutCheck(void);
int  RunNoSplineQ2CutCheck(void);
bool Check(bool ok, const char * msg);
bool NearlyEqual(double a, double b);
int  RunChild(const char * self, const char * mode, const std::string & xml);
std::string WritePhaseSpaceCutsXML(const std::string & xml_body);
std::string WriteSplineXML(const std::string & xml_body);
void RemoveTempXML(const std::string & xml);

//__________________________________________________________________________
int main(int argc, char ** argv)
{
  if(argc >= 2) {
    std::string mode = argv[1];
    if(mode == "--expect-missing-em-cut") _exit(RunMissingEMCutCheck());
    if(mode == "--expect-weak-config-cut") _exit(RunWeakConfigCutCheck());
    if(mode == "--expect-no-spline-q2-cut") _exit(RunNoSplineQ2CutCheck());
  }

  // -- get a DIS interaction object & access its kinematics

  int tgt         = kPdgTgtFe56;
  int hit_nucleon = kPdgProton;
  int neutrino    = kPdgNuMu;
  double Ev       = 3;

  Interaction * qelcc = Interaction::QELCC(tgt,hit_nucleon,neutrino,Ev);
  Interaction * rescc = Interaction::RESCC(tgt,hit_nucleon,neutrino,Ev);
  Interaction * discc = Interaction::DISCC(tgt,hit_nucleon,neutrino,Ev);

  PrintLimits(qelcc);
  PrintLimits(rescc);
  PrintLimits(discc);

  return RunDefaultChecks(argv[0]);
}
//__________________________________________________________________________
void PrintLimits(const Interaction * interaction)
{
  LOG("test", pNOTICE) << *interaction;

  const KPhaseSpace & phase_space = interaction->PhaseSpace();

  Range1D_t xl  = phase_space.Limits(kKVx);
  Range1D_t yl  = phase_space.Limits(kKVy);
  Range1D_t Q2l = phase_space.Limits(kKVQ2);
  Range1D_t Wl  = phase_space.Limits(kKVW);

  LOG("test", pNOTICE) << "x  e [" << xl.min  << ", " << xl.max  << "]";
  LOG("test", pNOTICE) << "y  e [" << yl.min  << ", " << yl.max  << "]";
  LOG("test", pNOTICE) << "Q2 e [" << Q2l.min << ", " << Q2l.max << "]";
  LOG("test", pNOTICE) << "W  e [" << Wl.min  << ", " << Wl.max  << "]";
}
//__________________________________________________________________________
int RunDefaultChecks(const char * self)
{
  int failures = 0;
  const double em_q2_min = 0.02;
  const double weak_default = 0.123;

  Interaction * em = Interaction::QELEM(
    kPdgTgtFe56, kPdgProton, kPdgElectron, 3.);
  Interaction * weak = Interaction::QELCC(
    kPdgTgtFe56, kPdgProton, kPdgNuMu, 3.);

  KPhaseSpaceCuts * cuts = KPhaseSpaceCuts::Instance();

  failures += !Check(cuts->HasQ2MinCut(em),
    "EM interactions have a configured Q2 cut");
  failures += !Check(NearlyEqual(cuts->Q2MinCut(em, 0.), em_q2_min),
    "EM Q2 cut is loaded from CommonPhaseSpaceCuts.xml");
  failures += !Check(em->PhaseSpace().Limits(kKVQ2).min >= em_q2_min,
    "KPhaseSpace applies the EM Q2 cut to Q2 limits");

  failures += !Check(!cuts->HasQ2MinCut(weak),
    "weak interactions do not use a configurable Q2 cut by default");
  failures += !Check(NearlyEqual(cuts->Q2MinCut(weak, weak_default), weak_default),
    "weak interactions keep the caller default without an explicit cut");
  failures += !Check(!cuts->HasSplineQ2MinCutForProbe(kPdgNuMu),
    "neutrino spline metadata ignores EM-only Q2 cuts");
  failures += !Check(cuts->HasSplineQ2MinCutForProbe(kPdgElectron),
    "electron spline metadata sees EM Q2 cuts");
  failures += !Check(NearlyEqual(cuts->SplineQ2MinCutForProbe(kPdgElectron), em_q2_min),
    "electron spline metadata uses EM-Q2-min");

  cuts->SetQ2MinOverride(0.25);
  failures += !Check(cuts->HasQ2MinCut(weak),
    "command-line Q2 override opts weak interactions into the cut");
  failures += !Check(NearlyEqual(cuts->Q2MinCut(weak, 0.), 0.25),
    "command-line Q2 override is applied to weak interactions");
  failures += !Check(cuts->HasSplineQ2MinCutForProbe(kPdgNuMu),
    "command-line Q2 override provides spline metadata");
  failures += !Check(NearlyEqual(cuts->SplineQ2MinCutForProbe(kPdgNuMu), 0.25),
    "spline metadata uses the command-line Q2 override");

  delete em;
  delete weak;

  std::string missing_xml = WritePhaseSpaceCutsXML(
    "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
    "<common_PhaseSpaceCuts_list>\n"
    "  <param_set name=\"Default\">\n"
    "  </param_set>\n"
    "</common_PhaseSpaceCuts_list>\n");
  if(missing_xml.empty()) {
    failures += !Check(false, "temporary XML for missing EM-Q2-min was created");
  } else {
    int missing_status = RunChild(self, "--expect-missing-em-cut", missing_xml);
    failures += !Check(WIFEXITED(missing_status) &&
      WEXITSTATUS(missing_status) == 78,
      "EM interactions fail loudly when EM-Q2-min is missing");
  }
  RemoveTempXML(missing_xml);

  std::string weak_xml = WritePhaseSpaceCutsXML(
    "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
    "<common_PhaseSpaceCuts_list>\n"
    "  <param_set name=\"Default\">\n"
    "    <param type=\"double\" name=\"EM-Q2-min\"> 0.02 </param>\n"
    "    <param type=\"double\" name=\"Weak-Q2-min\"> 0.25 </param>\n"
    "  </param_set>\n"
    "</common_PhaseSpaceCuts_list>\n");
  if(weak_xml.empty()) {
    failures += !Check(false, "temporary XML for Weak-Q2-min was created");
  } else {
    int weak_status = RunChild(self, "--expect-weak-config-cut", weak_xml);
    failures += !Check(WIFEXITED(weak_status) &&
      WEXITSTATUS(weak_status) == 0,
      "Weak-Q2-min from CommonPhaseSpaceCuts.xml opts weak interactions in");
  }
  RemoveTempXML(weak_xml);

  std::string empty_xml = WritePhaseSpaceCutsXML(
    "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
    "<common_PhaseSpaceCuts_list>\n"
    "  <param_set name=\"Default\">\n"
    "  </param_set>\n"
    "</common_PhaseSpaceCuts_list>\n");
  if(empty_xml.empty()) {
    failures += !Check(false, "temporary XML for empty Q2 metadata was created");
  } else {
    int empty_status = RunChild(self, "--expect-no-spline-q2-cut", empty_xml);
    failures += !Check(WIFEXITED(empty_status) &&
      WEXITSTATUS(empty_status) == 0,
      "empty CommonPhaseSpaceCuts.xml produces no spline Q2 metadata");
  }
  RemoveTempXML(empty_xml);

  std::string spline_xml = WriteSplineXML(
    "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
    "<genie_xsec_spline_list version=\"3.00\" uselog=\"1\">\n"
    "  <genie_tune name=\"UnitTestTune\">\n"
    "    <kinematics>\n"
    "      <q2 min=\"0.25\" unit=\"GeV^2\" source=\"unit-test\"/>\n"
    "    </kinematics>\n"
    "    <spline name=\"dummy\" nknots=\"3\">\n"
    "      <knot> <E> 0.10000 </E> <xsec> 1.0000000000e-38 </xsec> </knot>\n"
    "      <knot> <E> 1.00000 </E> <xsec> 2.0000000000e-38 </xsec> </knot>\n"
    "      <knot> <E> 2.00000 </E> <xsec> 3.0000000000e-38 </xsec> </knot>\n"
    "    </spline>\n"
    "  </genie_tune>\n"
    "</genie_xsec_spline_list>\n");
  if(spline_xml.empty()) {
    failures += !Check(false, "temporary spline XML was created");
  } else {
    XSecSplineList * xspl = XSecSplineList::Instance();
    failures += !Check(xspl->LoadFromXml(spline_xml) == kXmlOK,
      "spline XML with kinematics metadata loads");
    double q2min = -1.;
    std::string unit = "";
    std::string source = "";
    failures += !Check(xspl->GetTuneQ2MinKinematics(
      "UnitTestTune", q2min, &unit, &source),
      "loaded spline XML exposes Q2 metadata");
    failures += !Check(NearlyEqual(q2min, 0.25),
      "loaded spline XML preserves Q2 metadata value");

    std::string saved_spline_xml = spline_xml + ".roundtrip.xml";
    xspl->SaveAsXml(saved_spline_xml);
    failures += !Check(xspl->LoadFromXml(saved_spline_xml) == kXmlOK,
      "saved spline XML with kinematics metadata reloads");
    q2min = -1.;
    failures += !Check(xspl->GetTuneQ2MinKinematics(
      "UnitTestTune", q2min, &unit, &source),
      "round-tripped spline XML exposes Q2 metadata");
    failures += !Check(NearlyEqual(q2min, 0.25),
      "round-tripped spline XML preserves Q2 metadata value");
    unlink(saved_spline_xml.c_str());
  }
  RemoveTempXML(spline_xml);

  return failures == 0 ? 0 : 1;
}
//__________________________________________________________________________
int RunMissingEMCutCheck(void)
{
  Interaction * em = Interaction::QELEM(
    kPdgTgtFe56, kPdgProton, kPdgElectron, 3.);
  KPhaseSpaceCuts::Instance()->Q2MinCut(em, 0.);
  delete em;
  return 1;
}
//__________________________________________________________________________
int RunWeakConfigCutCheck(void)
{
  Interaction * weak = Interaction::QELCC(
    kPdgTgtFe56, kPdgProton, kPdgNuMu, 3.);
  KPhaseSpaceCuts * cuts = KPhaseSpaceCuts::Instance();
  int failures = 0;
  failures += !Check(cuts->HasQ2MinCut(weak),
    "weak interaction sees Weak-Q2-min from XML");
  failures += !Check(NearlyEqual(cuts->Q2MinCut(weak, 0.), 0.25),
    "weak interaction uses Weak-Q2-min from XML");
  failures += !Check(cuts->HasSplineQ2MinCut(),
    "Weak-Q2-min provides spline metadata");
  failures += !Check(NearlyEqual(cuts->SplineQ2MinCut(), 0.25),
    "spline metadata uses Weak-Q2-min from XML");
  failures += !Check(cuts->HasSplineQ2MinCutForProbe(kPdgNuMu),
    "Weak-Q2-min provides neutrino spline metadata");
  failures += !Check(NearlyEqual(cuts->SplineQ2MinCutForProbe(kPdgNuMu), 0.25),
    "neutrino spline metadata uses Weak-Q2-min from XML");
  delete weak;
  return failures == 0 ? 0 : 1;
}
//__________________________________________________________________________
int RunNoSplineQ2CutCheck(void)
{
  KPhaseSpaceCuts * cuts = KPhaseSpaceCuts::Instance();
  int failures = 0;
  failures += !Check(!cuts->HasSplineQ2MinCut(),
    "empty CommonPhaseSpaceCuts.xml has no spline Q2 metadata");
  failures += !Check(!cuts->HasSplineQ2MinCutForProbe(kPdgNuMu),
    "empty CommonPhaseSpaceCuts.xml has no neutrino spline Q2 metadata");
  failures += !Check(!cuts->HasSplineQ2MinCutForProbe(kPdgElectron),
    "empty CommonPhaseSpaceCuts.xml has no electron spline Q2 metadata");
  return failures == 0 ? 0 : 1;
}
//__________________________________________________________________________
bool Check(bool ok, const char * msg)
{
  if(ok) LOG("test", pNOTICE) << "PASS: " << msg;
  else   LOG("test", pERROR)  << "FAIL: " << msg;
  return ok;
}
//__________________________________________________________________________
bool NearlyEqual(double a, double b)
{
  return std::fabs(a - b) < 1e-12;
}
//__________________________________________________________________________
int RunChild(const char * self, const char * mode, const std::string & xml)
{
  size_t slash = xml.rfind('/');
  std::string dir = slash == std::string::npos ? "." : xml.substr(0, slash);

  pid_t pid = fork();
  if(pid == 0) {
    setenv("GXMLPATH", dir.c_str(), 1);
    execlp(self, self, mode, static_cast<char *>(0));
    _exit(127);
  }

  int status = 0;
  if(pid < 0 || waitpid(pid, &status, 0) < 0) return -1;
  return status;
}
//__________________________________________________________________________
std::string WritePhaseSpaceCutsXML(const std::string & xml_body)
{
  char tmp[] = "/tmp/genie-kphase-space-cuts-XXXXXX";
  char * dir = mkdtemp(tmp);
  if(!dir) return "";

  std::string xml = std::string(dir) + "/CommonPhaseSpaceCuts.xml";
  std::ofstream out(xml.c_str());
  out << xml_body;
  out.close();
  return xml;
}
//__________________________________________________________________________
std::string WriteSplineXML(const std::string & xml_body)
{
  char tmp[] = "/tmp/genie-xsec-spline-XXXXXX";
  char * dir = mkdtemp(tmp);
  if(!dir) return "";

  std::string xml = std::string(dir) + "/spline.xml";
  std::ofstream out(xml.c_str());
  out << xml_body;
  out.close();
  return xml;
}
//__________________________________________________________________________
void RemoveTempXML(const std::string & xml)
{
  if(xml.empty()) return;
  unlink(xml.c_str());
  size_t slash = xml.rfind('/');
  if(slash != std::string::npos) {
    std::string dir = xml.substr(0, slash);
    rmdir(dir.c_str());
  }
}
//__________________________________________________________________________
