#include "Tools/EvtLib/EvtLibPXSec.h"
#include "Tools/EvtLib/Utils.h"

#include "Framework/Conventions/Units.h"

#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"

#include "TFile.h"
#include "TGraph.h"

using namespace genie;
using namespace genie::evtlib;

//____________________________________________________________________________
EvtLibPXSec::EvtLibPXSec() :
XSecAlgorithmI("genie::evtlib::EvtLibPXSec")
{

}

//____________________________________________________________________________
EvtLibPXSec::EvtLibPXSec(string config) :
XSecAlgorithmI("genie::evtlib::EvtLibPXSec", config)
{

}

//____________________________________________________________________________
EvtLibPXSec::~EvtLibPXSec()
{
  ClearXSecs();
}

//____________________________________________________________________________
double EvtLibPXSec::XSec(const Interaction* /*in*/,
                         KinePhaseSpace_t /*kps*/) const
{
  LOG("ELI", pFATAL) << "This point should not be reached" ; 
  abort();
}

//____________________________________________________________________________
double EvtLibPXSec::Integral(const Interaction* in) const
{
  TGraph* g = GetXSec(in);
  if(!g) return 0; // Reason already printed

  const InitialState& init_state = in->InitState();
  const double E  = init_state.ProbeE(kRfLab);

  // Units of the cross-section graph are expected to be 10^-38 cm^2 / nucleus
  return g->Eval(E) * 1e-38 * genie::units::cm2;
}

//____________________________________________________________________________
bool EvtLibPXSec::ValidProcess(const Interaction* /*in*/) const
{
  return true;
}

//____________________________________________________________________________
void EvtLibPXSec::Configure(const Registry& config)
{
  Algorithm::Configure(config);
  LoadXSecs();
}

//____________________________________________________________________________
void EvtLibPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  LoadXSecs();
}

//____________________________________________________________________________
void EvtLibPXSec::ClearXSecs()
{
  for(auto it: fXSecs) delete it.second;
  fXSecs.clear();
}

//____________________________________________________________________________
void EvtLibPXSec::LoadXSecs()
{
  ClearXSecs();

  std::string libPath;
  GetParam("EventLibraryPath", libPath);
  Expand(libPath);

  PDGLibrary* pdglib = PDGLibrary::Instance();

  TFile fin(libPath.c_str());
  if(fin.IsZombie()) exit(1);

  TIter next(fin.GetListOfKeys());
  while(TObject* dir = next()){
    const std::string& tgtName = dir->GetName();
    const TParticlePDG* tgtPart = pdglib->DBase()->GetParticle(tgtName.c_str());
    if(!tgtPart){
      LOG("ELI", pWARN) << "Unknown nucleus " << tgtName
                        << " found in " << libPath
                        << " -- skipping";
      continue;
    }

    for(int pdg: {kPdgNuE,   kPdgAntiNuE,
                  kPdgNuMu,  kPdgAntiNuMu,
                  kPdgNuTau, kPdgAntiNuTau}){

      for(bool iscc: {true, false}){
        // NCs should be the same for all flavours. Use nu_mu as a convention
        // internal to this code to index into the xsecs map.
        if(!iscc && abs(pdg) != kPdgNuMu) continue;

        std::string nuName = pdglib->Find(pdg)->GetName();
        if(!iscc) nuName = pdg::IsAntiNeutrino(pdg) ? "nu_bar" : "nu";

        const std::string graphName =
          TString::Format("%s/%s/%s/xsec",
                          tgtName.c_str(),
                          iscc ? "cc" : "nc",
                          nuName.c_str()).Data();

        const Key key(tgtPart->PdgCode(), pdg, iscc);

        TGraph* g = (TGraph*)fin.Get(graphName.c_str());
        if(!g){
          LOG("ELI", pINFO) << graphName << " not found in "
                            << libPath << " -- skipping";
          continue;
        }

        fXSecs[key] = new TGraph(*g);
      } // end for iscc
    } // end for pdg
  } // end for dir
}

//____________________________________________________________________________
TGraph* EvtLibPXSec::GetXSec(const Interaction* in) const
{
  const InitialState& init_state = in->InitState();

  const ProcessInfo& proc = in->ProcInfo();
  if(!proc.IsWeakCC() && !proc.IsWeakNC()){
    LOG("ELI", pINFO) << "Skipping unknown process " << proc;
    return 0;
  }

  int pdg = init_state.ProbePdg();
  // Use nu_mu for NC as a convention internal to this code to index into the
  // xsec graphs map.
  if(proc.IsWeakNC()){
    if(pdg > 0) pdg = kPdgNuMu; else pdg = kPdgAntiNuMu;
  }

  const Key key(init_state.TgtPdg(), pdg, proc.IsWeakCC());

  auto it = fXSecs.find(key);
  if(it == fXSecs.end()){
    LOG("ELI", pINFO) << "Skipping process without xsec " << key;
    return 0;
  }

  return it->second;
}
