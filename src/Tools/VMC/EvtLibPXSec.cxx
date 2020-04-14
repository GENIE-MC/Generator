#include "Tools/VMC/EvtLibPXSec.h"
#include "Tools/VMC/Utils.h"

#include "Framework/Conventions/Units.h"

#include "Framework/ParticleData/PDGLibrary.h"

#include "TFile.h"
#include "TGraph.h"

using namespace genie;
using namespace genie::vmc;

//____________________________________________________________________________
EvtLibPXSec::EvtLibPXSec() :
XSecAlgorithmI("genie::vmc::EvtLibPXSec")
{

}
//____________________________________________________________________________
EvtLibPXSec::EvtLibPXSec(string config) :
XSecAlgorithmI("genie::vmc::EvtLibPXSec", config)
{

}
//____________________________________________________________________________
EvtLibPXSec::~EvtLibPXSec()
{

}
//____________________________________________________________________________
double EvtLibPXSec::XSec(const Interaction* /*in*/,
                         KinePhaseSpace_t /*kps*/) const
{
  std::cout << "EvtLibPXSec::XSec(): should not be reached" << std::endl;
  abort();
}

//____________________________________________________________________________
double EvtLibPXSec::Integral(const Interaction* in) const
{
  TGraph* g = GetXSec(in);
  if(!g) return 0; // Reason already printed

  const InitialState& init_state = in->InitState();
  const double E  = init_state.ProbeE(kRfLab);
  const Target& target = init_state.Tgt();

  return g->Eval(E) * target.N() * 1e-38 * genie::units::cm2;
}
//____________________________________________________________________________
bool EvtLibPXSec::ValidProcess(const Interaction* /*in*/) const
{
  return true;
}
//____________________________________________________________________________
void EvtLibPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
}
//____________________________________________________________________________
void EvtLibPXSec::Configure(string config)
{
  Algorithm::Configure(config);
}

//____________________________________________________________________________
void EvtLibPXSec::LoadXSecs() const
{
  std::string libPath;
  if(!GetParamDef( "LibraryPath", libPath, std::string() )){
    LOG("ELI", pFATAL) << "Must specify 'LibraryPath'";
    exit(1);
  }

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

    for(int sign: {+1, -1}){
      for(int pdg: {12, 14, 16}){
        for(bool iscc: {true, false}){
          // NCs should be the same for all flavours. Use numu by convention.
          if(!iscc && pdg != 14) continue;

          std::string nuName = pdglib->Find(sign*pdg)->GetName();
          if(!iscc) nuName = (sign > 0) ? "nu" : "nu_bar";

          const std::string graphName =
            TString::Format("%s/%s/%s/xsec",
                            tgtName.c_str(),
                            iscc ? "cc" : "nc",
                            nuName.c_str()).Data();

          const Key key(tgtPart->PdgCode(), sign*pdg, iscc);

          TGraph* g = (TGraph*)fin.Get(graphName.c_str());
          if(!g){
            LOG("ELI", pINFO) << graphName << " not found in "
                              << libPath << " -- skipping";
            continue;
          }

          fXSecs[key] = new TGraph(*g);
        } // end for nucleus
      } // end for pdg
    } // end for sign
  } // end for iscc
}

//____________________________________________________________________________
TGraph* EvtLibPXSec::GetXSec(const Interaction* in) const
{
  if(fXSecs.empty()) LoadXSecs();

  const InitialState& init_state = in->InitState();

  const ProcessInfo& proc = in->ProcInfo();
  if(!proc.IsWeakCC() && !proc.IsWeakNC()){
    LOG("ELI", pINFO) << "Skipping unknown process " << proc;
    return 0;
  }

  int pdg = init_state.ProbePdg();
  // Use nu_mu for NC by convention
  if(proc.IsWeakNC()){
    if(pdg > 0) pdg = +14; else pdg = -14;
  }

  const Key key(init_state.TgtPdg(), pdg, proc.IsWeakCC());

  auto it = fXSecs.find(key);
  if(it == fXSecs.end()){
    LOG("ELI", pINFO) << "Skipping process without xsec " << key;
    return 0;
  }

  return it->second;
}
