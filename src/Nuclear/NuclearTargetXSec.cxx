//____________________________________________________________________________
/*!

\class    genie::NuclearTargetXSec

\brief    Computes the cross section for a nuclear target from the free
          nucleon cross sections.

          Is a concrete implementation of the XSecAlgorithmI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  June 10, 2004

*/
//____________________________________________________________________________

#include "Nuclear/NuclearTargetXSec.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using namespace genie;

//____________________________________________________________________________
NuclearTargetXSec::NuclearTargetXSec() :
XSecAlgorithmI()
{
  fName = "genie::NuclearTargetXSec";
}
//____________________________________________________________________________
NuclearTargetXSec::NuclearTargetXSec(const char * param_set) :
XSecAlgorithmI(param_set)
{
  fName = "genie::NuclearTargetXSec";

  this->FindConfig();
}
//____________________________________________________________________________
NuclearTargetXSec::~NuclearTargetXSec()
{

}
//____________________________________________________________________________
double NuclearTargetXSec::XSec(const Interaction * interaction) const
{
  double xsec = 0;

  //--- get particle (nucleon, electron) cross section calculator
  const Algorithm * xsec_alg_base = this->SubAlg(
              "free-particle-xsec-alg-name", "free-particle-xsec-param-set");
  const XSecAlgorithmI * free_particle_xsec_alg =
                         dynamic_cast<const XSecAlgorithmI *>(xsec_alg_base);

  //--- get number of proton/neutrons/nucleons/e
  const Target & tgt = interaction->GetInitialState().GetTarget();
  int Z = tgt.Z();
  int N = tgt.N();
  int A = tgt.A();

  //--- check whether a struck nucleon was set
  int  pdgc = tgt.StruckNucleonPDGCode();
  bool struck_nuc_set = pdg::IsProton(pdgc) || pdg::IsNeutron(pdgc);

  //--- get the process type
  const ProcessInfo & proc_info = interaction->GetProcessInfo();

  //--- handle IMD:
  if(proc_info.IsInverseMuDecay()) {
    double xsec_e = free_particle_xsec_alg->XSec(interaction);
    xsec = Z*xsec_e;
    return xsec;
  }

  //--- handle QEL,DIS,RES,COH:

  //--- get a cloned interaction
  Interaction ci(*interaction);

  if(struck_nuc_set) {
    // get cross section for a free struck nucleon (p or n) and multiply
    // with the number of p or n in the target
    LOG("NuclXSec", pDEBUG)
         << "Computing xsec as a weighted average over nucl.target's "
                            << (pdg::IsProton(pdgc) ? "protons" : "neutrons");

    int nnuc = pdg::IsProton(pdgc) ? Z : N;
    double xsec_nuc = free_particle_xsec_alg->XSec(&ci);
    xsec = nnuc * xsec_nuc;

  } else {
    // get cross section for free nucleons (p,n) compute a weighted average
    LOG("NuclXSec", pDEBUG)
       << "Computing xsec as a weighted average over nucl.target's nucleons";

    Target * ctgt = ci.GetInitialStatePtr()->GetTargetPtr();

    ctgt->SetStruckNucleonPDGCode(kPdgProton);
    double xsec_p = free_particle_xsec_alg->XSec(&ci);

    ctgt->SetStruckNucleonPDGCode(kPdgNeutron);
    double xsec_n = free_particle_xsec_alg->XSec(&ci);

    xsec = (Z*xsec_p + N*xsec_n) / A;
  }
  LOG("NuclXSec", pDEBUG)
           << "xsec for a nucl.target with (Z,A) = ("
                                         << Z << ", " << A << ") = " << xsec;
  return xsec;
}
//____________________________________________________________________________
