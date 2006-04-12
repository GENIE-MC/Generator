//____________________________________________________________________________
/*!

\program testSAIDHadronXSec

\brief   test program used for testing the SAID hadron cross section data
         and corrections provided by S.Dytman.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created May 12, 2004
*/
//____________________________________________________________________________

#include "Numerical/Spline.h"
#include "Nuclear/SAIDHadronXSec.h"

using namespace genie;

int main(int argc, char ** argv)
{
  // load SAID hadron cross section data
  SAIDHadronXSec * sxs = SAIDHadronXSec::Instance();

  // save splines to a ROOT file
  sxs -> PiplusPElasXSecSpl()  -> SaveAsROOT("said.root", "pip_p_elas", true);
  sxs -> PiplusPInelXSecSpl()  -> SaveAsROOT("said.root", "pip_p_inel", false);
  sxs -> PiminusPElasXSecSpl() -> SaveAsROOT("said.root", "pim_p_elas", false);
  sxs -> PiminusPInelXSecSpl() -> SaveAsROOT("said.root", "pim_p_inel", false);
  sxs -> PiminusPCExXSecSpl()  -> SaveAsROOT("said.root", "pim_p_cex",  false);
  sxs -> PiplusAbsXSecSpl()    -> SaveAsROOT("said.root", "pip_abs",    false);
  sxs -> NPElasXSecSpl()       -> SaveAsROOT("said.root", "np_elas",    false);
  sxs -> NPInelXSecSpl()       -> SaveAsROOT("said.root", "np_inel",    false);
  sxs -> PPElasXSecSpl()       -> SaveAsROOT("said.root", "pp_elas",    false);
  sxs -> PPInelXSecSpl()       -> SaveAsROOT("said.root", "pp_inel",    false);
}
