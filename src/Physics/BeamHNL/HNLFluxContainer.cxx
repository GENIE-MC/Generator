//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 John Plows <komninos-john.plows@physics.ox.ac.yk>
 University of Oxford

 Robert Hatcher <rhatcher@fnal.gov>
 Fermi National Accelerator Laboratory

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool
*/
//____________________________________________________________________________

#include <cstdlib>
#include <fstream>
#include <vector>
#include <sstream>
#include <cassert>
#include <climits>

#include "libxml/xmlmemory.h"
#include "libxml/parser.h"

#include "Framework/Utils/XmlParserUtils.h"
#include "Framework/Utils/StringUtils.h"

#include <TFile.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TSystem.h>
#include <TStopwatch.h>

#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/GBuild.h"

#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Utils/UnitUtils.h"

using std::endl;

#include <vector>
#include <algorithm>
#include <iomanip>
#include "TRegexp.h"
#include "TString.h"

#include "Physics/BeamHNL/HNLFluxContainer.h"

using namespace genie;
using namespace genie::hnl;

//___________________________________________________________________________
FluxContainer::FluxContainer()
  : TObject()
{
  this->ResetCopy();
}

//___________________________________________________________________________
void FluxContainer::ResetCopy() const
{
  evtno = 0;

  pdg = 0;
  parPdg = 0;
  lepPdg = 0;
  nuPdg = 0;

  prodChan = 0;
  nuProdChan = 0;

  startPoint.SetXYZ(0.0, 0.0, 0.0);
  targetPoint.SetXYZ(0.0, 0.0, 0.0);
  startPointUser.SetXYZ(0.0, 0.0, 0.0);
  targetPointUser.SetXYZ(0.0, 0.0, 0.0);
  delay = 0.0;

  polz.SetXYZ(0.0, 0.0, 0.0);

  p4.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
  parp4.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
  p4User.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
  parp4User.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);

  Ecm = 0.0;
  nuEcm = 0.0;

  XYWgt = 0.0;
  boostCorr = 0.0;
  accCorr = 0.0;
  zetaMinus = 0.0;
  zetaPlus = 0.0;
  acceptance = 0.0;

  nimpwt = 0.0;
}

//___________________________________________________________________________
void FluxContainer::Print(const Option_t* /* opt */ ) const
{
  std::cout << *this << std::endl;
}
//___________________________________________________________________________

namespace genie {
namespace hnl  {
  ostream & operator << (
    ostream & stream, const FluxContainer & info)
    {
      // convert some stuff to string
      HNLProd_t hProdChan = static_cast<HNLProd_t>(info.prodChan);

      std::string sNuProdChan;
      int typeMod = (info.pdg > 0) ? 1 : -1;
      int switchChan = typeMod * info.nuProdChan;
      switch(switchChan){
      case 1:  sNuProdChan = std::string("K0L -> nue + pi- + e+"); break;
      case -1: sNuProdChan = std::string("K0L -> nuebar + pi+ + e-"); break;
      case 2:  sNuProdChan = std::string("K0L -> numu + pi- + mu"); break;
      case -2: sNuProdChan = std::string("K0L -> numubar + pi+ + mu-"); break;
      case 3:  sNuProdChan = std::string("K+ -> numu + mu+"); break;
      case -3: sNuProdChan = std::string("K- -> numubar + mu-"); break;
      case 4:  sNuProdChan = std::string("K+ -> nue + e+"); break;
      case -4: sNuProdChan = std::string("K- -> nuebar + e-"); break;
      case 5:  sNuProdChan = std::string("K+ -> numu + pi0 + mu+"); break;
      case -5: sNuProdChan = std::string("K- -> numubar + pi0 + mu-"); break;
      case 6:  sNuProdChan = std::string("K+ -> nue + pi0 + e+"); break;
      case -6: sNuProdChan = std::string("K- -> nuebar + pi0 + e-"); break;
      case 7:  sNuProdChan = std::string("pi+ -> numu + mu+"); break;
      case -7: sNuProdChan = std::string("pi- -> numubar + mu-"); break;
      case 8:  sNuProdChan = std::string("pi+ -> nue + e+"); break;
      case -8: sNuProdChan = std::string("pi- -> nuebar + e-"); break;
      case 9:  sNuProdChan = std::string("mu- -> numu + nuebar + e-"); break;
      case -9: sNuProdChan = std::string("mu+ -> numubar + nue + e+"); break;
      default: sNuProdChan = std::string("Unknown!"); break;
      }
      
      stream << "\nEvent number: " << info.evtno
	     << "\nHNL    PDG code: " << info.pdg
	     << "\nParent PDG code: " << info.parPdg
	     << "\nCo-produced lepton PDG code: " << info.lepPdg
	     << "\nParent weight: " << info.nimpwt
	     << "\nHNL polarisation vector [HNL rest frame, NEAR coords]: " << utils::print::Vec3AsString(&info.polz)
	     << "\nPDG code of equivalent SM neutrino: " << info.nuPdg
	     << "\nProduction channel: " << utils::hnl::ProdAsString(hProdChan)
	     << " ; code " << info.prodChan
	     << "\nEquivalent neutrino production channel: " << sNuProdChan
	     << " ; code " << info.nuProdChan
	     << "\nHNL parent rest-frame energy [GeV]: " << info.Ecm
	     << "\nEquivalent neutrino parent rest-frame energy [GeV]: " << info.nuEcm
	     << "\nStart point [NEAR, m]: " << utils::print::Vec3AsString(&info.startPoint)
	     << "\nStart point [USER, m]: " << utils::print::Vec3AsString(&info.startPointUser)
	     << "\nFlux passes through point [NEAR, m]: " << utils::print::Vec3AsString(&info.targetPoint)
	     << "\nFlux passes through point [USER, m]: " << utils::print::Vec3AsString(&info.targetPointUser)
	     << "\nHNL momentum [NEAR, GeV]: " << utils::print::P4AsString(&info.p4)
	     << "\nHNL momentum [USER, GeV]: " << utils::print::P4AsString(&info.p4User)
	     << "\nHNL delay wrt SM neutrino [ns]: " << info.delay
	     << "\nParent momentum [NEAR, GeV]: " << utils::print::P4AsString(&info.parp4)
	     << "\nParent momentum [USER, GeV]: " << utils::print::P4AsString(&info.parp4User)
	     << "\nDeviation angles zeta- = " << info.zetaMinus << ", zeta+ = " << info.zetaPlus
	     << "\nGeometric acceptance: " << info.XYWgt
	     << "\nBoost correction: " << info.boostCorr
	     << "\nAcceptance correction: " << info.accCorr
	     << "\nFull acceptance: " << info.acceptance;

     return stream;
  }
}//hnl
}//genie

//___________________________________________________________________________
