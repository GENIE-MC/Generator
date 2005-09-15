//____________________________________________________________________________
/*!

\class   genie::geant::GNtpPrimaryGenerator

\brief   A GENIE/GEANT4 PrimaryGenerator. The NtpRdPrimaryGenerator feeds
         GENIE events into GEANT by simply reading them from GENIE's ER Tree.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created Sepember 12, 2005

*/
//____________________________________________________________________________

#ifndef _G_NTP_PRIMARY_GENERATOR_H_
#define _G_NTP_PRIMARY_GENERATOR_H_

#include <string>

#include <G4VPrimaryGenerator.hh>

class TFile;
class TTree;
class G4Event;

using std::string;

namespace genie {

class NtpMCTreeHeader;
class NtpMCEventRecord;

namespace geant {

class GNtpPrimaryGenerator : public G4VPrimaryGenerator
{
public:
  GNtpPrimaryGenerator();
  virtual ~GNtpPrimaryGenerator();

  // overriden methods from GEANT's G4VPrimaryGenerator interface
  void GeneratePrimaryVertex(G4Event * event);

  // methods specific to this G4VPrimaryGenerator implementation
  void ReadFromFile(string filename);

private:

  void Initialize (void);
  void CleanUp    (void);

  EventRecord * ReadNextEvent(void);

  TFile *            fFile;
  TTree *            fTree;
  NtpMCTreeHeader *  fTreeHdr;
  NtpMCEventRecord * fNtpRec;
  Long64_t           fNEntries;
  Long64_t           fCurrentEvent;
};

} // geant namespace
} // genie namespace

#endif // _G_NTP_PRIMARY_GENERATOR_H_


