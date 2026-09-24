#include <cstdlib>
#include <string>
#include <utility>
#include <vector>

#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TObject.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVector3.h>
#include <TObjString.h>
#include <TString.h>

#include "Framework/EventGen/EventRecord.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/Ntuple/NtpMCFormat.h"
#include "Framework/Ntuple/NtpMCTreeHeader.h"
#include "Framework/Ntuple/NtpMCTreeMultiHeader.h"
#include "Framework/Ntuple/NtpWriter.h"
#include "Framework/Utils/CmdLnArgParser.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"

using namespace genie;



void GetCommandLineArgs(int argc, char **argv);
void MergeFiles(void);
void PrintSyntax(void);

NtpMCFormat_t kDefOptNtpFormat = kNFGHEP;
long int gOptRanSeed = -1;
Long_t gOptRunNu = 0;
std::string gOutFileName;
std::vector<std::string> inputfiles;
NtpMCTreeMultiHeader head;
bool ignoreMissmatch = false;
TString common_cvstag = TString("");
TString common_tune = TString("");



int main(int argc, char **argv) {
  GetCommandLineArgs(argc, argv);

  MergeFiles();

  return 0;
}

void GetCommandLineArgs(int argc, char **argv) {
  LOG("gEvMerge", pINFO) << "Parsing command line arguments";
  CmdLnArgParser parser(argc, argv);
  if (parser.OptionExists('h')) {
    PrintSyntax();
    std::exit(0);
  }

  if (parser.OptionExists('o')) {
    gOutFileName = parser.ArgAsString('o');
  } else {

    LOG("gEvMerge", pFATAL) << "No output file given";
    gAbortingInErr = true;
    std::exit(1);
  }

  if (parser.OptionExists('i')) {
    std::string files = parser.ArgAsString('i');
    inputfiles = utils::str::Split(files, ",");
    if (inputfiles.empty()) {
      LOG("gEvMerge", pFATAL) << "No input files given";
      gAbortingInErr = true;
      std::exit(1);
    }

  } else {

    LOG("gEvMerge", pFATAL) << "No input files given";

    gAbortingInErr = true;
    std::exit(1);
  }
  if (parser.OptionExists("ignore-missmatch")){
    ignoreMissmatch = true;
  }
  std::vector<std::string> valid_files;
  for (const std::string &file : inputfiles) {
    if (!file.empty()) {
      valid_files.push_back(file);
    }
  }
  inputfiles.swap(valid_files);
  if (inputfiles.empty()) {
    LOG("gEvMerge", pFATAL) << "No valid input files given";
    gAbortingInErr = true;
    std::exit(1);
  }
  LOG("gEvMerge", pNOTICE) << "Output file: " << gOutFileName;
  LOG("gEvMerge", pNOTICE) << "Number of input files: "
                           << inputfiles.size();
  for (const std::string &file : inputfiles) {
    LOG("gEvMerge", pINFO) << "  " << file;
  }
}

void MergeFiles(void) {

  NtpWriter ntpw(kDefOptNtpFormat, gOptRunNu, gOptRanSeed);
  NtpMCTreeHeader* thdrWriter = ntpw.EventTreeHeader();
  bool matchTune = true;
  bool matchTag = true;
  for (const std::string &filename : inputfiles) {
    TFile fin(filename.c_str(), "READ");
    NtpMCTreeHeader *thdr = dynamic_cast<NtpMCTreeHeader *>(fin.Get("header"));
    if (!thdr) {
      LOG("gEvMerge", pWARN)
          << "Input file does not contain a valid GENIE header: " << filename
          << ", skipping file";
      continue;
    } else {
       LOG("gEvMerge", pNOTICE) << "Retrieving headder info for " << filename;
      if (common_cvstag == TString("") || common_tune == TString("")){
        common_cvstag = thdr->cvstag.GetString ();
        common_tune = thdr->tune.GetString ();
      } else {
        matchTune &= common_tune == (thdr->tune.GetString ());
        matchTag &=  common_cvstag == (thdr->cvstag).GetString ();
      }
    }
  }
  LOG("gEvMerge", pNOTICE) << "Finished analyzing headders";

  if ((!matchTune || !matchTag) && !ignoreMissmatch){
    LOG("gEvMerge", pERROR) << "Missmatch in tunings or cvstag in files. If you really want to merge them use --ignore-missmatch";
    gAbortingInErr = true;
    std::exit(1);
  }


  ntpw.CustomizeFilename(gOutFileName);
  ntpw.Initialize();

  Long64_t ievt = 0;
  for (const std::string &filename : inputfiles) {
    std::set<int> nu;
    LOG("gEvMerge", pNOTICE) << "Opening input file: " << filename;

    TFile fin(filename.c_str(), "READ");
    if (fin.IsZombie() || !fin.IsOpen()) {
      LOG("gEvMerge", pERROR) << "Could not open input file: " << filename;
      continue;
    }
    TTree *er_tree = dynamic_cast<TTree *>(fin.Get("gtree"));
    if (!er_tree) {
      LOG("gEvMerge", pERROR) << "Input file does not contain a valid "
                              << "GENIE event tree 'gtree': " << filename;
      continue;
    }
    NtpMCTreeHeader *thdr = dynamic_cast<NtpMCTreeHeader *>(fin.Get("header"));
    if (!thdr) {
      continue;

    } else {
      LOG("gEvMerge", pINFO) << "Input header for " << filename << ":\n"
                             << *thdr;
    }
    LOG("gEvMerge", pNOTICE) << "Retrieving Multihead"; 
    NtpMCTreeMultiHeader *multiCurr =  dynamic_cast<NtpMCTreeMultiHeader *>(fin.Get("MultiHead"));
    if (multiCurr){
      for (size_t j = 0; j < multiCurr->size(); j++){
        std::pair<Long64_t, Long64_t> indi = multiCurr->getIndices(j);
        head.insertHead(multiCurr->getHead(j), 
        multiCurr->getFileName(j), 
        multiCurr->getNeutrinos(j), 
        std::make_pair(indi.first + ievt,  indi.second + ievt)
        );
      }
    }
    NtpMCEventRecord *mcrec = nullptr;
    if (er_tree->SetBranchAddress("gmcrec", &mcrec) < 0) {
      LOG("gEvMerge", pERROR)
          << "Could not set branch address for 'gmcrec' in " << filename;
      continue;
    }
    const Long64_t n = er_tree->GetEntries();
    LOG("gEvMerge", pNOTICE) << "Input file contains " << n << " events";
    if (n == 0) {
      LOG("gEvMerge", pWARN) << "Input file contains no events: " << filename;

      continue;
    }
    const Long64_t first_output_event = ievt;
    for (Long64_t iev = 0; iev < n; ++iev) {
      const Long64_t bytes = er_tree->GetEntry(iev);
      if (bytes <= 0) {
        LOG("gEvMerge", pERROR)
            << "Failed to read event " << iev << " from " << filename;

        continue;
      }

      if (!mcrec) {

        LOG("gEvMerge", pERROR)
            << "Null NtpMCEventRecord for event " << iev << " in " << filename;

        continue;
      }

      if (!mcrec->event) {

        LOG("gEvMerge", pERROR)
            << "Null EventRecord for event " << iev << " in " << filename;

        continue;
      }
      nu.insert(mcrec->event->Probe()->Pdg());
      ntpw.AddEventRecord(static_cast<int>(ievt), mcrec->event);

      ++ievt;
    }

    if (ievt > first_output_event) {
      head.insertHead(*thdr, 
        filename, 
        nu, 
        std::make_pair(first_output_event, ievt - 1)
        );
      LOG("gEvMerge", pNOTICE)
          << "Merged " << (ievt - first_output_event) << " events from "
          << filename << " -> output events [" << first_output_event << ", "
          << (ievt - 1) << "]";
    }
    fin.Close();
  }
  LOG("gEvMerge", pNOTICE) << "Writing " << ievt << " total events to "
                           << gOutFileName;
  std::string tune = matchTune ? common_tune.Data()    : "N.N.";
  std::string csvt = matchTag  ? common_cvstag.Data()  : "N.N.";
  LOG("gEvMerge", pNOTICE) << "Changing Head, old one " << thdrWriter;
  thdrWriter->tune.SetString(tune.c_str());
  thdrWriter->cvstag.SetString(csvt.c_str());
  ntpw.Save();
  TFile * fOutFile = TFile::Open(gOutFileName.c_str(),"UPDATE");
  std::cout << head << std::endl;
  head.Write("MultiHead");
  fOutFile->Write();
  fOutFile->Close();
  
  LOG("gEvMerge", pNOTICE) << "Merge completed successfully";
}

void PrintSyntax(void) {
  LOG("gEvMerge", pNOTICE)
      << "\n\n"
      << "Syntax:"
      << "\n"
      << "\n      gevmerge [-h]"
      << "\n"
      << "               -i input_file1.root,input_file2.root,..."
      << "\n"
      << "               -o outfile_name.root\n"
      << "               --ignore-missmatch"
      << "\n"
      << RunOpt::RunOptSyntaxString(true) << "\n";
}

