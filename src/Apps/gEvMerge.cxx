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

#include "Framework/EventGen/EventRecord.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/Ntuple/NtpMCFormat.h"
#include "Framework/Ntuple/NtpMCTreeHeader.h"
#include "Framework/Ntuple/NtpWriter.h"
#include "Framework/Utils/CmdLnArgParser.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Utils/StringUtils.h"

using namespace genie;

struct MultiHead : public TObject { // TODO: move to framework

  std::vector<std::string> inputfiles;

  // [first event, last event] in the merged output corresponding
  // to each input file.
  std::vector<std::pair<Long64_t, Long64_t>> indices;

  // generated neutrinos each run.
  std::vector<std::vector<int>> neutrinos;

  // Store the headers of all input files.
  std::vector<NtpMCTreeHeader> heads;

};

void GetCommandLineArgs(int argc, char **argv);
void MergeFiles(void);
void PrintSyntax(void);

NtpMCFormat_t kDefOptNtpFormat = kNFGHEP;
long int gOptRanSeed = -1;
Long_t gOptRunNu = 0;
std::string gOutFileName;
MultiHead head;

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

  if (parser.OptionExists('f')) {
    std::string files = parser.ArgAsString('f');
    head.inputfiles = utils::str::Split(files, ",");
    if (head.inputfiles.empty()) {
      LOG("gEvMerge", pFATAL) << "No input files given";
      gAbortingInErr = true;
      std::exit(1);
    }

  } else {

    LOG("gEvMerge", pFATAL) << "No input files given";

    gAbortingInErr = true;
    std::exit(1);
  }

  std::vector<std::string> valid_files;
  for (const std::string &file : head.inputfiles) {
    if (!file.empty()) {
      valid_files.push_back(file);
    }
  }
  head.inputfiles.swap(valid_files);
  if (head.inputfiles.empty()) {
    LOG("gEvMerge", pFATAL) << "No valid input files given";
    gAbortingInErr = true;
    std::exit(1);
  }
  LOG("gEvMerge", pNOTICE) << "Output file: " << gOutFileName;
  LOG("gEvMerge", pNOTICE) << "Number of input files: "
                           << head.inputfiles.size();
  for (const std::string &file : head.inputfiles) {
    LOG("gEvMerge", pINFO) << "  " << file;
  }
}

void MergeFiles(void) {

  NtpWriter ntpw(kDefOptNtpFormat, gOptRunNu, gOptRanSeed);
  ntpw.CustomizeFilename(gOutFileName);
  ntpw.Initialize();
  Long64_t ievt = 0;
  for (const std::string &filename : head.inputfiles) {
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
      LOG("gEvMerge", pWARN)
          << "Input file does not contain a valid GENIE header: " << filename
          << ", skipping file";
      continue;

    } else {
      LOG("gEvMerge", pINFO) << "Input header for " << filename << ":\n"
                             << *thdr;
      head.heads.push_back(*thdr);
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
      ntpw.AddEventRecord(static_cast<int>(ievt), mcrec->event);

      ++ievt;
    }

    if (ievt > first_output_event) {
      head.indices.push_back(std::make_pair(first_output_event, ievt - 1));
      LOG("gEvMerge", pNOTICE)
          << "Merged " << (ievt - first_output_event) << " events from "
          << filename << " -> output events [" << first_output_event << ", "
          << (ievt - 1) << "]";
    }
    fin.Close();
  }
  LOG("gEvMerge", pNOTICE) << "Writing " << ievt << " total events to "
                           << gOutFileName;
  ntpw.Save();
  LOG("gEvMerge", pNOTICE) << "Merge completed successfully";
}

void PrintSyntax(void) {
  LOG("gEvMerge", pNOTICE)
      << "\n\n"
      << "Syntax:"
      << "\n"
      << "\n      gevmerge [-h]"
      << "\n"
      << "               -f input_file1.root,input_file2.root,..."
      << "\n"
      << "               -o outfile_name.root"
      << "\n"
      << RunOpt::RunOptSyntaxString(true) << "\n";
}

