//_________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 For the class documentation see the corresponding header file.
*/
//_________________________________________________________________________

#include <cstdlib>
#include <string>
#include <fstream>
#include <cassert>
#include <algorithm>

#include "Framework/Conventions/Constants.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Numerical/BLI2D.h"
#include "Physics/Multinucleon/XSection/MECHadronTensor.h"

#include <TSystem.h>
#include <TFile.h>

using std::ostringstream;
using std::istream;
using std::ios;

using namespace genie;
using namespace genie::constants;

//_________________________________________________________________________
MECHadronTensor * MECHadronTensor::fgInstance = 0;
//_________________________________________________________________________
MECHadronTensor::MECHadronTensor()
{
  // This list should be abstracted to a configuration file
  // Along with the known target method.
  
  // Three targets have explicit tensor tables.
  //
  fKnownTensors.push_back(kPdgTgtC12);  // C12 
  fKnownTensors.push_back(kPdgTgtO16);  // O16
  fKnownTensors.push_back(1000140280);  // Si28
  fKnownTensors.push_back(kPdgTgtCa40); // Ca40, also for Ar40
  fKnownTensors.push_back(1000280560);  // Ni56,  actually it is pseudo-Fe56
  fKnownTensors.push_back(1000561120);  // Ba112, actually it is pseudo-Cd112
  fKnownTensors.push_back(1001042080);  // Rf208, actually it is pseudo-Pb208
  // why pseudo ?  The had tensor calculation requires the nuclear density function
  // simple two-parameter hartree fock Fermi function for heavy nuclei .
  // and in principle modified harmonic oscillator function for light nuclei
  // we are never really going to want isoscalar Ni56, so I used density for Fe56.
  // likewise never want Rf208, I used the density for Pb208
  // likewise never want Ba112, I used the density for Cd112

  // this one loads all known targets when instantiated, regardless of what the user wants.
  // is there a point to being efficient, check if the requested one is loaded ?
  // TODO  there is a fail if the spline file is missing one of the known targets
  // and the user requests a mix of two targets.
  vector<int>::const_iterator it=fKnownTensors.begin();
  for( ; it!=fKnownTensors.end(); ++it) {
    this->LoadTensorTables(*it);
    // this method will fail silently if the target list here does not match later.   
  }  
  fgInstance = 0;
}
//_________________________________________________________________________
MECHadronTensor::~MECHadronTensor()
{
	vector<int>::const_iterator it=fKnownTensors.begin();
	for( ; it!=fKnownTensors.end(); ++it) {
		MECHadronTensor::MECHadronTensorTable table = fTargetTensorTables[*it];
		for(int tensorType = 0; tensorType <= MECHadronTensor::kMHTValenciaDeltapn; ++tensorType) {
			vector<genie::BLI2DNonUnifGrid*> Grids = table.Table[(MECHadronTensor::MECHadronTensorType_t)tensorType];
			for (vector<genie::BLI2DNonUnifGrid*>::iterator gridit = Grids.begin() ; gridit != Grids.end(); ++gridit){
				genie::BLI2DNonUnifGrid *hadTensorGrid = *gridit;
				if(hadTensorGrid) {
				  delete hadTensorGrid;
				  hadTensorGrid=0;
				}
			}
		}
	}  
}
//_________________________________________________________________________
MECHadronTensor * MECHadronTensor::Instance()
{
  if(fgInstance == 0) {
    LOG("MECHadronTensor", pDEBUG) << "Late initialization";
    static MECHadronTensor::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fgInstance = new MECHadronTensor();
  }  
  return fgInstance;
}
//_________________________________________________________________________
bool MECHadronTensor::KnownTarget(int targetpdg)
{
  if(std::count(fKnownTensors.begin(), fKnownTensors.end(), targetpdg)!=0) {
    return true;
  }

  // abstract these limits.  They are used also in XSecFour.
  // hard coding things in two locations isn't a good idea.
  // also, future-person might give a had tensor and want to not allow scaling.

  int Arequest = pdg::IonPdgCodeToA(targetpdg);
  int Zrequest = pdg::IonPdgCodeToZ(targetpdg);
  if(Arequest >= 9 || (Arequest == 4 && Zrequest == 2)){
    return true;
  }
  return false;  
} 
//_________________________________________________________________________
bool MECHadronTensor::KnownTensor(int targetpdg)
{
  return std::count(fKnownTensors.begin(), fKnownTensors.end(), targetpdg)!=0;  
}
//_________________________________________________________________________
const vector<genie::BLI2DNonUnifGrid *> &
   MECHadronTensor::TensorTable(int targetpdg, MECHadronTensorType_t type)
{
  return fTargetTensorTables[targetpdg].Table[type];
}
//_________________________________________________________________________
void MECHadronTensor::LoadTensorTables(int targetpdg)
{
// Load the hadron tensor tables.
// For the Nieves model they are in ${GENIE}/data/evgen/mectensor/nieves/
  
  // define dimensions of data in the hadron tensor files
  int nwpoints = 5;
  const int nq0points = 120;  
  const int nqzpoints = 120;  
  const int nq0qzpoints = nq0points*nqzpoints;
  double arraystep = 0.01; // GeV
  // if later we use tables that are not 120x120
  // then extract these constants to the config file with the input table specs
  // or find a way for the input tables to be self-descriptive.

  string tensorLocation  = string("/data/evgen/mectensor/nieves");
  string tensorFileStart = "HadTensor120-Nieves-";
  string tensorFileEnd   = "-20150210.dat";

  if(!KnownTarget(targetpdg)){
    LOG("MECHadronTensor", pERROR) 
      << "No MEC tensor table for target with PDG code: " 
      << targetpdg;
    return;
  }

  // define directory of hadron tensor files
  // Ideally, the xml configuration can override the default location
  string data_dir = string(gSystem->Getenv("GENIE")) + tensorLocation;

  // define arrays to fill from data files
  double hadtensor_q0_array[nq0points];
  double hadtensor_qz_array[nqzpoints];
  double hadtensor_w_array[nwpoints][nq0qzpoints];

  // fill q0 array (in GeV)
  // 240 5MeV bins or 120 10 MeV bins
  for (int a = 0; a < nq0points; a++){
    hadtensor_q0_array[a]=double(a+1)*arraystep;
  }

  // fill qz array (in GeV)
  for (int a = 0; a < nqzpoints; a++){
    hadtensor_qz_array[a]=double(a+1)*arraystep;
  }
  
  // possible future feature, allow a model to not deliver Delta tensors.
  std::map<MECHadronTensor::MECHadronTensorType_t, std::string> tensorTypeNames;
  tensorTypeNames[MECHadronTensor::kMHTValenciaFullAll]  = "FullAll";
  tensorTypeNames[MECHadronTensor::kMHTValenciaFullpn]   = "Fullpn";
  tensorTypeNames[MECHadronTensor::kMHTValenciaDeltaAll] = "DeltaAll";
  tensorTypeNames[MECHadronTensor::kMHTValenciaDeltapn]  = "Deltapn";

  // iterate over all four hadron tensor types in the map above.
  for(int tensorType = 0; 
          tensorType <= MECHadronTensor::kMHTValenciaDeltapn; ++tensorType) {

    // build filenames from the bits of string
    ostringstream datafile;
    datafile << data_dir << "/" << tensorFileStart << targetpdg << "-" 
	     << tensorTypeNames[(MECHadronTensor::MECHadronTensorType_t)tensorType]
	     << tensorFileEnd;

    // make sure data files are available
    LOG("MECHadronTensor", pDEBUG) 
       << "Asserting that file " << datafile.str().c_str() << " exists...";      
    assert (! gSystem->AccessPathName(datafile.str().c_str()));
  
    // read data file
    ReadHadTensorqzq0File(
      datafile.str(), nwpoints, nqzpoints, nq0points, hadtensor_w_array
    );
  
    //loop over all 5 tensors 
    for (int i = 0; i < nwpoints; i++){
   
      // create a non uniform grid from tensor data
      genie::BLI2DNonUnifGrid *hadTensorGrid = 
           new genie::BLI2DNonUnifGrid(
               nqzpoints, 
               nq0points, 
               hadtensor_qz_array, 
               hadtensor_q0_array, 
               hadtensor_w_array[i]
           );

      // and store in a map using the target PDG as a key
      fTargetTensorTables[targetpdg].Table[
         (MECHadronTensor::MECHadronTensorType_t)tensorType].push_back(hadTensorGrid);
    }
  }
}
//_________________________________________________________________________
void MECHadronTensor::ReadHadTensorqzq0File( 
  string filename, int nwpoints, int nqzpoints, int nq0points, 
  double hadtensor_w_array[][14400])
{
  // open file
  std::ifstream tensor_stream(filename.c_str(), ios::in);

  // check file exists
  if(!tensor_stream.good()){
    LOG("MECHadronTensor", pERROR) << "Bad file name: " << filename;
    return;
  }

  double temp;  
  for (int ij = 0; ij < (nqzpoints*nq0points); ij++){  
    for (int k = 0; k < nwpoints; k++) {
      tensor_stream >> temp;
      hadtensor_w_array[k][ij]=temp;
    }
  }
}
//_________________________________________________________________________
