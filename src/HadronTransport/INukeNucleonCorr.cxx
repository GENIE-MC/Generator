#include "HadronTransport/INukeNucleonCorr.h"
#include "HadronTransport/INukeUtils2015.h"
#include "HadronTransport/INukeHadroData2015.h"
#include "PDG/PDGLibrary.h"
#include "Conventions/Units.h"
#include "Numerical/RandomGen.h"
#include "Messenger/Messenger.h"
#include <TSystem.h>
using namespace genie;

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <TGraph.h>
using namespace std;


INukeNucleonCorr* INukeNucleonCorr::fInstance = NULL; // initialize instance with NULL


const int NColumns  = 18;
const int NRows = 205;
       
string genie_dir = string(gSystem->Getenv("GENIE"));
//string directory = genie_dir + string("/data/evgen/nncorr/");
string dir = genie_dir + string("/data/evgen/nncorr/");
string infile;
vector<vector<double> > infile_values;
vector<string> comments;



// ----- STATIC CONSTANTS ----- //

const unsigned int INukeNucleonCorr::fRepeat = 1000; // how many times kinematics should be generated to get avg corr

const double INukeNucleonCorr::fRho0  = 0.16; // fm^-3

const double INukeNucleonCorr::fBeta1   = -116.00 / fRho0 / 1000.0;          // converted to GeV
const double INukeNucleonCorr::fLambda0 =    3.29 / (units::fermi);          // converted to GeV
const double INukeNucleonCorr::fLambda1 =  -0.373 / (units::fermi) / fRho0;  // converted to GeV



// ----- CALCULATIONS ----- //

//! \f$m^* (k,\rhp) = m \frac{(\Lambda^2 + k^2)^2}{\Lambda^2 + k^2)^2 - 2\Lambda^2\beta m}\f$
double INukeNucleonCorr::mstar (const double rho, const double k2)
{
  // density [fm^-3], momentum square [GeV^2]
  
  static const double m = (PDGLibrary::Instance()->Find(kPdgProton)->Mass() + 
                           PDGLibrary::Instance()->Find(kPdgNeutron)->Mass()) / 2.0;
  
  const double L = lambda (rho); // potential coefficient lambda
  const double B =   beta (rho); // potential coefficient beta
    
  const double L2 = L * L; // lambda^2 used twice in equation
  
  const double num = (L2 + k2) * (L2 + k2); // numerator
    
  return m * num / (num - 2.0 * L2 * B * m);
}

//! \f$k_F = (\frac{3}{2}\pi^2\rho)^{1/3}\f$
double INukeNucleonCorr :: localFermiMom (const double rho, const int A, const int Z, const int pdg)
{
  static double factor = 3.0 * M_PI * M_PI / 2.0;
    
  return pdg == kPdgProton ? pow (factor * rho * Z / A, 1.0 / 3.0) / (units::fermi) :
                             pow (factor * rho * (A - Z) / A, 1.0 / 3.0) / (units::fermi);
}

//! generate random momentum direction and return 4-momentum of target nucleon
TLorentzVector INukeNucleonCorr :: generateTargetNucleon (const double mass, const double fermiMom)
{
  RandomGen * rnd = RandomGen::Instance();
    
  // get random momentum direction
  const double costheta = 2.0 * rnd->RndGen().Rndm() - 1.0;  // random cos (theta)
  const double sintheta = sqrt (1.0 - costheta * costheta);  // sin (theta)  
  const double      phi = 2.0 * M_PI * rnd->RndGen().Rndm(); // random phi

  // set nucleon 4-momentum
  const double p = rnd->RndGen().Rndm() * fermiMom; // random nucleon momentum up to Fermi level

  const TVector3   p3 = TVector3 (p * sintheta * cos (phi), p * sintheta * sin (phi), p * costheta); // 3-momentum
  const double energy = sqrt (p3.Mag2() + mass * mass); // energy

  return TLorentzVector (p3, energy);
}

//! calculate correction given by eq. 2.9
double INukeNucleonCorr :: getCorrection (const double mass, const double rho,
                                          const TVector3 &k1, const TVector3 &k2,
                                          const TVector3 &k3, const TVector3 &k4)
{
  const double num = (k1 - k2).Mag() * mstar (rho, (k3.Mag2() + k4.Mag2()) / 2.0) / mass / mass;
  const double den = (k1 * (1.0 / mstar (rho, k1.Mag2())) - k2 * (1.0 / mstar (rho, k2.Mag2()))).Mag();
  
  return num / den;
}

//! generate kinematics fRepeat times to calculate average correction
double INukeNucleonCorr :: AvgCorrection (const double rho, const int A, const int Z, const int pdg, const double Ek)
{

  RandomGen * rnd = RandomGen::Instance();

  setFermiLevel (rho, A, Z); // set Fermi momenta for protons and neutrons
  
  const double mass   = PDGLibrary::Instance()->Find(pdg)->Mass(); // mass of incoming nucleon
  const double energy = Ek + mass;

  TLorentzVector p (0.0, 0.0, sqrt (energy * energy - mass * mass), energy); // incoming particle 4-momentum
  GHepParticle incomingParticle (pdg, kIStUndefined, -1,-1,-1,-1, p, TLorentzVector ()); // for IntBounce
  
  double corrPauliBlocking = 0.0; // correction coming from Pauli blocking
  double     corrPotential = 0.0; // correction coming from potential
  
  for (unsigned int i = 0; i < fRepeat; i++) // generate kinematics fRepeat times to get avg corrections
  {
    // get proton vs neutron randomly based on Z/A
    const int targetPdg = rnd->RndGen().Rndm() < (double) Z / A ? kPdgProton : kPdgNeutron;
    
    const double targetMass = PDGLibrary::Instance()->Find(targetPdg)->Mass(); // set nucleon mass
    
    const TLorentzVector target = generateTargetNucleon (targetMass, fermiMomentum (targetPdg)); // generate target nucl
    
    TLorentzVector outNucl1, outNucl2, RemnP4; // final 4-momenta
    
    // random scattering angle
    double C3CM = INukeHadroData2015::Instance()->IntBounce (&incomingParticle, targetPdg, pdg, kIHNFtElas);
    
    // generate kinematics
    utils::intranuke2015::TwoBodyKinematics (mass, targetMass, p, target, outNucl1, outNucl2, C3CM, RemnP4);

    // update Pauli blocking correction
    corrPauliBlocking += (outNucl1.Vect().Mag() > fermiMomentum (pdg) and outNucl2.Vect().Mag() > fermiMomentum (targetPdg));
    
    // update potential-based correction
    corrPotential += getCorrection (mass, rho, p.Vect(), target.Vect(), outNucl1.Vect(), outNucl2.Vect());
  }
  
  corrPauliBlocking /= fRepeat;
      corrPotential /= fRepeat;
      
      return corrPotential * corrPauliBlocking;
     
}


//This function reads the correction files that will be used to interpolate new correction values for some target//
void read_file(string rfilename) 
{  
  ifstream file;
  file.open((char*)rfilename.c_str(), ios::in);
  
  if (file.is_open()) 
    {
    string line;
    int cur_line = 0;            
    while (getline(file,line)) 
      { 
      if (line[0]=='#') 
	{ 
        comments.push_back(line); 
	} 
      else {
        vector<double> temp_vector;
        istringstream iss(line);            
        string s; 
        for (int i=0; i<18; i++) 
	  {
          iss >> s;
          temp_vector.push_back(atof(s.c_str()));
        }
        infile_values.push_back(temp_vector);                  
        cur_line++;
      }
    }
    LOG("INukeNucleonCorr",pNOTICE) << "Successful open file" << rfilename << "\n";

  } 
  else {
    LOG("INukeNucleonCorr",pNOTICE) << "Could not open " << rfilename << "\n";

  }
  file.close();
}


//This function interpolates and returns correction values//
double INukeNucleonCorr :: getAvgCorrection(double rho, double A, double ke)                                                                               
{                                                         
  //Read in energy and density to determine the row and column of the correction table//                                                                                                           
  int Column = round(rho*100) + 1;
  int Row = -1;
  for(int e = 0; e < (ke*1000); e++)
    {
      if(e > 100){e+=4;}
      if(e > 500){e+=15;}
      Row++;
    }
  //If the table of correction values has already been created, return a value. Else, interpolate the needed correction table//
  static double cache[NRows][NColumns] = {{-1}};
  static bool ReadFile;
  if(ReadFile == true) {
    LOG("INukeNucleonCorr",pDEBUG)  "Nucleon Corr interpolated value for correction factor = "<< cache[Row][Column] << " for rho, KE, A= "<<  rho << "  " << ke << "   " << A << "\n";
    return cache[Row][Column];}
  else{
    //Reading in correction files//
    //    string dir = genie_dir + string("/data/evgen/nncorr/");
    vector<vector<double> > HeliumValues;
    vector<vector<double> > CarbonValues;    
    vector<vector<double> > CalciumValues;
    vector<vector<double> > IronValues;
    vector<vector<double> > TinValues;
    vector<vector<double> > UraniumValues;
    vector<vector<double> > clear;

    read_file(dir+"NNCorrection_2_4.txt");
    HeliumValues = infile_values;
    infile_values = clear;

    read_file(dir+"NNCorrection_6_12.txt");
    CarbonValues = infile_values;
    infile_values = clear;

    read_file(dir+"NNCorrection_20_40.txt");
    CalciumValues = infile_values;
    infile_values = clear;

    read_file(dir+"NNCorrection_26_56.txt");
    IronValues = infile_values;
    infile_values = clear;

    read_file(dir+"NNCorrection_50_120.txt");
    TinValues = infile_values;
    infile_values = clear;

    read_file(dir+"NNCorrection_92_238.txt");
    UraniumValues = infile_values;
    infile_values = clear;

    LOG("INukeNucleonCorr",pNOTICE)
      "Nucleon Corr interpolation files read in successfully";
    double Interpolated[205][18];
    for(int i = 0; i < 205; i++){Interpolated[i][0] = (i*0.001);}
    for(int i = 0; i < 18; i++){Interpolated[0][i] = (i*0.01);}
    
    const int Npoints = 6;
    //Interpolate correction values at every energy and density on the correction table//
    for(int e = 0; e < 201; e++){
      for(int r = 1; r < 19; r++){
	TGraph * Interp = new TGraph(Npoints);
	Interp->SetPoint(0,4,HeliumValues[e][r]);
	Interp->SetPoint(1,12,CarbonValues[e][r]);
	Interp->SetPoint(2,40,CalciumValues[e][r]);
	Interp->SetPoint(3,56,IronValues[e][r]);
	Interp->SetPoint(4,120,TinValues[e][r]);
	Interp->SetPoint(5,238,UraniumValues[e][r]);


	Interpolated[e][r] = Interp->Eval(A);
	delete Interp;
      }
    }

    //Save the interpolated values and return the needed correction//
    for(int e = 0; e < 205; e++){
      for(int r = 1; r < 18; r++){
        cache[e][r] = Interpolated[e][r]; 
      }
    }
    ReadFile = true;
 
    return cache[Row][Column];
  }
}

//This function outputs new correction files a new target if needed//
void  INukeNucleonCorr :: OutputFiles(int A, int Z)
{
  string outputdir = genie_dir + string("/data/evgen/nncorr/");
  double pdgc;
  string file;
  string header;

 
    file = outputdir+"NNCorrection.txt";
    header = "##Correction values for protons (density(horizontal) and energy(vertical))";
    pdgc = 2212;


  double output[1002][18];
  //label densities (columns)//	  
  for(int r = 1; r < 18; r++)
    { double rho = (r-1)*0.01; 
      output[0][r] = rho;}
 
  //label energies (rows) //
  for(int e = 1; e < 1002; e++)
    { double energy = (e-1)*0.001;
      output[e][0] = energy;}
  
  //loop over each energy and density to get corrections and build the correction table//
  for(int e = 1; e < 1002; e++){
    for(int r = 1; r < 18; r++){
      double energy  = (e-1)*0.001;
      double density = (r-1)*0.01;
      double correction = INukeNucleonCorr::getInstance()-> AvgCorrection (density, A, Z, pdgc, energy);
      output[e][r] = correction;
    }
  }
  //output the new correction table //
  ofstream outfile;
  outfile.open((char*)file.c_str(), ios::trunc);
  if (outfile.is_open()) {
    
    outfile << header << endl;
    outfile << left << setw(12) << "##";
    for (int e = 0; e < 1002; e++) {    
      for (int r = 0; r < 18; r++) {  
	if((e == 0) && (r == 0)) {}
	else{outfile << left << setw(12) << output[e][r];}
      } 
      outfile << endl;
      
    }
  }
  outfile.close();
}
    
  


