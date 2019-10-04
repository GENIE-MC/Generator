//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Tomek Golan <tomasz.golan@uwr.edu.pl>, FNAL/Rochester
         Steve Dytman <dytman+@pitt.edu>, Pittsburgh Univ.
         Josh Kleckner <jok84@pitt.edu>, Pittsburgh Univ.

         Auguest 20, 2016

 Calculate the cross section attenuation factor due to medium
corrections for NA FSI from Pandharipande, Pieper (Phys Rev (2009).
Paper gives prescription for Pauli blocking and average nuclear
potential in (e,e'p).  GENIE code was adapted from NuWro implementation.

 Important revisions after version 2.12.0 :
 @ Aug, 2016 - TK
   adapted to GENIE from NuWro.  Use free NN xs.
 @ Aug, 2017 - SD, JK
   Original code stores values in rotating buffer.  This won't be accurate
   for many problems, esp. heavy targets and multiple nuclei.  New code has
   lookup tables in probe KE and nuclear density (rho) stored in text files
   for He4, C12, Ca40, Fe56, Sn120, and U238.  Use values from the text
   files for KE and rho, interpolation in A.
*/
//____________________________________________________________________________
#include "Physics/HadronTransport/INukeNucleonCorr.h"
#include "Physics/HadronTransport/INukeUtils2018.h"
#include "Physics/HadronTransport/INukeHadroData2018.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Messenger/Messenger.h"
using namespace genie;

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <cstdlib>

#include <TGraph.h>
using namespace std;


INukeNucleonCorr* INukeNucleonCorr::fInstance = NULL; // initialize instance with NULL


const int NRows     = 200;
const int NColumns  =  17;

string genie_dir(std::getenv("GENIE"));
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

//! \f$m^* (k,\rho) = m \frac{(\Lambda^2 + k^2)^2}{\Lambda^2 + k^2)^2 - 2\Lambda^2\beta m}\f$
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
    vector<vector<double> > HeliumValues;
    vector<vector<double> > CarbonValues;
    vector<vector<double> > CalciumValues;
    vector<vector<double> > IronValues;
    vector<vector<double> > TinValues;
    vector<vector<double> > UraniumValues;
    vector<vector<double> > clear;

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
    double C3CM = INukeHadroData2018::Instance()->IntBounce (&incomingParticle, targetPdg, pdg, kIHNFtElas);

    // generate kinematics
    utils::intranuke2018::TwoBodyKinematics (mass, targetMass, p, target, outNucl1, outNucl2, C3CM, RemnP4);

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


// This function interpolates and returns correction values
//
double INukeNucleonCorr :: getAvgCorrection(double rho, double A, double ke)
{
  //Read in energy and density to determine the row and column of the correction table - adjust for variable binning - throws away some of the accuracy
   int Column = round(rho*100);
   if(rho<.01) Column = 1;
   if (Column>=NColumns) Column = NColumns-1;
   int Row = 0;
   if(ke<=.002) Row = 1;
   if(ke>.002&&ke<=.1) Row = round(ke*1000.);
   if(ke>.1&&ke<=.5) Row = round(.1*1000.+(ke-.1)*200);
   if(ke>.5&&ke<=1) Row = round(.1*1000.+(.5-.1)*200+(ke-.5)*40);
   if(ke>1) Row = NRows-1;
   //LOG ("INukeNucleonCorr",pNOTICE)
   //  << "row, column = " << Row << "   " << Column;
  //If the table of correction values has already been created
  // return a value. Else, interpolate the needed correction table//
//  static double cache[NRows][NColumns] = {{-1}};
  static bool ReadFile = false;
  if( ReadFile == true ) {
   int Npoints = 6;
   TGraph * Interp = new TGraph(Npoints);
   //LOG("INukeNucleonCorr",pNOTICE)
   //  << HeliumValues[Row][Column];
   Interp->SetPoint(0,4,HeliumValues[Row][Column]);
   Interp->SetPoint(1,12,CarbonValues[Row][Column]);
   Interp->SetPoint(2,40,CalciumValues[Row][Column]);
   Interp->SetPoint(3,56,IronValues[Row][Column]);
   Interp->SetPoint(4,120,TinValues[Row][Column]);
   Interp->SetPoint(5,238,UraniumValues[Row][Column]);

   //   Interpolated[e][r] = Interp->Eval(A);
	//	LOG("INukeNucleonCorr",pNOTICE)
	//	  << "e,r,value= " << e << "   " << r << "   " << Interpolated[e][r];
   double returnval = Interp->Eval(A);
   delete Interp;
   LOG("INukeNucleonCorr",pINFO)
      << "Nucleon Corr interpolated correction factor = "
      << returnval  //cache[Row][Column]
      << " for rho, KE, A= "<<  rho << "  " << ke << "   " << A;
    //    return cache[Row][Column];
   return returnval;
  } else {
    //Reading in correction files//
    //    string dir = genie_dir + string("/data/evgen/nncorr/");

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
      << "Nucleon Corr interpolation files read in successfully";
    //get interpolated value for first event.
    int Npoints = 6;
    TGraph * Interp = new TGraph(Npoints);
    //LOG("INukeNucleonCorr",pNOTICE)
    //  << HeliumValues[Row][Column];
    Interp->SetPoint(0,4,HeliumValues[Row][Column]);
    Interp->SetPoint(1,12,CarbonValues[Row][Column]);
    Interp->SetPoint(2,40,CalciumValues[Row][Column]);
    Interp->SetPoint(3,56,IronValues[Row][Column]);
    Interp->SetPoint(4,120,TinValues[Row][Column]);
    Interp->SetPoint(5,238,UraniumValues[Row][Column]);

    //	LOG("INukeNucleonCorr",pNOTICE)
    //	  << "Row,Column,value= " << Row << "   " << Column << "   " << Interp->Eval(A);
    double returnval = Interp->Eval(A);
    delete Interp;
    ReadFile = true;
    LOG("INukeNucleonCorr",pINFO)
      << "Nucleon Corr interpolated correction factor = "
      << returnval
      << " for rho, KE, A= "<<  rho << "  " << ke << "   " << A;
    return returnval;
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
