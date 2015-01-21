//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Pauli Kehayias (Tufts Univ.)

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Jan 13, 2009 - PK
   This file was added in v2.5.1

*/
//____________________________________________________________________________

#include <iostream>

#include "validation/Merenyi/MerenyiAnalysis.h"

using namespace genie;
using namespace genie::vld_merenyi;

//____________________________________________________________________________
MerenyiAnalysis::MerenyiAnalysis() : 
binSizeE   (0.2), 
binSizeAng (0.1), 
binSizeMom (0.1)
{
  pEventTypes.resize(6);     
  nEventTypes.resize(8);
  nuECounts.resize(30);
  muAngCounts.resize(20);
  piMomCounts.resize(45);
  numPions = 0;
}
//____________________________________________________________________________
void MerenyiAnalysis::clearAll()
{
  numEvents = 0;
  numTossed = 0;
    
  pEventTypes.assign(6,0);
  nEventTypes.assign(8,0);

  ill1EventTypes = 0;
  ill2EventTypes = 0;
}
//____________________________________________________________________________
void MerenyiAnalysis::readData()
{
  clearAll();

  int numHadrons;
  double iHadron, hadPx, hadPy, hadPz, hadE;
  MerenyiNuEvent inEvent;

  in.open (myFile.c_str(), ifstream::in);
  in.ignore(256,'\n');

  while(in.good())
    {
      numEvents++;
  
      in.ignore(256,'\n');      //skip lines
      in.ignore(256,'\n');
  
      readNuEnergy(inEvent);            //read in neutrino energy
      in.ignore(256,'\n');      //skip the rest of the line

      readMuMom(inEvent);               //read in muon momentum
      in.ignore(256,'\n');      //skip the rest of the line

      in.ignore(256,'\n');      //skip one lines
 
      in >> numHadrons;
  
      for (int i = 1; i <= numHadrons; i++)
        {
          in >> iHadron >> hadPx >> hadPy >> hadPz >> hadE;
          inEvent.addParticle(iHadron, hadPx, hadPy, hadPz, hadE);
          in.ignore(256,'\n');
        }
      in.ignore(256,'\n'); 
 
      in >> numHadrons;
  
      for (int i = 1; i <= numHadrons; i++)
        {
          in >> iHadron >> hadPx >> hadPy >> hadPz >> hadE;
          inEvent.addParticle(iHadron, hadPx, hadPy, hadPz, hadE);
          in.ignore(256,'\n');
        }
      in.ignore(256,'\n');
 
      processEvent(inEvent);            //process event once all params are collected
                                                                        
      inEvent.zeroAll();
    }

  in.close();
  cout << "Done processing file" << endl << endl;
}
//____________________________________________________________________________
void MerenyiAnalysis::processEvent(MerenyiNuEvent myEvent)
{
  int myBin;
  double muAng;

  myBin = (int) (myEvent.getNuE() / binSizeE);
  nuECounts[myBin]++;           //store neutrino energy

  muAng = myEvent.getLepAngCos();
  myBin = (int) (muAng / binSizeAng + 10);
  muAngCounts[myBin]++; //store muon angle cosine

  //classifyHadrons is by reference, though processEvent is
  //by value. processEvent doesn't need to be by reference
  //because the above loop (readData) is finished with
  //the event and processEvent then gets its reaction

  classifyHadrons(myEvent);
  checkCandP(myEvent);                  //candidate proton algorithm

  if (myEvent.getNumGam() == 0)
    getReaction(myEvent);               //classify the event
  else
    numTossed++;
}
//____________________________________________________________________________
void MerenyiAnalysis::getReaction(MerenyiNuEvent myEvent)
{
  int hadCharge = myEvent.getNumPiP() + myEvent.getNumStrP() + myEvent.getNumP() - myEvent.getNumPiM() - myEvent.getNumStrM();

  switch (hadCharge)
    {
    case 2:             //proton target
      pTarget(myEvent);
      break;

    case 1:             //neutron target
      nTarget(myEvent);
      break;

    case 0:
    case 3:
    case 4:
      illTarget1(myEvent);
      break;

    case -1:
    case 5:
    case -2:
      illTarget2();
      break;
    default:                                                            
      cout << "Unknown charge " << hadCharge << " event: #" << numEvents << endl;
    }
}
//____________________________________________________________________________
void MerenyiAnalysis::pTarget(MerenyiNuEvent myEvent)
{
  int nPi0 = myEvent.getNumPi0(), nPiP = myEvent.getNumPiP(), nPiM = myEvent.getNumPiM();
  int nN = myEvent.getNumN(), nP = myEvent.getNumP();
  int nStr0 = myEvent.getNumStr0(), nStrP = myEvent.getNumStrP(), nStrM = myEvent.getNumStrM();

  //event type 0: nu p --> mu- pi+ p
  if (nPi0 == 0 && nStr0 == 0 && nPiP == 1 && nPiM == 0 && nStrP == 0 && nStrM == 0 && nN == 0 && nP == 1)
    pEventTypes[0]++;

  //event type 1: nu p --> mu- pi+ p l*pi0              l >= 1
  else if (nPi0 >= 1 && nStr0 == 0 && nPiP == 1 && nPiM == 0 && nStrP == 0 && nStrM == 0 && nN == 0 && nP == 1)
    pEventTypes[1]++;
    
  //event type 2: nu p --> mu- 2*pi+ m*pi0 p            m >= 0
  else if (nPi0 >= 0 && nStr0 == 0 && nPiP == 2 && nPiM == 0 && nStrP == 0 && nStrM == 0 && nN == 1 && nP == 0)
    pEventTypes[2]++;
 
  //event type 3: nu p --> mu- (p or n) >= 3 pi's
  else if ((nPi0 + nPiP + nPiM) >= 1 && nStr0 == 0 && nStrP == 0 && nStrM == 0 && (nN + nP) == 1)
    pEventTypes[3]++;
     
  //event type 4: nu p --> mu- + any strange particles (K+, K-, K0L, K0S, Lambda0)
  else if ((nStr0 + nStrP + nStrM) > 0)
    pEventTypes[4]++;
 
  else
    {cout << "Unknown proton target event: #" << numEvents << endl;
      for (unsigned int i = 0; i<myEvent.getHadList().size(); i++)
        cout << myEvent.getHadList()[i].GetPdgCode() << endl;

      cout << nP << " = num protons" << endl;
    }
}
//____________________________________________________________________________
void MerenyiAnalysis::nTarget(MerenyiNuEvent myEvent)
{
  int nPi0 = myEvent.getNumPi0(), nPiP = myEvent.getNumPiP(), nPiM = myEvent.getNumPiM();
  int nN = myEvent.getNumN(), nP = myEvent.getNumP();
  int nStr0 = myEvent.getNumStr0(), nStrP = myEvent.getNumStrP(), nStrM = myEvent.getNumStrM();

  //event type 0: nu n --> mu- p
  if (nPi0 == 0 && nStr0 == 0 && nPiP == 0 && nPiM == 0 && nStrP == 0 && nStrM == 0 && nN == 0 && nP == 1)
    nEventTypes[0]++;

  //event type 1: nu n --> mu- p pi0
  else if (nPi0 == 1 && nStr0 == 0 && nPiP == 0 && nPiM == 0 && nStrP == 0 && nStrM == 0 && nN == 0 && nP == 1)
    nEventTypes[1]++;

  //event type 2: nu n --> mu- pi+ n
  else if (nPi0 == 0 && nStr0 == 0 && nPiP == 1 && nPiM == 0 && nStrP == 0 && nStrM == 0 && nN == 1 && nP == 0)
    nEventTypes[2]++;

  //event type 3: nu n --> mu- pi- pi+ p
  else if (nPi0 == 0 && nStr0 == 0 && nPiP == 1 && nPiM == 1 && nStrP == 0 && nStrM == 0 && nN == 0 && nP == 1)
    nEventTypes[3]++;

  //event type 4: nu n --> mu- pi0 pi0 p
  else if (nPi0 == 2 && nStr0 == 0 && nPiP == 0 && nPiM == 0 && nStrP == 0 && nStrM == 0 && nN == 0 && nP == 1)
    nEventTypes[4]++;

  //event type 5: nu n --> mu- pi+ n l*pi0              l >= 1
  else if (nPi0 >= 1 && nStr0 == 0 && nPiP == 1 && nPiM == 0 && nStrP == 0 && nStrM == 0 && nN == 1 && nP == 0)
    nEventTypes[5]++;

  //event type 6: nu n --> mu- (p or n) >= 3 pi's
  else if ((nPi0 + nPiP + nPiM) >= 1 && nStr0 == 0 && nStrP == 0 && nStrM == 0 && (nN + nP) == 1)
    nEventTypes[6]++;

  //event type 7: nu n --> mu- + any strange particles (K+, K-, K0L, K0S, Lambda0)
  else if ((nStr0 + nStrP + nStrM) > 0)
    nEventTypes[7]++;

  else
    cout << "Unknown neutron target event: #" << numEvents << endl;
}
//____________________________________________________________________________
void MerenyiAnalysis::illTarget1(MerenyiNuEvent myEvent)
{
  int nPi0 = myEvent.getNumPi0(), nPiP = myEvent.getNumPiP(), nPiM = myEvent.getNumPiM();
  int nN = myEvent.getNumN(), nP = myEvent.getNumP();
  int nStr0 = myEvent.getNumStr0(), nStrP = myEvent.getNumStrP(), nStrM = myEvent.getNumStrM();

  //test if this is a proton event type 5
  //event type 5: nu p --> mu- pi- p
  if (nPi0 == 0 && nStr0 == 0 && nPiP == 0 && nPiM == 1 && nStrP == 0 && nStrM == 0 && nN == 0 && nP == 1)
    pEventTypes[5]++;
  else
    ill1EventTypes++;
}
//____________________________________________________________________________
void MerenyiAnalysis::illTarget2()
{
  ill2EventTypes++;
}
//____________________________________________________________________________
void MerenyiAnalysis::readNuEnergy(MerenyiNuEvent& myEvent)
{
  double nuE;

  for (int i = 1; i <= 4; i++)
    in >> nuE;

  myEvent.setNuE(nuE);
}
//____________________________________________________________________________
void MerenyiAnalysis::readMuMom(MerenyiNuEvent& myEvent)
{
  double px, py, pz, E;

  in >> px;             //ignore the first entry (13)
  in >> px >> py >> pz >> E;
  myEvent.setLepP(px, py, pz, E);
}
//____________________________________________________________________________
void MerenyiAnalysis::findPiMom(TParticle pion)
{
  double piMom;;
  int myBin;

  numPions++;
  piMom = pion.P();
  myBin = (int) (piMom / binSizeMom);
  piMomCounts[myBin]++;
}
//____________________________________________________________________________
void MerenyiAnalysis::dispPops()
{
  cout << "Fractional populations:" << endl;
  int numCC = numEvents - numTossed;

  cout << (double) nEventTypes[0] / numCC << endl;
  cout << (double) pEventTypes[0] / numCC << endl;
  cout << (double) nEventTypes[1] / numCC << endl;
  cout << (double) nEventTypes[2] / numCC << endl;
  cout << (double) nEventTypes[3] / numCC << endl;
  cout << (double) nEventTypes[4] / numCC << endl;
  cout << (double) pEventTypes[1] / numCC << endl;
  cout << (double) nEventTypes[5] / numCC << endl;
  cout << (double) pEventTypes[2] / numCC << endl;
  cout << (double) (nEventTypes[6] + pEventTypes[3]) / numCC << endl;
  cout << (double) (nEventTypes[7] + pEventTypes[4]) / numCC << endl;
  cout << (double) pEventTypes[5] / numCC << endl;
  cout << (double) ill1EventTypes / numCC << endl;
  cout << (double) ill2EventTypes / numCC << endl;
}
//____________________________________________________________________________
void MerenyiAnalysis::dispNuEHist()
{
  cout << "Neutrino energy histogram: " << endl;
  for (unsigned int i = 0; i < nuECounts.size(); i++)
    cout << (double) nuECounts[i] / numEvents << endl;
}
//____________________________________________________________________________
void MerenyiAnalysis::dispMuAngHist()
{
  cout << "Muon angle cosine histogram: " << endl;
  for (unsigned int i = 0; i < muAngCounts.size(); i++)
    cout << (double) muAngCounts[i] / numEvents << endl;
}
//____________________________________________________________________________
void MerenyiAnalysis::dispPiMomHist()
{
  cout << "Pi+- momentum histogram: " << endl;
  for (unsigned int i = 0; i < piMomCounts.size(); i++)
    cout << (double) piMomCounts[i] / numPions << endl;
}
//____________________________________________________________________________
void MerenyiAnalysis::printPops()
{
  int numCC = numEvents - numTossed;
  ofstream pops("populations.txt");

  pops << numEvents << " " << numCC << endl;

  pops << (double) nEventTypes[0] / numCC << endl;
  pops << (double) pEventTypes[0] / numCC << endl;
  pops << (double) nEventTypes[1] / numCC << endl;
  pops << (double) nEventTypes[2] / numCC << endl;
  pops << (double) nEventTypes[3] / numCC << endl;
  pops << (double) nEventTypes[4] / numCC << endl;
  pops << (double) pEventTypes[1] / numCC << endl;
  pops << (double) nEventTypes[5] / numCC << endl;
  pops << (double) pEventTypes[2] / numCC << endl;
  pops << (double) (nEventTypes[6] + pEventTypes[3]) / numCC << endl;
  pops << (double) (nEventTypes[7] + pEventTypes[4]) / numCC << endl;
  pops << (double) pEventTypes[5] / numCC << endl;
  pops << (double) ill1EventTypes / numCC << endl;
  pops << (double) ill2EventTypes / numCC << endl;
  pops.close();
}
//____________________________________________________________________________
void MerenyiAnalysis::printNuEHist()
{
  ofstream nuEs("nu_energies.txt");

  nuEs << numEvents << endl;
  for (unsigned int i = 0; i < nuECounts.size(); i++)
    nuEs << (double) nuECounts[i] / numEvents << endl;

  nuEs.close();
}
//____________________________________________________________________________
void MerenyiAnalysis::printMuAngHist()
{
  ofstream angs("muon_angles.txt");

  angs << numEvents << endl;
  for (unsigned int i = 0; i < muAngCounts.size(); i++)
    angs << (double) muAngCounts[i] / numEvents << endl;

  angs.close();
}
//____________________________________________________________________________
void MerenyiAnalysis::printPiMomHist()
{
  ofstream piMoms("pion_momenta.txt");

  piMoms << numPions << endl;
  for (unsigned int i = 0; i < piMomCounts.size(); i++)
    piMoms << (double) piMomCounts[i] / numPions << endl;

  piMoms.close();
}
//____________________________________________________________________________
vector <double> MerenyiAnalysis::getPops()
{
  int numCC = numEvents - numTossed;
  vector <double> pops;

  pops.push_back((double) nEventTypes[0] / numCC);
  pops.push_back((double) pEventTypes[0] / numCC);
  pops.push_back((double) nEventTypes[1] / numCC);
  pops.push_back((double) nEventTypes[2] / numCC);
  pops.push_back((double) nEventTypes[3] / numCC);
  pops.push_back((double) nEventTypes[4] / numCC);
  pops.push_back((double) pEventTypes[1] / numCC);
  pops.push_back((double) nEventTypes[5] / numCC);
  pops.push_back((double) pEventTypes[2] / numCC);
  pops.push_back((double) (nEventTypes[6] + pEventTypes[3]) / numCC);
  pops.push_back((double) (nEventTypes[7] + pEventTypes[4]) / numCC);
  pops.push_back((double) pEventTypes[5] / numCC);
  pops.push_back((double) ill1EventTypes / numCC);
  pops.push_back((double) ill2EventTypes / numCC);

  return pops;
}
//____________________________________________________________________________
vector <double> MerenyiAnalysis::getNuEHist()
{
  vector <double> nuE;

  for (unsigned int i = 0; i < 30; i++)
    nuE.push_back((double) nuECounts[i] / numEvents);

  return nuE;
}
//____________________________________________________________________________
vector <double> MerenyiAnalysis::getMuAngHist()
{
  vector <double> muAng;

  for (int i = 0; i < 20; i++)
    muAng.push_back((double) muAngCounts[i] / numEvents);

  return muAng;
}
//____________________________________________________________________________
vector <double> MerenyiAnalysis::getPiMomHist()
{
  vector <double> piMom;

  for (int i = 0; i < 40; i++)
    piMom.push_back((double) piMomCounts[i] / numPions);

  return piMom;
}
//____________________________________________________________________________
double MerenyiAnalysis::getNuEi(int i)
{
  vector <double> nuE = getNuEHist();
  return nuE[i];
}
//____________________________________________________________________________
double MerenyiAnalysis::getMuAngi(int i)
{
  vector<double> muAng = getMuAngHist();
  return muAng[i];
}
//____________________________________________________________________________
double MerenyiAnalysis::getPiMomi(int i)
{
  vector<double> piMom = getPiMomHist();
  return piMom[i];
}
//____________________________________________________________________________

