//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab
*/
//____________________________________________________________________________

#include <sstream>
#include <string>
#include <vector>

#include <TRootEmbeddedCanvas.h>
#include <TLorentzVector.h>
#include <TLine.h>
#include <TEllipse.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TStyle.h>

#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Messenger/Messenger.h"
#include "Tools/Masterclass/FastSimCherenkov.h"

using namespace genie;
using namespace genie::masterclass;

//______________________________________________________________________________
FastSimCherenkov::FastSimCherenkov() 
{

}
//______________________________________________________________________________
FastSimCherenkov::~FastSimCherenkov()
{

}
//______________________________________________________________________________
void FastSimCherenkov::SetEmbeddedCanvas(TRootEmbeddedCanvas * ec)
{
  fEmbeddedCanvas = ec;
}
//______________________________________________________________________________
void FastSimCherenkov::Draw(EventRecord * /*event*/) 
{
   LOG("MasterClass", pINFO) << "Drawing input event";

   fEmbeddedCanvas->GetCanvas()->cd();
   fEmbeddedCanvas->GetCanvas()->Clear();

   fEmbeddedCanvas->GetCanvas()->Range(0,0,200,100);
   fEmbeddedCanvas->GetCanvas()->SetFillColor(0);
   fEmbeddedCanvas->GetCanvas()->SetBorderMode(0);

   //
   // ...
   //

   fEmbeddedCanvas->GetCanvas()->Update();
}
//______________________________________________________________________________

