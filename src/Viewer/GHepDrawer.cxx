//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>, Rutherford Lab.
         November 30, 2007

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Nov 30, 2007 - CA
   Was added in 2.0.1, intended to replace an older set of GHEP rendering
   classes. No GHEP drawing at that time.

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

#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Messenger/Messenger.h"
#include "Viewer/GHepDrawer.h"

using namespace genie;

//______________________________________________________________________________
GHepDrawer::GHepDrawer() 
{

}
//______________________________________________________________________________
GHepDrawer::~GHepDrawer()
{

}
//______________________________________________________________________________
void GHepDrawer::SetEmbeddedCanvas(TRootEmbeddedCanvas * ec)
{
  fEmbeddedCanvas = ec;
}
//______________________________________________________________________________
void GHepDrawer::Draw(EventRecord * /*ev_rec*/) 
{
   LOG("gviewer", pINFO) << "Drawing input event";

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

