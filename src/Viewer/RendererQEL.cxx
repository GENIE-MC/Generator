//____________________________________________________________________________
/*!

\class    genie::RendererQEL

\brief    Draws the Feynman diagram for a QEL event

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  October 07, 2004

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
#include "Viewer/RendererQEL.h"

using namespace genie;

//______________________________________________________________________________
RendererQEL::RendererQEL() :
Renderer()
{

}
//______________________________________________________________________________
RendererQEL::~RendererQEL()
{

}
//______________________________________________________________________________
void RendererQEL::DrawDiagram(EventRecord * ev_rec) 
{
   LOG("gviewer", pINFO) << "Drawing Feynman diagram for a QEL event";

   fEmbeddedCanvas->GetCanvas()->cd();
   fEmbeddedCanvas->GetCanvas()->Clear();

   fEmbeddedCanvas->GetCanvas()->Range(0,0,200,100);
   fEmbeddedCanvas->GetCanvas()->SetFillColor(0);
   fEmbeddedCanvas->GetCanvas()->SetBorderMode(0);

   gStyle->SetLineWidth(2);
   TLatex t;
   TLatex tn;
   t.SetTextAlign(12);
   t.SetTextSize(0.04);
   tn.SetTextAlign(12);
   tn.SetTextSize(0.03);
   TLine * l, * ln;

   TEllipse * nucleus = 0;

   //Interaction * interaction = ev_rec->GetInteraction();

   // <-- leptons

   // incoming neutrino
   GHepParticle * nu = ev_rec->Probe();

   l = new TLine(15, 85, 50, 68);
   l->SetLineColor(2);
   l->Draw();
   t.DrawLatex(10,95, nu->Name().c_str());
   tn.DrawLatex(10,90, P4AsString(nu->Energy(),nu->Px(),nu->Py(),nu->Pz()));

   // final state primary lepton
   GHepParticle * fsl = ev_rec->FinalStatePrimaryLepton();

   l = new TLine(50, 68, 85, 85);
   l->SetLineColor(1);
   l->Draw();
   t.DrawLatex(80,95,fsl->Name().c_str());
   tn.DrawLatex(80,90, P4AsString(fsl->Energy(),fsl->Px(),fsl->Py(),fsl->Pz()));

   GHepParticle * p = 0;
   TIter rec_iter(ev_rec);

   // <-- nucleons

   while( (p = (GHepParticle *) rec_iter.Next()) ) {

      //-- if particle = nucleon
      if( p->PdgCode() == 2112 || p->PdgCode() == 2212 ) {

         // initial state
         if( p->Status() == 0 || p->Status() == 11) {

               l = new TLine(15, 15, 50, 40);
               l->SetLineColor(2);
               l->Draw();
               t.DrawLatex(15,40, p->Name().c_str());
               tn.DrawLatex(15,35, P4AsString(p->Energy(),p->Px(),p->Py(),p->Pz()));
         }
         // final state
         else if( p->Status() == 1 ) {

               l = new TLine(50, 40, 85, 15);
               l->SetLineColor(1);
               l->Draw();
               t.DrawLatex(60,40,p->Name().c_str());
               tn.DrawLatex(60,35, P4AsString(p->Energy(),p->Px(),p->Py(),p->Pz()));
         }
      }
   }

   rec_iter.Reset();

   // <-- nuclei

   while( (p = (GHepParticle *) rec_iter.Next()) ) {

      //-- if particle = nucleus
      if( p->PdgCode() == 0 ) {

         // initial state
         if( p->Status() == 0 ) {
                nucleus = new TEllipse(15, 15, 12, 12, 0, 360, 0);
                nucleus->SetLineColor(2);
                nucleus->Draw();
                t.DrawLatex(8, 8, p->Name().c_str());
         }
         // final state
         else if( p->Status() == 1 ) {
                nucleus = new TEllipse(85, 15, 12, 12, 0, 360, 0);
                nucleus->SetLineColor(1);
                nucleus->Draw();
                t.DrawLatex(78, 8, p->Name().c_str());
         }
      }
   }

   ln = new TLine(50, 68, 50, 40); ln->Draw(); // W boson
   ln->SetLineStyle(2);
   t.DrawLatex(58,56,"W^{+}");

   fEmbeddedCanvas->GetCanvas()->Update();
}
//______________________________________________________________________________

