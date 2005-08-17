//____________________________________________________________________________
/*!

\class   genie::EventRecord

\brief   Generated Event Record. It is a GHepRecord object that can accept /
         be visited by EventRecordVisitorI objects (event generation modules).
         All the other important container manipulation methods are defined
         at the base GHepRecord record.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 1, 2004

*/
//____________________________________________________________________________

#include <iomanip>

#include "EVGCore/EventRecord.h"
#include "EVGCore/EventRecordVisitorI.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"

using namespace genie;

using std::endl;
using std::setw;
using std::setprecision;
using std::setfill;
using std::ios;

ClassImp(EventRecord)

//____________________________________________________________________________
namespace genie {
 ostream & operator<< (ostream& stream, const EventRecord & event_record)
 {
   event_record.Print(stream);

   return stream;
 }
}
//___________________________________________________________________________
EventRecord::EventRecord() :
GHepRecord()
{
  fInteraction = 0;

  this->SetOwner(true);
}
//___________________________________________________________________________
EventRecord::EventRecord(int size) :
GHepRecord(size)
{
  this->SetOwner(true);
}
//___________________________________________________________________________
EventRecord::EventRecord(const EventRecord & record) :
GHepRecord(record)
{
  this->SetOwner(true);
}
//___________________________________________________________________________
EventRecord::~EventRecord()
{
  this->Delete();
}
//___________________________________________________________________________
void EventRecord::AcceptVisitor(EventRecordVisitorI * visitor)
{
  visitor->ProcessEventRecord(this);
}
//___________________________________________________________________________
void EventRecord::Print(ostream & stream) const
{
  stream << "\n\n |";
  stream << setfill('-') << setw(104) << "|";

  stream << "\n |";
  stream << setfill(' ') << setw(6)  << "Idx | "
         << setfill(' ') << setw(11) << "Name | "
         << setfill(' ') << setw(6)  << "Ist | "
         << setfill(' ') << setw(13) << "PDG | "
         << setfill(' ') << setw(12) << "Mother  | "
         << setfill(' ') << setw(12) << "Daughter  | "
         << setfill(' ') << setw(9)  << "Px | "
         << setfill(' ') << setw(9)  << "Py | "
         << setfill(' ') << setw(9)  << "Pz | "
         << setfill(' ') << setw(9)  << "E  | "
         << setfill(' ') << setw(9)  << "m  | ";

  stream << "\n |";
  stream << setfill('-') << setw(104) << "|";

  GHepParticle * p = 0;

  TObjArrayIter piter(this);

  unsigned int idx = 0;

  double sum_E  = 0;
  double sum_px = 0;
  double sum_py = 0;
  double sum_pz = 0;

  while( (p = (GHepParticle *) piter.Next()) ) {

     stream << "\n |";
     stream << setfill(' ') << setw(3)  << idx++               << " | ";
     stream << setfill(' ') << setw(8)  << p->Name()           << " | ";
     stream << setfill(' ') << setw(3)  << p->Status()         << " | ";
     stream << setfill(' ') << setw(10) << p->PdgCode()        << " | ";
     stream << setfill(' ') << setw(3)  << p->FirstMother()    << " | ";
     stream << setfill(' ') << setw(3)  << p->LastMother()     << " | ";
     stream << setfill(' ') << setw(3)  << p->FirstDaughter()  << " | ";
     stream << setfill(' ') << setw(3)  << p->LastDaughter()   << " | ";
     stream << setiosflags(ios::fixed) << setprecision(3);
     stream << setfill(' ') << setw(6)  << p->Px()             << " | ";
     stream << setfill(' ') << setw(6)  << p->Py()             << " | ";
     stream << setfill(' ') << setw(6)  << p->Pz()             << " | ";
     stream << setfill(' ') << setw(6)  << p->E()              << " | ";

     if( p->IsOnMassShell() )
        stream << setfill(' ') << setw(6)  << p->Mass()        << " | ";
     else
        stream << setfill('*') << setw(6)  << p->Mass()        << " | " << p->GetP4()->M();

     // compute P4Final - P4Initial
     //
     // Take into account real particles and fake (generator-specific)
     // particles (rootino, bindino, ...) used to record non-fake physics.
     // Ignore initial & final state ions (if any).

     if( p->IsParticle() || p->IsFake() ) {

       if(p->Status() == kIStStableFinalState) {

          sum_E  += p->E();
          sum_px += p->Px();
          sum_py += p->Py();
          sum_pz += p->Pz();
       }
       else if(p->Status() == kIStInitialState || p->Status() == kIstNucleonTarget) {

          sum_E  -= p->E();
          sum_px -= p->Px();
          sum_py -= p->Py();
          sum_pz -= p->Pz();
       }
     }// !nucleus

  } // loop over particles

  stream << "\n |";
  stream << setfill('-') << setw(104) << "|";

  // Print SUMS

  stream << "\n |";
  stream << setfill(' ') << setw(6)  << "    | "
         << setfill(' ') << setw(11) << "Fin-Init:| "
         << setfill(' ') << setw(6)  << "    | "
         << setfill(' ') << setw(13) << "    | "
         << setfill(' ') << setw(12) << "        | "
         << setfill(' ') << setw(12) << "          | ";
  stream << setiosflags(ios::fixed) << setprecision(3);
  stream << setfill(' ') << setw(6)  << sum_px  << " | ";
  stream << setfill(' ') << setw(6)  << sum_py  << " | ";
  stream << setfill(' ') << setw(6)  << sum_pz  << " | ";
  stream << setfill(' ') << setw(6)  << sum_E   << " | ";
  stream << setfill(' ') << setw(9)  << "   | ";

  stream << "\n |";
  stream << setfill('-') << setw(104) << "|";
  stream << "\n";
}
//___________________________________________________________________________
