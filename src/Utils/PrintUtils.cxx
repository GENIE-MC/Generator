//____________________________________________________________________________
/*!

\namespace  genie::print_utils

\brief      Simple printing utilities
          
\author     Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
            CCLRC, Rutherford Appleton Laboratory

\created    May 06, 2004

*/ 
//____________________________________________________________________________

#include <sstream>

#include "Utils/PrintUtils.h"

using std::ostringstream;

//____________________________________________________________________________
string genie::print_utils::P4AsString(const TLorentzVector * p)
{
  ostringstream fmt;

  fmt << " ( E = " << p->Energy() 
      << ", px = " << p->Px() 
      << ", py = " << p->Py() 
      << ", pz = " << p->Pz() << " )"
      << " / M = " << TMath::Sqrt( TMath::Max(0., p->Mag2()) )
      << " / P = " << p->P();

  return fmt.str();
}
//____________________________________________________________________________
string genie::print_utils::X4AsString(const TLorentzVector * vec4)
{
  ostringstream fmt;

  fmt << " ( t = " << vec4->T()
      << ", x = " << vec4->X()
      << ", y = " << vec4->Y()
      << ", z = " << vec4->Z() << " )";

  return fmt.str();
}
//____________________________________________________________________________
string genie::print_utils::Vec3AsString(const TVector3 * vec)
{
  ostringstream fmt;

  fmt << "( x = " << vec->X() 
      << ", y = " << vec->Y() 
      << ", z = " << vec->Z() << " )";

  return fmt.str();
}
//____________________________________________________________________________
string genie::print_utils::BoolAsString(bool tf)
{
  if(tf) return "[true]";
  else   return "[false]";
}
//____________________________________________________________________________

