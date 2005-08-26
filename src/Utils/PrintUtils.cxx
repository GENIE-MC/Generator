//____________________________________________________________________________
/*!

\namespace  genie::print_utils

\brief      Simple printing utilities

\author     Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
            CCLRC, Rutherford Appleton Laboratory

\created    May 06, 2004

*/
//____________________________________________________________________________

#include <iostream>
#include <fstream>
#include <sstream>

#include <TSystem.h>

#include "Utils/PrintUtils.h"

using std::ostringstream;
using std::cout;
using std::endl;
using std::ios;
using std::ifstream;

//____________________________________________________________________________
string genie::print_utils::P4AsString(const TLorentzVector * p)
{
  ostringstream fmt;

  fmt << "(E = " << p->Energy()
      << ", px = " << p->Px()
      << ", py = " << p->Py()
      << ", pz = " << p->Pz() << ")"
      << " / M = " << TMath::Sqrt( TMath::Max(0., p->Mag2()) )
      << " / P = " << p->P();

  return fmt.str();
}
//____________________________________________________________________________
string genie::print_utils::P4AsShortString(const TLorentzVector * p)
{
  ostringstream fmt;

  fmt << "(E = " << p->Energy()
      << ", px = " << p->Px()
      << ", py = " << p->Py()
      << ", pz = " << p->Pz() << ")";

  return fmt.str();
}
//____________________________________________________________________________
string genie::print_utils::X4AsString(const TLorentzVector * vec4)
{
  ostringstream fmt;

  fmt << "(t = "  << vec4->T()
      << ", x = " << vec4->X()
      << ", y = " << vec4->Y()
      << ", z = " << vec4->Z() << ")";

  return fmt.str();
}
//____________________________________________________________________________
string genie::print_utils::Vec3AsString(const TVector3 * vec)
{
  ostringstream fmt;

  fmt << "( x = " << vec->X()
      << ", y = " << vec->Y()
      << ", z = " << vec->Z() << ")";

  return fmt.str();
}
//____________________________________________________________________________
string genie::print_utils::BoolAsString(bool b)
{
  return BoolAsTFString(b);
}
//____________________________________________________________________________
string genie::print_utils::BoolAsTFString(bool b)
{
  if(b) return "[true]";
  else  return "[false]";
}
//____________________________________________________________________________
string genie::print_utils::BoolAsIOString(bool b)
{
  if(b) return "[ON]";
  else  return "[OFF]";
}
//____________________________________________________________________________
string genie::print_utils::BoolAsYNString(bool b)
{
  if(b) return "[YES]";
  else  return "[NO]";
}
//____________________________________________________________________________
void genie::print_utils::PrintBanner(void)
{
// loads & prints the GENIE banner

  string base_dir = string(gSystem->Getenv("GENIE"));
  string filename = base_dir + string("/data/banner/BANNER.txt");

  ifstream banner(filename.c_str(), ios::in);

  if( banner.is_open() ) {
      banner.seekg(0, ios::end);

      int    length = banner.tellg();
      char * buffer = new char[length];

      banner.seekg(0, ios::beg);
      banner.read(buffer, length);

      cout << "\n\n" << buffer << "\n" << endl;
      delete [] buffer;

      gSystem->Sleep(1000); // watch the banner for 1 sec
  }
}
//___________________________________________________________________________
