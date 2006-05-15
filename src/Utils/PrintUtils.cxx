//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - May 06, 2004

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

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
string genie::utils::print::P4AsString(const TLorentzVector * p)
{
  ostringstream fmt;

  fmt << "(E = " << p->Energy()
      << ", px = " << p->Px()
      << ", py = " << p->Py()
      << ", pz = " << p->Pz() << ")";

  double m2 = p->Mag2();
  if(m2>0.) fmt << " / M = " << TMath::Sqrt(m2);
  else      fmt << " / M^2 = " << m2;

  fmt << " / P = " << p->P();

  return fmt.str();
}
//____________________________________________________________________________
string genie::utils::print::P4AsShortString(const TLorentzVector * p)
{
  ostringstream fmt;

  fmt << "(E = " << p->Energy()
      << ", px = " << p->Px()
      << ", py = " << p->Py()
      << ", pz = " << p->Pz() << ")";

  return fmt.str();
}
//____________________________________________________________________________
string genie::utils::print::X4AsString(const TLorentzVector * vec4)
{
  ostringstream fmt;

  fmt << "(t = "  << vec4->T()
      << ", x = " << vec4->X()
      << ", y = " << vec4->Y()
      << ", z = " << vec4->Z() << ")";

  return fmt.str();
}
//____________________________________________________________________________
string genie::utils::print::P3AsString(const TVector3 * vec)
{
  ostringstream fmt;

  fmt << "(px = "  << vec->X()
      << ", py = " << vec->Y()
      << ", pz = " << vec->Z() << ")";

  return fmt.str();
}
//____________________________________________________________________________
string genie::utils::print::Vec3AsString(const TVector3 * vec)
{
  ostringstream fmt;

  fmt << "(x = "  << vec->X()
      << ", y = " << vec->Y()
      << ", z = " << vec->Z() << ")";

  return fmt.str();
}
//____________________________________________________________________________
string genie::utils::print::BoolAsString(bool b)
{
  return BoolAsTFString(b);
}
//____________________________________________________________________________
string genie::utils::print::BoolAsTFString(bool b)
{
  if(b) return "[true]";
  else  return "[false]";
}
//____________________________________________________________________________
string genie::utils::print::BoolAsIOString(bool b)
{
  if(b) return "[ON]";
  else  return "[OFF]";
}
//____________________________________________________________________________
string genie::utils::print::BoolAsYNString(bool b)
{
  if(b) return "[YES]";
  else  return "[NO]";
}
//____________________________________________________________________________
void genie::utils::print::PrintBanner(void)
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
string genie::utils::print::PrintFramedMesg(string mesg)
{
  string frame(4+mesg.size(), '*');

  string framed_mesg = string("\n") + frame + string("\n* ") +
                               mesg + string(" *\n") + frame + string("\n");
  return framed_mesg;
}
//___________________________________________________________________________
