//_____________________________________________________________________________
/*!

\class    genie::nuvld::GuiTablePrinter

\brief    Responds to GUI events for printing DBTable<T>'s

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  August 26, 2004
*/
//_____________________________________________________________________________

#include <iomanip>
#include <sstream>

#include <TGProgressBar.h>

#include "NuVldGUI/GuiTablePrinter.h"
#include "NuVldGUI/SysLogSingleton.h"
#include "NuVldGUI/BrowserSingleton.h"
#include "NuVldGUI/NuVldUserData.h"

using std::ostringstream;
using std::setfill;
using std::setw;
using std::setprecision;

using namespace genie;
using namespace genie::nuvld;

ClassImp(GuiTablePrinter)

//______________________________________________________________________________
GuiTablePrinter::GuiTablePrinter()
{
  fScaleXSecWithEnergy = false;
}
//______________________________________________________________________________
GuiTablePrinter::~GuiTablePrinter()
{

}
//______________________________________________________________________________
void GuiTablePrinter::PrintXSecTable(DBTable<vXSecTableRow> * table) const
{
  SysLogSingleton *  syslog  = SysLogSingleton::Instance();
  BrowserSingleton * browser = BrowserSingleton::Instance();

  NuVldUserData * user_data = NuVldUserData::Instance();

  // if >0 entries plot them, else display warning
  if( table->NRows() > 0) {

     syslog->StatusBar()->SetText( "Data are printed in 'Data Viewer' tab", 0 );
     syslog->StatusBar()->SetText( " ",    1 );

     syslog->ProgressBar()->SetPosition(100);
     syslog->ProgressBar()->SetPosition(0);

     // print table header & separator
     string header = this->PrintNuXSecTableHeader();

     string separator = this->PrintNuXSecTableSeparator();

     browser->TextBrowser()->AddLine( header.c_str()    );
     browser->TextBrowser()->AddLine( separator.c_str() );

     double dprogress = 100. / ( user_data->NuXSec()->NRows() );

     for(int i = 0; i< table->NRows(); i++) {

        const vXSecTableRow * row =
                          dynamic_cast<const vXSecTableRow * > (table->Row(i));

        string strrow = this->PrintXSecTableRowAsString(row);

        browser->TextBrowser()  -> AddLine     ( strrow.c_str()  );
        syslog->ProgressBar() -> SetPosition ( (i+1)*dprogress );
     }

     syslog->ProgressBar()->SetPosition(0);

  } else {

     syslog->StatusBar() -> SetText( "The xsec table has no data", 1 );
     syslog->Log()       -> AddLine( "The xsec table has no data"    );

     syslog->ProgressBar()->SetPosition(0);
  }
}
//______________________________________________________________________________
void GuiTablePrinter::PrintXSecTable(DBTable<eDiffXSecTableRow> * table) const
{
  SysLogSingleton *  syslog  = SysLogSingleton::Instance();
  BrowserSingleton * browser = BrowserSingleton::Instance();
  
  // if >0 entries plot them, else display warning
  if( table->NRows() > 0) {

     syslog->StatusBar()->SetText( "Data are printed in 'Data Viewer' tab", 0 );
     syslog->StatusBar()->SetText( " ",    1 );

     syslog->ProgressBar()->SetPosition(100);
     syslog->ProgressBar()->SetPosition(0);

     // print table header & separator

     string header    = this->PrintElDiffXSecTableHeader();
     string units     = this->PrintElDiffXSecTableHeaderUnits();
     string separator = this->PrintElDiffXSecTableSeparator();

     browser->TextBrowser()->AddLine( header.c_str()    );
     browser->TextBrowser()->AddLine( units.c_str()     );
     browser->TextBrowser()->AddLine( separator.c_str() );

     double dprogress = 100. / ( table->NRows() );

     for(int i = 0; i< table->NRows(); i++) {

        const eDiffXSecTableRow * row =
                     dynamic_cast<const eDiffXSecTableRow *>(table->Row(i));

        string strrow = this->PrintXSecTableRowAsString(row);

        browser -> TextBrowser() -> AddLine     ( strrow.c_str()  );
        syslog  -> ProgressBar() -> SetPosition ( (i+1)*dprogress );
     }

     syslog->ProgressBar()->SetPosition(0);

  } else {

     syslog->StatusBar() -> SetText( "The xsec table has no data", 1 );
     syslog->Log()       -> AddLine( "The xsec table has no data"    );

     syslog->ProgressBar()->SetPosition(0);
  }
}
//______________________________________________________________________________
string GuiTablePrinter::PrintNuXSecTableHeader(void) const
{
  ostringstream strrow;

  strrow << setfill(' ') << setw(7)   << "E ";
  strrow << setfill(' ') << setw(10)  << "Emin ";
  strrow << setfill(' ') << setw(10)  << "Emax ";

  if(fScaleXSecWithEnergy) strrow << setfill(' ') << setw(10)  << "xsec/E ";
  else                     strrow << setfill(' ') << setw(10)  << "xsec ";

  strrow << setfill(' ') << setw(14)  << "+dstat ";
  strrow << setfill(' ') << setw(14)  << "-dstat ";
  strrow << setfill(' ') << setw(10)  << "+dsyst ";
  strrow << setfill(' ') << setw(14)  << "-dsyst ";

  return strrow.str();
}
//______________________________________________________________________________
string GuiTablePrinter::PrintElDiffXSecTableHeader(void) const
{
  ostringstream strrow;

  strrow << "| ";
  strrow << setfill(' ') << "dxsec/dEdOmega ";
  strrow << setfill(' ') << setw(19)  << "Uncertainty ";
  strrow << setfill(' ') << setw(10)  << "E ";
  strrow << setfill(' ') << setw(10)  << "Ep ";
  strrow << setfill(' ') << setw(13)  << "Theta ";
  strrow << setfill(' ') << setw(10)  << "W^2 ";
  strrow << setfill(' ') << setw(12)  << "Q^2 ";
  strrow << setfill(' ') << setw(10)  << "v ";
  strrow << setfill(' ') << setw(10)  << "x ";
  strrow << setfill(' ') << setw(14)  << "Epsilon ";
  strrow << setfill(' ') << setw(10)  << "Gamma ";

  return strrow.str();
}
//______________________________________________________________________________
string GuiTablePrinter::PrintElDiffXSecTableHeaderUnits(void) const
{
  ostringstream strrow;

  strrow << "| ";
  strrow << setfill(' ') << "nb/GeV*sr      ";
  strrow << setfill(' ') << setw(19)  << "nb/GeV*sr   ";
  strrow << setfill(' ') << setw(10)  << "GeV";
  strrow << setfill(' ') << setw(10)  << "GeV";
  strrow << setfill(' ') << setw(13)  << "Deg   ";
  strrow << setfill(' ') << setw(8)   << "GeV^2 ";
  strrow << setfill(' ') << setw(10)  << "GeV^2 ";
  strrow << setfill(' ') << setw(8)   << "GeV ";
  strrow << setfill(' ') << setw(10)  << "  ";
  strrow << setfill(' ') << setw(14)  << "        ";
  strrow << setfill(' ') << setw(10)  << "      ";

  return strrow.str();
}
//______________________________________________________________________________
string GuiTablePrinter::PrintNuXSecTableSeparator(void) const
{
  ostringstream strrow;
  strrow << setfill('-') << setw(93) << " ";
  return strrow.str();
}
//______________________________________________________________________________
string GuiTablePrinter::PrintElDiffXSecTableSeparator(void) const
{
  ostringstream strrow;
  strrow << setfill('-') << setw(139) << " ";
  return strrow.str();
}
//______________________________________________________________________________
string GuiTablePrinter::PrintXSecTableRowAsString(
                                                const vXSecTableRow * row) const
{
  ostringstream strrow;

  strrow << "|";
  strrow << setfill(' ') << setw(7)  << setprecision(2) << row->E() << " ";
  strrow << "|" << setfill(' ') << setw(7)  << setprecision(1) << row->Emin() << " ";
  strrow << "|" << setfill(' ') << setw(7)  << setprecision(1) << row->Emax() << " ";
  strrow << "|" << setfill(' ') << setw(10) << setprecision(2) << row->XSec() << " ";
  strrow << "|" << " + " << setfill(' ') << setw(8)  << setprecision(2) << row->StatErrP() << " ";
  strrow << "|" << " - " << setfill(' ') << setw(8)  << setprecision(2) << row->StatErrM() << " ";
  strrow << "|" << " + " << setfill(' ') << setw(8)  << setprecision(2) << row->SystErrP() << " ";
  strrow << "|" << " - " << setfill(' ') << setw(8)  << setprecision(2) << row->SystErrM() << " ";
  strrow << "|";

  return strrow.str();
}
//______________________________________________________________________________
string GuiTablePrinter::PrintXSecTableRowAsString(
                                            const eDiffXSecTableRow * row) const
{
  ostringstream strrow;

  strrow << "|" << setfill(' ') << setw(21) << setprecision(6) << row->Sigma()   << " ";
  strrow << "|" << setfill(' ') << setw(13) << setprecision(6) << row->dSigma()  << " ";
  strrow << "|" << setfill(' ') << setw(9)  << setprecision(6) << row->E()       << " ";
  strrow << "|" << setfill(' ') << setw(9)  << setprecision(6) << row->EP()      << " ";
  strrow << "|" << setfill(' ') << setw(9)  << setprecision(2) << row->Theta()   << " ";
  strrow << "|" << setfill(' ') << setw(9)  << setprecision(6) << row->W2()      << " ";
  strrow << "|" << setfill(' ') << setw(9)  << setprecision(6) << row->Q2()      << " ";
  strrow << "|" << setfill(' ') << setw(9)  << setprecision(6) << row->Nu()      << " ";
  strrow << "|" << setfill(' ') << setw(9)  << setprecision(6) << row->x()       << " ";
  strrow << "|" << setfill(' ') << setw(9)  << setprecision(6) << row->Epsilon() << " ";
  strrow << "|" << setfill(' ') << setw(9)  << setprecision(6) << row->Gamma()   << " ";
  strrow << "|";

  return strrow.str();
}
//______________________________________________________________________________
