//_____________________________________________________________________________
/*!

\class    genie::nuvld::GuiTablePrinter

\brief    Responds to GUI events for printing DBTable<T>'s

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

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
void GuiTablePrinter::PrintTable(DBTable<vXSecTableRow> * table) const
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
     string separator = this->PrintTableSeparator(112);

     browser->TextBrowser()->AddLine( header.c_str()    );
     browser->TextBrowser()->AddLine( separator.c_str() );

     double dprogress = 100. / ( user_data->NuXSec()->NRows() );

     for(int i = 0; i< table->NRows(); i++) {

        const vXSecTableRow * row =
                          dynamic_cast<const vXSecTableRow * > (table->Row(i));

        string strrow = this->PrintTableRowAsString(row);

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
void GuiTablePrinter::PrintTable(DBTable<eDiffXSecTableRow> * table) const
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
     string separator = this->PrintTableSeparator(158);

     browser->TextBrowser()->AddLine( header.c_str()    );
     browser->TextBrowser()->AddLine( units.c_str()     );
     browser->TextBrowser()->AddLine( separator.c_str() );

     double dprogress = 100. / ( table->NRows() );

     for(int i = 0; i< table->NRows(); i++) {

        const eDiffXSecTableRow * row =
                     dynamic_cast<const eDiffXSecTableRow *>(table->Row(i));

        string strrow = this->PrintTableRowAsString(row);

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
void GuiTablePrinter::PrintTable(DBTable<SFTableRow> * table) const
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

     string header    = this->PrintSFTableHeader();
     string separator = this->PrintTableSeparator(130);

     browser->TextBrowser()->AddLine( header.c_str()    );
     //browser->TextBrowser()->AddLine( units.c_str()     );
     browser->TextBrowser()->AddLine( separator.c_str() );

     double dprogress = 100. / ( table->NRows() );

     for(int i = 0; i< table->NRows(); i++) {

        const SFTableRow * row = dynamic_cast<const SFTableRow *>(table->Row(i));        
        string strrow = this->PrintTableRowAsString(row);
        browser -> TextBrowser() -> AddLine     ( strrow.c_str()  );
        syslog  -> ProgressBar() -> SetPosition ( (i+1)*dprogress );
     }
     syslog->ProgressBar()->SetPosition(0);

  } else {
     syslog->StatusBar() -> SetText( "The table contains no data", 1 );
     syslog->Log()       -> AddLine( "The table contains no data"    );

     syslog->ProgressBar()->SetPosition(0);
  }
}
//______________________________________________________________________________
string GuiTablePrinter::PrintNuXSecTableHeader(void) const
{
  ostringstream strrow;

  strrow << "| ";
  strrow << setfill(' ') << "Experiment ";
  strrow << setfill(' ') << setw(6)   << "Tag ";

  strrow << setfill(' ') << setw(7)   << "E ";
  strrow << setfill(' ') << setw(10)  << "Emin ";
  strrow << setfill(' ') << setw(10)  << "Emax ";

  if(fScaleXSecWithEnergy) strrow << setfill(' ') << setw(10)  << "xsec/E ";
  else                     strrow << setfill(' ') << setw(10)  << "xsec ";

  strrow << setfill(' ') << setw(13)  << "+dstat ";
  strrow << setfill(' ') << setw(13)  << "-dstat ";
  strrow << setfill(' ') << setw(13)  << "+dsyst ";
  strrow << setfill(' ') << setw(13)  << "-dsyst ";

  return strrow.str();
}
//______________________________________________________________________________
string GuiTablePrinter::PrintElDiffXSecTableHeader(void) const
{
  ostringstream strrow;

  strrow << "| ";
  strrow << setfill(' ') << "Experiment ";
  strrow << setfill(' ') << setw(6)   << "Tag ";
  strrow << setfill(' ') << setw(19)  << "dxsec/dEdOmega ";
  strrow << setfill(' ') << setw(18)  << "Uncertainty ";
  strrow << setfill(' ') << setw(9)   << "E ";
  strrow << setfill(' ') << setw(11)  << "Ep ";
  strrow << setfill(' ') << setw(13)  << "Theta ";
  strrow << setfill(' ') << setw(10)  << "W^2 ";
  strrow << setfill(' ') << setw(10)  << "Q^2 ";
  strrow << setfill(' ') << setw(11)  << "v ";
  strrow << setfill(' ') << setw(10)  << "x ";
  strrow << setfill(' ') << setw(13)  << "Epsilon ";
  strrow << setfill(' ') << setw(10)  << "Gamma ";

  return strrow.str();
}
//______________________________________________________________________________
string GuiTablePrinter::PrintElDiffXSecTableHeaderUnits(void) const
{
  ostringstream strrow;

  strrow << "| ";
  strrow << setfill(' ') << setw(34)  << "nb/GeV*sr ";
  strrow << setfill(' ') << setw(19)  << "nb/GeV*sr ";
  strrow << setfill(' ') << setw(10)  << "GeV";
  strrow << setfill(' ') << setw(10)  << "GeV";
  strrow << setfill(' ') << setw(15)  << "Deg   ";
  strrow << setfill(' ') << setw(10)  << "GeV^2 ";
  strrow << setfill(' ') << setw(11)  << "GeV^2 ";
  strrow << setfill(' ') << setw(10)  << "GeV ";
  strrow << setfill(' ') << setw(10)  << "  ";
  strrow << setfill(' ') << setw(14)  << "        ";
  strrow << setfill(' ') << setw(10)  << "      ";

  return strrow.str();
}
//______________________________________________________________________________
string GuiTablePrinter::PrintSFTableHeader(void) const
{
  ostringstream strrow;

  strrow << "| ";
  strrow << setfill(' ') << "Experiment ";
  strrow << setfill(' ') << setw(6)   << "Tag ";
  strrow << setfill(' ') << setw(9)   << "p ";
  strrow << setfill(' ') << setw(11)  << "R ";
  strrow << setfill(' ') << setw(11)  << "x ";
  strrow << setfill(' ') << setw(12)  << "Q^2 ";
  strrow << setfill(' ') << setw(12)  << "S/F ";
  strrow << setfill(' ') << setw(13)  << "+dstat ";
  strrow << setfill(' ') << setw(13)  << "-dstat ";
  strrow << setfill(' ') << setw(13)  << "+dsyst ";
  strrow << setfill(' ') << setw(13)  << "-dsyst ";

  return strrow.str();
}
//______________________________________________________________________________
string GuiTablePrinter::PrintTableSeparator(int n) const
{
  ostringstream strrow;
  strrow << setfill('-') << setw(n) << " ";
  return strrow.str();
}
//______________________________________________________________________________
string GuiTablePrinter::PrintTableRowAsString(const vXSecTableRow * r) const
{
  ostringstream strrow;

  strrow << "|" << setfill(' ') << setw(11) << r->Experiment()      << " ";
  strrow << "|" << setfill(' ') << setw(4)  << r->MeasurementTag()  << " ";

  strrow << "|" << setfill(' ') << setw(7)
                                << setprecision(2) << r->E() << " ";
  strrow << "|" << setfill(' ') << setw(7)
                                << setprecision(1) << r->Emin() << " ";
  strrow << "|" << setfill(' ') << setw(7)
                                << setprecision(1) << r->Emax() << " ";
  strrow << "|" << setfill(' ') << setw(10)
                                << setprecision(2) << r->XSec() << " ";

  strrow << "|" << " + " << setfill(' ') << setw(8)
                         << setprecision(2) << r->StatErrP() << " ";
  strrow << "|" << " - " << setfill(' ') << setw(8)
                         << setprecision(2) << r->StatErrM() << " ";
  strrow << "|" << " + " << setfill(' ') << setw(8)
                         << setprecision(2) << r->SystErrP() << " ";
  strrow << "|" << " - " << setfill(' ') << setw(8)
                         << setprecision(2) << r->SystErrM() << " ";
  strrow << "|";

  return strrow.str();
}
//______________________________________________________________________________
string GuiTablePrinter::PrintTableRowAsString(const eDiffXSecTableRow * r) const
{
  ostringstream strrow;

  strrow << "|" << setfill(' ') << setw(11) << r->Experiment()     << " ";
  strrow << "|" << setfill(' ') << setw(4)  << r->MeasurementTag() << " ";
  strrow << "|" << setfill(' ') << setw(21)
                                << setprecision(6) << r->Sigma()   << " ";
  strrow << "|" << setfill(' ') << setw(13)
                                << setprecision(6) << r->dSigma()  << " ";
  strrow << "|" << setfill(' ') << setw(9)
                                << setprecision(6) << r->E()       << " ";
  strrow << "|" << setfill(' ') << setw(9)
                                << setprecision(6) << r->EP()      << " ";
  strrow << "|" << setfill(' ') << setw(9)
                                << setprecision(2) << r->Theta()   << " ";
  strrow << "|" << setfill(' ') << setw(9)
                                << setprecision(6) << r->W2()      << " ";
  strrow << "|" << setfill(' ') << setw(9)
                                << setprecision(6) << r->Q2()      << " ";
  strrow << "|" << setfill(' ') << setw(9)
                                << setprecision(6) << r->Nu()      << " ";
  strrow << "|" << setfill(' ') << setw(9)
                                << setprecision(6) << r->x()       << " ";
  strrow << "|" << setfill(' ') << setw(9)
                                << setprecision(6) << r->Epsilon() << " ";
  strrow << "|" << setfill(' ') << setw(9)
                                << setprecision(6) << r->Gamma()   << " ";
  strrow << "|";

  return strrow.str();
}
//______________________________________________________________________________
string GuiTablePrinter::PrintTableRowAsString(const SFTableRow * r) const
{
  ostringstream strrow;

  strrow << "|" << setfill(' ') << setw(11) << r->Experiment()      << " ";
  strrow << "|" << setfill(' ') << setw(4)  << r->MeasurementTag()  << " ";
  strrow << "|" << setfill(' ') << setw(10) << r->p()               << " ";
  strrow << "|" << setfill(' ') << setw(10) << r->R()               << " ";
  strrow << "|" << setfill(' ') << setw(9)
                                << setprecision(6) << r->x()        << " ";
  strrow << "|" << setfill(' ') << setw(9)
                                << setprecision(6) << r->Q2()       << " ";
  strrow << "|" << setfill(' ') << setw(9)
                                << setprecision(6) << r->SF()       << " ";
  strrow << "|" << setfill(' ') << setw(11)
                                << setprecision(6) << r->StatErrP() << " ";
  strrow << "|" << setfill(' ') << setw(11)
                                << setprecision(6) << r->StatErrM() << " ";
  strrow << "|" << setfill(' ') << setw(11)
                                << setprecision(6) << r->SystErrP() << " ";
  strrow << "|" << setfill(' ') << setw(11)
                                << setprecision(6) << r->SystErrM() << " ";
  strrow << "|";

  return strrow.str();
}
//______________________________________________________________________________
