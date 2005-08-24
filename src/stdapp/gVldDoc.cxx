//____________________________________________________________________________
/*!

\program gVldDoc

\brief   Creates GENIE's Validity Document

         Syntax : gVldDoc [-c document-config]

         where document-config is a string specifying one of the named
         VldDocument algorithm configurations [as they are found in the
         algorithm's XML config file in $GENIE/config/vld_document.xml]

         Some of the Validity Document configurations are the:
            >> Show-All-Pages [default]
            >> Test

         Example: gVldDoc -c Show-All-Pages

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created July 14, 2004
*/
//____________________________________________________________________________

#include "Messenger/Messenger.h"
#include "ValidityDoc/VldDocument.h"

using namespace genie;

int main(int argc, char ** argv)
{
  //------ get user's choice on validity document configuration
  //       by parsing the command line arguments

  char validity_doc_config[256] = "Show-All-Pages"; // default

  while( argc>1 && (argv[1][0] == '-'))
  {
    if (argv[1][1] == 'c') {

       if ( strlen(&argv[1][2]) ) strcpy(validity_doc_config, &argv[1][2]);
       else if( (argc>2) && (argv[2][0] != '-') ) {
          argc--;
          argv++;
          strcpy(validity_doc_config, &argv[1][0]);
       } else {
          LOG("Main", pWARN) << "Bad option ignored: " << argv[1];
       }
    }
    argc--;
    argv++;
  }

  //------ create GENIE's validity document

  LOG("Main", pINFO) 
          << "Creating ValidityDocument with configuration=" 
                                               << validity_doc_config;

  VldDocument vld_doc(validity_doc_config);
  vld_doc.CreateDocument();

  return 0;
}

