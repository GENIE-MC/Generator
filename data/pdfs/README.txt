The PDF set included here was producted by Graeme Watt (Durham, UK) to
correct a kink in the PDFs for the default version of GRV98lo that
is available in LHAPDF 5.9.x. Our thanks to Graeme and LHAPDF for the
patched file.

To install the PDF, run the included installer script. A successful run
might look like:

$ ./install_patched_lhapdf.sh 
LHAPATH is /mypath/GENIE/support/lhapdf
 Using this as the PDF install directory.
Installing the PDF set to /mypath/GENIE/support/lhapdf
`GRV98lo_pdflib.LHgrid' -> `/mypath/GENIE/support/lhapdf/GRV98lo.LHgrid' 
  (backup: `/mypath/GENIE/support/lhapdf/GRV98lo.LHgrid~')

Here `cp` created a backup to preserve the existing file with that
name (if, for example, you have already fetched `GRV98lo.LHgrid`
from LHAPDF.

G. N. Perdue, 2014/08/28
