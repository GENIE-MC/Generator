The important point: Set this directory to be your LHAPATH environment
variable! If you do that, you needn't worry much about the remarks below.

The PDF set included here was producted by Graeme Watt (Durham, UK) to
correct a kink in the PDFs for the default version of GRV98lo that
is available in LHAPDF 5.9.x. The patched version recreates:
  http://hepdata.cedar.ac.uk//hepdata/pdflib/grv/grv98/grv98.f 

Our thanks to Graeme and LHAPDF for the patched file. Please note that 
this file will work _only_ with LHAPDF 5.x and is not properly formatted
for LHAPDF 6.

There are three files here:

* GRV98lo.LHgrid
* GRV98lo_pdflib.LHgrid
* GRV98nlo.LHgrid

`GRV98nlo.LHgrid` is identical to the file that one would fetch from
the LHAPDF installation. `GRV98lo_pdflib.LHgrid` is the patched file
and `GRV98lo.LHgrid` is a copy of the patched file. We're keeping 
two with different names but the same content to help make it
clear what the content here _is_, but we're using the copy
with the name expected for the LHAPDF standard distrubution to make
things work more simply for the user.


G. N. Perdue, 2014/08/29
