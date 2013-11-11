Documentation for making Intranuke validation plots.
S. Dytman  10 November, 2013

- ntuple files are created with gevgen_hadron and gntpc -f ginuke
- there are many external data files
- these are separate from other validation programs.  
- these are built on top of existing programs written previously,
but this is transparent to the user.

CREATE good directory structure.
- make a directory for root files, one for data files, and one for plot files.
- suggestion for top directory=genie-out.  Put runfast.pl, intranukeplotter.pl.
genieStyle.C, and rootgINukeVal.C in that directory.
- make subdirectories Sim, ExtData, png_files for 3 kinds of file.
- move external data files into ExtData directory.

CREATE root ntuple files with hadron-nucleus data.
- run runfast.pl with perl to produce the root files you need.
Name convention is Month_day_yr_probe_target_beamenergy_version_model.ginuke.root,
 e.g. Mar_08_13_p_C_197_v280_hA.ginuke.root for proton carbon at 197 MeV using hA model 
and GENIE version 2.8.0 on March 8, 2013. This would be produced with
perl runfast.pl --type root --a hautala --m hA.  The --a option specifies the lead
author of the data publication and the program knows what to create.  Default is 100k events.
Program needs to run a setup program for each version, e.g. /usr/GENIE/setup_genie for 2.8.0.
- run intranukeplotter.pl with perl to produce the plots you need.
Program creates format files for the plots which interface to the plotting
program rootgINukeVal.C which uses style file genieStyle.C.  If setup is correct, 
this is transparent.
perl intranukeplotter.pl --type nrg --a hautala --dorf Mar_08_13 --v 280 --m hA
if no directory information needs to be specified.
--datadir specifies location of external data files (default = ./)
--rootdir specifies location of GENIE root files (default = ./)
--pngdir specifies location of png plot files (default = png_files)

input details are given when input specifications are in error.

