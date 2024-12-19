* DK2NU FLATTENING SCRIPT FOR GENIEv3 LONG-LIVED HNL MODULE

  This directory contains all the files necessary to convert dk2nu flux files into flat ROOT trees
  that contain all the necessary information for a dynamic HNL flux prediction from hadron decay
  in GENIEv3.

\author  John Plows <komninos-john.plows \at physics.ox.ac.uk>
         University of Oxford

\cpright Copyright (c) 2003-2025, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org

** Notice:

   ${LSTR} is taken to be the "baseline" for your flux files. For example, using the MINERvA publicly available fluxes:

   ${LSTR}="g4numiv6_dk2nu_minervame[bar]_me000z[-]200i"
   (where 000z means target is at z = 0 and \pm 200i means horn current is \pm 200 kA)
   
   The final dk2nu file has the name
   ${LSTR}_${FNUM}_${SUFFIX}.

   Again, a MINERvA example:
   ${FNUM}="0000"
   ${SUFFIX}="0006.root"

   ==> ${INSTRING}="root://fndca1.fnal.gov:1095///pnfs/fnal.gov/usr/minerva/persistent/DataPreservation/flux/g4numiv6_dk2nu_minervame_me000z200i_0000_0006.root"
   (for unauthenticated access. There exist exactly 2 files, 1 FHC and 1 RHC. Please run the flattening script once with minervame_me000z200i and once with minervamebar_me000z-200i)

   If you have a series of files, say from ${FNUM}="0000" to ${FNUM}="9999", then leave
   ${FNUM}="REPLACE_INDEX"
   You can run doFlattenDk2nus.sh with your desired range instead.

** File description:

   // -- setup
   + setup_dk2nu.sh : sources proper ROOT / dk2nu versions.

   // -- workhorse files, no need to modify these
   + write_dk2nus.{C,h} : the ROOT macro itself. Do not call directly.
   + write_dk2nus_reduced.C : same as write_dk2nus.C but fewer branches copied to save size.
                              These branches *must* exist in your flux files
                              (see src/contrib/beamhnl/exampleFluxes/ for ready-made trees you can use)
   
   // -- interface files
   + interactive_write_dk2nus.sh : Define important environment variables.
                                   See explanation above.

   + doFlattenDk2nus.sh : Runs the actual event loop and handles bookkeeping for files.
                          Make sure you've created $TGTDIR and ./testoutdir/test/ before running.

   Each flat ROOT tree is about 58 MB large for about 150k entries in each dk2nu file.
