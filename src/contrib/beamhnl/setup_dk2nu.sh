#!/bin/bash

# Necessary setup script for the dk2nu flattening script
#
# \author  John Plows <komninos-john.plows \at physics.ox.ac.uk>
#          University of Oxford
#
# \cpright Copyright (c) 2003-2025, The GENIE Collaboration
#          For the full text of the license visit http://copyright.genie-mc.org

# root
source /cvmfs/larsoft.opensciencegrid.org/products/setup #moved over to larsoft external module for SL7 UPS products
setup root v6_18_04d -f Linux64bit+3.10-2.17 -q debug:e19:py2 -z /cvmfs/larsoft.opensciencegrid.org/products

# dk2nu
setup dk2nu v01_05_01b -f Linux64bit+3.10-2.17 -q debug:e15

