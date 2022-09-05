#!/bin/bash

# root
source /cvmfs/larsoft.opensciencegrid.org/products/setup #moved over to larsoft external module for SL7 UPS products
setup root v6_18_04d -f Linux64bit+3.10-2.17 -q debug:e19:py2 -z /cvmfs/larsoft.opensciencegrid.org/products

# dk2nu
setup dk2nu v01_05_01b -f Linux64bit+3.10-2.17 -q debug:e15

