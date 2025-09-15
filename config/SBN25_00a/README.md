# SBN tune

This configuration is based upon AR23_20i_00_000 but has modifications
requested in summer 2025 by the Short Baseline Neutrino program experiments.

Notable components of the physics model are the following:
 - Valencia model for 1p1h, using z-expansion, with RPA turned off
 - SuSAv2 model for 2p2h
 - Spectral function for select nuclei, including one for 40Ar based on JLab
   measurements (see https://doi.org/10.1103/PhysRevD.105.112002
   and https://doi.org/10.1103/PhysRevD.107.012005)
 - Other nuclei use the AR23_20i_00_000 "spectral-function-like approach" for LFG
 - The parameters related to pion production are taken from the G18_10a_02_11b
   tune in order to ensure a better starting point.
 - De-exctitation photons are enabled for 40Ar
