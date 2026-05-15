# SBN tune with CRPA 1p1h

This configuration is based upon AR23_20i_00_000 but has modifications
requested in summer 2025 by the Short Baseline Neutrino program experiments.

Notable components of the physics model are the following:
 - CRPA model for 1p1h
 - For free nucleons, the LwlynSmith model with Dipole is used
 - SuSAv2 model for 2p2h
 - Spectral function like approach for the Local Fermi Gas
    - Includes events with SRC-like missing momentum within the bulk of the $p_{miss}$ distribution, 
      but not extending pass the LFG cutoff in $p_{miss}$.
 - The parameters related to pion production are taken from the G18_10a_02_11b
   tune in order to ensure a better starting point.
 - De-exctitation photons are enabled for 40Ar