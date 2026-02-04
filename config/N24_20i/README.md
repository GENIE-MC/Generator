# NOvA configuration (2024 edition)

This configuration is based on the AR23_20i "liquid argon" tune (see below), with a few small tweaks:
 * Restores the "correlated tail" to the "correlated local Fermi gas", which was lost in updates broadening the (Emiss, pmiss) distribution
 * Increases the "correlated tail" strength to 20%, which is what's used in spectral function work such as https://arxiv.org/abs/2407.18226

The physics content is the following (cribbed from AR23 notes):
 - Valencia model for 1p1h, using z-expansion form factor
 - SuSAv2 model for 2p2h
 - Spectral function-like (Emiss, pmiss) distribution for the Local Fermi gas with "correlated tail"
 - The parameters related to pion production are taken from the G18_10a_02_11b tune in order to ensure a better starting point. 
 - De-exitation photons are enabled for Argon