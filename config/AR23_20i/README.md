# Liquid Argon Experiment tune

This configuration is the result of requests from SBN and DUNE collaborations.
It is designed to serve as baseline model for SBN analyses. 

The physics content is the following
 - Valencia model for 1p1h, using z-expansion
 - SuSAv2 model for 2p2h
 - Spectral function like approach for the Local Fermi Gas
    - Includes events with SRC-like missing momentum within the bulk of the $p_{miss}$ distribution, 
      but not extending pass the LFG cutoff in $p_{miss}$.
 - The parameters related to pion production are taken from the G18_10a_02_11b tune in order to ensure a better starting point. 
 - De-exitation photons are enabled for Argon

## Modifications

  Three new tunes have been added, taking different QEL z-expansion axial FF parameters following
  https://arxiv.org/abs/2512.14097 and T. Cai et al (MINERvA), Nature 614 (2023) 48-53.
  Note the arxiv paper uses the opposite sign convention to GENIE.

  AR23_20i_01_000: uses `ZExp_minerva`       parameters, with    the full sum rules, Eq. (31) CVs
  AR23_20i_01_001: uses `ZExp_minervaNature` parameters, from Nature paper Supp. Table 4
  AR23_20i_02_000: uses `ZExp_lqcd`          parameters, with    the full sum rules, Eq. (39) CVs

  Note the `Q4limit` parameter, which applies the sum rules: 2512.14097 paper's kmax=6 is
  equivalent to GENIE kmax=2 with Q4limit set to true, which adds four extra constraints.
  This also means that Nature paper kmax=8 maps to kmax=4 in GENIE with Q4limit true

  Note also the Tcut values: The arxiv paper uses Tcut = 9 (m_pi0)^2 not 9 (m_pi+)^2

  A summary table of the parameters and tunes is below:

  |-----------------------|------|---------|--------------|----------------|--------------------------|
  | Tune: AR23_20i_X0_000 | Kmax | FA0     | T0 (GeV/c)^2 | Tcut (GeV/c)^2 |       A1,A2...           |
  |-----------------------|------|---------|--------------|----------------|--------------------------|
  |  Def: AR23_20i_00_000 |    4 | -1.2670 |        -0.28 |       0.1764   |      2.30,-0.6,-3.8,2.4  |
  |-----------------------|------|---------|--------------|----------------|--------------------------|
  |  Mnv: AR23_20i_01_000 |    2 | -1.2754 |        -0.5  |       0.161604 |              1.65,-0.94  |
  |-----------------------|------|---------|--------------|----------------|--------------------------|
  |  Mnv: AR23_20i_01_001 |    4 | -1.2723 |        -0.75 |       0.1764   |      1.50,-1.2,-0.1,0.2  |
  |-----------------------|------|---------|--------------|----------------|--------------------------|
  | LQCD: AR23_20i_02_000 |    2 | -1.2754 |        -0.5  |       0.161604 |             1.721,-0.31  |
  |-----------------------|------|---------|--------------|----------------|--------------------------|