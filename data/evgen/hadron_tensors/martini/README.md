The default implemented tables contain the **sum of all 2p2h and 3p3h contributions**.

**Additional tables**, corresponding to specific components, are provided in the `martini/npnh/` subdirectory:

- `Martini_XXXXXXXXXX_MEC_NN2p2h.dat` → nucleon–nucleon correlations for 2p2h (denoted *NN 2p2h* in [1])  
- `Martini_XXXXXXXXXX_MEC_DD2p2h.dat` → Δ–MEC contribution for 2p2h (denoted *ΔΔ 2p2h* in [1])  
- `Martini_XXXXXXXXXX_MEC_ND2p2h.dat` → NN correlations–Δ–MEC interference 2p2h (denoted *NΔ 2p2h* in [1])  
- `Martini_XXXXXXXXXX_MEC_2p2h.dat` → sum of the three 2p2h contributions above  
- `Martini_XXXXXXXXXX_MEC_3p3h.dat` → Δ–MEC contribution for 3p3h (denoted *ΔΔ 3p3h* in [1])  
- `Martini_XXXXXXXXXX_MEC_npnh.dat` → sum of all 2p2h and 3p3h contributions  

By default, GENIE uses the last option (`_npnh.dat`), accessed via the filename  
`Martini_XXXXXXXXXX_MEC_FullAll.dat`.

---

## How to use a different table

Since GENIE always looks for the table with the filename:

data/evgen/hadron_tensors/martini/Martini_XXXXXXXXXX_MEC_FullAll.dat (see also `~/config/MartiniMECHadronTensorModel.xml`).

To use a different option, simply replace this file with the desired one.  
For example, to run GENIE with only the 2p2h contributions:

1. Rename `Martini_XXXXXXXXXX_MEC_2p2h.dat` → `Martini_XXXXXXXXXX_MEC_FullAll.dat`  
2. Copy it to the directory `data/evgen/hadron_tensors/martini/`, replacing the existing file.

---

## Reference

For a comprehensive description of these interaction channels, see:

[1] M. Martini, M. Ericson, and G. Chanfray, *Phys. Rev. C 80*, 065501 (2009), [arXiv:0910.2622](https://arxiv.org/abs/0910.2622) [nucl-th]
