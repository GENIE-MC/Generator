**Additional tables** are provided in separate subdirectories, one for each npnh contribution:

- `martini/NN2p2h/` → nucleon–nucleon correlations for 2p2h (denoted *NN 2p2h* in [1])
- `martini/DD2p2h/` → Δ–MEC contribution for 2p2h (denoted *ΔΔ 2p2h* in [1])
- `martini/ND2p2h/` → NN–Δ interference for 2p2h (denoted *NΔ 2p2h* in [1])
- `martini/2p2h/` → total 2p2h (sum of NN, ΔΔ, and NΔ)
- `martini/3p3h/` → Δ–MEC contribution for 3p3h (denoted *ΔΔ 3p3h* in [1])

The hadron tensor tables used by default contain the **sum of all 2p2h and 3p3h contributions** and are located in the top-level `martini/` directory.

All files inside these subdirectories share the same filename, which for each nuclear target has the form:
`Martini_XXXXXXXXXX_MEC_FullAll.dat`

The **only element that changes is the subdirectory name**, which specifies the npnh component.  
This design allows users to switch between tensor sets simply by redirecting GENIE to a different subdirectory, without renaming any files.

GENIE uses the following tables by default, depending on the selected interaction channel:

- MEC channel → `data/evgen/hadron_tensors/martini/Martini_XXXXXXXXXX_MEC_FullAll.dat`
- CCQE channel → `data/evgen/hadron_tensors/martini/Martini_XXXXXXXXXX_QE_Full.dat`

---

## How to use a different table

Since the filename is fixed, selecting a different npnh component requires changing the directory where GENIE searches for the hadron tensor table file.  
This path is controlled by the line:

```
<param type="string" name="DataPath"> data/evgen/hadron_tensors/martini </param>
```

inside `config/MartiniMECHadronTensorModel.xml`.

To select another component (e.g. only NN 2p2h or only 3p3h), modify the `DataPath` entry so that it points to:

```
data/evgen/hadron_tensors/<desired_subdirectory>
```

All available subdirectories are already listed in  
`config/MartiniMECHadronTensorModel.xml`.  
To switch to a different component, **simply comment out the default line and uncomment the desired one**.  

No renaming of the hadron tensor table files is required — only the directory path changes.  
No recompilation is needed after editing `config/MartiniMECHadronTensorModel.xml`.

---

## Reference

For a comprehensive description of these interaction channels, see:

[1] M. Martini, M. Ericson, and G. Chanfray, *Phys. Rev. C 80*, 065501 (2009), arXiv:0910.2622 [nucl-th]
