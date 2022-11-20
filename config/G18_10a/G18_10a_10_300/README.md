# Configuring Nuclear Partial tunes on CC0pi data
GENIE published a number of partial tunes to CC0pi data (https://arxiv.org/abs/2206.11050). Each tune was performed against a different set of data: 

      - G10a: MicroBooNE numu CCQELike data

      - G11a: MicroBooNE numubar CCQELike data

      - G20a: T2K numu CC0p0pi data           

      - G30a: MINERvA numu CC0pi data         

      - G31a: MINERvA numubar CC0p0pi data    

      - G35a: MINERvA numu CCNp0pi data  

## How to configure the partial tunes 
The configuration of these tunes is not trivial. This CMC (G18_10a_10_300), offers with a template on how to setup these tunes. The default tune corresponds to the G30a tune, obtained by tuning GENIE against MINERvA numu CC0pi data. In order to enable other partial tunes from the paper, you can simply edit the following files: 
    - CommonParam.xml
    - MECScaleVsW.xml
    - XSecLinearCombinations.xml
The parameters for the other partial tunes are commented out in these files. In order to enable them, comment the current parameter associated to the G30a tune and uncomment the desired tune's parameters.

**NOTE**: These tunes are not the final version - the G18_10a_10_300 CMC is a simple placeholder for the tunes and it will be outdated after the next GENIE release. 