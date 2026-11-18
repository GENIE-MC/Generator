# Configuring Nuclear Partial tunes on CC0pi data
GENIE published a number of partial tunes to CC0pi data (https://journals.aps.org/prd/abstract/10.1103/PhysRevD.106.112001). Each tune was performed against a different set of data: 

      - G10a: MiniBooNE numu CCQELike data

      - G11a: MiniBooNE numubar CCQELike data

      - G20a: T2K numu CC0p0pi data           

      - G30a: MINERvA numu CC0pi data         

      - G31a: MINERvA numubar CC0p0pi data    

      - G35a: MINERvA numu CCNp0pi data  
The tunes are preformed on top of the G18_02a_02_11b CMC tune to free nucleon data. 

## How to configure the partial tunes 
The configuration of these tunes is not trivial. This CMC (GPRD18_10a) offers with a template on how to setup these tunes. The default tune corresponds to the G30a tune, obtained by tuning GENIE against MINERvA numu CC0pi data. In order to enable other partial tunes from the paper, you can simply edit the following files: 
    - CommonParam.xml
    - MECScaleVsW.xml
    - XSecLinearCombinations.xml
The parameters for the other partial tunes are commented out in these files. In order to enable them, comment the current parameter associated to the G30a tune and uncomment the desired tune's parameters.
To run the predictions from the PhysRevD.106.112001, use the GPRD18_10a_02_11b tune option. 

**NOTE**: These tunes are not part of an official GENIE release - the GPRD18_10a CMC is a simple placeholder for the tunes and it will be outdated after the next GENIE release. 
