# MS2Tox
MS2Tox is a machine learning tool for predicting ecotoxicity of unidentified chemicals in water by nontarget LC/HRMS

Before using function *FishLC50Prediction(folderwithSIRIUSfiles, LC50mode)* fingerprints must be calculated with SIRIUS+CSI:FingerID application.

Folder "Data" contains train, test and validation sets with calculated results. 

Folder "MS2Tox" is R package folder containing all the needed files for using the package *MS2Tox* in R.

Folder "Example data" contains SIRIUS+CSI:FingerID calculated fingerprints for testing if the code works.


Step by step instructions from MS data to LC50 values

## Fingerprint calculation with SIRIUS+CSI:FingerID

Before calculating toxicity values, SIRIUS+CSI:FingerID must be used for calculations of fingerprints. Fingerprint calculations with SIRIUS+CSI:FingerID can be done all at once and in this case all the calculated fingerprints will be in this one folder. This folder  is used as path to directory in the input for function *FishLC50Prediction()*. For too many peaks, the SIRIUS+CSI:FingerID calculation may crash. In this case we suggest loading the peaks for calculations in the patch of 200 at a time. Fingerprints can be calculated using SIRIUS GUI, command line or patRoon. More detailed instructions for SIRIUS file can be found from (https://boecker-lab.github.io/docs.sirius.github.io/io/). Here short desription of workflow used by the Kruvelab group while conducting the experiments is given. For this workflow SIRIUS (4.9.5) GUI was used.  

For fingerprint calculation it is important that correct input files are dragged into SIRIUS+CSI:FingerID GUI. Example of the .ms file is shown in Figure 1. In order for the toxicity prediction function to give out correct table, SIRIUS+CSI:FingerID input file should start with unique ID. Function uses this ID and later subfolder number to differentiate one compound from another. The first inthe file row must start with **>compound** and then after space can contain either name of the compound or unique id for the analysed peak/spectrum. **>parentmass** and **>ionization** are not mandatory but helpful for better identification. Row **>formula** can be added if known. Information from MS2 spectra are added after row **>collision** *collision value*. As SIRIUS+CSI:FingerID adds all the different MS2 spectra together before calculations it is not extremely important to add exactly correct collision energy values or even separate the spectra  (RAMP can be used as well). In this case instead of **>collision** you can use **>ms2** and add the MS2 information same way below. For last part after **>ms1peaks** (or **>ms1**) isotope pattern from MS1 spectrum can be added; however, this is not mandatory for calculations but can improve the formula prediction. 

![image](https://user-images.githubusercontent.com/68953270/153868916-528a8127-22a6-41f9-99c8-30880f7d18e9.png) 

**Figure 1.** Example input file for SIRIUS+CSI:FingerID


For calculations .ms files can be dragged into the GUI. Before that new project directory must be established by clicking to icon New (see Fig 2A). This path to the directory will be used in final R function for calculations. After dragging .ms files into the application all the information about each file is listed on the left side in the application (Fig 2B). By clicking on Compute All button (Fig 2C) compute window (Fig 2D) opens. In this window parameters for calculations can be chosen. In this work for Orbitrap measurements mass deviation 5 ppm was chosen. Calculations may take time from few minutes to hours depending on the number of .ms files and molecular mass of the compounds. After calculations a subfolder per .ms file has been created and it contains calculated fingerprints, scores, fragmentation trees, etc. All output folders start with order number that is shown in the final table in column named “foldernumber”. See Fig  3. 

![image](https://user-images.githubusercontent.com/68953270/153869370-9aaa1fc3-4fdb-41eb-b504-53e41f391ee3.png)

**Figure 2.** SIRIUS+CSI:FingerID application with chosen parameters for calculations

**Table 1.** Selected elements and their maximum amount for SIRIUS+CSI:FingerID calculations

| Element  | Amount |   | Element | Amount |   | Element | Amount |   | Element | Amount |
| - | - | - | - | -- | - | -- | - | - | - | -- |
| H  | Inf |   | C | Inf |   | N | Inf |   | O | Inf |
| P  | 8 |   | B | 11 |   | Si | 9 |   | S | 12 |
| Cl  | 18 |   | Se | 2 |   | Br | 10 |   | F | Inf |
| I  | 6 |   | K | 1 |   | Na | 1 |   | As | 2 |



![image](https://user-images.githubusercontent.com/68953270/153868996-770a007f-4f06-4dc5-bc9c-30fd57fc89cd.png)

**Figure 3.** Folder after fingerprint calculations that directory must be put into the function


## *MS2Tox* package for predicting LC50 values from fingerprints 

R package *MS2Tox* can be then used for toxicity value predictions.  R package can be found in GitHub and before installing *MS2Tox* package *Rdisop* is needed to be installed. To install this package, start R (version "4.1") and enter:

    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

    BiocManager::install("Rdisop")


Addition to Rdisop CRAN packages *dplyr, readr, rlist, stringr, tibble, xgboost, tidyselect, magrittr* are needed before installing and using *MS2Tox*. 

    library(dplyr)
    library(readr)
    library(rlist)
    library(stringr)
    library(tibble)
    library(xgboost)
    library(tidyselect)
    library(magrittr)
    library(Rdisop)

For installing MS2Tox.zip file can be downloaded from GitHub folder which contains addition to code also pretrained models for fish static LC50 and fish flow-through LC50. For installing .zip file can be unpacked and .tar.gz file can be installed using Install -> Package Archieve file (.zip; .tar.gz). In order for code to be able to reach trained models, **Working directory must be set** to downloaded folder from GitHub  **“/MS2Tox”**. 

*FishLC50Prediction ()* function takes in the path to the directory containing all the calculated fingerprints and information about the LC50 mode. Right now there are two possible LC50 modes for which predictions can be made: static and flow through. If LC50mode is not chosen, calculations will be done automatically in static mode. Directory that will be input for prediction function is shown in Fig 3. It returns a summary table with calculated toxicities for all of the peaks for which fingerprints were calculated by SIRIUS+CSI:FingerID. The calculation is based on the fingerprints associated with the first ranked molecular formula. SIRIUS+CSI:FingerID predicted formula, prediced ion and folder number are also displayed in the summary table. In order for the code to work properly input files for SIRIUS+CSI:FingerID must contain the compounds id/name in the first part of the file name (e.g AU231458_*restofthename*.ms). This first part of the name is used as id column in the summary table and can be used to link a specific peak with a specific toxicity value.

If data sets are very big, in SIRIUS+CSI:FingerID calculations results could be saved so that one directory contains up to 1000 subfolders, otherwise calculation with *FishLC50Prediction()* might exceed the memory limits of R.

The main function links to other functions, which gather the fingerprint data (*FingerPrintTable()*), process it, filter out only first rank chemical formulas (*FpTableForPredictions()*) based on SiriusScore (*SiriusScoreRank1()*) and calculate the toxicity value. PS! If for spectrum, two formulas are calculated with same SiriusScore, both of them are considered as first rank and in this case two different toxicity values may be seen for one peak/id. In final table, possible chemical formulas are given as well. 

In this work SIRIUS+CSI:FingerID version 4.9.5 was used so the prediction method is limited to these fingerprints. If newer or older versions do not calculate some of the fingerprints used here, the *FishLC50Prediction()* might not work . In this case please notify the author for retraining the model. 

Code example:
    
    setwd(C:/Desktop/MS2TOX)  #Directory must be set to downloaded MS2Tox folder as trained model is located there

    folderwithSIRIUSfiles <- "C:/Desktop/Folder"  #Add your directory with SIRIUS calculated fingerprints

    chosen_mode <- “static” #or "flow"

    results  <- FishLC50Prediction(folderwithSIRIUSfiles, chosen_mode)
