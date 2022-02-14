# MS2Tox
MS2Tox is a machine learning tool for predicting ecotoxicity of unidentified chemicals in water by nontarget LC/HRMS

Before using function FishLC50Prediction(folderwithSIRIUSfiles, LC50mode) fingerprints must be calculated with SIRIUS+CSI:FingerID application.
Folder "Data" contains train, test and validation sets with calculated results. 
Folder "MS2Tox" is R package folder containing all the needed files for using the package MS2Tox in R


Step by step instructions from MS data to LC50 values

## Fingerprint calculation with SIRIUS+CSI:FingerID

Before calculating toxicity values, SIRIUS+CSI:FingerID must be used for calculations of fingerprints. Fingerprint calculations with SIRIUS+CSI:FingerID can be done all at once and in this case all the calculated fingerprints will be in this one folder. This folder  is used as path to directory in the input for function FishLC50Prediction(). For too many peaks, the SIRIUS+CSI:FingerID calculation may crash. In this case we suggest loading the peaks for calculations in the patch of 200 at a time.

For fingerprint calculation it is important that correct input files are dragged into SIRIUS+CSI:FingerID GUI (version 4.9.5). Example of the .ms file is shown in Figure 1. The first row must start with **>compound** and then after space can contain either name of the compound or unique id for the analysed peak/spectrum. **>parentmass** and **>ionization** are not mandatory but helpful for better identification. Row **>formula** can be added if known. Information from MS2 spectra are added after row **>collision** *collision value*. As SIRIUS+CSI:FingerID adds all the different MS2 spectra together before calculations it is not extremely important to add exactly correct collision energy values or even separate the spectra  (RAMP can be used as well) . For last part after **>ms1peaks** isotope pattern from MS1 spectrum can be added; however, this is not mandotory for calculations can improve the formula prediction. 

@pildid

**Figure 1.** On the left is example input file for SIRIUS+CSI:FingerID, on the right folder after fingerprint calculations that directory must be put into the function

For calculations .ms files can be dragged into the GUI. Before that new project directory must be established by clicking to icon New (see Fig 2A). This path to the directory will be used in final R function for calculations. After dragging .ms files into the application all the information about each file is listed on the left side in the application (Fig 2B). By clicking on Compute All button (Fig 2C) compute window (Fig 2D) opens. In this window parameters for calculations can be chosen. In this work for Orbitrap measurements mass deviation 5 ppm was chosen. Calculations may take time from few minutes to hours depending on the number of .ms files and molecular mass of the compounds. After calculations a subfolder per .ms file has been created and it contains calculated fingerprints, scores, fragmentation trees, etc. All output folders start with order number that is shown in the final table in column named “foldernumber”.

@pilt siriuse programmist

**Figure 2.** SIRIUS+CSI:FingerID application with chosen parameters for calculations

**Table 1.** Selected elements and their maximum amount for SIRIUS+CSI:FingerID calculations

| Element  | Amount |   | Element | Amount |   | Element | Amount |   | Element | Amount |
| - | - | - | - | -- | - | -- | - | - | - | -- |
| H  | Inf |   | C | Inf |   | N | Inf |   | O | Inf |
| - | - | - | - | -- | - | -- | - | - | - | -- |
| P  | 8 |   | B | 11 |   | Si | 9 |   | S | 12 |


## MS2Tox package for predicting LC50 values from fingerprints 

R package MS2Tox can be then used for toxicity value predictions.  R package can be found in GitHub and before installing MS2Tox package Rdisop is needed to be installed. To install this package, start R (version "4.1") and enter:

    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

    BiocManager::install("Rdisop")


Addition to Rdisop CRAN packages dplyr, readr, rlist, stringr, tibble, xgboost, tidyselect, magrittr are needed before installing and using MS2Tox. 

For installing MS2Tox .zip file can be downloaded from GitHub folder  which contains addition to code also pretrained models for fish static LC50 and fish flow-through LC50. For installing .zip file can be unpacked and .tar.gz file can be installed using Install -> Package Archieve file (.zip; .tar.gz). In order for code to be able to reach trained models, Working directory must be set to downloaded folder from GitHub  “/MS2Tox”. 

FishLC50Prediction () function takes in the path to the directory containing all the calculated fingerprints and information about the LC50 mode. Right now there are two possible LC50 modes for which predictions can be made: static and flow through. If LC50mode is not chosen, calculations will be done automatically in static mode. Directory that will be input for prediction function is shown in Fig 1. It returns a summary table with calculated toxicities for all of the peaks for which fingerprints were calculated by SIRIUS+CSI:FingerID. The calculation is based on the fingerprints associated with the first ranked molecular formula. SIRIUS+CSI:FingerID predicted formula, prediced ion and folder number are also displayed in the summary table. In order for the code to work properly input files for SIRIUS+CSI:FingerID must contain the compounds id/name in the first part of the file name (e.g AU231458  _restofthename.ms). This first part of the name is used as id column in the summary table and can be used to link a specific peak with a specific toxicity value.

If data sets are very big, in SIRIUS+CSI:FingerID calculations results could be saved so that one directory contains up to 1000 subfolders, otherwise calculation with FishLC50Prediction() might exceed the memory limits of R.

The main function links to other functions, which gather the fingerprint data (FingerPrintTable()), process it, filter out only first rank chemical formulas (FpTableForPredictions()) based on SiriusScore (SiriusScoreRank1()) and calculate the toxicity value. PS! If for spectrum, two formulas are calculated with same SiriusScore, both of them are considered as first rank and in this case two different toxicity values may be seen for one peak/id. In final table, possible chemical formulas are given as well. 

In this work SIRIUS+CSI:FingerID version 4.9.5 was used so the prediction method is limited to these fingerprints. If newer or older versions do not calculate some of the fingerprints used here, the FishLC50Prediction()might not work . In this case please notify the author for retraining the model. 

Code example:

    folderwithSIRIUSfiles <- "C:/Desktop/Folder"

    chosen_mode <- “static”

    results  <- FishLC50Prediction(folderwithSIRIUSfiles, chosen_mode)




