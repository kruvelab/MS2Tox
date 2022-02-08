\name{MS2Tox}
\alias{hello}
\title{MS2Tox Machine Learning Tool for Predicting Ecotoxicity of Unidentified Chemicals in Water by Nontarget LC/HRMS}
\usage{FishLC50Prediction(folderwithSIRIUSfiles, LC50mode)
}
\arguments{
\item{folderwithSIRIUSfiles}{path to the folder that has subfolders with SIRIUS fingerprints}
\item{LC50mode}{LC50mode can be either "static" or "flow", if not specified default is "static"}
}
\description{
This predicts LC50 values for fathead minnow from fingerprints that are obtained from MS2 fragmentation spectra. For fingerprints SIRIUS+CSI:fingerID is used and code input is a folder containing all the results from SIRIUS+CSI:FingerID calculations.
}

