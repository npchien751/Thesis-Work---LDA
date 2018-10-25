# Thesis-Work-LDA
This is MATLAB code developed by myself and Dr. Laura Lautz to fingerprint sources of chloride in shallow groundwater. It uses a form of the discriminant analysis classification algorithm. 

This repositary contains four scripts related to this work. All with the prefix SFP standing for Swift FingerPrinting.
  - SFP_parameters.m is a function that reads in an excel file in a standard template and outputs an array that can be easily read into the main algorithm file.
  - SFP_classification.m is the main script that requires an input array output from SFP_parameters.m
  - SFP_SensitivityAnalysis.m Code to test the sensitivity of the classification results to different combinations of input parameters
  - SFP_FinalFigures.m Final figures used for journal publication
It is a Matlab script and is published in the journal Science of the Total Environment: Chien, N.P., and L.K. Lautz. 2018. Discriminant analysis as an improved quantitative method for geochemical fingerprinting of groundwater salinity. Science of the Total Environment 618: 379-387. 
