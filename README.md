# MergedElectron_Analyzer
This respository contains analysis code for the analysis of electron properties. The `MergedElectron_ana.cc` is the analysis code and `ana.C` is the driver script. The analysis code can be boradly divided into two parts:
1. This part of the code is used to plot a bunch of histograms for properties of electrons, GSF tracks and some hybrid properties. This part is specific to user and they can plot and fill histogram with whatever properties is available in the Tree used for the analysis. (which i
2. The later part of the code writes some of these tree variables and some more (which are a combination of already existing variables) into another root file - `myTree_..`. This second root file is used for the training of the MergedElectron_NNClassifier. All the codes of training and testing the classifier is given as part of another github repository - `MergedElectron_NNclassifier`.

The `ana.C` driver script takes as input root files and run the `MergedElectron_ana.cc` on those root files. This .cc analysis code is built for a specific structure of root files that are made from AOD root files using the AOD analyzer code given in the github repository - `Flat_NtupleMaker_AOD`

Hence the chain of analysis is (in terms of repository to be followed) : Flat_NtupleMaker_AOD  (to be used in CMSSW environment inside lxplus) ---> MergedElectron_Analyzer --->MergedElectronNNclassifier

The final updated root file that is produced by `NNscore_branchadder.py` with the NNscore branch added to the tree in the root file can then be used for further analysis using NN scores.
