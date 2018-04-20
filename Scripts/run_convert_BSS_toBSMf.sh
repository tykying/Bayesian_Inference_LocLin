#!/bin/bash
#SIV=32
#echo \"${SIV}\"
#nohup matlab -nodisplay -nodesktop -r "SIv_in=${SIV}; Diffusivity_inference_from_MCMC" &> ./zkappa_$SIV.txt &

nohup matlab -nodisplay -nodesktop -r "SIv_in=[1]; Script_Convert_BSS_to_BayesSampleMeanFine" &> ./zkappa1.txt &
nohup matlab -nodisplay -nodesktop -r "SIv_in=[2]; Script_Convert_BSS_to_BayesSampleMeanFine" &> ./zkappa2.txt &
nohup matlab -nodisplay -nodesktop -r "SIv_in=[4]; Script_Convert_BSS_to_BayesSampleMeanFine" &> ./zkappa4.txt &
nohup matlab -nodisplay -nodesktop -r "SIv_in=[8]; Script_Convert_BSS_to_BayesSampleMeanFine" &> ./zkappa8.txt &
nohup matlab -nodisplay -nodesktop -r "SIv_in=[16]; Script_Convert_BSS_to_BayesSampleMeanFine" &> ./zkappa16.txt &
nohup matlab -nodisplay -nodesktop -r "SIv_in=[32]; Script_Convert_BSS_to_BayesSampleMeanFine" &> ./zkappa32.txt &
nohup matlab -nodisplay -nodesktop -r "SIv_in=[64]; Script_Convert_BSS_to_BayesSampleMeanFine" &> ./zkappa64.txt &
