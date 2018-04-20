#!/bin/bash
#SIV=32
#echo \"${SIV}\"
#nohup matlab -nodisplay -nodesktop -r "SIv_in=${SIV}; Diffusivity_inference_from_MCMC" &> ./zkappa_$SIV.txt &

#SIV=1
#echo \"${SIV}\"
#nohup matlab -nodisplay -nodesktop -r "SIv_in=${SIV}; Diffusivity_inference_from_MCMC" &> ./zkappa_$SIV.txt &

#SIV=64
#echo \"${SIV}\"
#nohup matlab -nodisplay -nodesktop -r "SIv_in=${SIV}; Diffusivity_inference_from_MCMC" &> ./zkappa_$SIV.txt &


#nohup matlab -nodisplay -nodesktop -r "SIv_in=[1 2 4 8]; Diffusivity_inference_from_MCMC" &> ./zkappa_1248.txt &
#nohup matlab -nodisplay -nodesktop -r "SIv_in=[16 32 64]; Diffusivity_inference_from_MCMC" &> ./zkappa_163264.txt &

#nohup matlab -nodisplay -nodesktop -r "SIv_in=[1 8 4 2]; Diffusivity_inference_from_MCMC" &> ./zkappa_1842.txt &
#nohup matlab -nodisplay -nodesktop -r "SIv_in=[16 32 64]; Diffusivity_inference_from_MCMC" &> ./zkappa_163264.txt &
#nohup matlab -nodisplay -nodesktop -r "SIv_in=[16]; Diffusivity_inference_from_MCMC" &> ./zkappa_TransMat16.txt &
#nohup matlab -nodisplay -nodesktop -r "SIv_in=[32]; Diffusivity_inference_from_MCMC" &> ./zkappa_TransMat32.txt &
#nohup matlab -nodisplay -nodesktop -r "SIv_in=[64]; Diffusivity_inference_from_MCMC" &> ./zkappa_TransMat64.txt &


nohup matlab -nodisplay -nodesktop -r "kappa_Profile = 'sinusoidal'; veloc_Profile = 'childress_soward'; Script_generate_particle_trajectories" &> ./sinCS.txt &
#nohup matlab -nodisplay -nodesktop -r "kappa_Profile = 'const'; veloc_Profile = 'taylor_green_noisy'; Script_generate_particle_trajectories" &> ./constTGn.txt &

