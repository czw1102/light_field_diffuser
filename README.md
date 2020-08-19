# Code
This MATLAB code is developed for 4D light field reconstruction from 2D raw images captured by a lensless camera with a diffuser as an encoder [1]. The optimization algorithm is the alternating direction method of multipliers (ADMM), which was implemented by Nick Antipa et al. in their project of DiffuserCam (https://waller-lab.github.io/DiffuserCam/). Interested readers can refer to their paper [2].


### Demo
Users can run 'LF_demo.m' with the provided pattern and image for an example. Note that a higher spatio-angular sampling requires significantly higher memory space and runtime.


### References
[1] Zewei Cai et al. "Lensless light-field imaging through diffuser encoding," Light Sci. Appl. 9(1), 143-151 (2020).

[2] Nick Antipa et al. "DiffuserCam: Lensless single-exposure 3D imaging," Optica 5(1), 1-9 (2018).

