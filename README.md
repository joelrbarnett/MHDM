# MHDM
This repository contains files to run MHDM decomposition on images affected by multaplicative gamma noise. It assumes Test Images is stored in the same directory. 
The scripts to run the tests are:
* testImagesAdditive.m
* testImagesAdditiveTight.m
* testImagesAdditiveRefined.m

They rely on the following functions
- Osher.m
- OsherTight.m
- OsherRefined.m
- metrics.m
- plotFigsOsher.m

Images are stored in the folders under the "Additive" main folder. The subfolders describe the contents. All remaining files and folders are for the AA model with TV(log(u)) penalty term, which is incomplete at the moment. 
