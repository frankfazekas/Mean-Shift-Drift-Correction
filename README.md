# Mean Shift Drift Correction

Supplementary MATLAB software for the paper *A Mean Shift Algorithm for Drift Correction in Localization Microscopy* (https://www.biorxiv.org/content/10.1101/2021.05.07.443176v1).

Estimates and corrects for the drift of a 2D or 3D localization microscopy dataset using a mean shift algorithm. "compute_drift_2D.m" and "compute_drift_3D.m" are the main body functions which accept the inputted list of localization coordinates for the drift correction. "example_NPC.m" and "example_Bcell.m" call these functions and perform sample drift corrections on experimental datasets from the paper. These examples may prove useful to users in acquainting themselves with the software and in extending this analysis to their own data.

## Using the Software
To run the software, first install a C compiler compatible to your version of MATLAB and apply the commands

	mex Fcrosspairs.c
	mex Fcrosspairs_3D.c

to compile the source code into MEX files appropriate for your platform. Sample MEX files, compiled by gcc and MinGW64, are provided in the distribution, but these may not be compatible with other compilers. The software has been tested on Red Hat Enterprise Linux, Version 8.2 (MATLAB 2016b and 2019b) and Windows 10 Enterprise Version 20H2 (MATLAB 2019a).

## Attributions

The distribution makes use of the supplementary software from Wang et al., (https://github.com/yinawang28/RCC), employing a redundant least-squares minimization algorithm for optimizing the drift correction. Their code has been included in the "RCC" folder and is incorporated into the "compute_drift_2D.m" and "compute_drift_3D." functions in the main folder. The "compute_drift_2D_RCC.m" and "compute_drift_3D_RCC.m" functions are directly based on "RCC.m" and were used in the paper to compare the mean shift and nonlinear fitting algorithms.

To evaluate the performance of the drift correction, we used the Fourier Ring Correlation method described in Nieuwenhuizen et al., implemented as supplementary software to the paper and available in the "FIREFunctions" folder. This is dependent on the DIPImage image processing MATLAB library, licensed under the Apache License, Version 2.0.   
Copyright 2014-2021 Cris Luengo and contributors  
Copyright 1995-2014 Delft University of Technology


## References
1. Wang, Y., J. Schnitzbauer, Z. Hu, X. Li, Y. Cheng, Z.-L. Huang, and B. Huang. 2014. Localization events-based sample drift correction for localization microscopy with redundant cross-correlation algorithm. Opt. Express. 22:15982.

2. Nieuwenhuizen, R.P.J., K.A. Lidke, M. Bates, D.L. Puig, D. Grünwald, S. Stallinga, and B. Rieger. 2013. Measuring image resolution in optical nanoscopy. Nat Methods. 10:557–562.