# evo_migration
Model simulating micro-evolution of inherited migratory orientation programs using an evolutionary algorithm incorporating long-term geomagnetic data
Written in MATLAB

Simulates migration of naive migrating birds according to inherited magnetic (and other, e.g., star) compass use, 
with geomagnetic signposts to switch between inherited headings en route, using mean field (IGRF) geomagnetic data between 1900-2024.
The evolutionary strategy algorithm implemented simaultes inheritance of headings and signposts as averages, 
including intrinsic (stochastic) variability.

(see McLaren et al, https://www.biorxiv.org/content/10.1101/2022.06.29.498190v2, under review at Movement Ecology)

To run model, type "run_evo_migr" in the MATLAB command window, and follow the prompts to specify compass and signpost use, 
and whether to run a full 124-year or shorter (default 10-year) model simulation.

Other choices can be edited in the script itself (run_evo_migr.m). 

The folder "plot_scripts" contains a few additional scripts to visualize model output.

In order to run the model and plotting routines, the following MATAB toolboxes are required:

	Mapping toolbox
	Statistics and Machine Learning Toolbox
	Parallel Computing Toolbox

	Chad Greene's Climate Data Toolbox, which requires downloading
	https://www.chadagreene.com/CDT/CDT_Contents.html

Additionally, the model uses already-included data and scripts from 

	(i) the IGRF https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html, 
	as implemented by https://www.mathworks.com/matlabcentral/fileexchange/34388-international-geomagnetic-reference-field-igrf-model,
	found in the geo_data folder
	(ii) Consensus landcover https://www.earthenv.org/landcover, downscaled to 1x1 degrees (also found in the folder geo_data)
	(iii) a package to sample von Mises random variables vmrand(fMu, fKappa, varargin)
	https://de.mathworks.com/matlabcentral/fileexchange/37241-vmrand-fmu-fkappa-varargin
	(iv) Phillipp Berens' Circular statistics package 
	https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics
	(v) a Brewer colormap toolbox https://de.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer-attractive-and-distinctive-colormaps

For more info, email james.mclaren@uni-oldenburg.de
