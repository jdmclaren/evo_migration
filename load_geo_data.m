function [coastlat,coastlon,igrfcoefs, hi_elevs, barrens, ...
    poor_stops, ok_veg_breed, ndvi, veg_stp, thr_veg_stp] = ...
    load_geo_data(geo_filnm, thrsh_snow_barrn, ...
    thrsh_elev,thrsh_ndvi_breed, thrsh_ndvi_stop_min, thrsh_trees_breed, ...
    thrsh_trees_stop,max_trees_stop,thrsh_conif_breed,thrsh_low_veg_breed, ...
    thrsh_low_veg_stop,thrsh_all_veg_stop,min_any_veg_stp)

% load geomagnetic, elevation and landcover data
    igrfcoefs = load([geo_filnm '/igrf_data/igrfcoefs']);

% folder name for 1 degree data
one_deg_fnm = [geo_filnm '/mdl_1_deg_data'];


% load topography: coast and elevation data
load coastlines
coastlat = coastlat*pi/180;
coastlon = coastlon*pi/180;

% 1x1 elevation for high elevation
load([one_deg_fnm '/elev'])
hi_elevs = elev_1_deg > thrsh_elev;

% load landcover (snow, barren land and NDVI) 
load([one_deg_fnm '/snow_1_deg'])

load([one_deg_fnm '/barrn_1_deg'])
% load mean fall NDVI
load([one_deg_fnm '/ndvi_mn_fall'])
% load([one_deg_fnm '/ndvi_Sep'])
ndvi = double(ndvi_mn_fall); % double(ndvi_Sep); %  
load([one_deg_fnm '/all_trees_1_deg'])
load([one_deg_fnm '/conifers_1_deg'])
load([one_deg_fnm '/low_veg_1_deg'])
% load([one_deg_fnm '/all_veg_1_deg'])
load([one_deg_fnm '/all_veg_rWtr'])
load([one_deg_fnm '/non_grass_1_deg'])

load([one_deg_fnm '/all_veg_paleo']); % '/all_veg_paleo_rel']);
load([one_deg_fnm '/GPP_paleo']);

% Apr 2022 changed poor stops to both barren and low ndvi (rather than
% either)
% poor_stops = barrens & (ndvi_mn_fall(1:size(barrens,1),:) < thrsh_ndvi_stop);


% barrens = ndvi_mn_fall(1:size(snow,1),:) <  thrsh_ndvi_stop & ...
%     (snow > thrsh_snow_barrn | barrn > thrsh_snow_barrn);

ok_veg_breed = all_trees >= thrsh_trees_breed & ...
    conifers >= thrsh_conif_breed & ...
    low_veg >= thrsh_low_veg_breed & ...
    ndvi_mn_fall(1:size(snow,1),:) >= thrsh_ndvi_breed;
%     ndvi_Sep(1:size(snow,1),:) >= thrsh_ndvi_breed;

% keyboard;

if thrsh_trees_stop > 0

     veg_stp = all_trees;
    thr_veg_stp = thrsh_trees_stop;

elseif thrsh_low_veg_stop > 0

    veg_stp = low_veg;
    thr_veg_stp = thrsh_low_veg_stop;

else % thrsh_all_veg_stop >= 0

    veg_stp = all_vegp{1}; % GPPp{20}; % %    all_veg; %        non_grass; %   
    thr_veg_stp = thrsh_all_veg_stop;


end

% poor_stops = ndvi_mn_fall(1:size(snow,1),:) < thrsh_ndvi_stop | ...
% poor_stops = (ndvi_Sep(1:size(snow,1),:) < thrsh_ndvi_stop_min) & ...
%     (all_trees(1:size(snow,1),:) < thrsh_trees_stop | ...
%     all_trees(1:size(snow,1),:) > max_trees_stop | ...
%     low_veg(1:size(snow,1),:) < thrsh_low_veg_stop | ...
%     all_veg(1:size(snow,1),:) < thrsh_all_veg_stop);

% poor_stops = (ndvi_Sep(1:size(snow,1),:) < thrsh_ndvi_stop_min) | ...
poor_stops = (ndvi(1:size(snow,1),:) <= thrsh_ndvi_stop_min) | ...
    (veg_stp <= min_any_veg_stp);

barrens = poor_stops & (snow > thrsh_snow_barrn | barrn > thrsh_snow_barrn);

