function [successful, succs_DepPs, succs_ArrPs, succs_DepArrPs, ...
    mn_succ_hds, std_succ_hds, mn_zk_hds,std_zk_hds, ...
    mn_zk_signps, std_zk_signps, ...
    mn_succ_hd_DepPs, mn_zk_hd_DepPs, mn_zk_signp_DepPs, ...
    std_succ_hd_DepPs, std_zk_hd_DepPs, std_zk_signp_DepPs, ...
    mn_intn_DepPs, mn_incl_DepPs, mn_dcl_DepPs, ...
    mn_Arr_lat_DepPs,mn_Arr_lon_DepPs,std_Arr_lat_DepPs,std_Arr_lon_DepPs, ...
    mn_Dep_lat_DepPs,mn_Dep_lon_DepPs,std_Dep_lat_DepPs,std_Dep_lon_DepPs, ...
    mn_lat_zk_DepPs,mn_lon_zk_DepPs,std_lat_zk_DepPs,std_lon_zk_DepPs, ...
    mn_Arr_lat_ArrPs,mn_Arr_lon_ArrPs,std_Arr_lat_ArrPs,std_Arr_lon_ArrPs, ...
    mn_Dep_lat_ArrPs,mn_Dep_lon_ArrPs,std_Dep_lat_ArrPs,std_Dep_lon_ArrPs, ...
    mn_lat_zk_ArrPs,mn_lon_zk_ArrPs,std_lat_zk_ArrPs,std_lon_zk_ArrPs, ...
    mn_Arr_lat_DepArrPs,mn_Arr_lon_DepArrPs,std_Arr_lat_DepArrPs,std_Arr_lon_DepArrPs, ...
    mn_Dep_lat_DepArrPs,mn_Dep_lon_DepArrPs,std_Dep_lat_DepArrPs,std_Dep_lon_DepArrPs, ...
    mn_lat_zk_DepArrPs,mn_lon_zk_DepArrPs,std_lat_zk_DepArrPs,std_lon_zk_DepArrPs, ...
    mn_std_inh_hds,std_std_inh_hds,mn_std_inh_sps,std_std_inh_sps, ...
    mn_disp_d_dep,std_disp_d_dep,mn_std_inh_disp_km,std_std_inh_disp_km, ...
    fr_zks,succ_heads,succ_zug_inh_hds, speed,blat_succs, blon_succs, ...
    coeff_succ, succ_zug_signps,all_zug_signps, ...
    all_inher_heads,all_zugkn_heads,all_geo_heads, all_mag_heads,all_sun_heads, ...
    all_coeffs, blats, blons,lat_bs_zugs,lon_bs_zugs, passed_Stops, ...
    arrived, survived, lost, stoppedAndArrived, ...
    all_ys_inher_heads, all_ys_init_geo_heads, all_ys_zugkn_heads, all_ys_zugkn_sps, ...
    all_ys_stoppedAndArrived, all_ys_surv_fuel, all_ys_surv_wtr, all_ys_surv_barrn_elev,  ...
    all_ys_lost, all_ys_late, n_fl_step, date_jul, ...
    std_inher_head,std_inher_steer,std_inher_signp,day_start_succ,all_day_starts, ...
    all_inclns,all_declns,all_tots,all_sun_az,all_ovr_wtr,all_barrn, ...
    succ_disp_d_dep,all_disp_d_dep,succ_disp_d_arr,all_disp_d_arr, ...
    all_dates,all_fl_hs,all_max_hrs_flt,all_stp_ovs,all_stp_nrs, ...
    all_zg_nrs,all_cum_hrs_flt,all_survs,all_fins,all_arrs,all_losts, ...
    incln_0s,decln_0s,intns_0s,idx_deps,idx_arrs,mn_lat_EU_DepPs, ...
    std_lat_EU_DepPs,mn_lon_EU_DepPs,std_lon_EU_DepPs, ...
    all_lat_arr_EU, all_lon_arr_EU,all_lat_zg, all_lon_zg, ...
    all_ys_lat_fin,all_ys_lon_fin] = ...
  simul_mag_sun_comp_rtes ...
        (geo_filnm,N_inds_sim,N_mig_Sims,N_wrmp_rnd, ...
        N_wrmp_bck,N_wrmp_1st,N_init_fact,start_date, ...
    zug_opt,zug_signp,req_land_signp,n_zugs,n_req_zugs, ...
    zgknk_req_in_stop, zgknk_req_stop_warm, require_stop_poly, ...
    require_stop_warmup,require_stopover_poly, ...
    require_stopover_zug,Dep_verts,Arr_verts,Stop_verts, ...
    magn_star_night,max_mig_dur,Va_mps, ...
    min_ns_fl_hs,max_ns_fl_hs,min_init_fl_hs,max_init_fl_hs, ...
    min_ns_sgnp_hs, max_ns_sgnp_hs, ... % May 2022 added distinct refuelling at signpost
    reOrWtrOpt,reOrWtrStr,n_wtr_stps, ...
    rand_wtr_fact,reOrReComp,dtrCstOpt,dtrCstStr,short_detrs, ...
    kmsAheadDtrCst,d_stop_max,dtrReComp, ...
    maxDtrSteps,stopCstOpt,fracStopCst,stpCstWtrFlt,stpovCstWtrFlt, ...
    endog_comp,calibr_comp,std_err_calibr,sun_head_reset,tc_comp_reset, ...
    tc_comp_wtr,sign_clockw_sun,geo_mag_refl_sun_opt,perf_init_calibr, ...
    magcl_reset_opt,magcl_opt,geomagn_model, ....
    max_magn_shift,min_std_inh_hd,max_std_inh_hd, ...
    min_std_inh_steer,max_std_inh_steer,std_err_incln_deg,std_rel_err_inten, ...
    min_std_inh_sp,max_std_inh_sp,opt_no_evo_std_yrs,wt_inher_speed, ...
    offset_dep_dusk,var_fl_hrs,stop_durn,stop_durn_sp,endur_flt, ...
    min_fl_stp_hs,max_fl_stp_hs,min_init_fl_stp_hs, ...
    req_land_init,req_land_stop,req_land_arr, ...
    include_wind,wind_str,switch_wind_opt,dep_wtr_TWs, ...
    min_disp_d_dep_km,max_disp_d_dep_km,min_disp_d_arr_km, ...
    max_disp_d_arr_km,min_disp_d_zug_km,max_disp_d_zug_km, ...
    min_std_inh_disp_km,max_std_inh_disp_km,  ...
    fwd_back_no_time,daily_surv,daily_surv_barrn, ...
    daily_surv_wtr,daily_surv_elev, ...
    thrsh_snow_barrn, thrsh_elev,thrsh_ndvi_breed, ...
    thrsh_ndvi_stop_min,thrsh_ndvi_stop_max, ...
    thrsh_trees_breed,thrsh_trees_stop,thrsh_conif_breed, ...
    thrsh_low_veg_breed,thrsh_low_veg_stop,thrsh_all_veg_stop, ...
    max_trees_stop,min_any_veg_stp,min_n_fls_stp,lost_Lat,n_par_pools, ...
    wt_natal_intns,wt_natal_incln,wt_natal_decln,opt_same_loc_brd,accel_wrmp,arr_EU_opt)


% determine number of initial (warmup or spinup) simulations to 
% create a viable population
N_init_Sims = N_wrmp_rnd +  N_wrmp_bck + N_wrmp_1st;

% thrshs are for defining hi elev or barren (incl. snow, can't stopover, may die)
% and extended stopover from ndvi

% Optional extra weighting natal dispersal based on closeness to geomag signature
wt_natl_geomg = wt_natal_intns || wt_natal_incln || wt_natal_decln;

% coast for interpolation coastal behaviour              
Coast = load('land_mask.mat','land','lat','lon');
  
ndvi_rng = thrsh_ndvi_stop_max - thrsh_ndvi_stop_min;

% load geomagnetic, elevation and landcover data
[coastlat,coastlon, igrfcoefs, hi_elevs, barrens, ...
    poor_stops, ok_veg_breed, ndvi, veg_stp, thrsh_veg_stop_max] = ...
    load_geo_data(geo_filnm, thrsh_snow_barrn, ...
    thrsh_elev,thrsh_ndvi_breed, thrsh_ndvi_stop_min, thrsh_trees_breed, ...
    thrsh_trees_stop,max_trees_stop,thrsh_conif_breed,thrsh_low_veg_breed, ...
    thrsh_low_veg_stop,thrsh_all_veg_stop,min_any_veg_stp);

% raneg of thrsh veg for stop
veg_rng = thrsh_veg_stop_max - min_any_veg_stp; 

% job = createJob('confiuration','local');
pool_locs = parpool('local',n_par_pools);

stop_flag = false;

% For option to "look ahead", use visual cues (equivalent to star comp) 
star_cmp_flag = 2;

% Handy function, diff between two angles
normalizeDeg=@(x)(mod(x,360));
absDiffDeg=@(a,b)abs(normalizeDeg(normalizeDeg(a)-normalizeDeg(b)));  
% normalizeRad=@(x)(mod(x,2*pi));
absDiffRad=@(a,b)abs(angdiff(a,b));  

% Options to add offsets to (sun and magn) compasses and also inher heading based
% on one of these two (or if == 3, assume geo compass is imprinted onto 
% magn and/or sun compasses - presumably better able to calibrate on natal
% grounds)
incl_sun = calibr_comp(2) ~=0 || endog_comp == 2;

% flag to add declination to compass and inheritancefor magnetic compass
% use (1st element of calibration compass vector)
add_decl_comp = calibr_comp(1) ~= 0;
add_decl_inher = endog_comp==1;

% Options for sun compass opposite to azimuth (counter-clockwise)
refl_sun_mag = incl_sun && sign_clockw_sun== -1 && geo_mag_refl_sun_opt == 2;
% Flag to offset for sun compass 
add_az_comp = calibr_comp(2) ~= 0;
add_az_inher = endog_comp==2;

% initialize times
yr_start = start_date(1);
month_start = start_date(2);
month_start_strs = {'jan','feb','mar','apr','may', ...
    'jun','jul','aug','sep','oct','nov','dec'};
month_start_str = month_start_strs{month_start};
day_start = start_date(3);
std_day_start = start_date(4);
max_dev_day_st = start_date(5);

rad2deg = 180/pi; 
deg2rad = pi/180;
R_Earth_km = 6371;

% Conversion factos to scale mean dispersal to the gaussian term used in
% script update_successful_population

    % arc_1_deg = R_Earth_km*pi/180;
    % distances from natal and wintering locns computed in degs below
    % R_natal = R_natal_km/R_Earth_km*rad2deg;
    % R_winter = R_winter_km/R_Earth_km*rad2deg;
    % From May 2022 we should input mean disp radius 
    % but convert to std for Gaussian distr - 
    % so add extra factor sqrt(pi/2) and also extra sqrt(2)max_std_inh_disp_km
    % rather than use 2*sig^2 in update popn parallel loop
    % I.e., the sqrt(2)'s cancel
     mn_2_sig = R_Earth_km/sqrt(pi);
    
     min_disp_d_dep =  min_disp_d_dep_km/mn_2_sig;
     max_disp_d_dep =  max_disp_d_dep_km/mn_2_sig;
     min_disp_d_arr =  min_disp_d_arr_km/mn_2_sig;
     max_disp_d_arr =  max_disp_d_arr_km/mn_2_sig;
     min_disp_d_zug =  min_disp_d_zug_km/mn_2_sig;
     max_disp_d_zug =  max_disp_d_zug_km/mn_2_sig;
     rng_disp_d_dep = max_disp_d_dep - min_disp_d_dep;
    
     % Aug 2022 same as with hd and sp, use std disp at least in warmup
     % Idea being it helps convergence to 'optimum' ground state
     % before the actual simulations (march thru years)
     min_std_inh_disp =  min_std_inh_disp_km/mn_2_sig;
     max_std_inh_disp =  max_std_inh_disp_km/mn_2_sig;
     rng_std_inh_disp = max_std_inh_disp - min_std_inh_disp;

% Convert air speed to km/h (via m/s)
Va = Va_mps*3.6; % 0; % 
% calcuLate hourly change Lat if Southbound
% Lon will have 1/sin(theta) factor (below)

hourly_del_Lat = Va/R_Earth_km;    
% max distance aim for land leaving coast
radAheadDtrCst = kmsAheadDtrCst/R_Earth_km;
% and max dist stopping on coast when over water at dawn
max_rad_stop_dawn = d_stop_max/R_Earth_km;

min_std_inh_hd = min_std_inh_hd*deg2rad;
max_std_inh_hd = max_std_inh_hd*deg2rad;
rng_std_inh_hd = max_std_inh_hd - min_std_inh_hd;

if zug_opt && zug_signp <= 2 % convert angles degs to rads
    min_std_inh_sp = min_std_inh_sp*deg2rad;
    max_std_inh_sp = max_std_inh_sp*deg2rad;
    rng_std_inh_sp = max_std_inh_sp - min_std_inh_sp;
elseif zug_opt % percentages (keep)
    rng_std_inh_sp = max_std_inh_sp - min_std_inh_sp;
else
    rng_std_inh_sp = NaN; % shouldn't be used
end

is_magcl = (calibr_comp(1) == 2 || calibr_comp(1) == 3);

min_std_inh_steer = min_std_inh_steer;
max_std_inh_steer = max_std_inh_steer;

compass_mag_std = std_err_calibr(1)*deg2rad;
% Kappa for von Mises angular distr
compass_cal_kappa(1) = 1/compass_mag_std^2;

compass_sun_std = std_err_calibr(2)*deg2rad;
% Kappa for von Mises angular distr
compass_cal_kappa(2) = 1/compass_sun_std^2;

compass_geo_std = std_err_calibr(3)*deg2rad;
% Kappa for von Mises angular distr
compass_cal_kappa(3) = 1/compass_geo_std^2;

% pointer to non-zero compass std errs
cal_comp_nonzero = std_err_calibr ~=0;
% idx_cal_comp_errs = calibr_comp ~= 0 & cal_comp_nonzero;

% inher_steer_kappa = 1/std_inher_steer^2;
std_err_incl = std_err_incln_deg*deg2rad;
compass_kappa_incl = 1/std_err_incl^2;

n_sgnps = max(n_zugs,zug_opt);

% create arrival and stopover polygons from cell array of vertices
nDep_polys = numel(Dep_verts);
nArr_polys = numel(Arr_verts);
nStop_polys = numel(Stop_verts);

[Lon_dep_poly, Lat_dep_poly] = create_poly(Dep_verts);
[Lon_arr_poly, Lat_arr_poly] = create_poly(Arr_verts);

% separate poly for each stopover polygon
% which themselves could consist of several disjoint polys
for iSt = 1:nStop_polys
    [Lon_stop_poly{iSt}, Lat_stop_poly{iSt}] = create_poly(Stop_verts{iSt});
end

% calculate mean Lon in 1st stopover for determination of initial
% headings (needed when more than one choice in stopover area)
nSt_Ply_1 = numel(Stop_verts{1});
for iSt = 1:nSt_Ply_1
    mn_Lon_St_1(iSt) =  ...
        180/pi*mean(Lon_stop_poly{1}([(iSt-1)*6+2 (iSt-1)*6+3]));
    mn_Lat_St_1(iSt) =  ... % 180/pi*
        180/pi*mean(Lat_stop_poly{1}([(iSt-1)*6+1 (iSt-1)*6+2]));
end

% calculate subPolys for option to control miratory connectivity 
% e.g. during warm-up of polulation
[~, Lon_Arr_Subpoly, Lat_Arr_Subpoly] = ...
det_subPolys(Lon_arr_poly,Lat_arr_poly); % 180/pi* 180/pi*

for iArr = 1:nArr_polys 
    
    mn_Lon_Arr_deg(iArr) = 180/pi* ...
        mean(Lon_arr_poly([(iArr-1)*6+2 (iArr-1)*6+3]));
    
    % for allowing expanxion into neighbouring arr poly from
    % previous generation
       idx_close = max(iArr-1,1):min(iArr+1,nArr_polys);
       Lon_nbr_subP{iArr} = [];
       Lat_nbr_subP{iArr} = [];
       for idx_sub = idx_close
           Lon_nbr_subP{iArr} = [Lon_nbr_subP{iArr} Lon_Arr_Subpoly{idx_sub} NaN];
           Lat_nbr_subP{iArr} = [Lat_nbr_subP{iArr} Lat_Arr_Subpoly{idx_sub} NaN];
       end
%     mn_Lat_Arr(iArr) =  ... % 180/pi*
%         mean(Lat_arr_poly{1}([(iArr-1)*6+1 (iArr-1)*6+2]));
end

% for determining if flew past arr poly (used in update subroutine)
% arrLat_width2 = Arr_verts{end}(4) - Arr_verts{1}(3);
% arrLon_width2 = Arr_verts{end}(2) - Arr_verts{1}(1); % mod(Arr_verts(2),2*pi) - mod(Arr_verts(1),2*pi);
% mean_arr_Lat = (Arr_verts(4) + Arr_verts(3))/2;

min_Lat_Arr = nanmin(Lat_arr_poly);
max_Lat_Arr = nanmax(Lat_arr_poly);
arrLat_width2 = (max_Lat_Arr - min_Lat_Arr)/2;
min_Lon_Arr = nanmin(Lon_arr_poly);
max_Lon_Arr = nanmax(Lon_arr_poly);
arrLon_width2 = (max_Lon_Arr - min_Lon_Arr)/2;

mean_arr_Lat = min_Lat_Arr + arrLat_width2;
mean_arr_Lon = min_Lon_Arr + arrLon_width2;

% in case polys are > 180 use simple addition
min_Lon_Dep = nanmin(Lon_dep_poly);
max_Lon_Dep = nanmax(Lon_dep_poly);
depLon_width2 = (max_Lon_Dep - min_Lon_Dep)/2;
mean_dep_Lon = min_Lon_Dep + depLon_width2;
% mean_dep_Lon = (Dep_verts(1) + Dep_verts(2))/2;

% for overshoots
min_Lat_Dep = nanmin(Lat_dep_poly);
max_Lat_Dep = nanmax(Lat_dep_poly);
depLat_width2 = (max_Lat_Dep - min_Lat_Dep)/2;

% significant E-W route to rule out birds flying (E/W) past their arrival areas
% mean_arr_Lon = circ_mean([mod(Arr_verts(1),2*pi) mod(Arr_verts(2),2*pi)]');

% in case polys are > 180 use simple addition for 
% mean_dep_Lon = (Dep_verts(1) + Dep_verts(2))/2;
% circ_mean([mod(Dep_verts(1),2*pi) mod(Dep_verts(2),2*pi)]');
% 
wide_EW_route = abs(mean_arr_Lon-mean_dep_Lon) > arrLon_width2;
sign_mean_arr_dep_Lon = sign(mean_arr_Lon-mean_dep_Lon);

% transform to -pi to pi
% mean_arr_Lon = shiftAnglesFromMinus180To180(mean_arr_Lon*rad2deg)*deg2rad;
% mean_dep_Lon = shiftAnglesFromMinus180To180(mean_dep_Lon*rad2deg)*deg2rad;

req_stps = find(require_stop_poly);
req_stps_warmup = find(require_stop_warmup);
n_req_stp = sum(require_stop_poly);
n_req_stp_warmup = sum(require_stop_warmup);

if zug_opt ==1
    
    reqs_zug = 1:n_req_zugs; % find(zgknk_req_in_stop==1);
%     reqs_zug(end+1) = nStop_polys+1;   
    
else
    
    reqs_zug = []; % 1:nStop_polys;  % 1:max(n_zugs+1,nStop_polys);
    
end

mn_std_inh_hds = NaN*ones(N_mig_Sims,1);
std_std_inh_hds = mn_std_inh_hds;

mn_std_inh_sps = mn_std_inh_hds;
std_std_inh_sps = mn_std_inh_hds;

mn_disp_d_dep = mn_std_inh_hds;
std_disp_d_dep = mn_std_inh_hds;

decln_0s = NaN*ones(N_inds_sim,N_mig_Sims);
incln_0s = NaN*ones(N_inds_sim,N_mig_Sims); 
intns_0s = NaN*ones(N_inds_sim,N_mig_Sims);

all_ys_init_geo_heads = NaN*ones(N_inds_sim,N_mig_Sims);
all_ys_inher_heads = NaN*ones(N_inds_sim,N_mig_Sims);
all_ys_zugkn_heads = NaN*ones(N_inds_sim,N_mig_Sims);
all_ys_zugkn_sps = NaN*ones(N_inds_sim,N_mig_Sims);

all_ys_stoppedAndArrived = false(N_inds_sim,N_mig_Sims);

all_ys_surv_fuel = false(N_inds_sim,N_mig_Sims);
all_ys_surv_wtr = false(N_inds_sim,N_mig_Sims);
all_ys_surv_barrn_elev = false(N_inds_sim,N_mig_Sims);
all_ys_lost = false(N_inds_sim,N_mig_Sims);
all_ys_late = false(N_inds_sim,N_mig_Sims);

all_lat_arr_EU = NaN*ones(N_inds_sim,N_mig_Sims);
all_lon_arr_EU = NaN*ones(N_inds_sim,N_mig_Sims);
all_ys_lat_fin = NaN*ones(N_inds_sim,N_mig_Sims);
all_ys_lon_fin = NaN*ones(N_inds_sim,N_mig_Sims);
all_lat_zg = NaN*ones(N_inds_sim,N_mig_Sims);
all_lon_zg = NaN*ones(N_inds_sim,N_mig_Sims);


% factor shortening detoured flights
dtrShortFact = 1 + (short_detrs == 1);
       
    successful = [];
    mn_zk_hds = [];
    std_zk_hds = [];
    mn_zk_signps = []; 
    std_zk_signps = [];

    mn_zk_hd_DepPs = [];
    mn_zk_signp_DepPs = [];
    std_zk_hd_DepPs = [];
    std_zk_signp_DepPs = [];
    is_succ_zks = [];
    fr_zks = [];
    mn_lat_zks = [];
    mn_lon_zks = [];
    std_lat_zks = [];
    std_lon_zks = [];
    
    mn_lat_zk_DepPs = []; 
    mn_lon_zk_DepPs = []; 
    std_lat_zk_DepPs = []; 
    std_lon_zk_DepPs = []; 
    mn_lat_zk_ArrPs = []; 
    mn_lon_zk_ArrPs = []; 
    std_lat_zk_ArrPs = []; 
    std_lon_zk_ArrPs = []; 
    mn_lat_zk_DepArrPs = [];
    mn_lon_zk_DepArrPs = [];
    std_lat_zk_DepArrPs = [];
    std_lon_zk_DepArrPs = [];
    mn_Dep_lat_zk_DepArrPs = []; 
    mn_Dep_lon_zk_DepArrPs = []; 
    std_Dep_lat_zk_DepArrPs = []; 
    std_Dep_lon_zk_DepArrPs = []; 
    mn_Arr_lat_zk_DepArrPs = []; 
    mn_Arr_lon_zk_DepArrPs = []; 
    std_Arr_lat_zk_DepArrPs = []; 
    std_Arr_lon_zk_DepArrPs = [];
    std_Dep_lat_DepArrPs = [];
    std_Dep_lan_DepArrPs = [];   
    std_Arr_lat_DepArrPs = [];
    std_Arr_lan_DepArrPs = [];  

    mn_lat_EU_DepPs = [];
    std_lat_EU_DepPs = [];   
    mn_lon_EU_DepPs = [];
    std_lon_EU_DepPs = [];   

%       prev_succ_heads = [];
%     prev_succ_zug_inh_hds = [];
%     prev_speed = [];
% 
%     blat_prev_succs = [];
%     blon_prev_succs = []; 
% 
%     prev_coeff_succ = []; 
     day_start_succ = [];
     mig_dur = [];
     speed = [];
     succ_heads = [];
     lat_deps_succ = [];
     lon_deps_succ = [];   
     lat_arrs_succ = [];
     lon_arrs_succ = [];   
     blat_succs = [];
     blon_succs = [];
     succ_zug_inh_hds = [];
    %                  plot_lat_lon_succ_fail_itern
     coeff_succ = []; 
     succ_std_inher_head = []; 
     succ_std_inher_steer = []; 
     std_inher_disp = [];
    %                 coeff_succ = all_coeffs(stoppedAndArrived,:); 
    succ_zug_signps = [];                
     if zug_opt == 1
         succ_std_inher_signp = []; 
     end

     succ_disp_d_dep = [];   
     succ_disp_d_arr = [];
     succ_disp_d_zug = [];
    
    iYear = 1;
    iWarmup = min(1,N_init_Sims);
    
    N_succ = 1;

    %% main loop starts here
    while iYear <= N_mig_Sims && N_succ > 0   
        
        N_inds = N_inds_sim*(1 + (iWarmup <= N_init_Sims)*(N_init_fact-1));

        % all simluations finished for this (warmup/) year
        Sim_nr = iWarmup + iYear -1;
                     
        % First, determine if any required stopover (staging)
        % areas to visit
        if Sim_nr <= N_init_Sims % iWarmup % n_req_stp % && zug_opt ~= 1 % ~isempty(req_stps) % 

           n_zg_i = n_req_zugs; % zug_opt*n_zugs;
            
           zg_rq_stp =  zgknk_req_stop_warm;
           
           if isempty(req_stps_warmup)
               
               lst_req_stp = 0; 
               
           else
               
               lst_req_stp = min(Sim_nr,req_stps_warmup(end)); 
           
           end

           if fwd_back_no_time == 0
               curr_year = yr_start;
           elseif Sim_nr <= N_wrmp_rnd
                curr_year = yr_start + floor(rand*(N_mig_Sims)); % 80 + floor(rand*(43)); %  /2 % 
           elseif fwd_back_no_time ==1  && Sim_nr <= N_wrmp_rnd+N_wrmp_bck 
               curr_year = yr_start + N_wrmp_rnd+N_wrmp_bck - Sim_nr; % 2000-51 = 1949
           elseif fwd_back_no_time ==2  && Sim_nr <= N_wrmp_rnd+N_wrmp_bck
               curr_year = yr_start + N_mig_Sims -1 + Sim_nr - N_init_Sims; %  1966 % 
           elseif fwd_back_no_time ==1 
                curr_year = yr_start; %  1900
           else
               curr_year = yr_start + N_mig_Sims -1; % 2015
           end
                       
        else

            n_zg_i = n_req_zugs;
                       
            zg_rq_stp =  zgknk_req_in_stop;
            
           if isempty(req_stps)
               
               lst_req_stp = 0; 
               
           else
               
                lst_req_stp = min(Sim_nr,req_stps(end)); 
           
           end

%            if curr_year < yr_start + N_mig_Sims -1

               curr_year = (fwd_back_no_time == 0)*yr_start + ...
                (fwd_back_no_time == 1)*(yr_start + iYear-1) + ...   
                (fwd_back_no_time == 2)*(yr_start + N_mig_Sims - iYear);

%            else
%             
%                curr_year = curr_year - 1; % (100+ yr_start + iYear-1);
% 
%            end
           
        end
                   
        zero_vec = zeros(N_inds,1);
        zero_ten_vec = zeros(N_inds,10);
        one_vec = ones(N_inds,1);
        false_vec = false(N_inds,1);
        true_vec = true(N_inds,1);

        % day of year will be determined according to 
        % std in dep date (and inheritance in subsequent years)

%         if Sim_nr <= true %

%         curr_date.year = curr_year*one_vec;
%         curr_date.month = month_start*one_vec;
        
        % date nr jan 1st used to update julian date using flight
        % durn using suncycle below
%         date_jan_1 = datenum(curr_year,1,1);  

        % now derive arrivals and initial headings if first iteration (year)
        if ~isinf(compass_kappa_incl)
             rand_magn_0 = vmrand(0, compass_kappa_incl, [N_inds 1]); %      
        else
            rand_magn_0 = zero_vec;
        end

        if iWarmup <= 1 &&  iWarmup <= N_init_Sims && iYear < 2  
            
              create_init_population    

        else % "iYear" or "iWarmup" > 1, use 'inherited' values init heads and any 'steering'
            
           update_successful_population

        end
        
        location.longitude = llamda_0*rad2deg; 
        location.latitude = theta_0*rad2deg;  

        %% 
                
    % initial or updated headings now defined
    % except for magnetoclinic, defined below

    if n_zugs > 0
        lat_bs_zugs = NaN*ones(N_inds,n_zugs);
        lon_bs_zugs = NaN*ones(N_inds,n_zugs);
    else
        lat_bs_zugs = NaN*ones(N_inds,1);
        lon_bs_zugs = NaN*ones(N_inds,1);          
    end


% compass & steering options

    % compass imprecision (std in degrees)

    % max 'shift' in orientation by any change in magn compass
    % and exponential 'rate' of steering vs. relative change in magnetic parameter
    if calibr_comp(1) == 0

        d_alpha_max = 0;

    else

        d_alpha_max = max_magn_shift; % deg2rad*[90 180];

    end

    % check if starts over land
         over_water_0 = ~island(location.latitude,location.longitude); % ,cst_accuracy);   
         % land_or_ocean(location.latitude,location.longitude,cst_accuracy);    

    % initialize model

    % initially does not reorient (reorientation applies to some 
    % cases when over ocean) nor detour (reorient in advance when water
    % ahead)
    reorients_0 = 0;
    detours_0 = 0;

     % update magcl headings
    if is_magcl

        if iWarmup == 1 
            
            % sign_inher_heads used in magncl strats as bailout orientation
            % when angle can't be realized (so inherited angle is relative
            % to magn axis and sign detern=mines whether E or W is flown
            % with 'transverse' magncl strat)
            % sign_inher_heads = sign(fn_heads);
            sign_inher_heads = sign(inher_heads + (endog_comp == 1)*decln_0);

            if magcl_opt == 1 % as per Kiep use projn of incl along hz axis
                % check if magnetoclinic and set endogenous projected headings 

            %     magncl_0 = tan(incln_0)./sin(fn_heads);
            %      magncl_parll_0 = tan(incln_0)./cos(fn_heads);
                magncl_0 = tan(incln_0)./sin(inher_heads);
                 magncl_parll_0 = tan(incln_0)./cos(inher_heads);

            elseif magcl_opt == 2 % use "sfc" of incl along geomag axis

                magncl_0 = sin(inher_heads).*tan(incln_0);
                magncl_parll_0 = cos(inher_heads).*tan(incln_0);

            else % use all hz compt and hz projn of vert compt

                % to implement

            end

            % for magncl strats we still need to calc init heads from
            % (inherited) projections
            if calibr_comp(1) == 2

                coeff_magn_steer(:,1) = magncl_0;

            elseif calibr_comp(1) == 3

                coeff_magn_steer(:,1) = magncl_parll_0;

            end

        else % "iYear" or iWarmup >  1 

            update_popn_endog_magcl_heads

        end
        
%         if mod(iYear,20) == 0
%             
%             save('model_20_yrs')
%             
%         end
        % for magnetoclinic strats, coeffs already inherited
        % (otherwise not used, but all magn parameters need to exist 
        % for setting of headings)
        if calibr_comp(1) == 2
            magncl_0 = coeff_magn_steer(:,1); % (N_inds,1); %
            magncl_parll_0 = NaN*one_vec;
        else
            magncl_0 = NaN*one_vec;
            magncl_parll_0 = coeff_magn_steer(:,1);
        end
        
    else
       
        magncl_0 = NaN*one_vec; % (N_inds,1); % 
        magncl_parll_0 = NaN*one_vec;
        
        
    end

        % keep magn parameter at constant (1) if no magn steering
%         magn_param_0 = (calibr_comp(1) == 0) + (calibr_comp(1) == 1)*incln_0 + ...
%             (calibr_comp(1) == 2)*vt_inten_0 + (calibr_comp(1) == 2)*magncl_0 + ...
%             (calibr_comp(1) == 3)*magncl_parll_0 + (calibr_comp(1) == 5)*hz_inten_0 ...
%             + (calibr_comp(1) == 6)*tot_inten_0;
        if is_magcl
            switch calibr_comp(1)
                case  2
                    magn_param_0 = magncl_0; %
                case 3
                     magn_param_0 = magncl_parll_0;
            end
        else
                magn_param_0 = NaN*one_vec;
        end

        % initialize loop variables with dep values
        theta = theta_0;
        llamda = llamda_0;
        last_lon_land = llamda_0;
%         incln_low = incln_0;
        decln = decln_0;
               
        del_llamda = zero_vec;
        del_ll_last_land = zero_vec; 
        del_ll_ref_lon_tc = zero_vec;
        
        del_alph_dtr = zero_vec;
        prev_del_alph_dtr = zero_vec;
        
        stops_over = false_vec;
        cumul_flts = zero_vec;
        over_water = over_water_0;
        stop_nxt_land = false_vec;
        reorients = zero_vec;
        detours = zero_vec;
        recmp_dtr= zero_vec;
        
        d_time_zone = zero_vec;
        
        surv_wtr = true_vec;
        surv_elev = true_vec;
        surv_barrn = true_vec;
        surv_land = true_vec;

        % alpha's are geographical headings, adjusted for magn decl (needed
        % for movemnet, wind etc)
%         alpha_base_0 = inher_heads +  (calibr_comp(2) ==3)*decln_0;

        % set initial heading (alpha) error per sun or magn compass
        for i_comp = 1:3
  
                if calibr_comp(i_comp) > 0 || endog_comp == i_comp 
                    if cal_comp_nonzero(i_comp) > 0 && perf_init_calibr(i_comp) == 0
                        d_alpha_0(:,i_comp) = vmrand(0, compass_cal_kappa(i_comp), [N_inds 1]);
                    else
                        d_alpha_0(:,i_comp) = zero_vec;
                    end
                % add magn varn to sun compass if refl about mag S axis
                elseif i_comp == 1 && refl_sun_mag && ...
                        (calibr_comp(2) > 0 || endog_comp == 2) 
                    if cal_comp_nonzero(i_comp) > 0 && perf_init_calibr(i_comp) == 0
                        d_alpha_0(:,i_comp) = vmrand(0, compass_cal_kappa(i_comp), [N_inds 1]);
                    else
                        d_alpha_0(:,i_comp) = zero_vec;
                    end
                else
                    d_alpha_0(:,i_comp) = zero_vec;
                end
                             
        end
                               
        % need to define base compass headings dependeing on which is
        % endogenous and whether they are initially perfectly assessed
        % and what any transfer error is made for the non-endogenous
        % compasses
        
        % for sun azimuth compass, account for reflection about geomag axis
        % option
        
        date_jul = date_jul_start;
        
        decln_0_hat = decln_0 + d_alpha_0(:,1);
        sun_az0_hat = sun_az_0 + d_alpha_0(:,2);    
        sun_az_1st_nt = sun_az0_hat;
        
        % add decln and sun zenith if appropriate and random initial calibration
        % possibly deflected Sun Az
        sun_defl_0_hat = incl_sun*(sun_az0_hat*sign_clockw_sun + ...
            2*(refl_sun_mag)*decln_0_hat);

        % determine actual heading taken
        alpha_0 = mod(inher_heads + (endog_comp == 1)*decln_0_hat + ...
            (endog_comp == 2)*sun_defl_0_hat + ...
            (endog_comp == 3)*d_alpha_0(:,3) + pi,2*pi) - pi; 
        
        % define as SE or SW based on alpha_0
        % Due South counts as both (for overshoots in either direction)
        SE = alpha_0 <= 0;
        SW = alpha_0 >= 0;

% Identify appropriate arrival longitude limit given initial flight
% direction (helps flag lost migrants to speed up computation)

         arr_lon_lim = SW.*mod(min_Lon_Arr,2*pi) + ...
                            SE.*mod(max_Lon_Arr,2*pi);

% flag sign of original difference with lon limit
         pos_d_0_lon_lim = arr_lon_lim - mod(llamda_0,2*pi) > 0;

        % add error in perceived geographic compass relative to others 
        % (for case of endog mag or sun comp transfer to geogr)
        alpha_0_hat = mod(alpha_0 - (endog_comp ~= 3)*d_alpha_0(:,3) + pi,2*pi) - pi;
        
        % determine "base" headings for each compass, whether inherited
        % (endog) or transferred on (this) first departure
        
        switch endog_comp
            
            case 1 % magn
  
                alpha_bases{1} = inher_heads;
                absDiff_sun = absDiffDeg(alpha_0*rad2deg,sun_defl_0_hat*rad2deg)*pi/180;       
                alpha_bases{2} =  mod(absDiff_sun + pi,2*pi) - pi;     
                alpha_bases{3} =  alpha_0_hat;
                
            case 2 % sun
                
                alpha_bases{2} = inher_heads; %
                absDiff_mag = absDiffDeg(alpha_0*rad2deg,decln_0_hat*rad2deg)*pi/180;
                alpha_bases{1} =  mod(absDiff_mag + pi,2*pi) - pi;            
                alpha_bases{3} =   alpha_0_hat;
                
            case 3 % geogr
                
                alpha_bases{3} =  inher_heads;  
                absDiff_mag = absDiffDeg(alpha_0_hat*rad2deg,decln_0_hat*rad2deg)*pi/180;
                alpha_bases{1} =  mod(absDiff_mag + pi,2*pi) - pi;     
                absDiff_sun = absDiffDeg(alpha_0_hat*rad2deg,sun_defl_0_hat*rad2deg)*pi/180;       
                alpha_bases{2} =  mod(absDiff_sun + pi,2*pi) - pi;    
                 
        end

        alpha = alpha_0;
        can_dtr = false_vec; % dtrCstOpt*(dist2coast(theta_0*rad2deg,llamda_0*rad2deg)<kmsAheadDtrCst); % dtrCstOpt*one_vec;
        
    % keep track of flight hours since 'last over land' 
%         cum_ns_flt_dur = zero_vec;
        cum_hrs_flt = zero_vec; % fl_hrs;
        toCmpDtr = zero_vec;
        mn_Dtr = zero_vec;
        cum_steps_dtr = zero_vec;
        n_toCmpReOr = zero_vec;
        mn_ReOr = zero_vec;
        
        over_wtr_fracs = zero_ten_vec;
%         lat_will_land_fr = zero_ten_vec;
%         lon_will_land_fr = zero_ten_vec;        
        
        [u_wind, v_wind] = det_stoch_wind(theta,include_wind,wind_str,NaN,NaN,0); % d
                
        determine_detours_barrier_crossings
        
        % nonzero calibr_comp(1) means compass is calibrated
        % using magnetic field. magn_star_night used thru night 
        % (in update_loc below) differentiates between magnetic and star compass
%         alpha = alpha_base_0;
        alpha_lst_lnd = alpha;
        
%         passed_zug_nr = zeros(N_inds,n_zugs);
%         sun_az_1st_nt = ref_azi;
%         set_rise(1:N_inds,:) = [civ_dusk_0 civ_dawn_0];
        magn_param = magn_param_0;
%         date_jul = curr_date.day; % one_vec;
%         date_jul_0 = date_jul;
%         fl_hrs = fl_hrs_0; % (~endur_flt)* + 12*endur_flt*one_vec;
        stop_num = zero_vec;
        n_stops = zero_vec;
        skip_stop_over = false_vec;
        % initially no deviation from solar or geomagnetic fields
        d_alpha_sun = zero_vec; % (calibr_comp(2) ~= 0)*sun_az0_hat;
        % offset for time compensation (or great circle)
        d_time_lag = zero_vec;
        d_time = zero_vec;
        d_alpha_magn = zero_vec;
        d_alpha_rnd = zero_vec;
%         decln = zero_vec;
        
        prev_fwd_flt = 1;
%         prev_cv_dk = civ_dusk_0;
%         prev_cv_dw = civ_dusk_0;
%         civ_dusk = civ_dusk_0;
%         civ_dawn = civ_dawn_0;
        
        % u_wind, v_wind is wind leading to initial locn
        % needed for case of switch_wind_opt 

        % assume heading was same leading to initial locn (not crucial)
        
        if include_wind == 1
            tail_w = calc_wind_supprt(alpha, u_wind, v_wind);
        else
            tail_w = zero_vec;
        end
 
%         try
        % check if needs tailwinds before water crossing
        if dep_wtr_TWs == 1 && any(tail_w < 0)
            
            idx_hw = find(tail_w < 0);
            
            % location depends on geogr heading incl. possibly change in decln
             [ahead_theta(idx_hw), ahead_llamda(idx_hw)] = update_loc(...
                1,radAheadDtrCst,theta(idx_hw), ...
                llamda(idx_hw),alpha(idx_hw),0,0,0, ...
                star_cmp_flag,decln(idx_hw),date_jul(idx_hw),geomagn_model);
            
%              hws = idx_hw(land_or_ocean(ahead_theta(idx_hw)*rad2deg, ...
%                 ahead_llamda(idx_hw)*rad2deg,cst_accuracy) == 1);

             hws = idx_hw(island(ahead_theta(idx_hw)*rad2deg, ...
                ahead_llamda(idx_hw)*rad2deg) == 0);
            
               while ~isempty(hws)
                   
                    [u_wind(hws), v_wind(hws)] = ...
                        det_stoch_wind(theta(hws),include_wind,wind_str, ...
                        u_wind(hws),v_wind(hws),switch_wind_opt); % 
                    tail_w(hws) = calc_wind_supprt(alpha(hws), ...
                    u_wind(hws), v_wind(hws));

                    date_jul(hws) = date_jul(hws) + 1;
%                     date_nrs(hws) = date_nrs(hws) + 1;
                    
                    tail_w(hws) = calc_wind_supprt(alpha(hws), ...
                        u_wind(hws), v_wind(hws));
                    
                    hws = hws(tail_w(hws) <0);
                                
               end
            
        end
%         catch
%             keyboard
%         end
            
        % n_fl_step is flght step (date_step varies with flight and stopover durations)
        n_fl_step = one_vec;

        arrived = false_vec;
        survived = true_vec;

        % set max hours over water from max_init_fls_wtr
        % Can assume they have better potential for endurance flight 
        % at natal site than subsequently
        
        % Determine min and max nonstop fl hrs per bird
        % will depend on departure polygon (eg for Atlantic wheratears will
        % be higher than for others)
        % Initial flt hours may differ from further along
        % Typically will be less (builid) on route
        % but starting near barrier (Nearct. Whtr) will be higher
        if nDep_polys > 1
            max_init_hs_bs = max_init_fl_hs(idx_deps)'; 
            min_init_hs_bs = min_init_fl_hs(idx_deps)'; 
            max_ns_hs_bs = max_ns_fl_hs(idx_deps)'; 
            min_ns_hs_bs = min_ns_fl_hs(idx_deps)';  
            max_sgnp_hs_bs = max_ns_sgnp_hs(idx_deps)'; 
            min_sgnp_hs_bs = min_ns_sgnp_hs(idx_deps)'; 
        else
            max_init_hs_bs = max_init_fl_hs*one_vec; 
            min_init_hs_bs = min_init_fl_hs*one_vec;  
            max_ns_hs_bs = max_ns_fl_hs*one_vec; 
            min_ns_hs_bs = min_ns_fl_hs*one_vec; 
            max_sgnp_hs_bs = max_ns_sgnp_hs*one_vec; 
            min_sgnp_hs_bs = min_ns_sgnp_hs*one_vec; 
        end

        % Determine individual (initial) max num (n.s.) flight hrs
        curr_max_ns_fl_hs = ...
            rand([N_inds 1]).*(max_init_hs_bs - min_init_hs_bs) + min_init_hs_bs; 

        % passed Stops used when passing stopover polygons is obligatory 
        passed_Stops = false(N_inds,max(nStop_polys,1));
        passed_Stop_nr = zero_vec;
        % counter for which (of any) Zugknicks have been passed
        % currently using signposts rather than stopover region
        zug_nr = zero_vec;
        
        % initially not lost :)theta_0
        lost = false_vec;
        
        finished = false_vec;
        is_late = false_vec;
        
        incln = incln_0;
        decln = decln_0;
        tot_inten = tot_inten_0;
        vt_inten = vt_inten_0;
        hz_inten = hz_inten_0;
        
        % for last Sim
        blats = theta_0; % lat_bs_deps(1:N_inds,1);
        blons = llamda_0; % lon_bs_deps(1:N_inds,1);;
        % set ref lon for time compensation
        ref_lon_tc = llamda_0;
        all_inher_heads = inher_heads;
        all_geo_heads = alpha_0*180/pi + 180; % alpha_0*180/pi + 180;
        all_mag_heads = all_geo_heads - decln_0*180/pi;
        all_sun_heads =  (alpha_0 - sun_az_0)*180/pi;
        all_inclns = incln_0*180/pi;
        all_declns = decln_0*180/pi;
        all_tots =  tot_inten_0;
        all_sun_az = sun_az_0*180/pi+180;
        all_dates = date_jul_start;
        all_fl_hs = fl_hrs;
%         all_cv_dsks = civ_dusk_0;
%         all_cv_dwns = []; % civ_dawn_0;
        all_ovr_wtr = over_water;
        % begin not barren
        all_barrn = zero_vec;
        
        all_stp_ovs = zero_vec;
        all_stp_nrs = zero_vec;
        all_zg_nrs = zero_vec;
        all_cum_hrs_flt = zero_vec;
        all_max_hrs_flt = curr_max_ns_fl_hs;
        
        all_survs = one_vec;
        all_fins = zero_vec;
        all_arrs = zero_vec;
        all_losts = zero_vec;
        
        all_zug_signps = [];

        past_sgnp = zero_vec;
        is_new_sgnp = false_vec;

        all_disp_d_dep = [];
        all_disp_d_arr = [];
        all_coeffs = [];
        all_day_starts = [];
%         date_nrs = date_jul; 

        % for opt arr EU (whtrs)
%         lat_arr_EU = NaN*one_vec;
%         lon_arr_EU = NaN*one_vec;
        past_EU = false_vec;
%         arr_EU = false_vec;
        
        % march through flight steps until target Lat or maximum time
        % reached
        while sum(~finished) > 0 % 0.001*N_inds


            land_frac = 10*one_vec;  
            curr_dep_Lat = theta; %  cellfun(@(v)v(end),date_jul); %
            curr_dep_Lon = llamda;
            curr_over_wtr = over_water;
            curr_arr = arrived;
            curr_passed_Stops = passed_Stops;
            curr_zug_nr = zug_nr;
            curr_past_sgnp = past_sgnp;
            
            % alpha is current heading relative to geogr South
%             curr_head = alpha;
            
            % at this point the heading is determined
            % we now need to assess where and for how long to stop
         
            % now update next location
        
            % calculate candidate arrival location
            % may be adjusted below based on topography
            [theta(~finished), llamda(~finished)] = update_loc(...
                fl_hrs(~finished),hourly_del_Lat, ...
                curr_dep_Lat(~finished),curr_dep_Lon(~finished), ...
                alpha(~finished),include_wind,u_wind(~finished), ...
                v_wind(~finished),magn_star_night,decln(~finished), ...
                date_jul(~finished),geomagn_model);
            
            lat_deg_idx = 91 + floor(-180/pi*theta);
            lon_deg_idx = min(360,floor(180/pi*llamda) + 181);
            
            % remove polar cases for land or ocean compatibility
             lost = lost | (abs(theta) > lost_Lat*deg2rad);
             finished = finished | lost;                         

            % now update this candidate next location if lands at coast 
            
            % first check if candidate flight terminates over land

            over_water(~finished) = ~island(theta(~finished)*rad2deg, ...
                llamda(~finished)*rad2deg);
% try
            is_high_elev = hi_elevs(sub2ind([180 360],lat_deg_idx,lon_deg_idx));

    
            is_barrn = barrens(sub2ind(size(barrens),lat_deg_idx,lon_deg_idx));
% catch
%     keyboard
% end
            % next check options for stopping (current flight) at coast
            % If was over water at flight beginning, stop immediately
            stop_nxt_land(~finished) = curr_over_wtr(~finished) == 1 & ...
            stpCstWtrFlt == 1; % & cum_hrs_flt(~finished) > 2*max_fl_stp_hs
            stop_nxt_land(finished) = false;
            
            fracStop = stop_nxt_land + ~stop_nxt_land*fracStopCst;
            
%             if sum(stop_nxt_land > 0)
%                 keyboard
%             end

            thetas_fracs = ...
                curr_dep_Lat + (1:10).*(theta-curr_dep_Lat)/10;

            del_llamda(~finished) = absDiffDeg(llamda(~finished)*rad2deg,curr_dep_Lon(~finished)*rad2deg)*pi/180;
            del_llamda(~finished) = shiftAnglesFromMinus180To180(del_llamda(~finished)*rad2deg)*deg2rad;

            for ii = 1:10

                llamdas_fracs(~finished,ii) =  mod(curr_dep_Lon(~finished),2*pi) + del_llamda(~finished)*ii/10;

            end

            llamdas_fracs = shiftAnglesFromMinus180To180(llamdas_fracs*rad2deg)*deg2rad;

            % check if over land after full flight (for computational effeciency, 
            % we here assume th ebird cananticipate (see ahead) if will
            % land after full flight)

            % will land early if last idx inhospiible (with Cst option)
            % or if started over such (stop_nxt_land)
            will_land_early = ~finished & (stop_nxt_land | ...
                (over_water | is_high_elev | is_barrn));

            over_wtr_fracs = ~island(thetas_fracs*rad2deg, ...
                llamdas_fracs*rad2deg); %
    
            if sum(will_land_early) > 0 
%                 try
                lat_will_land_fr = thetas_fracs(will_land_early,:)*rad2deg;
                lon_will_land_fr = llamdas_fracs(will_land_early,:)*rad2deg;
%                 over_wtr_fracs
%                 over_wtr_fracs(will_land_early,:) = ~island(lat_will_land_fr, ...

                lat_fr_idx = 91 + floor(-lat_will_land_fr);
                lon_fr_idx = floor(lon_will_land_fr) + 181;
                is_high_elev_fr = hi_elevs(sub2ind([180 360],lat_fr_idx,lon_fr_idx));
                
                is_barrn_fr = barrens(sub2ind(size(barrens),lat_fr_idx,lon_fr_idx));
      
            end
%              artificially set water until last hour to choose last hour
%              when no land
%              over_wtr_fracs(~will_land_early,1:9) = 1;
%              over_wtr_fracs(~will_land_early,10) = 0;
             
            idx_not_fin = find(~finished);
            idx_wants_lnd = find(will_land_early);
                
             if stopCstOpt || sum(will_land_early) > 0
                                
                true_10th = [false*ones(1,9) true];
                
                lnd_fr_lands = land_frac(idx_wants_lnd);
                ovr_wtr_fr_lands = over_wtr_fracs(idx_wants_lnd,:);
                
%                 tic
%                 try
                   
                    % here we assume bird will see land below after
                    % user-input fraction of night but will not
                    % land at user-defined high elev or barren locs's

                parfor ii = 1:numel(idx_wants_lnd) %  
                    
                    ii_nf = idx_wants_lnd(ii);
%                     over_wtr_fracs(ii_nf,:) = ~landmask(thetas_fracs(ii_nf,:)*rad2deg,llamdas_fracs(ii_nf,:)*rad2deg,cst_accuracy);

%                     over_wtr_fracs(ii_nf,:) = land_or_ocean(thetas_fracs(ii_nf,:)*rad2deg,llamdas_fracs(ii_nf,:)*rad2deg,cst_accuracy);
                   
                        not_wtr_hi_brrn = ...
                            ~ovr_wtr_fr_lands(ii,fracStop(ii_nf):10) & ...
                            ~is_high_elev_fr(ii,fracStop(ii_nf):10) & ...
                            ~is_barrn_fr(ii,fracStop(ii_nf):10);
                        lnd_fr_lands(ii) = fracStop(ii_nf) + ...
                        find(not_wtr_hi_brrn | ...
                        true_10th(fracStop(ii_nf):10),1,'first') -1;
                    
%                     if land_frac(ii_nf,1) < 10
%                         
%                         keyboard
%                         
%                     end
                    
                end
                
%                 catch
%                     
%                     parfor_will_land
                    
%                 end
%                 toc
                
                land_frac(idx_wants_lnd) = lnd_fr_lands;
                
                fl_hrs = fl_hrs.*land_frac/10;
                
                lin_ld_frc = sub2ind(size(over_wtr_fracs), 1:size(over_wtr_fracs,1), land_frac');        

                theta(~finished) = thetas_fracs(lin_ld_frc(~finished))';
                llamda(~finished) = llamdas_fracs(lin_ld_frc(~finished))';    
                
%                 if (any(abs(llamda) == pi))
%                     keyboard
%                 end

             
                over_water(~finished) = over_wtr_fracs(lin_ld_frc(~finished))';
%                 over_water(~finished) = over_water(~finished) & land_frac(~finished) == 10;
                
                idx_over_wtr = find(over_water & ~finished);
                % if (still) not landed, stop at nearest coast from
                % within 200 km of last location;
                % keep flight hours same (for simplicity)

                if ~isempty(idx_over_wtr)
 
                    ows = over_water(idx_over_wtr);
%                     ow_10s = ones(size(idx_over_wtr));
                    
                    th_ws = theta(idx_over_wtr);
                    ll_ws = llamda(idx_over_wtr);   
                    
%                     th_10s = th_ws; % thetas_fracs(idx_over_wtr,10);
%                     ll_10s = ll_ws; % llamdas_fracs(idx_over_wtr,10); 
                    
%                     del_ll = del_llamda(idx_over_wtr);
%                     tic
                   %% 
                   parfor i_ow = 1:numel(idx_over_wtr) % par

                       % here we assume bird at water at dawn has linear
                       % decreasing chance of finding (any) land with
                       % increasing distance
                        
%                         [d_Stops1] = distance(th_ws(i_ow), ...
%                             ll_ws(i_ow),coastlat,coastlon,'radians'); % , ang_Stops

                        % replace by Cartesian (local planar) distance 
                        % since within bird's eye view spherical effects
                        % are neglible
                         d_Stops = sqrt((th_ws(i_ow)-coastlat).^2 + ...
                            (cos(th_ws(i_ow))*(ll_ws(i_ow) -coastlon)).^2); %
                        
                        % linear decrease in prob seeing coast ahead 
                        % from half the max dist
                        p_dtc = ...
                            min(2*max((1-d_Stops/max_rad_stop_dawn),0),1); 

%                         sees_cst = find(p_dtc >= rand(size(d_Stops)));
%                          sees_cst = p_dtc >= rand(size(d_Stops));
                        if sum(p_dtc) > 0
                        
                            nonz = find(p_dtc ~=0);
                            sees_cst = nonz(p_dtc(nonz) >= rand(size(d_Stops(nonz))));

                        else
                        
                            sees_cst = [];
                        
                        end

%                          p_dtc = zeros(size(d_Stops));
%                          p_dtc(d_Stops<max_rad_stop_dawn)

%                 try
                         if numel(sees_cst) > 0 
                             
                             ows(i_ow) = 0;
%                              ow_10s(i_ow) = 0;
                             idx_closest = sees_cst(d_Stops(sees_cst) == min(d_Stops(sees_cst)));

%                              t = find(sees_cst, n_closest(1), 'first');
%                              idx_closest = t(end);

%                              idx_closest = find(d_Stops == min(d_Stops),1,'first');
        %                                  alpha(n_fl_step+1) = ang_Stops(idx_closest)*deg2rad-pi;
                             % update locn to nearest coast
                             th_ws(i_ow) = coastlat(idx_closest(1));
                             ll_ws(i_ow) = coastlon(idx_closest(1));
%                              llamda(idx_over_wtr(i_ow)) = shiftAnglesFromMinus180To180(nxt_llamda_deg)*deg2rad;
%                              ll_10s(i_ow) = ll_ws(i_ow);
%                              th_10s(i_ow) = th_ws(i_ow);
                                                       
                            %    fl_hrs(idx_over_wtr(i_ow)) = -99;
                                                          
                         end

%                 catch
% 
%                     keyboard
% 
%                 end
                                                  
                    end
%                     toc
                     over_water(idx_over_wtr) = ows;
                     over_wtr_fracs(idx_over_wtr,10) = ows;
                     theta(idx_over_wtr) = th_ws;
                     llamda(idx_over_wtr) = shiftAnglesFromMinus180To180(ll_ws*rad2deg)*deg2rad;
                     idx_lands = sub2ind(size(thetas_fracs), ...
                         idx_over_wtr,10*ones(size(idx_over_wtr)));
                     thetas_fracs(idx_over_wtr,10) = th_ws;
                     llamdas_fracs(idx_over_wtr,10) = llamda(idx_over_wtr);
                     
                   % re-calculate change in longitude
                    del_ll = ...
                        absDiffDeg(llamda(idx_over_wtr)*rad2deg,curr_dep_Lon(idx_over_wtr)*rad2deg)*pi/180;
                    del_llamda(idx_over_wtr) = shiftAnglesFromMinus180To180(del_ll*rad2deg)*deg2rad;
                     
                end
                
             else
                 
                 % fractions of dateline crossers all set to "10"
%                  lin_ld_frc = sub2ind(size(over_wtr_fracs), 1:size(over_wtr_fracs,1), 10*one_vec');   
                  land_frac = 10*one_vec;
                  
                end

                       
            location.longitude = llamda*rad2deg; 
            location.latitude = theta*rad2deg; 
            
%             if sum(fl_hrs == 1) > 0
%                 keyboard
%             end            
            
            update_arrive_survive_lost_finished
                    
            
            % update indices for ndvi and elev
            lat_deg_idx = 91 + floor(-180/pi*theta);
            lon_deg_idx = min(360,floor(180/pi*llamda) + 181);
%             lat_deg_idx = round(180/pi*theta + 90.5 - eps_grid_fact);
%             lon_deg_idx = round(180/pi*llamda + 180.5 - eps_grid_fact);
            
%             try
%             over wtr updated above
            is_high_elev = hi_elevs(sub2ind([180 360],lat_deg_idx,lon_deg_idx));
            is_barrn = ~over_water & ...
                barrens(sub2ind(size(barrens),lat_deg_idx,lon_deg_idx));  
            
            % March 2022 add no replensihing max "water" flight hours
            % when stopped at high elevation or barren location
%             no_replenish = over_water(~finished) | is_high_elev(~finished) | ...
%                 is_barrn(~finished);
            
%             catch
%                 keyboard
%             end

            cum_hrs_flt(~finished) = cum_hrs_flt(~finished) + fl_hrs(~finished);% 
            
%             cum_hrs_flt(~finished) = cum_hrs_flt(~finished) + ...
%                 no_replenish.*fl_hrs(~finished);% 
            
%             cum_hrs_flt(~finished) = ...
%                 no_replenish.*cum_hrs_flt(~finished) +  ...
%                 over_water(~finished).*fl_hrs(~finished);% 
            
%             if any(~over_water(~finished) & ~island(theta(~finished)*rad2deg, ...
%                 llamda(~finished)*rad2deg))
%             
%                 keyboard
%             
%             end
                                         
            % next location now determined.
            % now, update date for next flight step 
            % If over land, this occurs either the following night 
            % or after stopover according to schedule 
            % (e.g., 3 flight days followed by 8 stopover days) 
            % Over water, occurs in remainder of
            % calendar day or next night as appropriate               
            % We also check case at ocean barrier below
            

%             try

%             catch
%                 keyboard
%             end


        % will not make extended stop over on water, high elev or barren land
        is_poor_stop = poor_stops(sub2ind(size(poor_stops),lat_deg_idx,lon_deg_idx));
        skip_stop_over(~finished) = over_water(~finished) | is_high_elev(~finished) | ...
            is_poor_stop(~finished);
        
%         wants_continue_fly = (cumul_flts(~finished) <= n_flts_seq -1);
        % Apr 2022 added possibility to tank after 1/2 of max fl hrs used
        % up. Also new option stpovCstWtrFlt to stop over
        % after water (helps after Mediterranean before
        % Sahara)
        wants_continue_fly =  ...
            (curr_max_ns_fl_hs(~finished) - cum_hrs_flt(~finished) > ...
            min_n_fls_stp*min_fl_stp_hs) & ... 
            ~is_new_sgnp(~finished) & ...
            ~(stop_nxt_land(~finished) & stpovCstWtrFlt); 
        % (cumul_flts(~finished) <= n_flts_seq -1) &
%         /2 ...
%             & curr_max_ns_fl_hs(~finished)

%         

        % curr_max_ns_fl_hs(~finished)/2); <  (cum_hrs_flt(~finished) + min_ns_fl_hs
%         max(curr_max_ns_fl_hs(~finished)/2, max_ns_hs_bs(~finished)/2)
        
        % update alpha_base (zugkniks) based on  signpost (zug_opt == 1)

        if zug_opt == 1 

             update_zugs

        end

        % check if water ahead (use same orientn i.e. alpha for simplicity)
        % assume can look "halfway" ahead (radAheadDtrCst is upper limit)
%         try
       if ~isempty(idx_not_fin)
           [ahead_theta, ahead_llamda] = update_loc(...
            ones(size(idx_not_fin)),radAheadDtrCst/2,theta(~finished), ...
            llamda(~finished),alpha(~finished),0,zeros(size(idx_not_fin)),zeros(size(idx_not_fin)), ...
            star_cmp_flag,decln(~finished),date_jul(~finished),geomagn_model);
        
           ahead_wtr = island(ahead_theta*rad2deg, ahead_llamda*rad2deg) == 0;

       else

           ahead_wtr = [];

        end
%         catch
% 
%             keyboard
% 
%         end

        req_stpovr_poly = any(~curr_passed_Stops(~finished,:) ...
            & passed_Stops(~finished,:),2) & require_stopover_poly;
        req_stpovr_sgnp = is_new_sgnp(~finished) & require_stopover_zug;

        % determine whether (can and "programmed" to initiate extended stopover
        stops_over(~finished) = ~skip_stop_over(~finished) & (~wants_continue_fly | ...
            req_stpovr_poly | ahead_wtr | req_stpovr_sgnp);
     
        stops_over(finished) = false;
           
        idx_stopovs = find(stops_over);
      
        N_stop_overs = sum(stops_over);

% If not scaling to ndvi use these lines
%         curr_max_ns_fl_hs(stops_over) = ...
%         max(rand([N_stop_overs 1]).* ...
%        (max_ns_hs_bs(stops_over) - min_ns_hs_bs(stops_over)) + ...
%        min_ns_hs_bs(stops_over), ...
%         curr_max_ns_fl_hs(stops_over)-cum_hrs_flt(stops_over)); 
%        cum_hrs_flt(stops_over) = 0;

        % Version scaled to ndvi:
        % stopovers have ndvi > 1/2 thresh (otherwise are classified as
        % barren) and satisfy other thresh (veg, trees)
        % Further separate good > thresh from so-so (1/2 to thresh
        
%        if sum(is_new_sgnp(~finished) & ~stops_over(~finished)) > 0
%            keyboard
%        end

        % next scale fuel gain to stopover quality (ndvi)
        if N_stop_overs > 0
       
           stop_num(~finished) = stop_num(~finished) + stops_over(~finished);
           cumul_flts(~finished) = ~stops_over(~finished).*(cumul_flts(~finished) +1);
     
% try
            % distinguish between regular and signposted stopovers
            reg_stop = stops_over & ~is_new_sgnp;
            N_reg_stp = sum(reg_stop);
            sgnp_stop = stops_over & is_new_sgnp;

            ndvi_stop_regs = ...
             ndvi(sub2ind([180 360],lat_deg_idx(reg_stop), ...
             lon_deg_idx(reg_stop)));  

             ndvi_stop_sgnps = ...
             ndvi(sub2ind([180 360],lat_deg_idx(sgnp_stop), ...
             lon_deg_idx(sgnp_stop))); 

             veg_stop_regs = ...
             veg_stp(sub2ind(size(veg_stp),lat_deg_idx(reg_stop), ...
             lon_deg_idx(reg_stop)));  

             veg_stop_sgnps = ...
             veg_stp(sub2ind(size(veg_stp),lat_deg_idx(sgnp_stop), ...
             lon_deg_idx(sgnp_stop))); 

% catch
% 
%     keyboard
% 
% end

%            good_stops = ndvi_stops >= thrsh_ndvi_stop_max;
%            so_so_stops =  ndvi_stops >= thrsh_ndvi_stop_min; % ~good_stops;

            frac_ndvi_reg = min((ndvi_stop_regs - thrsh_ndvi_stop_min)/ndvi_rng,1); 
            frac_ndvi_sgnp = min((ndvi_stop_sgnps - thrsh_ndvi_stop_min)/ndvi_rng,1); 
            frac_veg_reg = min((veg_stop_regs - min_any_veg_stp)/veg_rng,1); 
            frac_veg_sgnp = min((veg_stop_sgnps - min_any_veg_stp)/veg_rng,1); 
            frac_DFL_reg = frac_ndvi_reg; %.*frac_veg_reg; %  
            frac_DFL_sgnp = frac_ndvi_sgnp; % .*frac_veg_sgnp; % 

%             curr_max_ns_fl_hs(stops_over) = ...
%             max(frac_DFL.*rand([N_stop_overs 1]).* ...
%            (max_ns_hs_bs(stops_over) - min_ns_hs_bs(stops_over)) + ...
%            min_ns_hs_bs(stops_over), ...
%             curr_max_ns_fl_hs(stops_over)-cum_hrs_flt(stops_over)); 

           curr_max_ns_fl_hs(reg_stop) = ...
            max(frac_DFL_reg.*rand([N_reg_stp 1]).* ...
           (max_ns_hs_bs(reg_stop) - min_ns_hs_bs(reg_stop)) + ...
           min_ns_hs_bs(reg_stop), ...
            curr_max_ns_fl_hs(reg_stop)-cum_hrs_flt(reg_stop)); 

           curr_max_ns_fl_hs(sgnp_stop) = ...
            max(frac_DFL_sgnp.*rand([N_stop_overs-N_reg_stp 1]).* ...
           (max_sgnp_hs_bs(sgnp_stop) - min_sgnp_hs_bs(sgnp_stop)) + ...
           min_sgnp_hs_bs(sgnp_stop), ...
            curr_max_ns_fl_hs(sgnp_stop)-cum_hrs_flt(sgnp_stop)); 

           cum_hrs_flt(stops_over) = 0;

        end
        
        % now that locn and state (survived, stopped, arrived..) 
        % are updated, update new date, time lag (for TC Sun comp),
        % sun angle and flight hours
     
        % calculate change in time due to longitudinal movement 
        % (this is separate from the time lag in TC Sun compass
        % and relates to keeping track of the julian day for any sun
        % compass and also for the day of year afetr longer over water
        % flights (fof migration duration) 

        % only add time when no longer over water i.e., during 
        % transient or extended stopover on land
        % del_ll_last accounts for time zone change (in solar time!)
        
        del_ll_last_land(~finished,1) = ...
            absDiffDeg(llamda(~finished)*rad2deg, ...
            last_lon_land(~finished)*rad2deg)*pi/180;       
        del_ll_last_land(~finished) = ...
            shiftAnglesFromMinus180To180(del_ll_last_land(~finished)*rad2deg)*deg2rad;     

        if tc_comp_wtr == 1
            
            d_time_zone = ~finished.*(cum_hrs_flt/24 + ...
                0.5*del_ll_last_land/pi);

            d_time = sign(d_time_zone).*floor(abs(d_time_zone))  + ...
                ~finished.*(~stops_over + reg_stop*stop_durn + sgnp_stop*stop_durn_sp);
        
        else
            
            d_time_zone(~over_water,1) = ~finished(~over_water).*(cum_hrs_flt(~over_water)/24 + ...
                0.5*del_ll_last_land(~over_water)/pi);
            d_time_zone(over_water) = 0;
            d_time = sign(d_time_zone).*floor(abs(d_time_zone))  + ...
                ~finished.*(~stops_over + reg_stop*stop_durn + sgnp_stop*stop_durn_sp);

        end

        last_lon_land(~over_water) = llamda(~over_water);

%             end

% calculate previous day's sun az for 'fixing' heading
% on 1st day of stopover with TCSC

% Apr 2022 add barren for longer flights
        sun_az_1st_nt(~finished) = calc_sun_az_fl_hs( ...
           theta(~finished),doys(~finished)+1,endur_flt, ...
            is_barrn(~finished) | over_water(~finished), ...
            min_init_fl_stp_hs,max_fl_stp_hs,offset_dep_dusk);

        doys(~finished) = doys(~finished) + d_time(~finished);
        date_jul(~finished) = date_jul(~finished)  + d_time(~finished); % .*next_jul_dt(~finished); % 

%         fly_on = ~finished & (over_water | endur_flt);

        % then, calculate flight hours next time step
%             curr_date.day = date_jul;
% try
       [sun_az(~finished), fl_hrs(~finished)] = calc_sun_az_fl_hs( ...
           theta(~finished),doys(~finished),endur_flt, ...
            is_barrn(~finished) | over_water(~finished), ...
            min_init_fl_stp_hs,max_fl_stp_hs,offset_dep_dusk);  

       hs_left = floor(curr_max_ns_fl_hs(~finished) - cum_hrs_flt(~finished));
       fl_hrs(~finished) = max(min(fl_hrs(~finished), ...
           hs_left/2),1);
        
        % catch
%     keyboard
% end
% work out flight hours and sun azimuth for next step (top of loop)

   % option to reset sun compass for those stopping over
   % (avoids big jumps in sun compass heading)
           
    % first account for change d_alpha_sun to migrants' sun compass through
    % optional (sun_head_reset == 2) transfer on arrval at extended stopover  
    % from their nocturnal (geographic / magnetic) compass direction 
    % 
    % all "previous" time lag will have affected the current heading
    % so are also accounted for (next section removes lag under  temporal
    % reset option
            if calibr_comp(2)~=0 ...
                   && (sun_head_reset == 3 || ... 
                   (sun_head_reset == 2 && sum(stops_over) > 0))

                % correct this later
                keyboard

               d_alpha_sun(stops_over) = absDiffDeg(alpha(stops_over)*rad2deg, ...
                   (sun_az_1st_nt(stops_over) + alpha_bases{2}(stops_over) + ...
                   d_time_lag(stops_over))*rad2deg)*deg2rad;

               d_alpha_sun = ...
                     mod(d_alpha_sun + pi, 2*pi) -pi;
                            
            end
           
          % compute time shift from 'jet lag' for TC sun comp
          % since last stopover (or initial departure)

        del_ll_ref_lon_tc(~finished,1) = ...
            absDiffDeg(llamda(~finished)*rad2deg, ...
            ref_lon_tc(~finished)*rad2deg)*pi/180;    

        del_ll_ref_lon_tc(~finished) = ...
            shiftAnglesFromMinus180To180(del_ll_ref_lon_tc(~finished)*rad2deg)*deg2rad;

        d_time_lag = (~finished).*(calibr_comp(3) == 2 ...
             | calibr_comp(2) == 2).* ...
             del_ll_ref_lon_tc.*sin(theta); 
                 

        % then account for any loss of "time lag" from last extended stopover or
        % initial departure location
          if calibr_comp(2) == 2 ...
               && (tc_comp_reset==1 && sum(stops_over) > 0)

%               adj_alph = stops_over; % & ~past_sgnp; %

%      update the offset to the original sun compass (offset to azimuth)
%       heading such that the preferred TCSC heading remains the same as
%       on the 1st night of stopover, accounting for the just updated 
%       TCSC clock shift (which is subsequently set to yero as bird resets
%       its inner clock
%                    
               d_alpha_sun(stops_over) = d_alpha_sun(stops_over) + ...
                    sun_az_1st_nt(stops_over) - sun_az(stops_over) + ...
                   d_time_lag(stops_over) -  del_alph_dtr(stops_over);

               d_alpha_sun = ...
                     mod(d_alpha_sun + pi, 2*pi) -pi;
                 
              d_time_lag(stops_over) = 0; 
              ref_lon_tc(stops_over) = llamda(stops_over);

          end

           % check for reoriented flight over water: if so, we
           % increase the "rand_fact" uncertainty in orientation for below
               
           reorients = reOrWtrOpt == 1 & over_water; % cum_hrs_flt > 12; %  & cum_hrs_flt < 24;
           
           % double error if reorienting
           rand_fact = 1 + rand_wtr_fact*(over_water);
           
           % determine random orientation component (we increased
            % uncertainty when reoriented over water)
            
            % this will be sum of all averaged compass types                
                
%             d_alpha_rnd = (~finished).*zero_vec;            
%             
            for i_comp = 1:3
                if (calibr_comp(i_comp) > 0 || endog_comp == i_comp) ...
                        && cal_comp_nonzero(i_comp)
                     d_alpha_errs{i_comp} = ...
                        rand_fact.*(~finished).*vmrand(0, ...
                        compass_cal_kappa(i_comp), [N_inds 1]);
                % add magn varn to sun compass if refl about mag S axis
                elseif i_comp == 1 && refl_sun_mag && ...
                        calibr_comp(2) > 0 && cal_comp_nonzero(i_comp) > 0 
                        d_alpha_errs{i_comp} = vmrand(0, compass_cal_kappa(i_comp), [N_inds 1]);
                else
                    d_alpha_errs{i_comp} =  zero_vec;
                end
            end
%             
%              d_alpha_rnd = rand_fact.*d_alpha_rnd/sum(cal_comp_nonzero);
            
            % rand_fact.*compass_std.*randn(N_inds,1);
            
           % now compute relevant magnetic parameter at present location
            % (here vertical Intensity ~ sin(Lat) and Inclination angle ~ atan(tan(2*Lat))
            % which can be used as steering variable directly (calibr_comp(1) = 1) or 
            % indirectly (calibr_comp(1) = 2, so-called magnetoclinic compass)
            % Option 3 is novel parallel magnetoclinic compass (only works
            % well across longitudes as opposed to magnetoclinic across
            % latitudes).
            % Note for calibr_comp(1) == 0 we keep magn parameter constant 
            % (==1) i.e. no magn steering.
            
            prev_decln = decln;
            
          if sum(~finished) > 0 
                     
                if geomagn_model == 1

                    clear Bx By Bz
                    d_alpha_incl = vmrand(0, compass_kappa_incl, [sum(~finished) 1]);
                      udy = round(mean(date_jul(~finished))); % unique(date_jul(~finished));
%                       dy_no_fin = date_jul(~finished);
                      Lat_no_fin = location.latitude(~finished);
                      Lon_no_fin = location.longitude(~finished);
                      
                      try
                          
                     [Bx, By, Bz] = igrf(udy, ...
                               Lat_no_fin, Lon_no_fin, 0);
%                       for idy = 1:numel(udy)
%                            dy_i = udy(idy);
%                            is_dyi =  dy_no_fin == dy_i;
%                            [Bx(is_dyi,1), By(is_dyi,1), Bz(is_dyi,1)] = igrf(dy_i, ...
%                                Lat_no_fin(is_dyi,1), Lon_no_fin(is_dyi,1), 0);
%                       end

                      catch
                          
                          keyboard
                          
                      end

                    decln(~finished,1) = atan2(By, Bx); % pi/2-
                    Bhz = hypot(Bx,By);
                    tot_B = hypot(Bhz,Bz);
                    incln(~finished,1) = min(atan(Bz./Bhz) + d_alpha_incl,pi/2);
                    vt_inten(~finished,1) = Bz.*(1 + randn(sum(~finished),1)*std_rel_err_inten); % ./tot_B; % ; % 
                    hz_inten(~finished,1) = Bhz.*(1 + randn(sum(~finished),1)*std_rel_err_inten); % ./tot_B; % ; % 

                    tan_incl = tan(incln(~finished,1));
                    tot_inten(~finished,1) = tot_B.*(1 + randn(sum(~finished),1)*std_rel_err_inten);

                    if is_magcl
                        magn_param(~finished) = tan_incl; %
                    else
                        magn_param_0 = NaN*one_vec;
                    end
        
                    
                else

                    % Also, we use that for dipole, tan(inclin angle) = 2*tan(theta)
                    tan_incl = 2*tan(theta + d_alpha_incl);
                    magn_param(~finished) = tan_incl;

                end

                % any steering will be an adjustment to relative value of 
                % initial magnetic paraneter 

                if calibr_comp(1) == 2      % magnetoclinic              

                   if magcl_opt == 1
                       
                        ok_magncl = abs(tan_incl./magn_param_0(~finished)) < 1;
                        ok_incl(~finished,1) = ok_magncl;
                        ok_incl(finished,1) = false;
                        ok_n_fin = ok_incl & ~finished;
                        % d_alpha_magn maintains sign of magnetic parameter 
                        % magn_param_0 for correct E/W orientation
                        d_alpha_magn(ok_n_fin) = asin(tan_incl(ok_magncl)./magn_param_0(ok_n_fin)) ...
                            - alpha_bases{1}(ok_n_fin); %...
                        d_alpha_magn(~ok_incl & ~finished) = ...
                            -alpha_bases{1}(~ok_incl & ~finished) + pi/2*sign_inher_heads(~ok_incl & ~finished); 
        %                     min(abs(asin(tan_incl/magn_param_0) - alpha(1)),pi); 
        
                        if magcl_reset_opt == 1 && sum((~ok_incl)) > 0

                            magn_param_0(~ok_incl & ~finished) = ...
                                tan_incl(~ok_magncl)./sin(inher_heads(~ok_incl & ~finished));

                        end
                    
                   elseif magcl_opt == 2
                       
                        ok_magncl = abs(magn_param_0(~finished)./tan_incl) < 1;
                        ok_incl(~finished,1) = ok_magncl;
                        ok_incl(finished,1) = false;
                        ok_n_fin = ok_incl & ~finished;
                        % d_alpha_magn maintains sign of magnetic parameter 
                        % magn_param_0 for correct E/W orientation
                        d_alpha_magn(ok_n_fin) = asin(magn_param_0(ok_n_fin)./tan_incl(ok_magncl)) ...
                            - alpha_bases{1}(ok_n_fin); %...
                        d_alpha_magn(~ok_incl & ~finished) = ...
                            -alpha_bases{1}(~ok_incl & ~finished) + pi/2*sign_inher_heads(~ok_incl & ~finished); 
        %                     min(abs(asin(tan_incl/magn_param_0) - alpha(1)),pi); 
                       
                        if magcl_reset_opt == 1 && sum((~ok_incl)) > 0

                            magn_param_0(~ok_incl & ~finished) = ...
                                tan_incl(~ok_magncl).*sin(inher_heads(~ok_incl & ~finished));

                        end
                        
                   else % magcl_opt == 3 to implement
                       
                       
                   end


                elseif calibr_comp(1) == 3

                  if magcl_opt == 1
                                           
                        ok_magncl = abs(magn_param(~finished)./magn_param_0(~finished)) <1;                    
                        ok_incl(~finished,1) = ok_magncl;
                        ok_incl(finished,1) = false;
                        ok_n_fin = ok_incl & ~finished;    
                        d_alpha_magn(ok_n_fin) = ....
                            sign_inher_heads(ok_n_fin).*acos(magn_param(ok_n_fin)./ ...
                            magn_param_0(ok_n_fin)) - alpha_bases{1}(ok_n_fin); % ...

                        d_alpha_magn(~ok_incl & ~finished) = -alpha_bases{1}(~ok_incl & ~finished);
                        % + pi/4*sign_inher_heads(~ok_incl); % - pi/2*sign_inher_heads(~ok_incl); 
        %                     - min(acos(magn_param_0/magn_param(n_fl_step+1)) - abs(alpha(1)),pi);
 
                        if magcl_reset_opt == 1 && sum((~ok_incl)) > 0

                              magn_param_0(~ok_incl & ~finished) = ...
                                tan_incl(~ok_magncl)./cos(inher_heads(~ok_incl & ~finished));

                        end        
                        
                  elseif magcl_opt == 2
                        
                      ok_magncl = abs(magn_param_0(~finished)./magn_param(~finished)) <1;                    
                        ok_incl(~finished,1) = ok_magncl;
                        ok_incl(finished,1) = false;
                        ok_n_fin = ok_incl & ~finished;    

                        d_alpha_magn(ok_n_fin) = ....
                            sign_inher_heads(ok_n_fin).*acos(magn_param_0(ok_n_fin)./magn_param(ok_n_fin)) ...
                            - alpha_bases{1}(ok_n_fin); % ...
                        d_alpha_magn(~ok_incl & ~finished) = -alpha_bases{1}(~ok_incl & ~finished); 
                        % + pi/4*sign_inher_heads(~ok_incl); % - pi/2*sign_inher_heads(~ok_incl); 
        %                     - min(acos(magn_param_0/magn_param(n_fl_step+1)) - abs(alpha(1)),pi);
 
                        if magcl_reset_opt == 1 && sum((~ok_incl)) > 0

                              magn_param_0(~ok_incl & ~finished) = ...
                                tan_incl(~ok_magncl).*cos(inher_heads(~ok_incl & ~finished));

                        end   
                        
                  end
                  
                else 
                        
                        d_alpha_magn = 0;      
                
                end

                % base geo shift always zero (exclding calibr error)
%                 d_alpha_geo = (~finished).*zero_vec;

                alpha_perfects{1} = alpha_bases{1} + d_alpha_magn + decln;
                alpha_perfects{2} = alpha_bases{2} + sun_az + ...
                    d_alpha_sun + (calibr_comp(3) ~= 2)*d_time_lag;
                alpha_perfects{3} = alpha_bases{3} + ...
                    (calibr_comp(3) == 2)*d_time_lag;

                alpha_all = [];

                for i_comp = 1:3
                    if  calibr_comp(i_comp) > 0 
                       alpha_all = [alpha_all alpha_perfects{i_comp}+d_alpha_errs{i_comp}] ;
                    end
                end
%             ...
%                             && cal_comp_nonzero(i_comp)
                 % can add compensation for prev detour (when on land)
                 alpha_mean = circ_mean(alpha_all')';

               % define new (non-detoured) headings
                canCmpDtr = dtrReComp*(~over_water & ~reorients);
                
                alpha(~reorients & ~finished) = mod(alpha_mean(~reorients & ~finished)  ...
                   - reOrReComp.*((n_toCmpReOr(~reorients & ~finished)) > 0).*mn_ReOr(~reorients & ~finished) + ...
                   - canCmpDtr(~reorients & ~finished).*(toCmpDtr(~reorients & ~finished) > 0).*mn_Dtr(~reorients & ~finished) + ...
                       +pi,2*pi)-pi;  
                   
%                
%             end
                
               n_toCmpReOr = max(n_toCmpReOr +reorients - ~reorients,0);
               mn_ReOr = (n_toCmpReOr>0).*(mn_ReOr.*(n_toCmpReOr) + reorients.*alpha) ...
                   ./(max(n_toCmpReOr,1));

               toCmpDtr = max(toCmpDtr - canCmpDtr,0);
               mn_Dtr = (toCmpDtr>0).*mn_Dtr;

               % reorient heading where applicable
               alpha(reorients & ~finished) = mod( ...
                   (reOrWtrStr==1).*(-sign(alpha_mean(reorients & ~finished))*pi/2) + ...
                   (reOrWtrStr==2).*(-alpha_mean(reorients & ~finished)) + ...n_fl_step
                   (reOrWtrStr==3).*(-alpha_lst_lnd(reorients & ~finished)) + ...
                   (reOrWtrStr==4).*(alpha_lst_lnd(reorients & ~finished)+pi/2) + ...
                   (reOrWtrStr==5).*(alpha_lst_lnd(reorients & ~finished)-pi/2) + ...
                   (reOrWtrStr==6).*(alpha_lst_lnd(reorients & ~finished) ...
                        -sign(alpha_lst_lnd(reorients & ~finished))*pi/2) + ...
                   (reOrWtrStr==7).*(alpha_lst_lnd(reorients & ~finished) ... 
                    -sign(alpha_0(reorients & ~finished))*pi/2) + ...
                    +pi,2*pi)-pi; % d_alpha_tot(reorients & ~finished)

               % update last head over land if non-recompensating
                  alpha_lst_lnd = (over_water | recmp_dtr | del_alph_dtr ~= 0).*alpha_lst_lnd ...
                      + ~(over_water | recmp_dtr | del_alph_dtr ~= 0).*alpha;
               
                % next, check and adjust heading and arrival 
                % for any detoured coastal cases

                % this may depend on winds; update for first day; this may be
                % re-updated below if bird waits for TWs

                [u_wind(~finished), v_wind(~finished)] = ...
                det_stoch_wind(theta(~finished),include_wind,wind_str,u_wind(~finished), ...
                v_wind(~finished),switch_wind_opt); % 

            if include_wind == 1
                tail_w(~finished) = calc_wind_supprt(alpha(~finished),u_wind(~finished), v_wind(~finished));
            end
               % work out detours for those who did not begin but will be over
               % water 

               if dtrCstOpt==1
                   
                   d2cst = dist2coast(theta*rad2deg,llamda*rad2deg);
                   can_dtr = ~finished & ~over_water & ...
                       (dtrCstOpt==1)*(d2cst<kmsAheadDtrCst) & ...
                       cum_steps_dtr < maxDtrSteps;
                   
               else
                   
                   can_dtr = false_vec;
                   
               end

               if n_fl_step(1) >= 1 && ...
                       sum(can_dtr) > 0 || (~endur_flt && dep_wtr_TWs && sum(~over_water) >0)

                   
                    determine_detours_barrier_crossings % determine_detours_barrier_parl
                    
               else
                   
                   del_alph_dtr = zero_vec;

               end

               recmp_dtr = reorients | cum_steps_dtr  > 0;

               n_fl_step = n_fl_step + ~finished;

               arr_EU = ((llamda >-13*pi/180 & theta>35*pi/180 & ...
                   theta<66.6*pi/180) | ...
                 (llamda >-20*pi/180 & theta<35*pi/180 & theta>0 ))  & ...
                  ~over_water(:) & ~past_EU;

               idx_arr_EU = find(arr_EU);

              all_lat_arr_EU(idx_arr_EU,iYear) = theta(arr_EU)*180/pi;
              all_lon_arr_EU(idx_arr_EU,iYear) = llamda(arr_EU)*180/pi;
              past_EU = past_EU | arr_EU;
             
          end


          % if last year, annotate fields for post analysis
          if iYear == N_mig_Sims % - 2 % 1 % 
              
            annotate_fields
 
          end
      

        end 
        
        
        % First, needs to have arrived and any obligatory zugknicks
        stoppedAndArrived = arrived & zug_nr >= min(Sim_nr,n_zg_i);
 
        if Sim_nr <=  N_init_Sims && zgknk_req_stop_warm
           % For half of the warm-up, 
           % enforce migratory connectivity between 
           % departure and arrival areas                  
            if Sim_nr <=  N_init_Sims/2  %   Inf %  

               for iArr = 1:nArr_polys

                   is_iA = idx_arrs == iArr;
                   stoppedAndArrived(is_iA) = stoppedAndArrived(is_iA) & ...
                        inpolygon(llamda(is_iA),theta(is_iA,1), ...
                             Lon_Arr_Subpoly{iArr},Lat_Arr_Subpoly{iArr});

               end

            else % afterwards, allow arrival to neighbouring arrival zones 
                % note this can lead to drift away from zones

               for iArr = 1:nArr_polys 

                   is_iA = idx_arrs == iArr;
                   stoppedAndArrived(is_iA) = stoppedAndArrived(is_iA) & ...
                        inpolygon(llamda(is_iA),theta(is_iA,1), ...
                             Lon_nbr_subP{iArr},Lat_nbr_subP{iArr});
    %               mn_lat_Deps(iA,iYear) = circ_mean(theta_0(is_iA))*180/pi;
    %               mn_lon_Deps(iA,iYear) = circ_mean(llamda_0(is_iA))*180/pi;                    
               end          

            end
        end
                            
        % only successful if arrived and stopped in obligatory staging
        % region

        if lst_req_stp > 0
            
            stoppedAndArrived = stoppedAndArrived & ...
                 all(passed_Stops(:,1:lst_req_stp),2);
                        
        end

        if zug_opt == 1 % 

            if zg_rq_stp

                for iz = 1:min(Sim_nr,n_zugs) %  
                    
                    is_zug_i = zug_nr >= iz;
                                       
                    if sum(is_zug_i) > 0 
                        stoppedAndArrived(is_zug_i) = ...
                            stoppedAndArrived(is_zug_i)  & ...
                         inpolygon(lon_bs_zugs(is_zug_i,1), ...
                         lat_bs_zugs(is_zug_i,1), ...
                         Lon_stop_poly{1},Lat_stop_poly{1}); 
                    end
                                        
                end
                
            end

        end
                                    
        N_succ = sum(stoppedAndArrived);
       
       %% accumulate yearly stats for post analysis
       successful(iYear) = N_succ/N_inds; 
       mn_succ_hds(iYear) = circ_mean(inher_heads(stoppedAndArrived));
       std_succ_hds(iYear) = circ_std(inher_heads(stoppedAndArrived));
       all_ys_stoppedAndArrived(:,iYear) = stoppedAndArrived;
       all_ys_surv_fuel(:,iYear) = surv_fuel;
       all_ys_surv_wtr(:,iYear) = surv_wtr;
       all_ys_surv_barrn_elev(:,iYear) = surv_barrn & surv_elev;
       all_ys_lost(:,iYear) = lost;
       all_ys_late(:,iYear) = is_late;

%        try
       all_ys_init_geo_heads(:,iYear) = all_geo_heads(:,1);
       all_ys_inher_heads(:,iYear) = inher_heads*180/pi + 180;
       
%        catch
%            keyboard
%        end
       % First, summarize Zug and Dep and Arr Poly
       for iDep = 1:nDep_polys 
           
           is_iDep = idx_deps == iDep;
           is_succ_iDep = is_iDep  & stoppedAndArrived;
           succs_DepPs(iDep,iYear) = sum(is_succ_iDep)/N_inds_DepPs(iDep);
                       
           % mean and std arrival per arr poly (indep of Arr poly defns)
            mn_Arr_lat_DepPs(iDep,iYear) = circ_mean(theta(is_succ_iDep))*180/pi;
            mn_Arr_lon_DepPs(iDep,iYear) = circ_mean(llamda(is_succ_iDep))*180/pi; 
            std_Arr_lat_DepPs(iDep,iYear) = circ_std(theta(is_succ_iDep))*180/pi;
            std_Arr_lon_DepPs(iDep,iYear) = circ_std(llamda(is_succ_iDep))*180/pi; 
            % mean dep locs per Dep poly
            mn_Dep_lat_DepPs(iDep,iYear) = circ_mean(theta_0(is_succ_iDep))*180/pi;
            mn_Dep_lon_DepPs(iDep,iYear) = circ_mean(llamda_0(is_succ_iDep))*180/pi; 
            std_Dep_lat_DepPs(iDep,iYear) = circ_std(theta_0(is_succ_iDep))*180/pi;
            std_Dep_lon_DepPs(iDep,iYear) = circ_std(llamda_0(is_succ_iDep))*180/pi;  

            % succ heads and signps per Dep 
            mn_succ_hd_DepPs(iDep,iYear) = circ_mean(inher_heads(is_succ_iDep))*180/pi;
             std_succ_hd_DepPs(iDep,iYear) = circ_std(inher_heads(is_succ_iDep))*180/pi;
            if zug_opt

                mn_zk_hd_DepPs(iDep,iYear) = circ_mean(zugkn_inher_heads(is_succ_iDep))*180/pi;
                std_zk_hd_DepPs(iDep,iYear) = circ_std(zugkn_inher_heads(is_succ_iDep))*180/pi; 
                if zug_signp <= 2
                    mn_zk_signp_DepPs(iDep,iYear) = circ_mean(zug_signps(is_succ_iDep))*180/pi;
                    std_zk_signp_DepPs(iDep,iYear) = circ_std(zug_signps(is_succ_iDep))*180/pi;
                else
                    mn_zk_signp_DepPs(iDep,iYear) = mean(zug_signps(is_succ_iDep));
                    std_zk_signp_DepPs(iDep,iYear) = std(zug_signps(is_succ_iDep));
                end

            end

            % track changes in magn field
            mn_dcl_DepPs(iDep,iYear) = circ_mean(decln_0(is_iDep))*180/pi;
            mn_incl_DepPs(iDep,iYear) = circ_mean(incln_0(is_iDep))*180/pi;
            mn_intn_DepPs(iDep,iYear) = mean(tot_inten_0(is_iDep));
            

            % mean and std zg per Dep polystd_lat_EU_DepPs
            for iz = 1:n_zugs
                
              is_succ_zk_DepPs = ~isnan(lat_bs_zugs(:,iz)) & is_succ_iDep;
              mn_lat_zk_DepPs(iz,iDep,iYear) = ...
                  circ_mean(lat_bs_zugs(is_succ_zk_DepPs,iz))*180/pi;
              mn_lon_zk_DepPs(iz,iDep,iYear) = ...
                  circ_mean(lon_bs_zugs(is_succ_zk_DepPs,iz))*180/pi;
              std_lat_zk_DepPs(iz,iDep,iYear) = ...
                  circ_std(lat_bs_zugs(is_succ_zk_DepPs,iz))*180/pi;
              std_lon_zk_DepPs(iz,iDep,iYear) = ...
                  circ_std(lon_bs_zugs(is_succ_zk_DepPs,iz))*180/pi;
              
            end

            if arr_EU_opt

                mn_lat_EU_DepPs(iDep,iYear) =  nanmean(all_lat_arr_EU(is_iDep,iYear));
                std_lat_EU_DepPs(iDep,iYear) = nanstd(all_lat_arr_EU(is_iDep,iYear));
                mn_lon_EU_DepPs(iDep,iYear) =  nanmean(all_lon_arr_EU(is_iDep,iYear));
                std_lon_EU_DepPs(iDep,iYear) = nanstd(all_lon_arr_EU(is_iDep,iYear));     

            end
                             
       end
           
       for jArr = 1:nArr_polys 
  
            % Next, summarize per Arrival Poly
            in_jArr = stoppedAndArrived & inpolygon(llamda,theta, ...
                                 Lon_Arr_Subpoly{jArr},Lat_Arr_Subpoly{jArr});

            succs_ArrPs(jArr,iYear) = sum(in_jArr)/N_inds;
            
            if sum(in_jArr) > 0
                
                % mean and std arrival per arr poly
                mn_Arr_lat_ArrPs(jArr,iYear) = circ_mean(theta(in_jArr))*180/pi;
                mn_Arr_lon_ArrPs(jArr,iYear) = circ_mean(llamda(in_jArr))*180/pi; 
                std_Arr_lat_ArrPs(jArr,iYear) = circ_std(theta(in_jArr))*180/pi;
                std_Arr_lon_ArrPs(jArr,iYear) = circ_std(llamda(in_jArr))*180/pi; 
                % mean dep locs per arrival poly
                % I.e., here indepednently of how the dep polys are
                % defined
                mn_Dep_lat_ArrPs(jArr,iYear) = circ_mean(theta_0(in_jArr))*180/pi;
                mn_Dep_lon_ArrPs(jArr,iYear) = circ_mean(llamda_0(in_jArr))*180/pi; 
                std_Dep_lat_ArrPs(jArr,iYear) = circ_std(theta_0(in_jArr))*180/pi;
                std_Dep_lon_ArrPs(jArr,iYear) = circ_std(llamda_0(in_jArr))*180/pi;    
                % mean and std zg per Dep poly
                for iz = 1:n_zugs
                  is_succ_zk_ArrPs = ~isnan(lat_bs_zugs(:,iz)) & in_jArr;
                  mn_lat_zk_ArrPs(iz,jArr,iYear) = ...
                      circ_mean(lat_bs_zugs(is_succ_zk_ArrPs,iz))*180/pi;
                  mn_lon_zk_ArrPs(iz,jArr,iYear) = ...
                      circ_mean(lon_bs_zugs(is_succ_zk_ArrPs,iz))*180/pi;
                  std_lat_zk_ArrPs(iz,jArr,iYear) = ...
                      circ_std(lat_bs_zugs(is_succ_zk_ArrPs,iz))*180/pi;
                  std_lon_zk_ArrPs(iz,jArr,iYear) = ...
                      circ_std(lon_bs_zugs(is_succ_zk_ArrPs,iz))*180/pi;
                end
            
            else
                
                mn_Arr_lat_ArrPs(jArr,iYear) = NaN;
                mn_Arr_lon_ArrPs(jArr,iYear) = NaN; 
                std_Arr_lat_ArrPs(jArr,iYear) = NaN;
                std_Arr_lon_ArrPs(jArr,iYear) = NaN; 
                mn_Dep_lat_ArrPs(jArr,iYear) = NaN;
                mn_Dep_lon_ArrPs(jArr,iYear) = NaN; 
                std_Dep_lat_ArrPs(jArr,iYear) = NaN;
                std_Dep_lon_ArrPs(jArr,iYear) = NaN;   
                mn_lat_zk_ArrPs(iz,jArr,iYear) = NaN;
                mn_lon_zk_ArrPs(iz,jArr,iYear) = NaN;
                std_lat_zk_ArrPs(iz,jArr,iYear) = NaN;
                std_lon_zk_ArrPs(iz,jArr,iYear) = NaN;
                  
            end
                
            % Lastly, summarize per Arr and (user-defined) Dep Poly               
            for iDep = 1:nDep_polys                  
               
                jArr_iDep = stoppedAndArrived & idx_deps == iDep &  ...
                    inpolygon(llamda,theta, ...
                             Lon_Arr_Subpoly{jArr},Lat_Arr_Subpoly{jArr});

                succs_DepArrPs(iDep,jArr,iYear) = sum(jArr_iDep)/N_inds_DepPs(iDep);                      

               if sum(jArr_iDep) > 0

                      % mean and std arrival per dep poly
                      mn_Arr_lat_DepArrPs(iDep,jArr,iYear) = circ_mean(theta(jArr_iDep))*180/pi;
                      mn_Arr_lon_DepArrPs(iDep,jArr,iYear) = circ_mean(llamda(jArr_iDep))*180/pi; 
                      std_Arr_lat_DepArrPs(iDep,jArr,iYear) = circ_std(theta(jArr_iDep))*180/pi;
                      std_Arr_lon_DepArrPs(iDep,jArr,iYear) = circ_std(llamda(jArr_iDep))*180/pi; 
                      % mean and std departure locs within Dep polys
                      mn_Dep_lat_DepArrPs(iDep,jArr,iYear) = circ_mean(theta_0(jArr_iDep))*180/pi;
                      mn_Dep_lon_DepArrPs(iDep,jArr,iYear) = circ_mean(llamda_0(jArr_iDep))*180/pi; 
                      std_Dep_lat_DepArrPs(iDep,jArr,iYear) = circ_std(theta_0(jArr_iDep))*180/pi;
                      std_Dep_lon_DepArrPs(iDep,jArr,iYear) = circ_std(llamda_0(jArr_iDep))*180/pi; 
                      % mean and std zg per Dep poly
                      for iz = 1:n_zugs
                          is_succ_zk_DepArrPs = ~isnan(lat_bs_zugs(:,iz)) & jArr_iDep;
                          mn_lat_zk_DepArrPs(iz,iDep,jArr,iYear) = ...
                              circ_mean(lat_bs_zugs(is_succ_zk_DepArrPs,iz))*180/pi;
                          mn_lon_zk_DepArrPs(iz,iDep,jArr,iYear) = ...
                              circ_mean(lon_bs_zugs(is_succ_zk_DepArrPs,iz))*180/pi;
                          std_lat_zk_DepArrPs(iz,iDep,jArr,iYear) = ...
                              circ_std(lat_bs_zugs(is_succ_zk_DepArrPs,iz))*180/pi;
                          std_lon_zk_DepArrPs(iz,iDep,jArr,iYear) = ...
                              circ_std(lon_bs_zugs(is_succ_zk_DepArrPs,iz))*180/pi;
                      end

                  else

                      mn_Arr_lat_DepArrPs(iDep,jArr,iYear) = NaN;
                      mn_Arr_lon_DepArrPs(iDep,jArr,iYear) = NaN; 
                      std_Arr_lat_DepArrPs(iDep,jArr,iYear) = NaN;
                      std_Arr_lon_DepArrPs(iDep,jArr,iYear) = NaN;
                      mn_Dep_lat_DepArrPs(iDep,jArr,iYear) = NaN;
                      mn_Dep_lon_DepArrPs(iDep,jArr,iYear) = NaN;
                      std_Dep_lat_DepArrPps(iDep,jArr,iYear) = NaN;
                      std_Dep_lon_DepArrPs(iDep,jArr,iYear) = NaN;
                      for iz = 1:n_zugs
                          mn_lat_zk_DepArrPs(iz,iDep,jArr,iYear) = NaN;
                          mn_lon_zk_DepArrPs(iz,iDep,jArr,iYear) = NaN;
                          std_lat_zk_DepArrPs(iz,iDep,jArr,iYear) = NaN;
                          std_lon_zk_DepArrPs(iz,iDep,jArr,iYear) = NaN;
                      end

                end
            
            end

      end
        
      for iz = 1:n_zugs
          
          mn_zk_hds(iz,iYear) = circ_mean(zugkn_inher_heads(stoppedAndArrived,iz));
          std_zk_hds(iz,iYear) = circ_std(zugkn_inher_heads(stoppedAndArrived,iz));

          is_succ_zks = ~isnan(lat_bs_zugs(:,iz)) & stoppedAndArrived;
          fr_zks(iz,iYear) = sum(is_succ_zks)/N_succ;

         if zug_opt == 1
             if zug_signp > 2
                mn_zk_signps(iz,iYear) = mean(zug_signps(stoppedAndArrived,iz)); 
                std_zk_signs(iz,iYear) = std(zug_signps(stoppedAndArrived,iz));
             else
                mn_zk_signps(iz,iYear) = circ_mean(zug_signps(stoppedAndArrived,iz));
                std_zk_signs(iz,iYear) = circ_std(zug_signps(stoppedAndArrived,iz));
             end
         end
         
      end   
           
      all_ys_lon_fin(:,iYear) = llamda*180/pi;
      all_ys_lat_fin(:,iYear) = theta*180/pi;

%% retain successful traits for next generation
      if N_succ > 0

             day_start_succ = day_starts(stoppedAndArrived);
             mig_dur = date_jul(stoppedAndArrived) - ...
                 date_jul_start(stoppedAndArrived);
             speed = (wt_inher_speed==1)*1./(mig_dur) ...
                 + (wt_inher_speed==0);
             succ_heads = inher_heads(stoppedAndArrived);
             lat_deps_succ = theta_0(stoppedAndArrived); % lat_bs_deps(stoppedAndArrived);
             lon_deps_succ = llamda_0(stoppedAndArrived); % lon_bs_deps(stoppedAndArrived);   
             lat_arrs_succ = theta(stoppedAndArrived)*rad2deg;
             lon_arrs_succ = llamda(stoppedAndArrived)*rad2deg;     
             
             succ_zug_inh_hds = zugkn_inher_heads(stoppedAndArrived,:);
%                  plot_lat_lon_succ_fail_itern
             coeff_succ = coeff_magn_steer(stoppedAndArrived,:); 

             succ_std_inher_head = std_inher_head(stoppedAndArrived,:); 
             mn_std_inh_hds(iYear) = mean(succ_std_inher_head); 
             std_std_inh_hds(iYear) = std(succ_std_inher_head); 

             succ_std_inher_steer = std_inher_steer(stoppedAndArrived,:); 
%                 coeff_succ = all_coeffs(stoppedAndArrived,:); 
            succ_zug_signps = zug_signps(stoppedAndArrived,:);                
             if zug_opt == 1
                 succ_std_inher_signp = std_inher_signp(stoppedAndArrived,:);
                 mn_std_inh_sps(iYear) = mean(succ_std_inher_signp); 
                 std_std_inh_sps(iYear) = std(succ_std_inher_signp); 
                 % only "track" 1st Zugknick for weighting probability
                 % of next-gen Zugknick signposts :-)
                 lat_zugs_succ = lat_bs_zugs(stoppedAndArrived,1);
                 lon_zugs_succ = lon_bs_zugs(stoppedAndArrived,1);
                 all_lat_zg(:,iYear) = lat_bs_zugs*180/pi;
                 all_lon_zg(:,iYear) = lon_bs_zugs*180/pi;
                 all_ys_zugkn_heads(:,iYear) = rad2deg*zugkn_inher_heads + 180;
                 all_ys_zugkn_sps(:,iYear) = zug_signps;
             end

             succ_disp_d_dep = disp_d_dep(stoppedAndArrived);  

             mn_disp_d_dep(iYear) = mean(succ_disp_d_dep)*mn_2_sig;
             std_disp_d_dep(iYear) = std(succ_disp_d_dep)*mn_2_sig;

             succ_disp_d_arr = disp_d_arr(stoppedAndArrived);
             succ_disp_d_zug = disp_d_zug(stoppedAndArrived);

             % same for (instrinsic / optimizing) std in inh of disp
             succ_std_inher_disp = std_inher_disp(stoppedAndArrived);  

             mn_std_inh_disp_km(iYear) = mean(succ_std_inher_disp)*mn_2_sig;
             std_std_inh_disp_km(iYear) = std(succ_std_inher_disp)*mn_2_sig;

            if opt_no_evo_std_yrs ~= 0  && ...
                 ((opt_no_evo_std_yrs ~= 3 && Sim_nr ==  N_wrmp_rnd) ...
                 || (opt_no_evo_std_yrs == 3 && Sim_nr == N_init_Sims)) %  N_init_Sims %
                 % make intrinsic std in inh hds and disp uniform within popn
                 % ie time scale of simulation not enough to change the genetics
                 % but inh can still vary according to this "evolved" variability
                
                 % = 0 keep optimizing with inherited std's 
                 % = 1 use mean std values; 
                 % now consider opt 2 just orientn variable
                 % option 3 none are variable (post warmup)

                std_inh_hd_0 = (opt_no_evo_std_yrs <= 2)*mn_std_inh_hds(iYear);
                succ_std_inher_head = std_inh_hd_0*ones(size(succ_std_inher_head));            
                rng_std_inh_hd = 0;
                 min_std_inh_hd = std_inh_hd_0; % mn_std_inh_hds(iYear);
                 max_std_inh_hd = std_inh_hd_0; % mn_std_inh_hds(iYear);

                 if opt_no_evo_std_yrs >= 2 % otherwise keep heterogeneity
                     succ_disp_d_dep = mean(succ_disp_d_dep)* ...
                         ones(size(succ_disp_d_dep));
                 end

                 % Both opt 2 and 3 have zero std disp
                 std_inh_dsp_0 = (opt_no_evo_std_yrs == 1)*mean(succ_std_inher_disp);
                 succ_std_inher_disp = std_inh_dsp_0* ... 
                     ones(size(succ_std_inher_disp));
                 min_std_inh_disp = std_inh_dsp_0; % mean(succ_std_inher_disp);
                 max_std_inh_disp = std_inh_dsp_0; % mean(succ_std_inher_disp);
                  
                  % no more need to vary inheritance regarding evoln disp dist
                  rng_std_inh_disp = 0;
 
                 if zug_opt
                     std_inh_sp_0 = (opt_no_evo_std_yrs <= 2)*mn_std_inh_sps(iYear);
                     min_std_inh_sp = std_inh_sp_0; % mean(succ_std_inher_signp);
                     max_std_inh_sp = std_inh_sp_0; % mean(succ_std_inher_signp);
                     rng_std_inh_sp =  0; %   
%                      succ_std_inher_signp = mean(succ_std_inher_signp)* ...
%                          ones(size(succ_std_inher_signp));
                     succ_std_inher_signp = std_inh_sp_0* ...
                         ones(size(succ_std_inher_signp));

                 end

             end

      end
         
      %% if last generation (year), save all field variables
          if iYear == N_mig_Sims && iWarmup > N_init_Sims % -1
              
             succ_heads = rad2deg*succ_heads + 180;
             all_inher_heads = rad2deg*inher_heads + 180;
             all_zugkn_heads = rad2deg*zugkn_inher_heads + 180;
             succ_zug_inh_hds = rad2deg*succ_zug_inh_hds + 180;             
             blat_succs = blats(stoppedAndArrived,:);
             blon_succs = blons(stoppedAndArrived,:);
             lat_deps_succ = lat_bs_deps(stoppedAndArrived);
             lon_deps_succ = lon_bs_deps(stoppedAndArrived);                
             all_day_starts = day_starts;

%              coeff_succ = coeff_magn_steer(stoppedAndArrived); 
               if calibr_comp(1) == 2
                     
                   if magcl_opt == 1
                       
                     all_coeffs = atan(coeff_magn_steer)*rad2deg;
                     
                   elseif magcl_opt == 2
                       
                    all_coeffs = atan(coeff_magn_steer)*rad2deg;
                       
                   else
                       
                       
                   end
                     
               elseif calibr_comp(1) == 3
                     
%                      all_coeffs = acos(coeff_magn_steer)*rad2deg; 
                    if magcl_opt == 1
                       
                         all_coeffs = atan(coeff_magn_steer)*rad2deg;
                         
                    elseif magcl_opt == 2
                        
                         all_coeffs = atan(coeff_magn_steer)*rad2deg;
                                                 
                    else
                        
                        
                    end
                     
               else
                     
                     all_coeffs = coeff_magn_steer; 
                     
               end
               
               all_zug_signps = zug_signps;
               all_disp_d_dep = disp_d_dep*mn_2_sig;   
               all_disp_d_arr = disp_d_arr*mn_2_sig;
                 
               succ_disp_d_dep = succ_disp_d_dep*mn_2_sig;
               succ_disp_d_arr = succ_disp_d_arr*mn_2_sig;
               
          end
                                
          % increment year or warmup counter
          if iWarmup <= N_init_Sims % || (N_init_Sims == 0 && iYear == 1)
              
              iWarmup = iWarmup +1;
              
          else
              
              iYear = iYear +1;
              
          end
    
    end   
     
% upadate last steps 
          annotate_fields

      end

