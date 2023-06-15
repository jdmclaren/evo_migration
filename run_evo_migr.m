disp('Choose an inherited compass');

comp_list = {'magnetic', 'geographic (star)'};

prompt_inh = {'Choose inherited compass'};
endg_cmp = listdlg('PromptString',prompt_inh,'ListString',comp_list);
endog_comp =  (endg_cmp==1)*1 + (endg_cmp==2)*3 % 3 %  

prompt_prim = {'Choose primary compass'};
prim_cmp = listdlg('PromptString',prompt_prim,'ListString',comp_list); 
calibr_comp =  [prim_cmp==1 0 prim_cmp==2];

prompt_star = {'Choose in-flight compass'};
star_cmp = listdlg('PromptString',prompt_star,'ListString',comp_list); 
star_comp=  star_cmp == 2;

prompt_sp = {'Choose signpost option'};
sp_list = {'none','inclination','intensitsy','declination'};
sp_opt = listdlg('PromptString',prompt_sp,'ListString',sp_list);
zug_opt = sp_opt > 1; 
zug_signp = (sp_opt==2)*1 + (sp_opt==3)*3 + (sp_opt==4)*2;

prompt_full = {'Full or shorter simulation'};
full_list = {'full (124 years)','short (10 years)'};
full_opt = listdlg('PromptString',prompt_full,'ListString',full_list);
full_run = full_opt == 1;

delete(gcp('nocreate'))

n_par_pools = 4; % 7; % 

dir_nm  = 'evo_mdl_out/'; %  
addpath 'circ_stats/'
addpath 'geo_data/igrf_data/'
addpath 'geo_data/mdl_1_deg_data/'

dir_nm  = 'evo_mdl_out/'; % 

% folder with geographic, geomagnetic and landcover data
geo_filnm = 'geo_data';

close_prev_figs =  0; % 1; %
if close_prev_figs == 1
    close all
end

% Is it a full or shortened run?
% % (124 year with 75 year spin-up) or partial run options
% Short run currently set (below) to 10-year simulation and 12-yr spin-up
% full_run = false % true % 

% plot and save options
plot_map_opt = 1 % 0; %  
save_opt =  1 % 0 % 

% option to weight prob of next gen (ie parental)
% locn by difference in incln or decln from
% natal site of (new) parents
wt_natal_decln = false; % true; %   
wt_natal_incln = false; % true; %  
wt_natal_intns = false; %  true; %  

% allow inds from same loc to breed
opt_same_loc_brd = false; % true; %   

% set option to "accelerate" evoln of warmup by selecting
% parents on their similarity regardig  disp, std_hds and std_inh_sp
% Currently with probability weight quadratic to these factors
accel_wrmp = true; % false; % 
if accel_wrmp
    accel_wrm_str = '_acc_wrm';
else
    accel_wrm_str = '_no_acc_wrm';
end

% file name markers (of parameter options)
wt_natl_str = [num2str(wt_natal_intns), num2str(wt_natal_incln), ...
    num2str(wt_natal_decln)];

% ---------------------------------------------------
% Input parameters:

% Choose migration system
% 1 = W Nearctic - S America,  
% 2 = W Europe - Afr, 
% 3  Palearctic wheatear distribution %  Europe-Alaska - Africa,
% 4 = Alaska NZ
% 5 = Nearctic (Can - Grld) plus Iceland wheatear distribution %
% 6 = Germany - Afr
% 7 = Grd - Icel Eva Kok red knot
% 8 - high Cdn Arctic to S Am (shorebirds)
% 9 - wheatear Alaska - Africa,
% 10-11 Blackpolls to SAm (East / West popns)
mig_sys = 5 % 2 % 1 %  4 %  10 % 3 % 10 %  11 %

% Option to return arrival at EU
arr_EU_opt = mig_sys == 5 % false % 

% flag for wheatears (for options below)
is_whtr = mig_sys >= 3 && mig_sys <= 5;

% Primary compass including
% 1st index: mag 0 none, 1 Lox (options 2 MgclT, 3 MgclP, removed)
% 2nd index: sun 0 none, 1 Fix Az, 2 time-comp
% 3rd index: geogr 0 none, 1 Loxodrome (star), 2 great circle
% calibr_comp =  [1, 0, 0]; %  [0, 2, 0]; %  [1, 0, 0]; %   [1, 0, 1]; %     %  (mig_sys == 2)*[1, 0, 0] + (mig_sys ~= 2)*
% Star or mag compass in flight
% star_comp =  (calibr_comp(1)~=1); % false; % true % 
% Flag for simulatons
magn_star_night = 1 + (star_comp == true); % 1 = mag lox (offset by decln), 2 = geogr lox
% Which geomagnetic model 
% 1 = IGRF data (mean core field)
% 2 is geomagnetic dipole
% (should check if Dipole option still works)
geomagn_model = 1; % 2; % 

% Specify which compass is inherited, i.e., endogenous
% 1 is inher magn comp, 2 is inher sun comp 
% 3 is inh geogr head
% is_mg_endg = calibr_comp(1)~=0 && calibr_comp(3)==0; % false; % true; % 
% Assume nonmagnetic is geographic inheritance
% endog_comp =  ~is_mg_endg*3 + is_mg_endg*1 % 3 %                   

% Zugknick no (0) or yes (1)
% zug_opt = 1 % 0 % *(sun_head_reset~=2 || endog_comp ~= 2) %  1 %            
n_zugs = zug_opt*1 %  0 %  (zug_opt ==1)*1 % 2; % 
% Is Zugknick needed to be successful? (default, no)
n_req_zugs = zug_opt*0 %  1 % 1*(is_whtr); %1*mig_sys~=2;
%

% Which signpost? 1 incl 2 decln 3 total inten 
% 4 = vert inten 5 = hz inten
% zug_signp = zug_opt*(is_whtr*3 + ~is_whtr*3) %
% Require to be on land when reacting to signpost threshold (or not)
req_land_signp = true % ~is_whtr % false % 
% Set flags for file name indicators
if zug_opt == 1 && ~req_land_signp 
    req_lnd_str = '_sgn_wtr';
else
    req_lnd_str = '';
end

% Stepwise error (from any cue detection, transfers, maintenance or drift)
base_err = 15 % 2.5 % 0 % 5 % 
std_err_calibr = [base_err base_err base_err]; % (mig_sys~=7)*
% On commencement  of migration, assume
% perfect initial calibration (i.e., as inherited), for each compass
perf_init_calibr = 1*[1 1 1] % 1 %  
% additive factor to error when over water
rand_wtr_fact = 0; % 0.5; % 1; % 

% Precision of gaiuging signpost components
% effective std error (inclination will be converted to von Mises kappa)
std_err_incl_deg = 1*5; % in inclination in degs
std_rel_err_inten = 1*0.02; % inintensity in fraction (i.e. 0.02 = 2%)

if full_run
    
    N_inds = is_whtr*50000 + ~is_whtr*20000; % + (~is_whtr && mig_sys~=5)*1000 % 20 %   10; %
    
    N_init_fact = 1; %  + (mig_sys==5 && calibr_comp(1) == 4)*9;
    N_mig_Sims =  124; % 4; % 247; % 124; % 44; %  3 % 116; %  10 %  10; %  50 % 20 % 
    N_wrmp_rnd =  50; % 0; %  75; % 2; %
    N_wrmp_bck = 0; % 25; %  2; % 10; % 15; % 2; %
    N_wrmp_1st = 24; % 50; %  2; % 10; % 
    
    N_init_Sims = N_wrmp_rnd +  N_wrmp_bck + N_wrmp_1st;
    % 150; %  3 % 0 % 20; %- (mig_sys==5 && calibr_comp(1) == 4)*3; 5; % 12 %
 
    N_plot = 350;
    
    transp_plot = 0.75;
    
else
    
    N_inds = (is_whtr && mig_sys==3)*20000 + (mig_sys~=3)*10000; % + (~is_whtr && mig_sys~=5)*1000 % 20 %   10; %
    
    N_init_fact = 1; %  + (mig_sys==5 && calibr_comp(1) == 4)*9;
    N_mig_Sims =  is_whtr*10 + ~is_whtr*20; % 10; %  
    N_wrmp_rnd = 10; %
    N_wrmp_bck = 0; % 5;
    N_wrmp_1st = 2; % 5;
    N_init_Sims = N_wrmp_rnd +  N_wrmp_bck + N_wrmp_1st; % is_whtr*20 + ~is_whtr*2; % 6; % 
    N_plot = 150; % (mig_sys==9)*3000 + (mig_sys~=9)*1000
    transp_plot = 0.75;
    
end

N_init_str = num2str([N_wrmp_rnd N_wrmp_bck N_wrmp_1st]);

% Year and date options
yr_start =  full_run*1900 + ~full_run*2002; % 1965; % % 1995; %  0 % 2015 %  1950 %
month_start = 8; % 9; %
day_start = 20; % 1; % 15; %
std_day_start = 5; %  5 % 0 
max_dev_day_st = 14; % 
start_date = [yr_start month_start day_start ...
    std_day_start max_dev_day_st];

% Forward or bacwards simulation through years of IGRF data
% 0 is always 1900, 1 is fwd 1900-2015, 2 bckwd 2015-1900
fwd_back_no_time = 1 % 0 %  2 %     
if fwd_back_no_time == 0
    fwd_bck_str = '_1_yr';
elseif fwd_back_no_time == 1
    fwd_bck_str = '_fwd';
elseif fwd_back_no_time == 2
    fwd_bck_str = '_bckw';
end    

% Constraints on migration (for wheatear paper, minimal)
% Maximum duration
max_mig_durs = [110 100 140 140 90 90 50 70 100, 70, 70];
max_mig_dur = max_mig_durs(mig_sys);
% Thresholds (%) landcover within 1-degs) for resources and hazards
thrsh_snow_barrn =  0.85; %  (~is_whtr)*0.85 + (is_whtr)*1;% 1; %
% Thresh for max fuelling down to zero for half of ndvi
thrsh_ndvi_stop_max =  0; % ~(is_whtr)*0.2; %; %(is_whtr)*0.05 + (mig_sys == 2)*0.25 + (mig_sys > 9)*0.3 + (is_whtr)*0.25; %0.05; % 0.3;  %  
thrsh_ndvi_stop_min = thrsh_ndvi_stop_max/5; %
thrsh_ndvi_breed = (~is_whtr)*0.25 + (is_whtr)*-1; %0.05; 
thrsh_low_veg_breed = (is_whtr)*0.005; % 
thrsh_low_veg_stop = 0; % (is_whtr)*0.1; % 0.05; %;0.1; % 
thrsh_trees_breed = 0; % ~(is_whtr)*0.05; %;
thrsh_trees_stop =  0; %  ~(is_whtr)*0.001; % (mig_sys>9)*0.05; % 
thrsh_conif_breed = (mig_sys>9)*0.2;
thrsh_all_veg_stop = ~(is_whtr)*0.1; % 0; % ~(is_whtr)*0.1; % *0.4; %  
max_trees_stop = 1; % ~(is_whtr)*1 + (is_whtr)*0.5;
% min frac any landcover var for stopover (will be prorated to relevant
% field, eg thrsh_low_veg_stop for whtrs
min_any_veg_stp = (is_whtr)*thrsh_low_veg_stop/5 + thrsh_all_veg_stop/5; % ~(is_whtr)*thrsh_trees_stop/5; % 
% thrsh_trees = (mig_sys > 9)*0.2; % ~(is_whtr)*0.2;
thrsh_elev =  Inf; % 2500; % 2000 % (is_whtr)*3500 + ~(is_whtr)*2500; % 2000; %  4000; %  (thresh metres)

% Lat above which bird assumed lost (won't survive)
lost_Lat = (mig_sys~=5)*82.5 + (mig_sys==5)*87.5; % 80; % 

% Options for daily survival rates
daily_surv = 1; % is_whtr*0.99 + ~is_whtr*1; %  0.98 %95; %  is_whtr* ~is_whtr*1 + 0.995; %  0.9975; 0.99; % 
daily_surv_elev = 1; % 0.99; %  0.99; % 0.5; %  0.99;
daily_surv_barrn = 1; %~is_whtr*0.99 + is_whtr*1; %; % 1; %  0.95; % 
daily_surv_wtr =  1; % 0.99; %

% Sun compass time-compensation options
% Setting these to sun_head_reset = 1 and tc_comp_reset = 2 is best
% and other options not complete or checked
% sun_head_reset opt for zugkn's == 0 is no new zugk 
% (ie reset endog heading; for use with endog == 2)
% == 1 means set new inherited (Zug) heading
% == 2 means sun heading resets at every stopover (zug should be zero)
% == 3 means heading reset every night (even during multiple flights)
sun_head_reset = 0 + 1*(calibr_comp(2)==2)%    1 %  1 % 1 %    
% time comp reset == 1 to reset any time-zone bias after nonstop flt seq's
tc_comp_reset = 1*(calibr_comp(2) == 2) % 0; %   1*(calibr_comp(2)~=0) % 1 %       0 %          
% Assume can tc comp when over water
tc_comp_wtr = 1; % 0; % 

% Dearture, stopover and arrival polygons per migration system
% [Lon Lon Lat Lat]
% Iql_Gr_Lons = [-77.5 -62.5]; % [-30 -15]; %[-50 -40]; % 
Dep_verts_mig_sys = {{[-110 -90 50 65]}, ...
    {[-11 1.5 50 63], [1.5 17.5 50 72], [17.5 30 50 72], [30 60 50 72]}, ... % -220 -140  190[-160 -140 60 71],[-190 -150 50 61]
    {[-80 -15 55 80],[-15 17.5 50 74],[17.5 60 50 74],[60 110 50 74], ...
    [110 140 55 72],[140 190 60 72],[190 230 60 72]}, ... {[-80 -15 58 80],[120 220 60 72]}, ...
    {[155 220 60 70]},...
    {[-80 -60 57.5 62.5],[-80 -60 62.5 70],[-80 -45 74 85],[-60 -45 59 74], ...
    [-45 -25 59 69],[-45 -10 69 75],[-45 -10 75 85],[-80 -65 69 74]} ... % ,[-25 -12.5 63 67]{[-80 -30 59 84]},  ... % {[-80 -13 55 84]},  ... % [-80 40 55 80],[40 230 55 80]
    {[5 15 45 55]},{[-70 -70 80 80]},{[-140 -110 65 70]}, ...
    {[190 210 54 70]}, {[-170 -85 55 70],[-85 -72.5 50 70],[-72.5 -50 42 60]}, {[-80 -50 47 60]}}; %,185 210 155 170 [ 60 70][-80 -60 62.5 68]  [-170 -160 60 65] [-80 -10 52.5 82.5]

Arr_verts_mig_sys = {{[-70 -50 -10 0]}, ...
    {[-17.5 45 5 15]}, ...  ...{[-17.5 -5 5 15],[-5 5 5 15],[5 15 5 15],[15 30 0 15],[30 45 5 15]}, ...
    {[-17.5 0 5 15],[0 25 5 15],[25 45 0 15]}, ... {[-17.5 45 5 15]}, ...{[-20 12.5 5 15],[12.5 42.5 5 15]}, ... {[-17.5 0 5 15],[0 25 5 15],[25 45 5 15]}, ... 
    {[30 45 5 15]},{[-20 10 5 15]},{[-20 10 10 20]}, ...
    {[0 15 50 55]},{[-70 -50 -10 0]}, ...
    {[30 40 -10 15]},{[-79 -57 5 13]},{[-79 -60 5 13]}}; %[-20 5 0 15][-20 0 60 70] [169 184 -40 -33]};[30 50 -10 15]

Stop_verts_mig_sys = {{{[-70 -50 -10 0]}}, ...
     {{[-10 10 30 45],[20 50 30 45]}}, ... % {{[-20 0 15 30],[20 35 30 40]}}, ... ,[10 25 30 45]
    {{[-17.5 15 15 60],[15 35 30 45], ...,[5 15 15 45]
    [35 50 35 45],[50 65 35 55],[65 80 45 60]}}, ... %  [90 110 50 72],[110 140 55 72],[140 190 60 72]}}, ... % {{[-10 10 30 65],[15 30 30 40],[30 50 40 45],[50 60 35 45],[70 120 50 65],[120 190 60 70]}},
    {{[45 70 35 50]}}, ... % [-30 15 30 40],[45 60 35 50]{[115 125 50 70],[70 80 40 65],[40 60 27 45]} {[-15 -4 30 36],[-25 -5 0 15]}
    {{[-20 10 40 60]}},{{}},{{}},{{}},{{}}, ... [-20 10 20 63]
    {{[-79 -75 21 35],[-77 -74 35 40.5],[-74 -69 40.5 42.5],[-69 -60 17 46]}}, ... [-75 -59 13 24]42.5 46],[-75 -61 42.5 46
    {{[-76 -60 38 46]}}};

% Start, stop or arrive on land options
% Start on land required = 1 (or not, 0)
req_land_init = 1; % *(mig_sys~=7);
% Require land on arrivals at stopovers and final destn
req_land_arrs = 1*[1 1 1 1 1 1 1 1 1 1 1]; % [1 1 1 0]; % 
req_land_stops = 1*[1 1 1 1 1 1 1 1 1 1 1]; % [1 1 1 0]; % 
req_land_arr = req_land_arrs(mig_sys);
req_land_stop = req_land_stops(mig_sys);

% Options to detour preemptively around coast within certain angle
% (for wheatear study, none. See main script for options)
dtrCstOpt = 0 % 1 %   0*(mig_sys < 3)*1 %  (mig_sys < 3)*1 %  (mig_sys < 3 || mig_sys == 5)*1 %% % 1 %  
dtrReComp =  0; % 1 % 
% cst strats 4 = closest pref first
% Else 5 = most numerous first,
% 6 = most numerous median (both pref side if tied)
% Else 8, 9 as 5, 6 but opp. side first
dtrCstStr = 5; % 4; %   6; % 9; % 8; % 8 %   
% hours ahead to "check" if over ocean (for detour option)
% hrsAheadDtrCst = 2; % 1; %
kmsAheadDtrCst = 300; % 200; % 
short_detrs =  0; %1; % shorter nights when detouring
maxDtrSteps = 6;

% Option to stop prematurely over coast when (would be) over water on
% arrival e.g. at dawn
stopCstOpt = 1;
% Proactive chacking of visible coast: which decile of (nightly) flight? 
fracStopCst = 7; % 5 %   mig_sys ~= 4*7 + mig_sys == 4*10; % 1; % 
% stop at coast when (daily/nightly) flight initiated over water
stpCstWtrFlt = 1; %  0 % (mig_sys ~= 4)*1 % 0; %
stpovCstWtrFlt = 1; % 0 %  (mig_sys ~= 4)*1 % 0; %
% Which max distance to find land in flight over water (depends on flt altitude)
d_stop_maxs = [100 300 300 300 300 100 200 100 400 400 400]; % [200 200 200 400 400 200];
d_stop_max = d_stop_maxs(mig_sys);
 
% Reorient water strategy i.e. when flight-step begins over water
% Tured off for wheatear study
% 1 = head E or W (opp. sense of pref dir)
% 2 = mirror image current pref dir
% 3 = mirror image most recent over-land pref Dir
% 4 = clockwise perpindic. to  most recent over-land pref Dir
% 5 = anti-clockwise perpindic. to  most recent over-land pref Dir
% 6 = perpindic. to most recent over-land pref Dir on 'opposite' side South
% 7 = perpindic. to initial pref Dir on 'opposite' side South
% 0 = South
reOrWtr_opts = 0; % 1 % 1*(mig_sys == 2); % *(mig_sys == 3) % 0*(mig_sys < 3)*1  % 1 % test_comp_type~=3; %
reOrWtrOpt = reOrWtr_opts; % this is always reOrient over water at dawn! 
reOrWtrStr = reOrWtr_opts*3; % * 0 % 5 % 
reOrReComp = 3; % 0; % 1; % reOrWtr_opts;
% this option limited the number of steps before reorienting 
% % - should check if still works
n_wtr_stps = 0; %  0; %   3; %  

% Options for requisite stopover areas and required Zugknicks
% Not used for wheatear study
% Here "stop" refers to requisite stopover areas (i.e., not any
% intermittent stopovers)
stop_opt = 0 % 1 % mig_sys > 9  || is_whtr % 1 % 0 % 1*(mig_sys == 2  |(zug_opt == 1)* %  0 % 1 %   calibr_comp(1) ~= 4 %  
% Does the stopover require a Zugknick (e.g., if we presume ecological
% constrainted to do so)
require_stop_polys = {[0], stop_opt*[1], stop_opt*[1], ... % stop_opt*[1 1]
    stop_opt*[1], stop_opt*[1],[0],[0],stop_opt*[1],stop_opt*[1 0 1], ...
    stop_opt*[1], stop_opt*[1]}; % stop_opt*[1 1 0]
require_stop_poly = require_stop_polys{mig_sys}; % false % 
require_stopover_polys = [0 0 0 0 0 0 0 0 0 0 0];
require_stopover_poly = require_stopover_polys(mig_sys);
zgknk_only_in_stops = {0; (zug_opt == 1)*[0]; ...
    (zug_opt == 1)*[0]; (zug_opt == 1)*0; (zug_opt == 1)*0; ...
    0; 0; 0; 0; (zug_opt == 1)*1; (zug_opt == 1)*1}; %  {0; [1 1];[1 0]; 1; 0}; % (mig_sys == 2 || mig_sys == 3)*1; % 1 %  
zgknk_only_in_stop = false; % true %  is_whtr %  mig_sys ~=2 % mig_sys > 9 %  true %  mig_sys ~=2 %false %  mig_sys~=2 %% false %  zgknk_only_in_stops{mig_sys}; % false % mig_sys > 9; %  | %stop_opt*
require_stopover_zug = true; % false; % mig_sys ~= 2 %  mig_sys < 10 % 

% Stopover warmup option requires both lands in stopover polygon and 
% (currently) Zugknick for success
% Not used for wheatear study
% These constraints are possibly released during subsequent (yearly)
% simulations
require_stop_warmup = require_stop_poly; % is_whtr; % (stop_opt==1)*(ones(size(require_stop_poly)));
req_zg_in_stpPly_wrm =  mig_sys == 3 || mig_sys == 4; %  true; %  ~is_whtr; %  require_stopover_zug; % zug_opt == 1 && stop_opt==1 % 

% Bird flight and stopover behaviour
% Options to add proxy for wind (autocorrelated but not airspeed related)
incl_wind =  0; % 1; %    
% flight speed m/s
Va_mps_mig_sys = [12.5 10.5 15 15 ...
    15 12.5 15 15 12.5+(incl_wind==0)*2.5 12 12]; % = +(incl_wind==0)*2.5
Va_mps = Va_mps_mig_sys(mig_sys);
% Specify number hours before/after dawn/dusk
% 1.5 is approx. civil dusk dawn +- 1 hour
offset_dep_dusk = 1.5; % & 2.5*is_whtr + 1.5*~is_whtr;
% max extra or fewer flight hours per ind each night
var_fl_hrs = 1;

% wind strength (if used) 
% Can be seen as effective drift after (partial) compensation
wind_str = 0.25 + (mig_sys==4)*0.1  + (mig_sys==5)*0.1 + (mig_sys>9)*0.15;
% option for stochastic winds
switch_wind_opt = (incl_wind==1); % 0; % 
dep_wtr_TWs = 0; % (incl_wind==1); % 

% old option number consecutive flight nights before extended stopover
% n_flts_seq = (mig_sys==1)*4 + (mig_sys==2)*4 + (is_whtr)*4 + ...
%    (mig_sys==6)*4 +(mig_sys==7)*10 + ...
%     (mig_sys==8)*5 + (mig_sys==9)*10 + (mig_sys>9)*5;

% stopover duration (days)
stop_durns = [6 5 5 5 5 6 7 2 5 7 7];
stop_durn = stop_durns(mig_sys);
% longer durns after Zug sp
stop_durn_sps = [6 14 14 14 14 6 7 2 5 7 7];
stop_durn = stop_durns(mig_sys);
stop_durn_sp = stop_durn_sps(mig_sys);
% max cumulative (half-day) flights (over barriers)

% Threhold min number of flights before wanting to stop over
min_n_fls_stp =~(is_whtr)*1.5 + (is_whtr)*3; % 1.5; % 2;

% May 2022 add signposted refuelling as in at Zugknicks or just pre barrier
% These refuellings are now only after signpotsed
hs_trSah_WW = is_whtr*72 + ~is_whtr*60; % 54; %  60; %
hs_xAsia = 60; % % 80; 36; % 
hs_trSah = 60; % 80;
hs_NA = 80; % 48;
% These are initial trans Atl fl hours
% (once wheatears reach signpost they have trans Sah potl fl hrs)
hs_tr_Atl = 72; %  96; % 100; %   120; %    

min_hs_fact = is_whtr*1.5 + ~is_whtr*2; % 1.5; % 1; %  2; % 1.2; %   

max_ns_sgnp_hs_s = [{[24]}, ...
    {[hs_trSah_WW hs_trSah_WW hs_trSah_WW hs_trSah_WW  ...
    hs_trSah_WW hs_trSah_WW]}, ...
    {[hs_trSah hs_trSah hs_trSah hs_trSah hs_xAsia hs_xAsia hs_xAsia]}, ...
    {[hs_xAsia]},{hs_trSah_WW*ones(size(Dep_verts_mig_sys{mig_sys}))},{[24]},{[24]},{[24]},{[24]}, ... hs_NA
    {[hs_NA hs_NA hs_NA]},{[hs_NA]}];
min_ns_sgnp_hs_s = cellfun(@(x) x/min_hs_fact,max_ns_sgnp_hs_s, ...
    'UniformOutput',false);

max_ns_sgnp_hs = max_ns_sgnp_hs_s{mig_sys}; % 7;
min_ns_sgnp_hs = min_ns_sgnp_hs_s{mig_sys}; % 7;

hs_reg_max = is_whtr*48 + ~is_whtr*36 + (mig_sys == 5)*24; % 12; % 
% hs_reg_min = is_whtr*32 + ~is_whtr*24 + (mig_sys == 5)*12; % 
hs_reg_min = hs_reg_max/min_hs_fact; % 

min_ns_fl_hs = hs_reg_min*ones(size(min_ns_sgnp_hs));
max_ns_fl_hs = hs_reg_max*ones(size(max_ns_sgnp_hs));

% For Atl wheatears & BLPWs, start with for pre barrier beginning
% in which case add even more flt hours

max_init_fl_hs = max_ns_fl_hs; %[6 4 5 40 5 5 7 7]; % 3 7 5; 64

if mig_sys == 3 
    max_init_fl_hs(1) = hs_tr_Atl; % 7;
    max_init_fl_hs(end) = hs_xAsia; % 7;
    max_ns_fl_hs(1) = max_ns_fl_hs(1) + 12;
    max_ns_fl_hs(end) = max_ns_fl_hs(end) + 12;
    min_ns_fl_hs(1) = min_ns_fl_hs(1) + 12;
    min_ns_fl_hs(end) = min_ns_fl_hs(end) + 12;
elseif mig_sys == 4
%     max_ns_fl_hs(1) = max_ns_fl_hs(1) + 12; 
%     min_ns_fl_hs(1) = min_ns_fl_hs(1) + 12;
elseif  mig_sys == 5
    max_init_fl_hs(1:end) = hs_tr_Atl; % 7;
    max_ns_fl_hs(1:end) = hs_tr_Atl; %
    min_ns_fl_hs(1:end) = hs_tr_Atl/min_hs_fact;
elseif mig_sys > 9
    max_init_fl_hs(3) = hs_tr_Atl; % 7;
end

min_init_fl_hs = max_init_fl_hs/min_hs_fact;

% max_init_fl_hs = max_init_fl_hs_s{mig_sys}; % 7;
% min_init_fl_hs = min_init_fl_hs_s{mig_sys}; % 7;

% endurance flight simulates 12 hour flights without stopping 
% until arrived (or bust)
endur_flts = [0 0 0 0 0 0 1 0 0 0 0]; % 3[1 1 1 1 1]; % 
endur_flt = endur_flts(mig_sys);
% if not endurance, use min / max fl hours (with extra flight hours for
% initial barrier crossings)
min_fl_stp_hs = 6; % 5;
max_fl_stp_hs = 12;
min_init_fl_stp_hs = min_fl_stp_hs; % min_fl_hs*(mig_sys ~= 4 & mig_sys ~= 5)  + ...
%     12*(mig_sys == 4) + 12*(mig_sys == 5);

% options for inheriance of headings and signposts
% and relative weights based on mutual proximity
% and speed of successful migration (all in prev generation)
% R_natal is max natal dispersal of returning migrants
% R_natal = 200; % 500; % 1000 % 1500; %  250; % 
d_min = is_whtr*(zug_opt*5 + ~zug_opt*5) + ...
    (mig_sys > 9)*100 + (mig_sys == 2)*50; % 50 % + is_whtr*200; % 200; % 20; %
d_max = 50; %  d_min; % 400; %50 % + is_whtr*200; % 2000; % 
min_disp_d_dep_km = 0.025; % 100; % 25; %  250; % 25; %  0.0025; % 25; % % 0.01; % 0.01; %50 % 0.4; % 2; % d_min; % 20;
max_disp_d_dep_km = 25; % 100; % 250; % 10; % 50;% 250; %   250; % 100; % 25; %  50 % 10; % 25; %0.4; %2; % d_max; % 400;
R_natal = round(100*(min_disp_d_dep_km + max_disp_d_dep_km)/2);

% Aug 2022 separate the intrinsic var in inh disp dist from allowable range
% (which will be mixed at least during warmup)
min_std_inh_disp_km = 0.0025; % 0.005; %  0 0.4; % 2; % d_min; % 20;
max_std_inh_disp_km =  1; %    5; % 0.4; %2; % d_max; % 400;
std_R_natal = round(100*(min_std_inh_disp_km + max_std_inh_disp_km)/2);
% if ~is_whtr % mig_sys > 9  % mig_sys ~= 2 % > 0 % 9 % mig_sys  
%     min_disp_d_arr_km = Inf; % d_min; % d_min; %20;
%     max_disp_d_arr_km = Inf; % d_max; %   d_max; %400
%     min_disp_d_zug_km =  Inf; % d_min; %  d_min; %20;
%     max_disp_d_zug_km = Inf; % d_max; %   d_max; %400
% else
    min_disp_d_arr_km = Inf; %5*d_min; %d_min; %   + ~is_whtr*2*d_min; % d_min; %20;
    max_disp_d_arr_km = Inf; %5*d_max; %d_max; %    d_max; %400
    min_disp_d_zug_km = Inf; % d_min; %  d_min; %20;
    max_disp_d_zug_km = Inf; % d_max; %   d_max; %400 
% end

% option to weight selectibility on speed (inverse  duration migration)
% seems reasonable but could confound known (programmed) migratory
% connectivity
wt_inher_speed = 0 % 1 %   

min_std_inh_sp = 1*(zug_opt == 1)*((zug_signp ~= 3)*0.1 + ... 0.05
    (zug_signp == 3)*0.001); % .005); %  + (calibr_comp(1) ==7)*[0.05 5];
max_std_inh_sp = 1*(zug_opt == 1)*((zug_signp ~= 3)*1 + ...
    (zug_signp == 3)*0.01); %0.005); % .05  + (calibr_comp(1) ==7)*[0.05 5];
min_std_inh_hd = 0.025; % 250; % .1; %2; % 0.5 %5 % 0 % 5 % 1*.5 % 0 % (~magcl || (calibr_comp(1) == 0 && zug_opt ~= 2))*5 % 0 %     5 %  
max_std_inh_hd = 5; % .1; % 0.1; % 5; % 2; %  5; % 5 %.1 % 50 % 1*5  % 0 % 10 %  (~magcl || (calibr_comp(1) == 0 && zug_opt ~= 2))*5 % 0 %     5 % 
% --------------------------------------

% Aug 2022 option to stop 'evolving' std's in inh hds and sgnps 
% after warmup ie during sequence of years (they will be held fixed)
% This arguably better represents actual inheritance at micro evolutionary
% time scales

% = 0 keep optimzing std's (in inh, sp & nat disp) after warmup
% = 1 set std's inh and disp to popn mean values
% = 2 set stds inh to mean vals but std disp to zero (nat disp to pop mean)
% 3 = set all std's to zero, nat disp to pop mean (Dec 2022: here warmup will evolve 
% until simul, i.e., std's inher used to create viable init popn)
opt_no_evo_std_yrs = 2; % 0; %  3; %  1; %    true; % false; % 

if opt_no_evo_std_yrs == 0
    opt_std_yr_str = '_var_all_evo_std';
elseif opt_no_evo_std_yrs == 1
    opt_std_yr_str = '_var_no_evo_std';
elseif opt_no_evo_std_yrs == 2
    opt_std_yr_str = '_mn_dsp_no_evo_std';
else % == 3
   opt_std_yr_str = '_mn_inh_dsp_no_evo_std';   
end
% map projections: Mercator is angle preserving (rhumblines are straight lines)
% ortho is azimuthal (projection from 'outer space')
% stereo is cylindrical (shows Earth curvature)
% gnomonic is strange but shows great circles as straight lines

% Option to plot changes in heading due to sun and magnetic effects
plot_heads_opt = 0; % 1; %     
% colorbar is time in days (plot_cb == 0) 
% or heading in degs (plot_cb_opt = 1)
plot_cb_opt = 1; % 0; %

%%%%%% -------------------------------------
% these are vestigial and should be removed
% pos_neg_coeffs = [1 2 2 3 1 3 3 3 2 3 3]; % [1 2 2 3 1 3];
% pos_neg_coeff = pos_neg_coeffs(mig_sys) %
% max_magn_coeff = pi*(calibr_comp(1)~= 0 && zug_opt~= 2 ...
%     && calibr_comp(1)~=3 && calibr_comp(1)~=4) %5 %
% mincoeff_magn_steer =  -(pos_neg_coeff>1)*max_magn_coeff;
% maxcoeff_magn_steer = (pos_neg_coeff~=2)*max_magn_coeff; % 0.5; %  
% % exponent for shifting compass (calibr_comp(2) == 1)
% magn_exp = 1; % (calibr_comp(1)<5)*1 - (calibr_comp(1)max_disp_d_arr_km>=5); % 0; %  
% rel_abs_magn_steer = (calibr_comp(1)<5 && calibr_comp(1)~=2)*1 %  % 0 %  0 = rel, 1 = abs

% Old options for reflecting about axes with sun compass
% if true, refl about (geo or Sun) South axis
refl_sun_mag_opt =  calibr_comp(2) == 1 && mig_sys == 5 % false %
clockwise_sun_opt =  refl_sun_mag_opt*(mig_sys == 1 || ...
    mig_sys == 5 || mig_sys == 8) % true % false % 
sign_clockw_sun = ~clockwise_sun_opt -(clockwise_sun_opt)
% fixed comp rel to geo South (opt 1) or magn S (opt 2, may
% be more reliable to calibrate but could result 
% in strange routes if decln changes sign)
geo_mag_refl_sun_opt = 1 % + (clockwise_sun_opt) % 1  %  

% old options magnetoclinic compass (can be removed)
magcl = (calibr_comp(1)==3 || calibr_comp(1)==4); 
min_std_inh_steer = (~magcl && calibr_comp(1) ~=0)*0.01; %*0.05; % 0 %   0 % *pi/90 % 
max_std_inh_steer = (~magcl && calibr_comp(1) ~=0)*0.1; %*0.05; % 0 %   0 % *pi/90 % 
% max mag shift (needed?)
max_magn_shift = 2*pi;
% magcl_opt specifies the reference for the magnetoclinic parameter 
% (only has effect if calibr_comp(1) = 3 or 4)
% 1 is as per Kiepenheuer (all vert compt, hz projn either transv as Kiep or pll
% 2 is along "surface" of inclination along geomagn (S) axis (either tran or pll to it)
% 3 is all hz compt and projn of vert compt along heading (angle still trans or pll)
magcl_opt = 3; % 2; %  3; % 
% if reset opt then magcl projn is reset to init head when fails
% should be as alternative i.e. not to combine with zug_opt 2
magcl_reset_opt = 0; %(calibr_comp(1) == 3 || calibr_comp(1) == 4)*1; %

%%%%%% ------------------------------------------

% strings for output and filenames
comp_strs = {'mag',['sun_' num2str(sun_head_reset) ...
    '_' num2str(tc_comp_reset)],'geo'};
endog_str = comp_strs{endog_comp};
calibr_str = '';
z_str = {'', 'zugkn '};
zug_str = {'', 'incl sgnp ','dec sgnp ','ttl sgnp ','vert sgnp ','hz sgnp '};
mag_str = {'','Lox ','trsv magncl ','pll magncl '}; 
sun_str = {'','fix sun ','tc sun '};
geo_str = {'','Lox ','gt crc '};
cst_str = {'', 'coast '};
reO_str = {'','reOr '};
magn_star_night_str = {'magN ','ggrN '};
any_zug = n_req_zugs > 0;
contn = mig_sys;
cont_str = {'NAm SAm ','N Eur Afr ', 'Palearctic Afr ', 'Alk NZ ', ...
    'NWAtlAfr ','Ger Afr ','BlSwan Grnl knot','Hi Arct SAm', ...
    'Blackpoll West','Blackpoll East'};
strat_string = [z_str{any_zug+1} zug_str{zug_signp+1} mag_str{calibr_comp(1)+1} ...
    sun_str{calibr_comp(2)+1} geo_str{calibr_comp(3)+1} ...
    magn_star_night_str{magn_star_night} cst_str{dtrCstOpt+1} ...
    reO_str{reOrWtr_opts+1}]
is_err_mag = std_err_calibr(1) > 0;

err_string = ['_err_'  num2str(std_err_calibr(1)) '_' ...
    num2str(10*std_err_incl_deg) '_' num2str(100*std_rel_err_inten)];

wind_strs = {'_no_wind','_wind'};
wind_string = wind_strs{incl_wind+1};

% now set arguments spcific to mig strat to call main routine
% Dep_verts = Dep_verts_mig_sys{mig_sys}*dg2rad;
dg2rad = pi/180;
Dep_verts = cellfun(@(x) x*dg2rad,Dep_verts_mig_sys{mig_sys}, ...
    'UniformOutput',false);
% Arr_verts = Arr_verts_mig_sys{mig_sys}*dg2rad;
Arr_verts = cellfun(@(x) x*dg2rad,Arr_verts_mig_sys{mig_sys}, ...
    'UniformOutput',false);
% need to convert each stop cell array)
for iSt = 1:numel(Stop_verts_mig_sys{mig_sys})
    Stop_verts{iSt} = cellfun(@(x) x*dg2rad,Stop_verts_mig_sys{mig_sys}{iSt}, ...
        'UniformOutput',false);
end
for i_comp = 1:3
    
    if calibr_comp(i_comp) ~= 0
        
        calibr_str = [calibr_str ' ' comp_strs{i_comp}];
        
    end
    
end

tic

    [successful, succs_DepPs, succs_ArrPs, succs_DepArrPs, ...
        mn_succ_hds,  std_succ_hds, mn_zk_hds,std_zk_hds, ...
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
        fr_zks,succ_heads, zugkns,speed,blat_succs, blon_succs, ...
       coeff_succ, succ_zug_inclns,all_zug_signps, ...
       all_init_heads, all_zugkn_heads, all_geo_heads, all_mag_heads, all_sun_heads, ...
       all_coeffs, blats, blons,lat_bs_zugs, ...
       lon_bs_zugs,passed_Stop, arrived, survived, lost, stoppedAndArrived, ...
       all_ys_inher_heads, all_ys_init_geo_heads, all_ys_zugkn_heads, all_ys_zugkn_sps, ...
       all_ys_stoppedAndArrived, all_ys_surv_fuel, all_ys_surv_wtr, all_ys_surv_barrn_elev,  ...
       all_ys_lost, all_ys_late, n_fl_step,date_jul,std_inher_head,std_inher_steer,std_inher_signp,day_start_succ, ...
       all_day_starts,all_inclns,all_declns,all_tots,all_sun_az,all_ovr_wtr, ...
       all_barrn,succ_disp_d_dep,all_disp_d_dep,succ_disp_d_arr,all_disp_d_arr, ...
       all_dates,all_fl_hs,all_max_hrs_flt,all_stp_ovs,all_stp_nrs, ...
       all_zg_nrs,all_cum_h_flt,all_survs,all_fins,all_arrs,all_losts, ...
       incln_0s,decln_0s,intns_0s,idx_deps,idx_arrs, ...
       mn_lat_EU_DepPs,std_lat_EU_DepPs,mn_lon_EU_DepPs,std_lon_EU_DepPs, ...
        all_ys_lat_arr_EU,all_ys_lon_arr_EU,all_ys_lat_zg,all_ys_lon_zg, ...
        all_ys_lat_fin,all_ys_lon_fin] = ...
     simul_mag_sun_comp_rtes(geo_filnm,N_inds,N_mig_Sims,N_wrmp_rnd, ...
        N_wrmp_bck,N_wrmp_1st,N_init_fact,start_date, ...
      zug_opt,zug_signp,req_land_signp,n_zugs,n_req_zugs, ...
      zgknk_only_in_stop,req_zg_in_stpPly_wrm,require_stop_poly,require_stop_warmup, ...
      require_stopover_poly,require_stopover_zug, ...
      Dep_verts,Arr_verts,Stop_verts,magn_star_night,max_mig_dur,Va_mps,...
      min_ns_fl_hs,max_ns_fl_hs,min_init_fl_hs,max_init_fl_hs, ...
      min_ns_sgnp_hs, max_ns_sgnp_hs, ...
      reOrWtrOpt,reOrWtrStr,n_wtr_stps,rand_wtr_fact, ...
      reOrReComp,dtrCstOpt,dtrCstStr, ...
      short_detrs,kmsAheadDtrCst,d_stop_max,dtrReComp,maxDtrSteps, ...
      stopCstOpt,fracStopCst,stpCstWtrFlt,stpovCstWtrFlt, ...
      endog_comp,calibr_comp,std_err_calibr,sun_head_reset, ...
      tc_comp_reset,tc_comp_wtr,sign_clockw_sun,geo_mag_refl_sun_opt, ...
      perf_init_calibr,magcl_reset_opt,magcl_opt,geomagn_model, ...
      max_magn_shift,min_std_inh_hd,max_std_inh_hd, ...
      min_std_inh_steer,max_std_inh_steer, ...
      std_err_incl_deg,std_rel_err_inten, ...
      min_std_inh_sp,max_std_inh_sp,opt_no_evo_std_yrs,wt_inher_speed, ...
      offset_dep_dusk,var_fl_hrs,stop_durn,stop_durn_sp,endur_flt, ...
      min_fl_stp_hs,max_fl_stp_hs,min_init_fl_stp_hs, ...
      req_land_init,req_land_stop,req_land_arr,incl_wind,wind_str,switch_wind_opt, ...
      dep_wtr_TWs,min_disp_d_dep_km,max_disp_d_dep_km,min_disp_d_arr_km,max_disp_d_arr_km, ...
      min_disp_d_zug_km,max_disp_d_zug_km,min_std_inh_disp_km,max_std_inh_disp_km, ...
      fwd_back_no_time, daily_surv,daily_surv_barrn,daily_surv_wtr, ...
      daily_surv_elev,thrsh_snow_barrn, thrsh_elev, ...
      thrsh_ndvi_breed,thrsh_ndvi_stop_min,thrsh_ndvi_stop_max, ...
      thrsh_trees_breed,thrsh_trees_stop,thrsh_conif_breed, ...
      thrsh_low_veg_breed,thrsh_low_veg_stop,thrsh_all_veg_stop, ...
      max_trees_stop,min_any_veg_stp,min_n_fls_stp,lost_Lat,n_par_pools, ...
      wt_natal_intns,wt_natal_incln,wt_natal_decln, ...
      opt_same_loc_brd,accel_wrmp,arr_EU_opt);
toc
                 
disp_summry_stats

if plot_map_opt == 1

    % plot options
    plot_zug = false;
    plot_zugs = false;
    plot_Zug_ellipse = false;
    plot_fails = true;
    plot_ends = true;
    plot_opt=1;

    plot_trajs_evo
    plot_heads_zugs_connect
    
end

if save_opt == 1
    
        [status, msg, msgID] = mkdir(dir_nm);

         allvars = whos;

        % Identify the variables that ARE NOT graphics handles. This uses a regular
        % expression on the class of each variable to check if it's a graphics object
        tosave = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));
%         
        zusp_str = '';      
        
        if max_std_inh_hd >= 1 && max_std_inh_sp >= 1
            max_std_str = [num2str(round(max_std_inh_hd)) '_'  num2str(round(max_std_inh_sp))];
        elseif max_std_inh_hd < 1 && max_std_inh_sp >= 1

            max_std_str = ['0' num2str(round(10*max_std_inh_hd)) '_' num2str(round(max_std_inh_sp))];

        elseif max_std_inh_hd >= 1 && max_std_inh_sp < 1

            max_std_str = [num2str(round(max_std_inh_hd)) '_00' num2str(round(100*max_std_inh_sp))];

        else

            max_std_str = ['0' num2str(round(10*max_std_inh_hd)) '_00' num2str(round(100*max_std_inh_sp))];

        end


        save([dir_nm cont_str{contn} strat_string err_string '_N_' ...
        num2str(N_inds/1000) '_natSig_' wt_natl_str '_n_zugs_'  ...
        '_' num2str(n_zugs) '_' num2str(n_req_zugs) '_Ninit_' N_init_str fwd_bck_str '_R_' ...
        num2str(round(R_natal)) '_' num2str(std_R_natal) '_inher_' endog_str  calibr_str ...
        '_sig_' num2str(round(min_std_inh_hd)) '_' max_std_str ... 
        opt_std_yr_str accel_wrm_str '_' num2str(10*min_hs_fact) '_' num2str(hs_tr_Atl) req_lnd_str], ...
        allvars(tosave).name)
    
end

delete(gcp('nocreate'))
