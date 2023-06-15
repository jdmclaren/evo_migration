blats = [blats theta];
blons = [blons llamda];
all_geo_heads = [all_geo_heads mod(alpha*rad2deg + 180,360)];
siz_ag = size(all_geo_heads,2);
% if  siz_ag > 1
% 
%                   if any(all_geo_heads(:,siz_ag) > 275 & ...
%                           all_geo_heads(:,siz_ag-1) < 180)
%                       keyboard
%                   end
% 
% end
all_mag_heads = [all_mag_heads mod((alpha - decln)*rad2deg + 180,360)];
% all_sun_heads = [all_sun_heads mod(alpha*rad2deg - ...
%   (circ_median(sun_az(1:min(10000,N_inds)))*sign_clockw_sun + ...
%    2*(refl_sun_mag)*decln)*rad2deg ...
%    + 180,360)];
all_sun_i = absDiffDeg(alpha*rad2deg,sun_az*rad2deg);       
all_sun_i =  mod(all_sun_i + 180,360) - 180; 
all_sun_heads = [all_sun_heads all_sun_i];
% all_sun_heads = [all_sun_heads mod(alpha*rad2deg - ...
%   (circ_median(sun_az(1:min(10000,N_inds)))*sign_clockw_sun + ...
%    2*(refl_sun_mag)*decln)*rad2deg ...
%    + 180,360)];
% incln_i = NaN*one_vec;
% try
% incln_i(~finished) = incln;
% catch
%    keyboard
% end
%                decln_i = NaN*one_vec;
%                decln_i(~finished) = decln;
% tot_i = NaN*one_vec;
% tot_i(~finished) = tot_inten;
all_inclns = [all_inclns incln*rad2deg];
all_declns = [all_declns decln*rad2deg];
all_tots = [all_tots tot_inten*rad2deg];
all_sun_az = [all_sun_az sun_az*180/pi+180];
all_dates = [all_dates date_jul]; % date_nrs];
all_fl_hs = [all_fl_hs fl_hrs];
% all_cv_dsks = [all_cv_dsks civ_dusk];
% all_cv_dwns = [all_cv_dwns civ_dawn];
all_ovr_wtr = [all_ovr_wtr over_water];
all_barrn = [all_barrn is_barrn];

all_stp_ovs = [all_stp_ovs stops_over];
all_stp_nrs = [all_stp_nrs stop_num];
all_zg_nrs = [all_zg_nrs zug_nr];
all_cum_hrs_flt = [all_cum_hrs_flt cum_hrs_flt];
all_max_hrs_flt = [all_max_hrs_flt curr_max_ns_fl_hs];

all_survs = [all_survs survived];
all_fins = [all_fins finished];
all_arrs = [all_arrs arrived];
all_losts = [all_losts lost];