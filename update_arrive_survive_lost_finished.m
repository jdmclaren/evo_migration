  % check if survived, arrived at next step
  
  % cumul steps_elev_fact are currently half-day steps_elev_fact
  
  try

surv_fuel = cum_hrs_flt <= curr_max_ns_fl_hs; 
surv_wtr = surv_wtr & ...
    (~over_water | (over_water.*(rand(N_inds,1) <=  daily_surv_wtr))); 

surv_elev = surv_elev & ...
    (~is_high_elev | is_high_elev.*rand(N_inds,1) <= daily_surv_elev.^max(d_time,1));

surv_barrn = surv_barrn & ...
   (~is_barrn | is_barrn.*rand(N_inds,1) <= daily_surv_barrn.^max(d_time,1));

is_std_land = ~(is_high_elev | over_water | is_barrn);

surv_land = surv_land & ... 
    (~is_std_land | is_std_land.*rand(N_inds,1) <= daily_surv.^max(d_time,1));

survived = survived & surv_fuel & surv_wtr & surv_land & surv_barrn & surv_elev; 

  % search depends on longitude (in case over dateline)
  
  % Assume here that migrat finds the coast of the topover for entire
  % flight (backtracking if need be), e.g., using signpost trigger when
  % over water

  catch
      keyboard
  end
  
 if nStop_polys > 0 && lst_req_stp > 0

     for iSt = 1:nStop_polys
         % check if stopped in obligatory staging region
%              passed_Stops(~finished,ist) =  passed_Stops(~finished,ist) | ...
%                  any(((req_land_stop & ~over_wtr_fracs(~finished,:)) | ~req_land_stop) &  ...
%                  llamdas_fracs(~finished,:) >= Stop_verts{ist}(1) & llamdas_fracs(~finished,:) <= Stop_verts{ist}(2) & ...
%                   thetas_fracs(~finished,:) >= Stop_verts{ist}(3) & thetas_fracs(~finished,:) <= Stop_verts{ist}(4),2);

            passed_Stops(~finished,iSt) =  passed_Stops(~finished,iSt) | ...
                 any(((req_land_stop & ~over_wtr_fracs(~finished,:)) | ~req_land_stop) &  ...
                inpolygon(llamdas_fracs(~finished,:),thetas_fracs(~finished,:), ...
                Lon_stop_poly{iSt}',Lat_stop_poly{iSt}'),2); 

             new_stop = find(passed_Stops(:,iSt) & ~curr_passed_Stops(:,iSt));

%               if numel(new_stop) > 0
%                   keyboard
%               end
     end 

     if ~isempty(new_stop) > 0
          for inw = 1:numel(new_stop)
              idx_new = find( ...
                 ((req_land_stop & ~over_wtr_fracs(new_stop(inw),1:9)) |  ~req_land_stop) & ...
                 inpolygon(llamdas_fracs(new_stop(inw),1:9),thetas_fracs(new_stop(inw),1:9), ...
                Lon_stop_poly{iSt}',Lat_stop_poly{iSt}'),1,'last');
            
            if isempty(idx_new)
                
                idx_new = 10;
                
            end
% try
              llamda(new_stop(inw)) = llamdas_fracs(new_stop(inw),idx_new);
              theta(new_stop(inw)) = thetas_fracs(new_stop(inw),idx_new);
% catch
%     keyboard
% end
              passed_Stop_nr(new_stop(inw)) = passed_Stop_nr(new_stop(inw)) +1;

          end
     end     

%      min_Lon = min(mod(Arr_verts(1),2*pi),mod(Stop_verts{end}(1),2*pi));
%      max_Lon = max(mod(Stop_verts{end}(2),2*pi),mod(Arr_verts(2),2*pi));

%  else
% 
%      new_arr = [];

%      min_Lon = mod(Arr_verts(1),2*pi);
%      max_Lon = mod(Arr_verts(2),2*pi);

 end

%      arrived(~finished) = any(((req_land_arr & ~over_wtr_fracs(~finished,:)) |  ~req_land_arr) & ...
%           llamdas_fracs(~finished,:) >= Arr_verts(1) & llamdas_fracs(~finished,:) <= Arr_verts(2) & ...
%           thetas_fracs(~finished,:) >= Arr_verts(3) & thetas_fracs(~finished,:) <= Arr_verts(4),2); % 
%  for iAr = 1:nArr_polys
     arrived(~finished) = ...
       any(((req_land_arr & ~over_wtr_fracs(~finished,:)) |  ~req_land_arr) & ...
        inpolygon(llamdas_fracs(~finished,:),thetas_fracs(~finished,:), ...
          Lon_arr_poly',Lat_arr_poly'),2); 
%  end  % arrived(~finished) | ... {iAr}
 new_arr = find(arrived & ~curr_arr);

 if ~isempty(new_arr) > 0
      for inw = 1:numel(new_arr)
          idx_new = find( ...
             ((req_land_arr & ~over_wtr_fracs(new_arr(inw),:)) | ~req_land_arr) & ...
             inpolygon(llamdas_fracs(new_arr(inw),:),thetas_fracs(new_arr(inw),:), ...
             Lon_arr_poly',Lat_arr_poly'),1,'last');
%          try
          llamda(new_arr(inw)) = llamdas_fracs(new_arr(inw),idx_new);
          theta(new_arr(inw)) = thetas_fracs(new_arr(inw),idx_new);
%          catch
%              keyboard
%          end
      end
 end

 % different "overshot" lons for SE and SW routes
%  diff_Lons_SW = arr_lon_lim-mod(llamda,2*pi); % 
%  diff_Lons_SE = max_Lon_Arr-mod(llamda,2*pi); % 
% 
%  diff_Lons_SE = mod(diff_Lons_SE + pi,2*pi) - pi;
%  diff_Lons_SW = mod(diff_Lons_SW + pi,2*pi) - pi;

%  overshot_Lon = (SE & diff_Lons_SE > arrLon_width2) | ...
%                 (SW & diff_Lons_SW > arrLon_width2);

 overshot_Lon = ...
     (pos_d_0_lon_lim & (mod(llamda,2*pi) - arr_lon_lim > arrLon_width2*2)) | ...
     (~pos_d_0_lon_lim & (arr_lon_lim - mod(llamda,2*pi) > arrLon_width2*2));
            
 overshot_Lat =  (min_Lat_Arr - theta) > arrLat_width2 | ... 
                (theta - max_Lat_Dep) > depLat_width2;
            
 lost = (abs(theta) > lost_Lat*deg2rad) | overshot_Lat; % | overshot_Lon; % 
 
%  diff_Arr1 = abs(Arr_verts(1)-llamda); % mod(Arr_verts(1)-llamda,pi);
%  diff_Arr2 = abs(Arr_verts(2)-llamda); % mod(Arr_verts(2)-llamda,pi);

% past_arr_area = (min(diff_Arr1,diff_Arr2) > arrLon_width2) ...
%     & (sign(mean_arr_Lon-mod(llamda,2*pi)) == -sign_mean_arr_dep_Lon);
% 
% lost = (abs(theta) > 88*deg2rad) | ... % | (theta-theta_0) > arrLat_width |  ...
%  (min_Lat_Arr - theta) > arrLat_width2 | (wide_EW_route & past_arr_area);

  is_late = is_late | ...
      (~finished & (date_jul - date_jul_start) > max_mig_dur);

  finished = lost | ~survived | arrived | is_late;

  idx_not_fin = find(~finished);
%   n_not_fin = numel(idx_not_fin);
  
%    n_late = sum((date_jul - date_jul_start) > max_mig_dur);
   
%    if n_late > 0
%        keyboard
%    end
   
%   fifi = 2;

%           if sum(arrived) > 1 % N_inds/2
%               keyboard
%           end

    % update indices for ndvi and elev
%     lat_deg_idx = 91 + floor(-180/pi*theta);
%     lon_deg_idx = min(360,floor(180/pi*llamda) + 181);