%% initialize dep. stop and Zug locations
[lat_bs_deps, lon_bs_deps, idx_deps] = ...
    initialize_locs(N_inds,rad2deg*Lon_dep_poly, ...
    rad2deg*Lat_dep_poly,req_land_init, ...
    hi_elevs, barrens, false(size(barrens)), ok_veg_breed); % 
%             lon_bs_deps = shiftAnglesFromMinus180To180(lon_bs_deps*rad2deg)*deg2rad;
theta_0 = shiftAnglesFromMinus180To180(lat_bs_deps)*deg2rad;
llamda_0 = lon_bs_deps*deg2rad;    

for iDep = 1:nDep_polys

    N_inds_DepPs(iDep) = sum(idx_deps == iDep);
    
end

% If multiple choices in Arr poly, 
% find closest poly for each departure Lon 

lat_bs_stops = NaN*one_vec;
lon_bs_stops = NaN*one_vec;
lat_bs_arrs = NaN*one_vec;
lon_bs_arrs = NaN*one_vec;

% if require_stopover_zug % ~isempty(req_stps)

%         Lon_stop_degs = cellfun(@(x) x*rad2deg,Lon_stop_poly{1}, ...
%             'UniformOutput',false);

       % If multiple choices in 1st stop, and zug required at stop,
       % find closest stopover option to mean of departure and arrival Lons 
       % for each individual each individual
% if nSt_Ply_1 > 1

   if zgknk_req_in_stop(1) == 1 || zgknk_req_stop_warm 

      % match dep locs to stop Polys
      dLon_Sts = abs(mn_Lon_St_1-lon_bs_deps);
      [~, idx_stops] = min(dLon_Sts,[],2);
      
      for iSt = 1:nSt_Ply_1
        idx_iSt = (iSt-1)*6+1:(iSt-1)*6+6;
        [lat_bs_stops(idx_stops==iSt), lon_bs_stops(idx_stops==iSt)] = ...
        initialize_locs(sum(idx_stops==iSt),rad2deg*Lon_stop_poly{1}(idx_iSt), ...
                rad2deg*Lat_stop_poly{1}(idx_iSt),req_land_stop, ...
        hi_elevs, barrens, poor_stops, true(size(ok_veg_breed))); %
      end
      
      % then match Arr polys to stops 
     dLon_Arrs = abs(mn_Lon_Arr_deg-lon_bs_stops);
%       dLon_Arrs = abs(mn_Lon_Arr_deg-mn_Lon_St_1(idx_stops)');
     [~, idx_arrs] = min(dLon_Arrs,[],2);

     for iArr = 1:nArr_polys
         idx_iA = (iArr-1)*6+1:(iArr-1)*6+6;
        [lat_bs_arrs(idx_arrs==iArr), lon_bs_arrs(idx_arrs==iArr)] = ...
            initialize_locs(sum(idx_arrs==iArr),rad2deg*Lon_arr_poly(idx_iA), ...
            rad2deg*Lat_arr_poly(idx_iA),req_land_arr, ...
            hi_elevs, barrens, poor_stops, true(size(ok_veg_breed))); % 
      end     
      
   else

     % match dep locs to Arr locs 
     dLon_Arrs = abs(mn_Lon_Arr_deg-lon_bs_deps);
     [~, idx_arrs] = min(dLon_Arrs,[],2);

     for iArr = 1:nArr_polys
         idx_iA = (iArr-1)*6+1:(iArr-1)*6+6;
        [lat_bs_arrs(idx_arrs==iArr), lon_bs_arrs(idx_arrs==iArr)] = ...
            initialize_locs(sum(idx_arrs==iArr),rad2deg*Lon_arr_poly(idx_iA), ...
            rad2deg*Lat_arr_poly(idx_iA),req_land_arr, ...
            hi_elevs, barrens, poor_stops, true(size(ok_veg_breed))); % 
      end

      mn_dep_arr_Lons =(lon_bs_deps + lon_bs_arrs)*pi/360; %c 
      dLon_Sts = abs(mn_Lon_St_1-mn_dep_arr_Lons); 
      [~, idx_stops] = min(dLon_Sts,[],2); 
      
      for iSt = 1:nSt_Ply_1
        idx_iSt = (iSt-1)*6+1:(iSt-1)*6+6;
        [lat_bs_stops(idx_stops==iSt), lon_bs_stops(idx_stops==iSt)] = ...
        initialize_locs(sum(idx_stops==iSt),rad2deg*Lon_stop_poly{1}(idx_iSt), ...
                rad2deg*Lat_stop_poly{1}(idx_iSt),req_land_stop, ...
        hi_elevs, barrens, poor_stops, true(size(ok_veg_breed))); %
      end

   end


%% initialize variation in inerited parameters

std_inher_disp = min_std_inh_disp + ...
    (max_std_inh_disp-min_std_inh_disp)*rand([N_inds,1]);

std_inher_head = min_std_inh_hd + ...
    (max_std_inh_hd-min_std_inh_hd)*rand([N_inds,1]);
inher_head_kappa = 1./std_inher_head.^2;

std_inher_steer = min_std_inh_steer + ...
(max_std_inh_steer - min_std_inh_steer)*rand([N_inds,1]);

% add variance to start date
dev_day = round(std_day_start*randn([N_inds 1]));
dev_day(dev_day > max_dev_day_st) = max_dev_day_st;
dev_day(dev_day < -max_dev_day_st) = -max_dev_day_st;
day_starts = day_start + dev_day; 
% curr_date.day = day_starts;
date_jul = datenum(curr_year*one_vec,month_start*one_vec,day_starts);
date_jul_start = date_jul;
doys = day(datetime(datevec(date_jul)),'dayofyear' );

% initialize sun az and geomagn field given locs and start (dep) dates
% initialize_sun_fl_hrs_geomag 
init_sun_geomag_fl_hrs

% first check if zug_opt == 1 to set zugknik magn signposts

if zug_opt == 1

    std_inher_signp = min_std_inh_sp + ...
    (max_std_inh_sp - min_std_inh_sp)*rand([N_inds,1]);
    inher_signp_kappa = 1./std_inher_signp.^2;

    % initial error magnetic compass parameters
    if ~isinf(compass_kappa_incl)
         rand_magn_a = vmrand(0, compass_kappa_incl, [N_inds 1]); %      
    else
         rand_magn_a = zero_vec;
    end

    % determine inclination at arrival (use median date, for
    % guesstimate)
    
    if  zgknk_req_in_stop(1) == 1 || zgknk_req_stop_warm % (zug_opt == 1 &&  n_zugs > 0) %  || 
        
         [Bxa, Bya, Bza] = igrf(round(median(date_jul)), ...
            lat_bs_stops, lon_bs_stops, 0); % lat_bs_arrs, lon_bs_arrs, 0); %  
         
    else
        
         [Bxa, Bya, Bza] = igrf(round(median(date_jul)), ...
              lat_bs_arrs, lon_bs_arrs, 0); % lat_bs_stops, lon_bs_stops, 0); % 
         
    end

    % decln_0 used to set up mean initial (first-year) offset 
    % in initial heads gicen the local declination field
    decln_a = atan2(Bya,Bxa);

    hz_inten_a = hypot(Bxa,Bya);
    incln_a = min(atan(Bza./hz_inten_a) + rand_magn_a,pi/2);

    if zug_signp == 0                   

        zug_signps =  NaN*ones(N_inds,n_zugs);
%                     zug_signps = sort(zug_signps,2,'descend');

    elseif zug_signp == 1

        if  zgknk_req_in_stop(1) == 1 || zgknk_req_stop_warm % (zug_opt == 1 &&  n_zugs > 0) %  || 

            zug_signps = incln_a; % + rand(N_inds,n_zugs).*(incln_0 - incln_a);
            
        else
            
            zug_signps = incln_a  + rand(N_inds,n_zugs).*(incln_0 - incln_a);
              
        end
        
        zug_signps = sort(zug_signps,2,'descend'); 

    elseif zug_signp == 2

%                     tot_B_a = hypot(hz_inten_a,Bza).*(1 + randn(N_inds,1)*std_rel_err_inten);
%                 zug_signps = rand(N_inds,n_zugs).*(1-tot_B_a./tot_inten_0); 

        if zgknk_req_in_stop(1) == 1 || zgknk_req_stop_warm % (zug_opt == 1 &&  n_zugs > 0) %  || 
            
             zug_signps = decln_a; 
             % Bza./tot_B_a + rand(N_inds,n_zugs).*(vt_inten_0 - Bza./tot_B_a);
             
        else
            
             zug_signps = decln_a + rand(N_inds,n_zugs).*(decln_0 - decln_a); 
             
        end 
        
        
        if mean(decln_0) > mean(decln_a)

            zug_signps = sort(zug_signps,2,'descend');   

        else

             zug_signps = sort(zug_signps,2,'ascend');                          

        end

    elseif zug_signp == 3
        
        tot_B_a = hypot(hz_inten_a,Bza).*(1 + randn(N_inds,1)*std_rel_err_inten);

        if  zgknk_req_in_stop(1) == 1 || zgknk_req_stop_warm % (zug_opt == 1 &&  n_zugs > 0) %  || 

            zug_signps = tot_B_a; %  + rand(N_inds,n_zugs).*(tot_inten_0 - tot_B_a); % Bza./tot_B_a + rand(N_inds,n_zugs).*(vt_inten_0 - Bza./tot_B_a);
            
        else
            
            zug_signps = tot_B_a  + rand(N_inds,n_zugs).*(tot_inten_0 - tot_B_a); % Bza./tot_B_a + rand(N_inds,n_zugs).*(vt_inten_0 - Bza./tot_B_a);
              
        end
        
        zug_signps = sort(zug_signps,2,'descend');


    elseif zug_signp == 4  

        vt_inten_a = Bza.*(1 + randn(N_inds,1)*std_rel_err_inten);
        if  zgknk_req_in_stop(1) == 1 || zgknk_req_stop_warm % (zug_opt == 1 &&  n_zugs > 0) %  || 
                zug_signps = vt_inten_a; % + rand(N_inds,n_zugs).*(vt_inten_0 - vt_inten_a);
        else
                zug_signps = vt_inten_a + rand(N_inds,n_zugs).*(vt_inten_0 - vt_inten_a);
        end
%                 zug_signps = Bza./tot_B_a + rand(N_inds,n_zugs).*(vt_inten_0 - Bza./tot_B_a);
        zug_signps = sort(zug_signps,2,'descend');
        
    elseif zug_signp == 5  

        hz_inten_a = hz_inten_a.*(1 + randn(N_inds,1)*std_rel_err_inten);
        
        if  zgknk_req_in_stop(1) == 1 || zgknk_req_stop_warm % (zug_opt == 1 &&  n_zugs > 0) %  || 
            
            zug_signps = hz_inten_a;
            
        else

            zug_signps = hz_inten_a + rand(N_inds,n_zugs).*(hz_inten_0 - hz_inten_a);
                        
        end
        
%                 zug_signps = Bza./tot_B_a + rand(N_inds,n_zugs).*(vt_inten_0 - Bza./tot_B_a);
        zug_signps = sort(zug_signps,2,'descend');


    end

else

    zug_signps = NaN*ones(N_inds,n_zugs);
    std_inher_signp = NaN*ones(N_inds,n_zugs);

end

% next set up the initial heads, given the departure arrival areas
% as well as any required stopover locations (and initial
% declination field which could be high e.g. in the Arctic)

              
%             if calibr_comp(1)  == 4 % use 90 to 270


if zgknk_req_in_stop(1) == 1  || zgknk_req_stop_warm %% (zug_opt == 1 &&  n_zugs > 0) %  || 
% nStop_polys > 0

    angs_dep_stop_rh = pi/180*azimuth('rh',lat_bs_deps,mod(lon_bs_deps,360), ...
          lat_bs_stops,mod(lon_bs_stops,360));
    angs_dep_stop_gc = pi/180*azimuth('gc',lat_bs_deps,mod(lon_bs_deps,360), ...
         lat_bs_stops,mod(lon_bs_stops,360));

%     mn_angs_arr = angs_dep_stop_rh; % (angs_dep_stop_rh ...
%         +  (calibr_comp(2)~=0)*angs_dep_stop_gc)/(1+(calibr_comp(2)~=0));
    
      mn_angs_arr =  circ_mean([angs_dep_stop_rh angs_dep_stop_gc]')'; %circ_mean([angs_dep_stop_rh' angs_dep_stop_gc']');
  
      sign_mns = -sign(lon_bs_stops-lon_bs_deps); 
     large_mns = abs(lon_bs_arrs-lon_bs_deps) > 90;      
     
      %     var_ang_arr = circ_var([angs_dep_stop_rh' angs_dep_stop_gc']');
% 
%     % use 10 times the variance for offsets to means
%     init_offs_lots = vmrand(0, 1/(10*var_ang_arr), [N_inds*10 1]); % rand(N_inds,1)*2*pi;   mn_angs_arr      

%     [~, std_ang_arr] = circ_std([angs_dep_stop_rh' angs_dep_stop_gc']);
    std_ang_arrs = max(absDiffRad(mn_angs_arr,angs_dep_stop_rh),absDiffRad(mn_angs_arr,angs_dep_stop_gc)); % circ_std([angs_dep_stop_rh' angs_dep_stop_gc']);
    std_ang_arr = quantile((std_ang_arrs),0.75);    % use 10 times the variance 

    init_offs_lots = vmrand(0, 1/(std_ang_arr^2), [N_inds*10 1]); % rand(N_inds,1)*2*pi;    mn_angs_arr  

else

    angs_dep_arr_rh = pi/180*azimuth('rh',lat_bs_deps,mod(lon_bs_deps,360), ...
      lat_bs_arrs,mod(lon_bs_arrs,360));
    angs_dep_arr_gc = pi/180*azimuth('gc',lat_bs_deps,mod(lon_bs_deps,360), ...
      lat_bs_arrs,mod(lon_bs_arrs,360));

%   mn_angs_arr = angs_dep_arr_rh; % 
   mn_angs_arr = circ_mean([angs_dep_arr_rh angs_dep_arr_gc]')';%  angs_dep_arr_gc; %  angs_dep_arr_rh; %

   sign_mns = -sign(lon_bs_arrs-lon_bs_deps);
   large_mns = abs(lon_bs_arrs-lon_bs_deps) > 90;
   
%     [~, std_ang_arr] = circ_std([angs_dep_arr_rh angs_dep_arr_gc]')';
    std_ang_arrs = max(absDiffRad(mn_angs_arr,angs_dep_arr_rh),absDiffRad(mn_angs_arr,angs_dep_arr_gc)); % circ_std([angs_dep_stop_rh' angs_dep_stop_gc']);
    std_ang_arr = quantile((std_ang_arrs),0.75);    % use 10 times the variance 

    % use 10 times the variance 
    init_offs_lots = vmrand(0, 1/(std_ang_arr^2), [N_inds*10 1]); % rand(N_inds,1)*2*pi;    mn_angs_arr           

end

%% adjust mean angles for cases where longitude changes sign (need to
% enter abs(lons) > 180 for this to work
% Don't forget mn_angs_arr is clockwise from N at this point
mn_angs_arr(large_mns & sign_mns == 1 & mn_angs_arr > 0) = ...
    -mn_angs_arr(large_mns & sign_mns == 1 & mn_angs_arr > 0);

mn_angs_arr(large_mns & sign_mns == -1 & mn_angs_arr < 0) =  ...
    -mn_angs_arr(large_mns & sign_mns == -1 & mn_angs_arr < 0);

% add offsets for compass strategy to assess if init heads are
% 'on track' within +/- 90 degs
inher_offs_lots = mod(init_offs_lots - ...
    (add_decl_inher)*median(decln_0) - ...
    (add_az_inher)*(circ_median(sun_az_0(1:min(10000,N_inds)))*sign(sign_clockw_sun) + ...
        2*(sign_clockw_sun==-1 && geo_mag_refl_sun_opt == 2)*median(decln_0))...
        + pi,2*pi) -pi;
    

%     absDiff = min(360-absDiffDeg(init_offs_lots*rad2deg,0), ...
%     absDiffDeg(init_offs_lots*rad2deg,0));
%     inh_offs_fun = inher_offs_lots; % will modulate and set vs. 180 = South below
   within_bounds = abs(init_offs_lots) <  pi; %(mig_sys~=7)*90 + (mig_sys==7)*15; %  > pi/2;    

% end

sign_mn_ang =  -sign(mn_angs_arr); % -pi     

 coeff_magn_steer = zeros(N_inds,n_zugs+1);      

%     inher_offs = randsample(inh_offs_fun,N_inds,true);
    inher_offs = randsample(init_offs_lots(within_bounds),N_inds,true);
% end

clear init_offs_lots inher_offs_lots within_bounds absDiff inh_hd_fun

%             inher_heads(inher_heads <0) = inher_heads(inher_heads <0) +2*pi;
inher_heads = mod(mn_angs_arr + inher_offs-pi+pi,2*pi)-pi;

if zgknk_req_in_stop(1) == 1  || zgknk_req_stop_warm %% (zug_opt == 1 &&  n_zugs > 0) %  || 
% nStop_polys > 0 % &&  zug_opt ~= 2
    for iz = 1:max(min(nStop_polys,n_zugs),1)

        % set last stopover as "previous"
        lon_bs_prev = lon_bs_stops;    
        lat_bs_prev = lat_bs_stops;            

        % if TC sun comp and no magn set zugkn head to endog head
        if ~(calibr_comp(2)==2 && tc_comp_reset == 0 && endog_comp == 2) && ...
                (require_stopover_zug(iz) == 1 || (zug_opt == 1 && iz <= n_zugs)) % ismember(iz,req_stps))

            if iz == nStop_polys

                angs_stop_arr_rh = pi/180*azimuth('rh',lat_bs_prev,mod(lon_bs_prev,360), ...
                      lat_bs_arrs,mod(lon_bs_arrs,360));
                angs_stop_arr_gc = pi/180*azimuth('gc',lat_bs_deps,mod(lon_bs_deps,360), ...
                      lat_bs_arrs,mod(lon_bs_arrs,360));   

            else 

                % update to next stopover
                [lat_bs_stops, lon_bs_stops] = ...
                    initialize_locs(N_inds,rad2deg*Lon_stop_poly{iz+1}, ...
                            rad2deg*Lat_stop_poly{iz+1},req_land_stop, ...
                    hi_elevs, barrens, poor_stops, true(size(ok_veg_breed))); %
%                [lat_bs_stops, lon_bs_stops] = ...
%                 initialize_locs(N_inds,rad2deg*Stop_verts{iz+1},req_land_stop, ...
%                 hi_elevs, barrens, poor_stops, true(size(ok_veg_breed))); %

                angs_stop_arr_rh = pi/180*azimuth('rh',lat_bs_prev,mod(lon_bs_prev,360), ...
                      lat_bs_stops,mod(lon_bs_stops,360));
                angs_stop_arr_gc = pi/180*azimuth('gc',lat_bs_prev,mod(lon_bs_prev,360), ...
                      lat_bs_stops,mod(lon_bs_stops,360));                    

            end
     
            mn_angs_arr =  circ_mean([angs_stop_arr_rh angs_stop_arr_gc]')';%angs_stop_arr_gc; % circ_mean([angs_dep_arr_rh'; angs_dep_arr_gc'])';%  

%             mn_angs_arr = (calibr_comp(2)~=2)*circ_mean(angs_stop_arr_rh) ...
%                 + (calibr_comp(2)==2)*circ_mean(angs_stop_arr_gc); % circ_mean([angs_dep_arr_rh' angs_dep_arr_gc']');
%            [~,var_ang_arr] =  circ_var([angs_stop_arr_rh' angs_stop_arr_gc']'); % pi; %
           
            std_ang_arrs = max(absDiffRad(mn_angs_arr,angs_stop_arr_rh),absDiffRad(mn_angs_arr,angs_stop_arr_gc)); % circ_std([angs_dep_stop_rh' angs_dep_stop_gc']);
            std_ang_arr = quantile((std_ang_arrs),0.75);  
         
%             [~, std_ang_arr] = circ_std([angs_stop_arr_rh' angs_stop_arr_gc']');

            % to get ballpark range for inherited
            % zugkn angle (i.e. rel to magn or sun az) 
            % guesstimate offset for heading from inherited
            % angle

            % first estimate the decln and sun az at
            % stopover. for simplicity keep dates same as first
            % departure (may be light bias but ang spread
            % should catch proper angles)
            
            udy = round(mean(date_jul)); 
            [Bx_stop, By_stop,~] = igrf(udy,lat_bs_stops,lon_bs_stops, ...
                zero_vec);
       
                decln_stop = atan2(By_stop,Bx_stop);
                
%             end

            if ~((calibr_comp(1) == 3 && magcl_opt == 1) || ...
                    (calibr_comp(1) == 2 && magcl_opt == 2))   
                
%                 zgkn_lots = vmrand(mn_angs_arr, 1/(10*var_ang_arr), [N_inds*10 1]); %
                zugkn_inher_heads(:,iz) = vmrand(mn_angs_arr, 1/(std_ang_arr^2)); % 
                % assume same sun azimuth as a ballpark guess
                if incl_sun
%                                 [sun_az_stop,~] = SolarAzEl(dates,lat_bs_stops,lon_bs_stops,zero_vec);
                    sun_az_stop = sun_az_0; % sun_az_stop*deg2rad -pi;
                else
                    sun_az_stop = zero_vec;
                end


                % add offsets for compass strategy to assess if init heads are
                % 'on track' within +/- 90 degs
                if add_az_inher
                    
                    zugkn_inher_heads(:,iz) = zugkn_inher_heads(:,iz) - ...
                    (add_decl_inher == 1)*median(decln_stop) - ...
                    (add_az_inher)*(circ_median(sun_az_stop)*sign(sign_clockw_sun) + ...
                        2*(sign_clockw_sun==-1 && geo_mag_refl_sun_opt == 2)*median(decln_stop));
                    
                else
                    
                    zugkn_inher_heads(:,iz) = zugkn_inher_heads(:,iz) - ...
                    (add_decl_inher == 1)*median(decln_stop);
                    
                end
                
            else

                % magcl transverse (Kiep) or pll inverted,
                % restricted init range works better
                zugkn_inher_heads(:,iz) = randinterval(pi,pi+sign_mn_ang*pi/12,N_inds) ...
                    - (endog_comp == 1)*median(decln_stop); %
%                             zgkn_fun = zgkn_lots;

            end

                zugkn_inher_heads(:,iz) = zugkn_inher_heads(:,iz)-pi;

        else

            if iz == 1

               zugkn_inher_heads(:,iz) = inher_heads;

            else % if ~ismember(iz,req_stps(1:end-1))

               zugkn_inher_heads(:,iz) = zugkn_inher_heads(:,iz-1);

            end

        end

    end

    for iz = iz+1:n_zugs

        zugkn_inher_heads(:,iz) = zugkn_inher_heads(:,iz-1);

    end

    clear zgkn_lots zugkn_init_lots  within_bounds absDiff zgkn_fun

else

    if n_zugs > 0
        for iz = 1:n_zugs
            zugkn_inher_heads(:,iz) = inher_heads(randperm(length(inher_heads)));
        end
    else
        zugkn_inher_heads = inher_heads(randperm(length(inher_heads)));
        
    end

end

zugkn_inher_heads = mod(zugkn_inher_heads-pi,2*pi)-pi;

% define half length selection on breeding grounds 
% from natal distance and winter connectivity
mn_d_dep = max_disp_d_dep/2;
std_d_dep = mn_d_dep/4;

if ~isinf(min_disp_d_dep)
%     disp_d_dep = min_disp_d_dep + rand(N_inds,1)*(max_disp_d_dep-min_disp_d_dep);
    disp_d_dep = min(max(mn_d_dep + randn(N_inds,1)*std_d_dep,min_disp_d_dep),max_disp_d_dep);
else
    disp_d_dep = Inf*one_vec;
end
if ~isinf(min_disp_d_arr)
    disp_d_arr = min_disp_d_arr + rand(N_inds,1)*(max_disp_d_arr-min_disp_d_arr);
else
       disp_d_arr = Inf*one_vec;
end 
if ~isinf(min_disp_d_zug)
    disp_d_zug = min_disp_d_zug + rand(N_inds,1)*(max_disp_d_zug-min_disp_d_zug);
else
       disp_d_zug = Inf*one_vec;
end 