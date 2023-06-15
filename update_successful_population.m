% keep track of field at dep locs
% date_jul_start is initialized in create_init_popn and 
% updated in update_succ_popn

% save last value for selection vs. decln etc. (note is same or 
% random yr in wrmup so better to do that during actual sim)
if  iYear > 1 % Sim_nr >= N_init_Sims/2 % +1 %

    prv_dcln_succ = decln_0s(stoppedAndArrived,max(iYear-1,1));
    prv_incln_succ = incln_0s(stoppedAndArrived,max(iYear-1,1));
    prv_intns_succ = intns_0s(stoppedAndArrived,max(iYear-1,1));

    % proxy day start based on previous year's depart date 
    % with curr_year updated
    date_jul_start = datenum(curr_year*one_vec,month_start*one_vec,day_starts); 

    [Bx0, By0, Bz0] = igrf(round(median(date_jul_start)), ...
      lat_bs_deps, lon_bs_deps, 0); % lat_bs_stops, lon_bs_stops, 0); % 
    % decln_0 used to set up mean initial (first-year) offset 
    % in initial heads gicen the local declination field
    dcln_curr = atan2(By0,Bx0)*180/pi;
    incln_curr = sqrt(Bx0.^2 + By0.^2 + Bz0.^2);
    intns_curr = atan(Bz0./hz_inten_0)*180/pi;  

    % we repopulate all locations, and (optionally) weight
    % % besides geographically also rel to geomagn signp
%     dcln_curr = decln_0s(:,iYear);
%     incln_curr = incln_0s(:,iYear);
%     intns_curr = intns_0s(:,iYear);
     
end

% prv_intns_suc = intns_0s(iYear,stoppedAndArrived);
% prv_incln_suc = incln_0s(iYear,stoppedAndArrived);  

%     prev_decln(stoppedAndArrived);

zugkn_inher_heads = NaN*ones(N_inds,n_zugs); % zugkn_inher_heads(1:N_inds,:);

% keep locations same but take random successful weighted by distance

% N_succ = sum(stoppedAndArrived);

i_trt_1 = NaN*ones(N_inds,1);
i_trt_2 = NaN*ones(N_inds,1);

% try 
    
% loop through departure locations and choose
%  parental traits i_trt_1 and i_trt_2 for each new locn 
parfor i_dep = 1:N_inds % par
%     try
    p_d2_dep = NaN*ones(N_succ,1);
    p_d2_dep_arr = NaN*ones(N_succ,1);
    p_d2_all = NaN*ones(N_succ,1);

    % compute distance from selected successful migrants
%     arclens_deps = distance('gc',[lat_bs_deps(i_dep),lon_bs_deps(i_dep)], ...
%         [lat_deps_succ,lon_deps_succ]);

    % For speed (10 x faster), replace by Cartesian (local planar) distance 
    % since within 1000 km spherical effects
    % are small
    arclens_deps =  sqrt((theta_0(i_dep)-lat_deps_succ).^2 + ...
        (cos(theta_0(i_dep))*(llamda_0(i_dep) -lon_deps_succ)).^2); %
    
    % expnential decay in selection beyond dispersal distance
    % modulated by (any) speed selectivity
    r_exps = arclens_deps./succ_disp_d_dep;
    r_exp_min = min(r_exps);
    rel_ps = speed.*exp(-r_exps.^2+r_exp_min^2)./succ_disp_d_dep; % .^2

    if  iYear > 1 && wt_natl_geomg % Sim_nr > N_init_Sims/2 && wt_natl_geomg % 
        if wt_natal_decln 
            d_decln_i = abs(dcln_curr(i_dep) - prv_dcln_succ);
            rel_ps = rel_ps./(d_decln_i/median(d_decln_i)).^2;
        end
        if wt_natal_incln
            d_incln_i = abs(incln_curr(i_dep) - prv_incln_succ);
            rel_ps = rel_ps./(d_incln_i/median(d_incln_i)).^2;
        end
        if wt_natal_intns
            d_intns_i = abs(intns_curr(i_dep) - prv_intns_succ);
            rel_ps = rel_ps./(d_intns_i/median(d_intns_i)).^2;
        end

    else
        d_dcln_i = NaN;
        d_incln_i = NaN;
        d_intns_i = NaN;
    end   

    pmax = max(rel_ps);
%     rmin = min(r_exps);
%     pmax = max(succ_disp_d_dep);
    p_d2_dep = rel_ps/pmax;
    
%     p_d2_dep(arclens_deps < succ_disp_d_dep) = ...
%         speed(arclens_deps < succ_disp_d_dep);
%     beyond_thr = arclens_deps >= succ_disp_d_dep;
%     
% 
%     p_d2_dep(beyond_thr) = speed(beyond_thr).* ...
%         exp(-(arclens_deps(beyond_thr)./succ_disp_d_dep(beyond_thr) -rmin).^2);
% 
%     % scale rel ps to avoid underflow & all zero probabilities
%     rel_dist_beynd= arclens_deps(beyond_thr)./succ_disp_d_dep(beyond_thr);
%     p_d2_dep(beyond_thr) = speed(beyond_thr).* ...
%         exp(-(rel_dist_beynd./min(rel_dist_beynd) -1).^2);
    
    % choose first "parent" traits based on natal proximity to new location
    i_trt_1(i_dep) = randsample(1:N_succ,1,true,p_d2_dep);

    % choose second based on distance as above and...
    % on distance of winter arrival locn from 1st chosen
            

    % expnential decay in selection with distance
    % modeulated by (any) speed selectivity

    if ~isinf(max_disp_d_arr) % any(reqs_zug)

        arclens_arrs = distance('gc',[lat_arrs_succ(i_trt_1(i_dep)), ...
          lon_arrs_succ(i_trt_1(i_dep))],[lat_arrs_succ,lon_arrs_succ]);       
                
        p_d2_dep_arr(arclens_arrs < succ_disp_d_arr) = ...
            p_d2_dep(arclens_arrs < succ_disp_d_arr);
        
        beyond_thr = arclens_arrs >= succ_disp_d_arr;
                
        % scale rel ps to avoid underflow & all zero probabilities
        rel_dist_beynd= arclens_arrs(beyond_thr)./succ_disp_d_arr(beyond_thr);        
        p_d2_dep_arr(beyond_thr) = p_d2_dep(beyond_thr).* ...
            exp(-(rel_dist_beynd./min(rel_dist_beynd) -1).^2);
        
    else
        
        p_d2_dep_arr = p_d2_dep;
        
    end
                
    % also weight by first zugkn locn if applicable
    if zug_opt == 1 && ~isinf(max_disp_d_zug) % any(reqs_zug)
        
        arclens_zugs = distance('gc',[lat_zugs_succ(i_trt_1(i_dep)), ...
        lon_zugs_succ(i_trt_1(i_dep))],[lat_zugs_succ,lon_zugs_succ]);  
    
        p_d2_all(arclens_zugs < succ_disp_d_zug) = ...
        p_d2_dep_arr(arclens_zugs < succ_disp_d_zug);
    
        beyond_thr = arclens_zugs >= succ_disp_d_zug;
               
        % scale rel ps to avoid underflow & all zero probabilities
        rel_dist_beynd= arclens_zugs(beyond_thr)./succ_disp_d_zug(beyond_thr);        
        p_d2_all(beyond_thr) = p_d2_dep_arr(beyond_thr).* ...
            exp(-(rel_dist_beynd./min(rel_dist_beynd) -1).^2);
        
    else
        
        arclens_zugs = zeros(size(arclens_deps));  
        p_d2_all = p_d2_dep_arr;
        
    end
    
     % add contingency for dep date if uses sun compass
     % note with "speed" we've already accounted for return timing
     % Here we inversely weight by fractional week in
     % diffrence of departure between parents
     if incl_sun

        ddays = abs(day_start_succ - ...
            day_start_succ(i_trt_1(i_dep)))/7 +1;
        p_d2_all = p_d2_all./ddays;

     end
     
    % during warmup 'select' sigma (headings) and thresh dists
    % Option (Aug 2022, now turned off) only for first third of warmup so that it doesn't overfit
    % the narrowness of "adaptive" std in inherited traits
    if accel_wrmp && Sim_nr <= N_wrmp_rnd +  N_wrmp_bck + 1 % N_init_Sims+1 % 
        
        if rng_std_inh_hd > 0
            d_stds = abs(succ_std_inher_head - succ_std_inher_head(i_trt_1(i_dep)));
            p_d2_all = p_d2_all.*(1-d_stds/rng_std_inh_hd).^2; % ...
        end

        if rng_disp_d_dep > 0
            d_deps = abs(succ_disp_d_dep - succ_disp_d_dep(i_trt_1(i_dep)));
            p_d2_all = p_d2_all.*(1-d_deps/rng_disp_d_dep).^2;
        end

        if zug_opt && rng_std_inh_sp > 0  
            d_stds = abs(succ_std_inher_signp - succ_std_inher_signp(i_trt_1(i_dep)));
            p_d2_all = p_d2_all.*(1-d_stds/rng_std_inh_sp).^2; % ...
        end

    end

     % If option not to use (parent) inds from same locn,
     % then exclude same parent to promote diversity and converegence
     if ~opt_same_loc_brd && sum(p_d2_all > 1e-10) > 1 && Sim_nr <=  N_wrmp_rnd +  N_wrmp_bck + 1 % N_init_Sims+1 % 
         p_d2_all(i_trt_1(i_dep)) = 0;
     end

     % now choose 2nd parent traits based on arrival, zug  & sched proximities
     i_trt_2(i_dep) = randsample(1:N_succ,1,true,p_d2_all);
  
 end

% catch
%     
%     keyboard
%     
% endinher_signp

if iWarmup == N_init_Sims +1 % max(n_req_stp+1,2)

    lat_bs_deps = lat_bs_deps(1:N_inds);
    lon_bs_deps = lon_bs_deps(1:N_inds);
    theta_0 = lat_bs_deps*deg2rad;
    llamda_0 = lon_bs_deps*deg2rad; 

    inher_head_kappa = [];
    std_inher_signp = [];
    inher_signp_kappa = [];   

     prev_succ_heads = NaN*one_vec;
     prev_succ_zug_inh_hds = NaN*one_vec;
     prev_speed = NaN*one_vec;

     blat_prev_succs = NaN*one_vec;
     blon_prev_succs = NaN*one_vec; 

     prev_coeff_succ = NaN*one_vec; 

end

% add variance to start date
day_starts = round(0.5*(day_start_succ(i_trt_1) + ...
   day_start_succ(i_trt_2)) + ...
   max(min(randn(N_inds,1)*std_day_start,max_dev_day_st), ...
   -max_dev_day_st));
curr_date.day = day_starts;     

% update variation in inerited headings 
% add "minimum" variability to inherited variability
% If range is zero it will remain unchanged for all inds
% if rng_std_inh_hd > 0
    std_inher_head = max(min(0.5*(succ_std_inher_head(i_trt_1) + ...
       succ_std_inher_head(i_trt_2)) + ...
       randn(N_inds,1)*min_std_inh_hd,max_std_inh_hd), ...
       min_std_inh_hd);

% inher_head_kappa(1:N_inds,1) = (1./succ_std_inher_head(i_trt_1).^2 + ...
%    1./succ_std_inher_head(i_trt_2).^2)/2;
inher_head_kappa(1:N_inds,1) = 1./std_inher_head.^2;
%            inher_head_kappa(1:N_inds,1) = kap_hds_bar.*(1 + randn(N_inds,1)/10);

% if rng_std_inh_disp > 0

    % update variation in inerited headings 
    % add "minimum" variability to inherited variability
    std_inher_disp = max(min(0.5*(succ_std_inher_disp(i_trt_1) + ...
    succ_std_inher_disp(i_trt_2)) + ...
    randn(N_inds,1)*min_std_inh_disp,max_std_inh_disp), ...
    min_std_inh_disp);

     % have now updated std inh disp dist in new pop
     % Use this value to obtain the inherited "mean disp" param
     % (i.e., the "behavioural" tendancy to disperse according to a
     % Gaussian distribution)

% end

    disp_d_dep = max(min(0.5*(succ_disp_d_dep(i_trt_1) + ...
        succ_disp_d_dep(i_trt_2)) + randn(N_inds,1).*std_inher_disp, ...
        max_disp_d_dep),min_disp_d_dep); % (1:N_inds,1)
    
    disp_d_arr = max(min(0.5*(succ_disp_d_arr(i_trt_1) + ...
        succ_disp_d_arr(i_trt_2)) + randn(N_inds,1).*std_inher_disp, ...
        max_disp_d_arr),min_disp_d_dep); % (1:N_inds,1)
    
    disp_d_zug = max(min(0.5*(succ_disp_d_zug(i_trt_1) + ...
        succ_disp_d_zug(i_trt_2)) + randn(N_inds,1).*std_inher_disp, ...
        max_disp_d_zug),min_disp_d_dep); % (1:N_inds,1)


if ~is_magcl || endog_comp ~= 1
    % init heads magn clinic depends on inherited projected angle
    % and will be calculated below

%                 try
    if ~isinf(inher_head_kappa(1))
        inher_heads = vmrand(circ_mean([succ_heads(i_trt_1,1)  succ_heads(i_trt_2,1)]')', ...
            inher_head_kappa,[N_inds,1]);
    else
        inher_heads = circ_mean([succ_heads(i_trt_1,1)  succ_heads(i_trt_2,1)]')';
    end


else % magnclinic, and non-geographic inheritance 
    % (I'm not currently allowing testing of ....
    % sun inheritance and magcl orientn)

    coeff_magn_steer = [];

    % vary inherited projection inclination angle                
    if calibr_comp(1) == 2

        init_incls = (coeff_succ(i_trt_1,:) + coeff_succ(i_trt_2,:))/2;
        proj_incl = atan(init_incls);
%                     try
        for ih = 1:n_zugs+1                                    
            inh_angles_ih = ...
                min(max(vmrand(proj_incl(:,ih), inher_head_kappa, ...
                [N_inds n_zugs+1]),-pi/2+pi/18000),pi/2-pi/18000);
            coeff_magn_steer(1:N_inds,ih) = tan(inh_angles_ih); 
        end


    else % calibr_comp(1) == 3

        proj_incl = atan((coeff_succ(i_trt_1,:) + coeff_succ(i_trt_2,:))/2);
        for ih = 1:n_zugs+1

%                         coeff_magn_steer(1:N_inds,ih) = ...
%                             cos(vmrand(proj_incl(:,ih), inher_head_kappa, [N_inds 1]));
            inh_angles_ih = ...
                min(max(vmrand(proj_incl(:,ih), inher_head_kappa, ...
                [N_inds n_zugs+1]),-pi/2+pi/18000),pi/2-pi/18000);
            coeff_magn_steer(1:N_inds,ih) = tan(inh_angles_ih);                          
        end
        % proj_incl+ randn(N_inds,1)*std_inher_steer);

    end   

    sign_inher_heads = sign(circ_mean([succ_heads(i_trt_1,1) ...
        succ_heads(i_trt_2,1)]')' ...
        - (endog_comp ~= 1)*(decln_0(i_trt_1,1) + decln_0(i_trt_2,1))/2 );

    if calibr_comp(1) == 3
      sign_inher_heads(sign_inher_heads == 0) = sign_mn_ang;
    end

%                 end

end

inher_heads = mod(inher_heads + pi,2*pi) - pi;

if zug_opt == 1 % need inclination or intensity signposts for zugknks

    
    std_inher_signp = max(min(0.5*(succ_std_inher_signp(i_trt_1) + ...
        succ_std_inher_signp(i_trt_2)) + randn(N_inds,1)*min_std_inh_sp, ...
        max_std_inh_sp),min_std_inh_sp); % 

%     try
    inher_signp_kappa(1:N_inds,1) = 1./std_inher_signp.^2;

%     catch
%         keyboard
%     end
   % update variation in inerited signposts
   % depends if directional or absolute (intensity)


    succ_zug_signps = zug_signps(stoppedAndArrived,:);

    if zug_signp == 1 || zug_signp ==  2 % incln and decln are circular var 

%         zug_signps = min(vmrand(0.5*(succ_zug_signps(i_trt_1,:) + ...
%             succ_zug_signps(i_trt_2,:)), ...
%            repmat(inher_signp_kappa,[1 n_zugs]), [N_inds n_zugs]),pi/2);

            for iz = 1:n_zugs
                   
                mn_zsp_iz = circ_mean([succ_zug_signps(i_trt_1,iz) ...
                    succ_zug_signps(i_trt_2,iz)]')';
                
                if ~isinf(inher_signp_kappa(1))
                    zug_signps(:,iz) = min(vmrand(mn_zsp_iz, ...
                       inher_signp_kappa, [N_inds n_zugs]),pi/2);  
                else
                    zug_signps(:,iz) = mn_zsp_iz;
                end

                  if zug_signp ~= 2  || mean(decln_0) <= mean(decln_a)

                        zug_signps = sort(zug_signps,2,'descend');   

                  else

                         zug_signps = sort(zug_signps,2,'ascend');                          

                  end
%                   zug_signps = sort(zug_signps,2,'descend'); 
       
            end
       
%            inher_heads = vmrand(circ_mean([succ_heads(i_trt_1,1)  succ_heads(i_trt_2,1)]')', ...
%         inher_head_kappa,[N_inds,1]);

%                 elseif calibr_comp(1) == 2 || calibr_comp(1) == 3 % same as above, currently
%                     
%                     zug_signps = min(vmrand(0.5*(succ_zug_signps(i_trt_1,:) + ...
%                         succ_zug_signps(i_trt_2,:)), ...
%                        repmat(inher_signp_kappa,[1 n_zugs]), [N_inds n_zugs]),pi/2);

    else % zug_signp == 3 - 5 : non-circ field intens

        zug_signps = 0.5*(succ_zug_signps(i_trt_1,:) + ...
            succ_zug_signps(i_trt_2,:)).* ...
            (1 + randn(N_inds,n_zugs).*repmat(std_inher_signp,[1 n_zugs])); 

    end

%                 .*(1 + randn(N_inds,2)*std_inher_head); 

else

    zug_signps = NaN*ones(N_inds,n_zugs);

end

if zug_signp ~= 2 ||  mean(decln_0) > mean(decln_a)
    
    zug_signps = sort(zug_signps,2,'descend');
    
else
    
    zug_signps = sort(zug_signps,2,'ascend');
    
end

% now zugknick heads 

if ~is_magcl

    % need to inherit zugkn dirns as well
    % magncl cases dealt with below after
    % defining inher_heads via inherited projections
    for iirz = 1:n_zugs % numel(reqs_zug)-1

        if ~isinf(inher_head_kappa(1))
            zugkn_inher_heads(:,iirz) = vmrand( ...
                circ_mean([succ_zug_inh_hds(i_trt_1,iirz)  ...
                succ_zug_inh_hds(i_trt_2,iirz)]')', ...
                inher_head_kappa); % , [N_inds 1]);
        else
           zugkn_inher_heads(:,iirz) = circ_mean([succ_zug_inh_hds(i_trt_1,iirz)  ...
                succ_zug_inh_hds(i_trt_2,iirz)]')'; % , [N_inds 1]);
        end

    end

end

% for magncl strats we still need to calc init heads from
% (inherited) projections
% initialize vector for angles where magncl heads are well defined
 ok_incl = false_vec;
 
date_jul_start = datenum(curr_year*one_vec,month_start*one_vec,day_starts); 
doys = day(datetime(datevec(date_jul_start)),'dayofyear' );

% initialize sun az and geomagn field given locs and start (dep) dates
% initialize_sun_fl_hrs_geomag  
init_sun_geomag_fl_hrs  