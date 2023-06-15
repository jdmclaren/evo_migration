        
% if incl_sun
 
[sun_az, fl_hrs] = calc_sun_az_fl_hs(lat_bs_deps*pi/180,doys, ...
    endur_flt,zero_vec,min_init_fl_stp_hs,max_fl_stp_hs,offset_dep_dusk);

% date_nrs(ib,1) = date_jan_1 + t_finish;

fl_hrs_0 = fl_hrs;

sun_az_0 = sun_az;   

% end


if geomagn_model == 1

  %                 tic
  % sort through dates among ind's
  
    udy = round(mean(date_jul)); 
    [Bx, By, Bz] = igrf(udy, ...
           lat_bs_deps, lon_bs_deps, 0);
%   udy = unique(date_jul);
%   for idy = 1:numel(udy)
%        dy_i = udy(idy);
%        is_dyi = date_jul == dy_i;
%        [Bx(is_dyi,1), By(is_dyi,1), Bz(is_dyi,1)] = igrf(dy_i, ...
%            lat_bs_deps(is_dyi,1), lon_bs_deps(is_dyi,1), 0);
%   end
    decln_0 = atan2(By,Bx);
    hz_inten_0 = hypot(Bx,By);
    incln_0 = atan(Bz./hz_inten_0); % min( % + rand_magn_0,pi/2);
    tot_B = hypot(hz_inten_0,Bz);
%             Bz = Bz./tot_B;
%             hz_inten_0 = hz_inten_0./tot_B;
%             tot_inten_0 = sqrt(1+3*hz_inten_0.^2);
    vt_inten_0 = Bz; % .*(1 + randn(N_inds,1)*std_rel_err_inten); % ./tot_B; % sin(theta_0);
    tot_inten_0 = tot_B; % .*(1 + randn(N_inds,1)*std_rel_err_inten);

    % initialiaze all year fields (pewrfect for now..)
    decln_0s(:,iYear) = decln_0*180/pi;
    intns_0s(:,iYear) = tot_B;
    incln_0s(:,iYear) = incln_0*180/pi;  


else

    % Also, we use that for dipole, tan(inclin angle) = 2*tan(theta)
    % initialize magnetic components assuming (non-tilted) dipole geomagnetic
    % field (depends on initial orientation init_head in case of magnetoclinic "compass")
    vt_inten_0 =  sin(theta_0 + rand_magn_0);
    incln_0 =  min(atan(2*tan(theta_0 + rand_magn_0)),pi/2);

%             magncl_parll_zs = cos(zugkn_inher_heads).*cos(incln_0);
    hz_inten_0 =  cos(theta_0 + rand_magn_0);
    tot_inten_0 = sqrt(1+3*vt_inten_0.^2);

    % don't include declination in dipole model
    decln_0 = zero_vec;

end
