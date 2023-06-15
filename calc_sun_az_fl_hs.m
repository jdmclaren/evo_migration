function  [sun_az, fl_hs] = calc_sun_az_fl_hs(lats,dates,endur_flt, ...
    over_barrier,fl_hs_min,fl_hs_max,dusk_offset)

% 360 - acosd(sind(-23.44* ...
%     cosd(360*255/360))/cosd(60))

c1 = 0.39779; % *pi/180;
c2 = 1.914*pi/180;
c3 = 0.98565*pi/180;

sin_decls = -c1*cos(c3*(dates+10) + ...
    c2*sin(c3*(dates-2)));

cos_lats = cos(lats);

is_sunset = abs(sin_decls) < cos_lats;

sun_az(is_sunset,1) = acos( ...
    -sin_decls(is_sunset)./cos_lats(is_sunset));

sun_az(~is_sunset,1) = pi;

if ~endur_flt
    
    % calc fl_hs
    tan_decls = sin_decls./sqrt(1-sin_decls.^2);

    % solar sunset
    H_sun = acos(-tan(lats(is_sunset & ~over_barrier)).*tan_decls(is_sunset & ~over_barrier));

    fl_hs(is_sunset & ~over_barrier,1) = max(min(24*(1-dusk_offset/12-H_sun/pi), ...
        fl_hs_max),fl_hs_min);
    fl_hs(is_sunset & over_barrier,1) = 12;
    fl_hs(~is_sunset) = fl_hs_min;
    
else
    
    fl_hs = 24*ones(size(lats));
    
end