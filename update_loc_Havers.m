function [Lat, Lon] = update_loc(fl_hrs,hourly_del_Lat,curr_Lat,curr_Lon, ...
    curr_head,incl_wind,u_wind,v_wind)

del_Lon = fl_hrs*hourly_del_Lat/cos(curr_Lat)* ...
   (-sin(curr_head)*(1 - incl_wind*sign(sin(curr_head))*u_wind)); % mod( ,2*pi);
del_Lat = fl_hrs*hourly_del_Lat*cos(curr_head)*(-1 + incl_wind*v_wind);

R_Earth = 6378.137*1000; % in meters
% TR = dt/R_Earth;
            d_tb = fl_hrs*hourly_del_Lat/R_Earth;
            sdtb = sin(d_tb);
            cdtb = cos(d_tb);
            c_pd = cos(curr_head);
            s_pd = sin(curr_head);
            
            s_lat = sin(lat_bs(idx_locs_pre));
            c_lat =  cos(lat_bs(idx_locs_pre));
            
            lat_bs(idx_locs) = ~Stopped.* ...
                asin(s_lat.*cdtb + c_lat.*sdtb.*c_pd) + ...
                Stopped.*Stop_poly(min_Stops,2);
            
           as_term =  asin(s_pd.*sdtb./cos(lat_bs(idx_locs)));
           lon_bs(idx_locs) = ~Stopped.*(mod(lon_bs(idx_locs_pre) ...
                + as_term+pi,2*pi) -pi) + Stopped.* ...
                Stop_poly(min_Stops,1);
            
% set Lon between -pi and pi conforming to function determining 
% if over water 
if abs(curr_Lon + del_Lon) > pi
   mod_Lon = mod(curr_Lon + del_Lon,2*pi);
   if mod_Lon > pi
       Lon = -2*pi + mod_Lon;
   else
        Lon = mod_Lon;
   end
else
   Lon = curr_Lon + del_Lon;
end 

if curr_Lat + del_Lat > pi/2      
    Lat = pi - (curr_Lat + del_Lat);
else
    Lat = curr_Lat + del_Lat;       
end