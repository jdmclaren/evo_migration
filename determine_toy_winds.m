function [u_wind, v_wind] = determine_toy_winds(curr_Lat, incl_wind, wind_str)            
is_mid_Lats = abs(curr_Lat)>pi/6;
% Westerlies vs. Easterlies 
u_wind = incl_wind*(is_mid_Lats - ~is_mid_Lats)*wind_str;
% winds from South (headwind) vs. from North (tailwind) assymetrical between Hemispheres
v_wind = incl_wind*(is_mid_Lats*sign(curr_Lat) ...
   - ~is_mid_Lats*sign(curr_Lat))*wind_str;