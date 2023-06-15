function [u_wind, v_wind] = det_stoch_wind(curr_Lat, incl_wind, wind_str,u_pr,v_pr,use_prev)  

is_mid_Lats = abs(curr_Lat)>=pi/6  & abs(curr_Lat)<=pi/3;

% new add random switch to tailwinds (also 2/3 chance)
p_EW = 0.75;
p_NS = 0.75;

% Westerlies vs. Easterlies 
EW_Ws = sign(p_EW - rand(size(curr_Lat)));
u_new = incl_wind*EW_Ws.*(is_mid_Lats - ~is_mid_Lats)*wind_str;
% winds from South (headwind) vs. from North (tailwind) assymetrical between Hemispheres

% Northeries vs. Southerlies 
NS_Ws = sign(p_NS - rand(size(curr_Lat)));
v_new = incl_wind*NS_Ws.*(is_mid_Lats.*sign(curr_Lat) ...
   - ~is_mid_Lats.*sign(curr_Lat))*wind_str;

if use_prev == 1
    
    % 2/3 chance same as last time
    p_steady = 0.75;
    steadys = rand(size(curr_Lat)) < p_steady;
    
    u_wind = steadys.*u_pr + ~steadys.*u_new;
    v_wind = steadys.*v_pr + ~steadys.*v_new;
    
else
    
    u_wind = u_new;
    v_wind = v_new;
    
end