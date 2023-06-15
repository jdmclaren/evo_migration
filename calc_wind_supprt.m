function tailw = calc_wind_supprt(alpha,u,v)

w_str = sqrt(u.^2 + v.^2);
w_dir = atan2(u,v);
tailw = w_str.*cos(pi+alpha-w_dir);