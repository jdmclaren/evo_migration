function [Lat, Lon] = update_loc(fl_hrs,hourly_del_Lat,curr_Lat,curr_Lon, ...
    curr_head,incl_wind,u_wind,v_wind,magn_star_night,decln,dateyear,magn_model)
% , sum_abs_decln,sum_decln,all_decs

u_wnd_fact = 1+incl_wind.*u_wind;
v_wnd_fact = 1+incl_wind.*v_wind;

   udy = round(mean(dateyear)); 
% initialize declination (will remain fixed for star compass)
    d_decln = zeros(size(decln));
%     sum_abs_decln = zeros(size(decln));
%     sum_decln = zeros(size(decln));
%     all_decs = decln;
%     decln_new = NaN*ones(size(decln));
% end

    flr_fl_hrs = floor(fl_hrs);
    
            
    t_res = 1; % 10;  %
        
    
    for iH = 1:max(flr_fl_hrs)*t_res
        
        flies = flr_fl_hrs >= iH/t_res;
      try  
        curr_head(flies) = mod(pi + curr_head(flies) + d_decln(flies),2*pi)-pi;
      catch
          keyboard
      end
        dLon = hourly_del_Lat.* ...
            (sin(pi+curr_head(flies)).*u_wnd_fact(flies))./ ...
            cos(curr_Lat(flies))/t_res;
        dLat = hourly_del_Lat.*(cos(pi+curr_head(flies)).*v_wnd_fact(flies))/t_res;
        
        curr_Lon(flies) = curr_Lon(flies) + dLon; % mod( ,2*pi);
        curr_Lat(flies) = curr_Lat(flies) + dLat;
    
%         if any(curr_Lat > pi/2)
%             keyboard
%         end
        
        % account for polar crossing
        curr_Lat(curr_Lat > pi/2) = pi - curr_Lat(curr_Lat > pi/2);
        curr_Lat(curr_Lat < -pi/2) = -pi - curr_Lat(curr_Lat < -pi/2);

        % account for datelne passing
        curr_Lon(flies) = shiftAnglesFromMinus180To180(curr_Lon(flies)*180/pi)*pi/180;
%         clear Bx By
        if magn_star_night == 1 && magn_model == 1
            % update declination hourly if magnetic loxodrome
       % unique(date_jul(~finished));
%                       dy_no_fin = date_jul(~finished);
       
%     try
                [Bx, By, ~] = igrf(udy, ...
                   curr_Lat(flies)*180/pi, curr_Lon(flies)*180/pi, 0);

%     catch
%         keyboard
%     end
    
%             dy_fly = dateyear(flies);
%             udy = unique(dy_fly);
%             currLat_fl = curr_Lat(flies)*180/pi;
%            currLon_fl = curr_Lon(flies)*180/pi;
%           for idy = 1:numel(udy)
%                dy_i = udy(idy);
%                is_dyi = dy_fly == dy_i;
% %                try
%                [Bx(is_dyi,1), By(is_dyi,1), ~] = igrf(dy_i, ...
%                    currLat_fl(is_dyi,1), currLon_fl(is_dyi,1), 0);
% %                catch
% %                    keyboard
% %                end
%           end
          
%             [Bx, By, ~] = igrf(dateyear, curr_Lat(flies)*180/pi, curr_Lon(flies)*180/pi, 0);
% try
            decln_new = atan2(By,Bx);
            d_decln(flies,1) = decln_new-decln(flies);
%             sum_abs_decln = sum_abs_decln + abs(d_decln);
%             sum_decln  = sum_decln + d_decln;
%             d_decln(~flies) = 0;
%             decln(flies,1) = decln_new;
%             all_decs = [all_decs decln];
% catch
%     keyboard
% end
            
%         else d_decln == 0 unchanged
%             
%             

        end
        
%         if any(curr_Lat < 10*pi/180 & curr_Lon > 0)
%             keyboard
%         end
        
    end
    
    % now complete last fractional flight hour (all records)
    fl_h_frac = rem(fl_hrs*t_res,1)/t_res;
    curr_head(flies) = mod(pi + curr_head(flies) + d_decln(flies),2*pi)-pi;
    dLon = fl_h_frac.*hourly_del_Lat./cos(curr_Lat).* ...
                (sin(pi+curr_head).*u_wnd_fact)/t_res;
    dLat = fl_h_frac.*hourly_del_Lat.*(cos(pi+curr_head).*v_wnd_fact)/t_res;

    % (final) increment location
    Lon = curr_Lon + dLon; % mod( ,2*pi);
    Lat = curr_Lat + dLat;

% end

% (final) accounting for polar crossing
Lat(Lat > pi/2) = pi - Lat(Lat > pi/2);
Lat(Lat < -pi/2) = -pi - Lat(Lat < -pi/2);

% and for datelne passing
Lon = shiftAnglesFromMinus180To180(Lon*180/pi)*pi/180;



