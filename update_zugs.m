% update to (possible) zugknik angle if passed Stopover
% (if so, first post stop record will have d_alpha_magn == 0)

if zug_opt == 1
    
        signp_vals = (zug_signp == 1)*incln(~finished) + ...
            (zug_signp == 2)*decln(~finished) + ...
            (zug_signp == 3)*tot_inten(~finished) +  ...
            (zug_signp == 4)*vt_inten(~finished) + ...
             (zug_signp == 5)*hz_inten(~finished);
        % May 2022 keep track of signps and zugs separately 
        % (for n_zugs = 0 case)
%         try
        if  req_land_signp
             poss_past =  ~over_water(~finished)  ...
                 & ~is_high_elev(~finished) & ~is_barrn(~finished) ...
                 & ~(skip_stop_over(~finished) & require_stopover_zug) ; % %
             past_sgnp(~finished) = rowwiseLast(repmat(poss_past,[1,n_sgnps]) & ...
                repmat(signp_vals ,[1,n_sgnps]) <= zug_signps(~finished,:)); 
        else
             past_sgnp(~finished) = rowwiseLast(repmat(find(~finished),[1,n_sgnps]) & ...
                repmat(signp_vals ,[1,n_sgnps]) <= zug_signps(~finished,:)); 
        end
%         catch
%             keyboard
%         end
%         if sum(past_sgnp) > 0
%             keyboard
%         end
            
         zug_nr(~finished) = (n_zugs > 0)*max(past_sgnp(~finished),zug_nr(~finished));   
         
        
    if size(zugkn_inher_heads,2) > 1
        zug_nr_lin = max(zug_nr,1);
        lin_zk = sub2ind(size(zugkn_inher_heads), 1:size(zugkn_inher_heads,1),zug_nr_lin');      
    else
        lin_zk = 1:size(zugkn_inher_heads); 
    end
           
    is_new_zug = false_vec;
    [new_z_i, new_z_j] = find(zug_nr > curr_zug_nr);
    is_new_zug(new_z_i) = true;
% try
    % new for signp for refuel only (n_zugs = 0)
     is_new_sgnp = false_vec;
    [new_sgnp_i, new_sgnp_j] = find(past_sgnp & ~curr_past_sgnp);
    is_new_sgnp(new_sgnp_i) = true;
% catch
%     keyboard
% end
    
%     if any(over_water(new_z_i))
%         keyboard
%     end
    
    if numel(new_z_i) > 0

         zug_heads = zugkn_inher_heads(lin_zk(new_z_i)');
                 
            % set initial heading (alpha) error per sun or magn compass
            for i_comp = 1:3

                    if calibr_comp(i_comp) > 0 || endog_comp == i_comp 
                        if perf_init_calibr(i_comp) == 0
                            d_alpha_z(:,i_comp) = vmrand(0, compass_cal_kappa(i_comp), [N_inds 1]);
                        else
                            d_alpha_z(:,i_comp) = zero_vec;
                        end
                    % add magn varn to sun compass if refl about mag S axis
                    elseif i_comp == 1 && refl_sun_mag && ...
                            calibr_comp(2) > 0 && cal_comp_nonzero(i_comp) > 0 
                            d_alpha_z(:,i_comp) = vmrand(0, compass_cal_kappa(i_comp), [N_inds 1]);
                    else
                        d_alpha_z(:,i_comp) = zero_vec;
                    end

            end

            % need to define base compass headings dependeing on which is
            % endogenous and whether they are initially perfectly assessed
            % and what any transfer error is made for the non-endogenous
            % compasses

            % for sun azimuth compass, account for reflection about geomag axis
            % option

            decln_z_hat = decln(new_z_i) + d_alpha_z(new_z_i,1);
            sun_defl_z = sun_az(new_z_i)*sign(sign_clockw_sun) + ...
                2*(refl_sun_mag)*decln(new_z_i);
            sun_defl_z_hat = sun_defl_z + d_alpha_z(new_z_i,2);
            
            % time-comp sun comp has no sun azimuth deviation
            
            if ismember(calibr_comp(2),[1 2])
              d_alpha_sun(new_z_i) = 0; %sun_az(new_z_i); % (calibr_comp(2) == 2)*zug_heads +
%               prev_sun_az(new_z_i) =  sun_az(new_z_i);
            end

            % determine actual heading taken
            alpha_z = mod(zug_heads + (endog_comp == 1)*decln_z_hat + ...
                (endog_comp == 2)*sun_defl_z_hat + ...
                (endog_comp == 3)*d_alpha_z(new_z_i,3) + pi,2*pi) - pi; 

            % add error in perceived geographic compass relative to others 
            alpha_z_hat = mod(alpha_z - (endog_comp ~= 3)*d_alpha_z(new_z_i,3) + pi,2*pi) - pi;
            
            % determine "base" headings for each compass, whether inherited
            % (endog) or transferred on (this) first departure

            switch endog_comp

                case 1 % magn

                    alpha_bases{1}(new_z_i) = zug_heads;
                    absDiff_sun = absDiffDeg(alpha_z*rad2deg,sun_defl_z_hat*rad2deg)*pi/180;       
                    alpha_bases{2}(new_z_i) =  mod(absDiff_sun + pi,2*pi) - pi;     
                    alpha_bases{3}(new_z_i) =  alpha_z_hat;

                case 2 % sun

                    alpha_bases{2}(new_z_i) = zug_heads; %
                    absDiff_mag = absDiffDeg(alpha_z*rad2deg,decln_z_hat*rad2deg)*pi/180;
                    alpha_bases{1}(new_z_i) =  mod(absDiff_mag + pi,2*pi) - pi;            
                    alpha_bases{3}(new_z_i) =   alpha_z_hat;

                case 3 % geo taken care of

                    alpha_bases{3}(new_z_i) =  zug_heads;  
                    absDiff_mag = absDiffDeg(alpha_z_hat*rad2deg,decln_z_hat*rad2deg)*pi/180;
                    alpha_bases{1}(new_z_i) =  mod(absDiff_mag + pi,2*pi) - pi;     
                    absDiff_sun = absDiffDeg(alpha_z_hat*rad2deg,sun_defl_z_hat*rad2deg)*pi/180;       
                    alpha_bases{2}(new_z_i) =  mod(absDiff_sun + pi,2*pi) - pi;    

            end

    
% update Jun 2020 base is just inherited angle (first or Zugkn)
%      alpha_base = (zug_nr==0).*inher_heads + ...
%          (zug_nr>0).*(zugkn_inher_heads(lin_zk'));    

% update if magnetocinic and new zug
% if numel(new_z_i) >0

        for ii = 1:numel(new_z_i)

            if  calibr_comp(1) == 2

                if iWarmup == 1 % set (inherited) inclination projection

                    if magcl_opt == 1

                         magn_param_0(new_z_i(ii)) = ...
                             magn_param(new_z_i(ii))./ ...
                             sin(zugkn_inher_heads(new_z_i(ii), reqs_zug(new_z_j(ii))));
                         % ,new_z_j(ii) lin_zk(new_z_stop)
                         coeff_magn_steer(new_z_i(ii),1+zug_nr(new_z_i(ii))) = magn_param_0(new_z_i(ii));

                    elseif magcl_opt == 2

                          magn_param_0(new_z_i(ii)) = ...
                             sin(zugkn_inher_heads(new_z_i(ii), reqs_zug(new_z_j(ii)))).* ...
                             magn_param(new_z_i(ii));
                         % ,new_z_j(ii) lin_zk(new_z_stop)
                         coeff_magn_steer(new_z_i(ii),1+zug_nr(new_z_i(ii))) = magn_param_0(new_z_i(ii));                   

                    else


                    end

                else % use successful zugknick projection

                    magn_param_0(new_z_i(ii)) = coeff_magn_steer(new_z_i(ii),1+zug_nr(new_z_i(ii)));

                end

            elseif calibr_comp(1) == 3

                if iWarmup == 1 % set (inherited) inclination projection

                   if magcl_opt == 1

        %                 magn_param_0(new_z_i(ii)) = magn_param(new_z_i(ii)).* ...
        %                  cos(zugkn_inher_heads(new_z_i(ii), reqs_zug(new_z_j(ii)))); % lin_zk'));
                      magn_param_0(new_z_i(ii)) = magn_param(new_z_i(ii))./ ...
                         cos(zugkn_inher_heads(new_z_i(ii), reqs_zug(new_z_j(ii)))); 
                      coeff_magn_steer(new_z_i(ii),1+zug_nr(new_z_i(ii))) = magn_param_0(new_z_i(ii));

                   elseif magcl_opt == 2

                         magn_param_0(new_z_i(ii)) = magn_param(new_z_i(ii)).* ...
                         cos(zugkn_inher_heads(new_z_i(ii), reqs_zug(new_z_j(ii)))); % lin_zk'));

                      coeff_magn_steer(new_z_i(ii),1+zug_nr(new_z_i(ii))) = magn_param_0(new_z_i(ii));

                   else

                   end

                  if any(imag(magn_param_0)~=0)
                      keyboard
                  end

                else

                     magn_param_0(new_z_i(ii)) = coeff_magn_steer(new_z_i(ii),1+zug_nr(new_z_i(ii)));               
                   if any(imag(magn_param_0)~=0)
                      keyboard
                   end

                end

            end

            
%             if ~island(theta(new_z_i(ii))*rad2deg, ...
%                 llamda(new_z_i(ii))*rad2deg)
%                [d_Stops] = distance(theta(new_z_i(ii)), ...
%                         llamda(new_z_i(ii)),coastlat,coastlon,'radians'); % , ang_Stops
%                 if nanmin(d_Stops) > 100/R_Earth_km
%                     keyboard
%                 end
%             end
            lat_bs_zugs(new_z_i(ii),zug_nr(new_z_i(ii))) = theta(new_z_i(ii));
            lon_bs_zugs(new_z_i(ii),zug_nr(new_z_i(ii))) = llamda(new_z_i(ii));
        

        end

    elseif numel(new_sgnp_i) > 0

            lat_bs_zugs(new_sgnp_i,1) = theta(new_sgnp_i);
            lon_bs_zugs(new_sgnp_i,1) = llamda(new_sgnp_i);

    end

end