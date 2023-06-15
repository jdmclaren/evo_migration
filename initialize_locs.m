function [lat_bs_deps, lon_bs_deps,idxSubPoly] = ...
    initialize_locs(N_inds,Lon_poly,Lat_poly, ...
    require_land,hi_elevs, barrens, poor_stops, ok_veg_breed)

viable_Loc = false(N_inds,1);

lond_i = NaN*ones(N_inds,1); % zeros(N_inds,1);
latd_i = NaN*ones(N_inds,1); %  zeros(N_inds,1);

idxSubPoly = NaN*ones(N_inds,1);

% convert lats to sin(lats) to obtain 
% uniform latitudinal ditribution
sin_lo = sin(nanmin(Lat_poly)*pi/180);
sin_hi = sin(nanmax(Lat_poly)*pi/180);

min_Lon = nanmin(Lon_poly);
max_Lon = nanmax(Lon_poly);
% min_Lat = nanmin(Lat_poly);
% max_Lat = nanmax(Lat_poly);

% if Poly_verts(1) > Poly_verts(2)
%     Poly_verts(1:2) = mod(Poly_verts(1:2),360);
% end

[nSubPolys, LonSub_poly, LatSub_poly] = ...
det_subPolys(Lon_poly,Lat_poly);

if require_land 
    
    while sum(~viable_Loc) > 0

        N_left = N_inds - sum(viable_Loc);

% determine whether randomly selected points are in one of subpolygons

         lond_i(~viable_Loc) = randinterval(min_Lon,max_Lon,N_left);
         
         ran_sin_lats = randinterval(sin_lo,sin_hi,N_left);
         latd_i(~viable_Loc) = asin(ran_sin_lats)*180/pi;

         % update indices for ndvi and elev
        lat_deg_idx = 91 + floor(-latd_i);
        % need shifted for index not for poly (if overlaps 180)
        shift_lond_i = shiftAnglesFromMinus180To180(lond_i);
        lon_deg_idx = min(360,floor(shift_lond_i) + 181);

        if nSubPolys > 1

            inPoly = viable_Loc;

             for iSub = 1:nSubPolys

                inSubPoly = false(N_inds,1);
                inSubPoly(~inPoly) = inpolygon(lond_i(~inPoly),latd_i(~inPoly), ...
                                LonSub_poly{iSub},LatSub_poly{iSub});

                idxSubPoly(~inPoly) = inSubPoly(~inPoly)*iSub;
                inPoly(~inPoly) = inSubPoly(~inPoly);

             end

         else

            inPoly = true(N_inds,1);
            inSubPoly = true(N_inds,1);
            idxSubPoly = ones(N_inds,1);

         end
     
%         try
        is_high_elev = hi_elevs(sub2ind([180 360],lat_deg_idx,lon_deg_idx));
        is_barrn = barrens(sub2ind(size(barrens),lat_deg_idx,lon_deg_idx));
        is_poor_stop = poor_stops(sub2ind(size(poor_stops),lat_deg_idx,lon_deg_idx)); 
        is_ok_veg_br = ok_veg_breed(sub2ind(size(ok_veg_breed),lat_deg_idx,lon_deg_idx)); 
        
%         catch
%             keyboard
%         end
        
         viable_Loc = island(latd_i,lond_i) & ~is_barrn & ~is_high_elev ...
             & ~is_poor_stop & is_ok_veg_br & inPoly;
         % landmask... ,cst_accuracy); % ~land_or_ocean_no_180_message(latd_i,lond_i,5);

    end
    
else
    
    while sum(~viable_Loc) > 0
        
         lond_i = viable_Loc.*lond_i + ...
             ~viable_Loc.*randinterval(min_Lon,max_Lon,N_inds);
         ran_sin_lats = ~viable_Loc.*randinterval(sin_lo,sin_hi,N_inds);
         latd_i = viable_Loc.*latd_i + ~viable_Loc.*asin(ran_sin_lats)*180/pi;

         viable_Loc = inpolygon(lond_i,latd_i,Lon_poly,Lat_poly);
         
    end
    
end

% lon_bs_deps(:,1)= shiftAnglesFromMinus180To180(lond_i);
% don't shift so that sign of headings will be on average consistent
% with how lon_deps are enetered (eg use > for Alaska to Africa)

lon_bs_deps(:,1)= lond_i;
lat_bs_deps(:,1) = latd_i;
