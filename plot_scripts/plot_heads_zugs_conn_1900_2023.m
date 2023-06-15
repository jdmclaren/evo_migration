addpath 'D:\generic_comp_rte_mdl\brewer'

yrs = [1 10 60 124];

hds_rng = [110 190];
arr_lons_rng = [-15 10];
zg_lats_rng = [10 70];

edge_sz = 2.5; % 
[Lon_dep_poly, Lat_dep_poly] = create_poly(Dep_verts);
[Lon_arr_poly, Lat_arr_poly] = create_poly(Arr_verts);
min_lon = nanmin(Lon_dep_poly)*180/pi-edge_sz;
max_lon = nanmax(Lon_dep_poly)*180/pi+edge_sz;
min_lat = nanmin(Lat_dep_poly)*180/pi-edge_sz;
max_lat = nanmax(Lat_dep_poly)*180/pi+edge_sz;


fig_sz = [200 200 300 300]; % [200 200 400 400];

if max_lon > 70
    max_lat = 90; 
end

for iy = 1:numel(yrs)

    yr = yrs(iy);
    
    thous = min(sum(all_ys_stoppedAndArrived(:,yr)),5000);
    thous_arr = find(all_ys_stoppedAndArrived(:,yr), thous,'first');
    thous_not = min(sum(~all_ys_stoppedAndArrived(:,yr)),5000);
    thous_not_arr = find(~all_ys_stoppedAndArrived(:,yr), thous_not,'first');
    
    n_Arrs_plt = numel(thous_arr);
    
    figure('Position',fig_sz)
    land = shaperead('landareas', 'UseGeoCoords', true);
    ax = worldmap([min_lat max_lat],[min_lon max_lon]);
    geoshow(land, 'FaceColor', [0.75 0.75 0.75]) %
    sc = scatterm(blats(thous_arr,1)*180/pi, ...
    blons(thous_arr,1)*180/pi,20, ...
    all_ys_inher_heads(thous_arr,yr),'fill'); % 'LineWidth',1)
    sc.Children.MarkerFaceAlpha = 0.375; % min(.025*N_inds/min(size(blat_succs,1), ...
    colormap(brewermap([],'*YlOrRd'));
    mlabel('off'); plabel('off'); gridm('off')
    clim(hds_rng)

    if iy == numel(yrs)
        figure('Position',fig_sz)
        cb = colorbar('FontSize',11);
        clim(hds_rng)
        title(cb,{'Magnetic', 'Heading (^o)'},'FontSize',11);
        colormap(brewermap([],'*YlOrRd'));
        box off
        axis off
    end
    
    figure('Position',fig_sz)
    lon_arrs = all_ys_lon_fin(thous_arr,yr);
    land = shaperead('landareas', 'UseGeoCoords', true);
    ax = worldmap([min_lat max_lat],[min_lon max_lon]);
    geoshow(land, 'FaceColor', [0.75 0.75 0.75]) %
    scatterm(blats(thous_arr,1)*180/pi, ...
    blons(thous_arr,1)*180/pi,20,lon_arrs,'o','fill')
    mlabel('off'); plabel('off'); gridm('off')
    colormap(brewermap([],'*RdYlBu')) 
    clim(arr_lons_rng)

    if iy == numel(yrs)
        figure('Position',fig_sz)
        cb = colorbar; % ('Location','SouthOutside');
        clim(arr_lons_rng)
        title(cb,'Arrival Longitude (^o)');
        a =  cb.Position; %gets the positon and size of the color bar
        set(cb,'Position',[a(1) a(2)-0.11 a(3) a(4)]);% To change size
        colormap(brewermap([],'*RdYlBu')) 
        box off
        axis off
    end
    
    if zug_opt == 1
        
        figure('Position',fig_sz)
%         lon_arrs = mod(blons(thous_arr,end)*180/pi+180,360)-180;
        land = shaperead('landareas', 'UseGeoCoords', true);
        ax = worldmap([min_lat max_lat],[min_lon max_lon]);
        geoshow(land, 'FaceColor', [0.75 0.75 0.75]) %
        scatterm(blats(thous_arr,1)*180/pi, ...
        blons(thous_arr,1)*180/pi,20, ...
        all_ys_lat_zg(thous_arr,yr),'o','fill')
        colormap(brewermap([],'RdYlBu')) 
        % caxis([min(lon_arrs)) max(abs(lon_arrs))])
        mlabel('off'); plabel('off'); gridm('off')
        clim(zg_lats_rng)
    %     title('Migratory connectivity')

        if iy == numel(yrs)
            figure('Position',fig_sz)
            cb = colorbar('FontSize',11); % ('Location','SouthOutside');
            clim(zg_lats_rng)
            title(cb,{'Zugknick', 'Latitude (^o)'},'FontSize',11);
            a =  cb.Position; %gets the positon and size of the color bar
            set(cb,'Position',[a(1) a(2)-0.11 a(3) a(4)]);% To change size
            colormap(brewermap([],'RdYlBu'))
            box off
            axis off
        end
    
    end

end