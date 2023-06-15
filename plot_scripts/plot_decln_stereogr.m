addpath D:\evo_mdl_clean\brewer\
addpath D:\evo_mdl_clean\circ_stats\

 load coastlines; 
 lnd_clr =  0.775; % 0.875; % 

% create colormaps
% 
cmaps_prd = colormap(brewermap(12,'Paired')); % 6,'Dark2'));
lt_gr_col = cmaps_prd(3,:);
dk_gr_col = cmaps_prd(4,:);
dk_bl_col =  cmaps_prd(2,:);
dk_red_col = cmaps_prd(6,:);
% lt_or_col = cmaps_prd(7,:);
% dk_or_col =  cmaps_prd(8,:);
% lil_col =  cmaps_prd(9,:);
% prp_col =  cmaps_prd(10,:);
% yl_col =  cmaps_prd(11,:);
% br_col = cmaps_prd(12,:);

cmaps_solid = colormap(brewermap(9,'Set1')); % 6,'Dark2'));
pnk_col =  cmaps_solid(8,:);
prp_col =  cmaps_solid(4,:);
br_col =  cmaps_solid(7,:);

% virtual track
NBI_col = pnk_col; 

br_rng_col = br_col; % prp_col; % dk_red_col; % 
BFI_col = br_col; % prp_col; % dk_red_col; %  dk_bl_col; % 
SWG_col = pnk_col; % dk_red_col; % 

arr_col = dk_gr_col; % prp_col; % lt_gr_col; % 'w'; % 

% detrmine data and tracks for leucoroha migr and records

% for dept & arr data load a random decent simulation
load 'dep_arr_locs_whtrs'

% virtual N Baff Is migr to Eire
[NBI_Lt_gc,NBI_Ln_gc] = track('gc',[71 54],[-67.5 -9]);
[NBI_Lt,NBI_Ln] = track('rh',[71 54],[-67.5 -9]);

% BaffIs and SWGrnlnd locs
SWG_lats = lat_bs_deps(idx_deps==4);
SWG_lons = lon_bs_deps(idx_deps==4);
Bf_lats = lat_bs_deps(idx_deps==2);
Bf_lons = lon_bs_deps(idx_deps==2);

% Ottosson wheatear data
SP_Lt = mean([41.5 42 40.75 43 44.5 44.75 46.25]);
SP_Ln = mean([-7 -7.5 0.25 -.5 -2 -1 -1.5]);

[Gd_SP_Lt,Gd_SP_Ln] = track('rh',[64.33 SP_Lt],[-52 SP_Ln]);
[Gd_SP_Lt_gc,Gd_SP_Ln_gc] = track('gc',[64.33 SP_Lt],[-52 SP_Ln]);

% Bairlein Iqauit track
[IqUK_Lt,IqUK_Ln] = track('rh',[63.7 55.5],[-68.5 -6]);
[UKAf_Lt,UKAf_Ln] = track('rh',[55.5 19],[-6 -15]); % 25 -6 

[IqUK_Lt_gc,IqUK_Ln_gc] = track('gc',[63.7 55.5],[-68.5 -6]);
[UKAf_Lt_gc,UKAf_Ln_gc] = track('gc',[55.5 19],[-6 -15]);

% line wdth trajs
LWtrj = 2.; % 1.5;

% alpha breeding (to see contours)
Alf = 0.5; 

% md_Zg_lat_SG = mean(mn_lat_EU_DepPs(4,:));
% md_Zg_lon_SG = mean(mn_lon_EU_DepPs(4,:));
% md_Arr_lat_SG = mean(mn_Arr_lat_DepPs(4,:));
% md_Arr_lon_SG = mean(mn_Arr_lon_DepPs(4,:));

% md_Zg_lat_IQ = mean(mn_lat_EU_DepPs(2,:));
% md_Zg_lon_IQ = mean(mn_lon_EU_DepPs(2,:));
% md_Arr_lat_IQ = mean(mn_Arr_lat_DepPs(2,:));
% md_Arr_lon_IQ = mean(mn_Arr_lon_DepPs(2,:));

% see astFormatContourLines  

% model_epochs = {'2020','2020'};
genrn = 13; 
decimal_years = [1900 1960 2023]; % 2010.75; %; % [1900 1960 2023]; %    
n_Yrs = numel(decimal_years);

params = [5 7 6]; % [7 6]; % 5; % 

para_nms = {{'Geomagnetic','Declination (^o)'},'Intensity (10^4 nT)','Inclination (^o)'};

lims_cnts{1} = [-90 90]; % [-110 110];
lims_cnts{2} = [3.5e4 6.e4];
lims_cnts{3} = [0 120];
% lims_cnts{1} = [3.5e4 6.e4];
% lims_cnts{2} = [0 120];

lbl_lns{1} = -90:10:90; %-110:10:110; % [-75 -60 -45 -30 -20 -10 0 10 20 30 75];
lbl_lns{2} = 3.e4:5e3:6.e4;
lbl_lns{3} = [0:10:90 120];
% lbl_lns{1} = [3.5e4:5e3:6.e4];
% lbl_lns{2} = [0:10:90 120];

cbls{1} = -90:20:90; %
cbls{2} = 3.:.5:6.;
cbls{3} = 0:20:120; %

d_cnts = [5 5e3 10];
lev_Steps = [2.5 2.5e3 5]; % [2.5e3 5]; % 

cmap_strs = {'*RdYlBu','YlOrRd','RdPu'}; % {'YlOrRd','RdPu'}; % 'Reds'}; % '*RdBu'};

% WMMFileName = 'astWMMResults_Epoch_2020_decyear_2020.mat';
% load(WMMFileName);
% landAreas = shaperead('landareas.shp','UseGeoCoords',true);
% plotWMM = load('astPlotWMM.mat');
% model_epoch = '2015v2';
% decimal_year = 2009;

map_proj = 'stereo'; % 'Mercator';  % 
if strcmp(map_proj,'stereo')
        min_lat = -5; % -5 % 
        max_lat =  80; %70 % 
        min_lon = -15; % + (magn_comp_nr == 4)*40;
        max_lon = 15; % 
        
        map_lat_cntr = 45; % 60; % (min_lat + max_lat)/2;
        map_lon_cntr = -11; % (min_lon + max_lon)/2;
        FLatLims = [43.5 43.5]; % [35 35]; %[0[30 30]; %
else
             min_lon = -85;
             max_lon = 25;

             min_lat = 0;
             max_lat = 85
end


land = shaperead('landareas','UseGeoCoords', true);

for iYr = 2:n_Yrs
    
    %     model_epoch = model_epochs{iYr};
        decimal_year = decimal_years(iYr);
    
    height = 0;
    % Geodetic Longitude value in degrees to use for latitude sweep.
    
    geod_lon = -180:1:180;
    % Geodetic Latitude values in degrees to sweep.
    
    geod_lat = -89.5:.5:89.5;
    % Loop through longitude values for each array of latitudes -89.5:89.5.
    
    for lonIdx = size(geod_lon,2):-1:1
        for latIdx = size(geod_lat,2):-1:1
    % To obtain magnetic parameters for each latitude and longitude value, use the wrldmagm function.
     [xyz, h, dec, dip, f]  = igrfmagm(height, geod_lat(latIdx),geod_lon(lonIdx), decimal_year, genrn);
    %     [xyz, h, dec, dip, f] = wrldmagm(height, geod_lat(latIdx),geod_lon(lonIdx), decimal_year, model_epoch);             
    % Store the results.
    
        IGRFResults(latIdx,1:7,lonIdx) = [xyz h dec dip f];
        
        end
    end
    
    for jVar = 3:numel(params)
    
        param = params(jVar); % declination
        levStep = lev_Steps(jVar);
        cmap_str = cmap_strs{jVar};
    
        H = figure('Position',[200 200 300 300]); % [left bottom width height],'Color','w');
        hold on
         if strcmp(map_proj,'stereo')
            ax = axesm ('stereo', 'Frame', 'on', 'Grid', 'off','Origin',[map_lat_cntr map_lon_cntr 0],'FlineWidth',1);
          setm(gca,'FLatLimit',FLatLims) % ,'FLonLimit',[-30 90])
         
         elseif strcmp(map_proj,'Mercator')

               ax = worldmap([min_lat max_lat],[min_lon max_lon]);
              setm(ax,'mapprojection',map_proj)

                map_lat_cntr = (min_lat + max_lat)/2;
                map_lon_cntr = (min_lon + max_lon)/2;

         end
            %            ax.PositionConstraint = 'innerposition';
               axis('off')
              
               tightmap
               
%           geoshow(land, 'FaceColor',[0.875 0.9 0.865 ],'EdgeColor',[0.65 0.65 0.65]) %
          geoshow(land, 'FaceColor', [lnd_clr lnd_clr lnd_clr],'LineWidth',0.2) %[0.85 0.85 0.85]
    
        plotm(coastlat,coastlon,'Color',[0.06 0.35 0.2],'LineWidth',0.25) % 'g')


%           if iYr ~= n_Yrs || iYr == 1 %  == 2 % 

%               idx_succ = find(stoppedAndArrived); % ,500,'first');

              % for ii = 1:500
            
              sc = scatterm(lat_bs_deps, lon_bs_deps,3, ...
                        br_rng_col,'fill');
              sc.Children.MarkerFaceAlpha = Alf;
               scatterm(lat_bs_arrs, lon_bs_arrs,3, ...
                     arr_col,'fill')  
               sc2 = scatterm(Bf_lats,Bf_lons,3, ...
                        BFI_col,'fill');
               sc2.Children.MarkerFaceAlpha = Alf;

%           elseif iYr == n_Yrs
% 
%                sc1 = scatterm(SWG_lats,SWG_lons,3, ...
%                         SWG_col,'fill');
%                sc1.Children.MarkerFaceAlpha = Alf;
%                sc2 = scatterm(Bf_lats,Bf_lons,3, ...
%                         BFI_col,'fill');
%                sc2.Children.MarkerFaceAlpha = Alf;
%                scatterm(lat_bs_arrs, lon_bs_arrs,3, ...
%                      arr_col,'fill')  
%         
% 
%           end

         % Data is stored in the matrix IGRFResults (length(lat),1:7,length(lon)):

        if jVar == 3
             decs_ys{iYr} =  squeeze(IGRFResults(:,param,:));
        end

        [ccc,hhh] = contourm(geod_lat,geod_lon,squeeze(IGRFResults(:,param,:)),...
            'levelStep',levStep,...
            'edgecolor',[.5 .5 .5],'linewidth',.015); % [.7 .7 .7]
            % Format major contour lines' color and width to match published reports.
        % http://www.ngdc.noaa.gov/geomag/WMM/
        labelLine = lbl_lns{jVar}; % lims_cnts{jVar}(1):d_cnts(jVar):lims_cnts{jVar}(2); % -180:10:180; % plotWMM.formatRange{param}(1):plotWMM.text_stepN{param}:plotWMM.formatRange{param}(2);
        
%         n_lns = numel(lbl_lns); % (lims_cnts{jVar}(2) - lims_cnts{jVar}(1))/lev_Steps(jVar);
        FormatContourLines_magmaps(ccc,hhh,labelLine,cmap_str); % ,levStep,n_lns); 


%         if iYr == 2 || n_Yrs == 1

                plotm(IqUK_Lt, IqUK_Ln,'Color',BFI_col,'LineStyle',':','LineWidth',LWtrj);
                plotm(UKAf_Lt, UKAf_Ln,'Color',BFI_col,'LineStyle',':','LineWidth',LWtrj);
                scatterm(UKAf_Lt(end-2), UKAf_Ln(end-2),80,BFI_col,'v','fill');

%         elseif iYr == 1
% 
% %                 plotm(NBI_Lt_gc,NBI_Ln_gc,'Color',NBI_col,'LineStyle','--','LineWidth',LWtrj);
%                 plotm(NBI_Lt,NBI_Ln,'Color',NBI_col,'LineStyle',':','LineWidth',LWtrj);
% 
%         elseif iYr == 3
% 
% 
%                 plotm(IqUK_Lt, IqUK_Ln,'Color',BFI_col,'LineStyle','--','LineWidth',LWtrj);
%                 plotm(UKAf_Lt, UKAf_Ln,'Color',BFI_col,'LineStyle','--','LineWidth',LWtrj);
%                 plotm(Gd_SP_Lt, Gd_SP_Ln,'Color',SWG_col,'LineStyle','--','LineWidth',LWtrj);
% 
%         end

    end

         if iYr > 1 % == 3 % 
 
            ddecs = decs_ys{iYr} - decs_ys{iYr-1};            
            plot_diff_dec

        end

end

% colorbars for each Var
for jVar = 1:numel(params)

    Hcb = figure('Position',[200 200 300 300]); % [left bottom width height],'Color','w');
%     cmap = colormap(brewermap(numel(lbl_lns{jVar}),cmap_strs{jVar}));
%     colormap(cmap(1:numel(lbl_lns{jVar}),:)) % (cmap(1:18,:))1:numel(lbl_lns{jVar}),:)) % (cmap(1:18,:))
    colormap(brewermap([],cmap_strs{jVar}));
    cb = colorbar;
    set(cb,'TickLength',0, ...'TickDirection','Out', ...
        'TickLabels',cbls{jVar},'FontSize',10); % -100:20:60); % ,'Ticks',-90:20:45); % ,'TickLabels',-90:15:47.5); % ('XTick',45:10:65,'FontSize',8); % 'NorthOutside',
%     set(cb,'YTick',[43 48 53 57],'YTickLabel',{'N Spain','N France','Ireland','N Scotland'})
    %     cb = colorbar('XTick',25:25:100,'FontSize',9); % 'NorthOutside',
    %     clim([25 100])
    clim(lims_cnts{jVar}) % [-65 45]) % clim([-115 65]) % })
    %     set(cb)
    title(cb,para_nms{jVar},'FontSize',10)

    axis off 
    box off

end
     