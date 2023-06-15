addpath 'D:\Oldenburg_models\generic_comp_mig_model\brewer'

% plot 1 std or 2 stds (~95% range)
n_stds = 1; %  2; % 

% create colormaps
cmaps_solid = colormap(brewermap(7,'Set1')); % 6,'Dark2'));
bl_col =  cmaps_solid(2,:);
gr_col = cmaps_solid(3,:);
brn_col = cmaps_solid(7,:);

plot_col =  gr_col; % bl_col; % brn_col; %

data_col = 'k'; 

cmap = colormap(brewermap(N_mig_Sims,'YlOrRd')); %'*Greys')); % RdYlBu')); %,'*YlGn')); %  *Reds')); %  '*YlGnBu'));
% cmap = colormap(brewermap(100,'*YlOrRd')); %'*Reds')); % *RdYlBu')); % '*YlGnBu'));

plot_stop_poly =  false % true % 
plot_fails = false; % true; % 
plot_order_Lon =  is_whtr | mig_sys == 2; % false; % 
ord_EW = 1; % 2; % 

% resolution of maps to look like trajectories
% rather than scatter plots
n_substp = 12; 
d_substp = 1/n_substp;

NptsEll = 500; % 100; % 

if ~exist('new_fig')
    % 1 geo 2 mag 3 sun
    new_fig = false %  true % 
end
if ~exist('plot_arrs')
    % 1 geo 2 mag 3 sun
    plot_arrs = false %  true % 
end
if ~exist('plot_deps')
    % 1 geo 2 mag 3 sun
    plot_deps =  false % true % 
end
if ~exist('plot_zugs')
    % 1 geo 2 mag 3 sun
    plot_zugs =  false %true %  
end

if ~exist('plot_fail')
    % 1 geo 2 mag 3 sun
    plot_fail = true % false %  
end

if ~exist('plot_poly_opt')
    % 1 geo 2 mag 3 sun
    plot_poly_opt = 0 %1   3 %  
end

if ~exist('heads_type')
    % 1 geo 2 mag 3 sun
    heads_type = 1 % 2 %  3 %  
end

if ~exist('plot_opt')
    % plot option 1 = init heading 2 = init Lon
    plot_opt = 2*(mig_sys ~= 2) + 1*(mig_sys == 2);
end

if heads_type == 1
    
    all_heads = all_geo_heads; % mod(all_geo_heads+180,360)-180;

elseif heads_type == 2
    
    all_heads = all_mag_heads;
    
else
    
    all_heads = all_sun_heads; 
    
end

if ~exist('transp_plot')
    
    transp_plot = 0.2; % 0.55
    
end

% mag_str = {'', 'incl signpost ','decln signpost ','total str signpost ','vert str ','trsv magncl ','pll magncl ','horiz str '};
strat_string = [z_str{any_zug+1} zug_str{zug_signp+1} mag_str{calibr_comp(1)+1} ...
    sun_str{calibr_comp(2)+1} geo_str{calibr_comp(3)+1} magn_star_night_str{magn_star_night} cst_str{dtrCstOpt+1} ...
    reO_str{reOrWtr_opts+1}];

% disp_summry_stats

% map ranges
if mig_sys == 1
    min_lat = -15;
    max_lat = 70;
    min_lon = -130; % + (magn_comp_nr == 4)*40;
    max_lon = -45; % 25; %  + (magn_comp_nr == 4)*40;
    map_proj = 'stereo'; % 'Mercator';      
elseif mig_sys == 2 % || mig_sys == 6
    min_lat = -10;
    max_lat = 72.5;
    min_lon = -22.5; % + (magn_comp_nr == 4)*40;
    max_lon = 65; % 25; %  + (magn_comp_nr == 4)*40;
    map_proj =  'Mercator'; % 'stereo'; %
    map_lat_cntr = 35; % (min_lat + max_lat)/2;
    map_lon_cntr = 0; % (min_lon + max_lon)/2;
    FLatLims = [90 90]; % [30 40]; % [0
elseif mig_sys == 3 % mig_sys = 3; [-190 -150 60 71], ...
    min_lat = 7.5; % -5 % ;
    max_lat = 90;
    min_lon = -180; % + (magn_comp_nr == 4)*40;
    max_lon = 180; % 25;
%     map_proj = 'stereo'; %  'eqdazim'; %  'Robinson'; % 'gnomonic'; % 'mollweid'; % 'ortho'; % 'stereo'; % 'Mercator'; %   % 'stereo'; % 'Mercator'; %  %    %
     map_proj = 'ortho'; % 'Mercator';  % 'stereo'; %  
    FLatLims = [40 70]; %[60 70]; %[0
    map_lat_cntr = 55; %  90; %(min_lat + max_lat)/2;
    map_lon_cntr = 60; % (min_lon + max_lon)/2;
elseif mig_sys == 4 % mig_sys = 3; [-190 -150 60 71], ...
    min_lat = -15;
    max_lat = 80;
    min_lon = 0; % + (magn_comp_nr == 4)*40;
    max_lon = 225; % 360; % 25;
    map_proj =  'Mercator'; % 'stereo'; %
        map_lat_cntr = 50; % (min_lat + max_lat)/2;
    map_lon_cntr = 95; % (min_lon + max_lon)/2;
    FLatLims = [80 80]; % [90 90]; % [0
elseif mig_sys == 6 % mig_sys == % mig_sys = 4; 
    min_lat = 2.5;
    max_lat = 57.5;
    min_lon = -25; % + (magn_comp_nr == 4)*40;
    max_lon = 20; % 25;
    map_proj =   'stereo'; % 'Mercator';  %
elseif mig_sys == 8
    min_lat = -15;
    max_lat = 75;
    min_lon = -150; % + (magn_comp_nr == 4)*40;
    max_lon = -45; % 25; %  + (magn_comp_nr == 4)*40;
    map_proj = 'stereo'; % 'Mercator';   
elseif mig_sys == 5 
    min_lat = 10; % -5 % 
    max_lat =  85; %70 % 
    min_lon = -85; % + (magn_comp_nr == 4)*40;
    max_lon = 40; % 
    map_proj =  'stereo'; % 'Mercator';  %
    map_lat_cntr = 55; % (min_lat + max_lat)/2;
    map_lon_cntr = -10; % (min_lon + max_lon)/2;
    FLatLims = [40 40]; % [35 35]; %[0[30 30]; %
elseif mig_sys == 9 
    min_lat = -5 % 10;
    max_lat =  80; % 70 %
    min_lon = 0; % + (magn_comp_nr == 4)*40;
    max_lon = 220; % 
    map_proj =  'stereo'; % 'Mercator';  %
elseif mig_sys > 9  % Blackpolls
    min_lat = 7.5 % 10;
    max_lat =  70; % 70 %
    min_lon = 170; % + (magn_comp_nr == 4)*40;
    max_lon = 310; % 
%     map_proj =  'Mercator'; % 'stereo'; %
    map_lat_cntr = 42.5; % (min_lat + max_lat)/2;
    map_lon_cntr = 265; % (min_lon + max_lon)/2;
    FLatLims = [52.5 50]; % [90 90]; % [0
    map_proj =  'stereo'; % 'Mercator';  %
else % == 7 just Greenland knot tracks
    min_lat = 35 % 10;
    max_lat =  85; % 70 %
    min_lon = -90; % + (magn_comp_nr == 4)*40;
    max_lon = 20; % 
    map_proj =  'stereo'; % 'Mercator';  %'stereo'; %'ortho'; %
end

zug_kn_cols = {'mx','wx','rx'};

bg_clr = 0.825;
set(gcf,'color',[bg_clr bg_clr bg_clr]); % [0.9 0.9 0.9]);
land = shaperead('landareas','UseGeoCoords', true);

if mig_sys > 9 
  MkSzElls = 1.5;
  MkSzMns = 15;
  N_cntr = 300;
  Alf = 1;% 0.5;
elseif mig_sys == 2
  MkSzElls = 1.5;
  MkSzMns = 15;
  N_cntr = 350;
  Alf = 1;% 0.5; 
else
  MkSzElls = 1.5;
  MkSzMns = 25;
  N_cntr = 350; % 600; % 
  Alf = 1;% 0.5;
end

if N_mig_Sims >= 116
    
     skip = 10;
     
else
    
     skip = 1;
     
end

mrkr = '>';
mrks = {'o','h','s','d'};

% num Dep and Arr Polys
nDs = size(succs_DepPs,1);
nAs = size(succs_DepPs,2);


max_Arr_Ds = round(100*max(succs_DepPs(:)));
    
for iD = 1:nDs 

%     succ_Ds = squeeze(max_Arr_Ds*succs_DepPs(iD,:,:)); 
%     min_Arr_Ds = min(succ_Ds(succ_Ds(:)>0));

%     range_Arr_Ds = max_Arr_Ds - min_Arr_Ds;
    
    figure
lnd_clr = 0.75;
geoshow(land, 'FaceColor', [lnd_clr lnd_clr lnd_clr],'LineWidth',0.2) %[0.85 0.85 0.85]
     ax = worldmap([min_lat max_lat],[min_lon max_lon]); % ,'Origin',[40 50]);
    
    if strcmp(map_proj,'Mercator')
        setm(ax,'mapprojection',map_proj)
        map_lat_cntr = (min_lat + max_lat)/2;
        map_lon_cntr = (min_lon + max_lon)/2;
        FLatLims = [100 100]; %[50 50]; %[0
    elseif strcmp(map_proj,'ortho')
        setm(ax,'mapprojection',map_proj,'Origin',[map_lat_cntr map_lon_cntr])
    else
        setm(ax,'mapprojection',map_proj)
       ax = axesm ('stereo', 'Frame', 'on', 'Grid', 'off','Origin',[map_lat_cntr map_lon_cntr 0],'FlineWidth',1);
       ax.PositionConstraint = 'innerposition';
       axis('off')
       setm(gca,'FLatLimit',FLatLims)
    end

    geoshow(land, 'FaceColor', [0.85 0.85 0.85],'LineWidth',0.2) %
    hold
    
% if plot_zugs
    
    for iz = 1:n_zugs
%         for jA = 1:nAs 
            for iY = 1:skip:N_mig_Sims
                
                % colour to Arr success
                pArr_ijY = round(100*succs_DepPs(iD,iY));
                
                if succs_DepPs(iD,iY) > 0
                        
                    stdLat = n_stds*std_lat_zk_DepPs(iz,iD,iY);
                    stdLon = n_stds*std_lon_zk_DepPs(iz,iD,iY);
                    if stdLat >= stdLon
                        ecci = axes2ecc(stdLat,stdLon);
                        [elat,elon] = ellipse1(mn_lat_zk_DepPs(iz,iD,iY), ...
                        squeeze(mn_lon_zk_DepPs(iz,iD,iY)),[stdLat ecci],0, ...
                        [],[],[],NptsEll);
                    else
                        ecci = axes2ecc(stdLon,stdLat); 
                        [elat,elon] = ellipse1(mn_lat_zk_DepPs(iz,iD,iY), ...
                        mn_lon_zk_DepPs(iz,iD,iY),[stdLon ecci],90, ...
                        [],[],[],NptsEll);
                    end

            %         [lati,loni] = scircle1(mn_lat_zks(iz,iA,1:skip:end)*180/pi, ...
            %             mn_lon_zks(iz,iA,1:skip:end)*180/pi,rad);
                    hh1 = scatterm(elat,elon,MkSzElls,cmap(iY,:),mrkr,'fill');
                    hh2 = scatterm(mn_lat_zk_DepPs(iz,iD,iY), ...
                      mn_lon_zk_DepPs(iz,iD,iY),MkSzMns, ...
                      cmap(iY,:),mrkr,'fill');
%                     set(get(hh2,'Children'),'MarkerEdgeColor','k','LineWIdth',0.01)

                end
            end
            
%         end
    end

    thous = min(sum(stoppedAndArrived),200);
    thous_arr = find(stoppedAndArrived, thous,'first');
    MkSz = 15; % 2.5;
        for ib = 1:thous
            id_b = thous_arr(ib);
            if all_zg_nrs(thous_arr(ib),end) == 1
                sc = scatterm(blats(id_b,1)*180/pi, ...
                blons(id_b,1)*180/pi,MkSz, ...
                'MarkerEdgeColor','w','MarkerFaceColor',plot_col);
    %             'bo','fill');
            
    %             sc.Children.MarkerEdgeColor = 'w'; 
    %             sc.Children.MarkerFaceAlpha = 0.75; 
                sc.Children.LineWidth = 0.05; 
            else
                sc = scatterm(blats(id_b,1)*180/pi, ...
                blons(id_b,1)*180/pi,MkSz, ...
                'MarkerEdgeColor',plot_col,'MarkerFaceColor','w');
    %             'wo','fill');
            
    %             sc.Children.MarkerEdgeColor = 'b'; 
    %             sc.Children.MarkerFaceAlpha = 0.75; 
                sc.Children.LineWidth = 0.05; 
            end
        end


% cc = colorbar;
% title(cc,'year','FontSize',10)
% xlabel('Longitude (^o)','FontSize',12)
% ylabel('Latitude (^o)','FontSize',12)
% colormap(brewermap([],'*YlGnBu')) % colormap winter

% now for deps and Arrs
% if plot_arrs ||  plot_deps  

%     for jA = 1:nAs 
        
%     end
    
    % if new_fig
    mlabel('off'); plabel('off'); gridm('off')

    cc = colorbar('FontSize',14);
    colormap(cmap) % cmap = colormap(brewermap(N_mig_Sims,'*RdYlBu')); % *YlGn')); % '*YlOrRd')); %Reds')); % 
    caxis([yr_start yr_start+N_mig_Sims-1])
%     tot_Arr_Ds = max(sum(succ_Ds,2));
%     cmap = colormap(brewermap(100,'*YlOrRd')); %'*Reds')); % *RdYlBu')); % 
%     caxis([0 max_Arr_Ds]) % min_Arr_Ds
%     title(cc,'Arrival (%)','FontSize',10)
    title(cc,'Year','FontSize',16)

    
end


% end
% xlabel('Longitude (^o)','FontSize',12)
% ylabel('Latitude (^o)','FontSize',12)
% colormap(brewermap([],'*YlGnBu')) % colormap winter

%% for wheatear plots

% keyboard % save Fig for efficiency with any revisions

figure;

thous = min(sum(stoppedAndArrived),5000);
thous_arr = find(stoppedAndArrived, thous,'first');
MkSz = 32.5;
    for ib = 1:250
        id_b = thous_arr(ib);
        if all_zg_nrs(thous_arr(ib),end) == 1
            sc = scatterm(blats(id_b,1)*180/pi, ...
            blons(id_b,1)*180/pi,MkSz, ...
            'MarkerEdgeColor','w','MarkerFaceColor',plot_col);
%             'bo','fill');
        
%             sc.Children.MarkerEdgeColor = 'w'; 
%             sc.Children.MarkerFaceAlpha = 0.75; 
            sc.Children.LineWidth = 0.05; 
        else
            sc = scatterm(blats(id_b,1)*180/pi, ...
            blons(id_b,1)*180/pi,MkSz, ...
            'MarkerEdgeColor',plot_col,'MarkerFaceColor','w');
%             'wo','fill');
        
%             sc.Children.MarkerEdgeColor = 'b'; 
%             sc.Children.MarkerFaceAlpha = 0.75; 
            sc.Children.LineWidth = 0.05; 
        end
    end

%     keyboard % save Fig for efficiency revisions

% Sp & Port from Thorup 2006 whtrs AUK
SP_Lt = [41.5 42 40.75 43 44.5 44.75 46.25];
SP_Ln = [-7 -7.5 0.25 -.5 -2 -1 -1.5];
MkSzRng = 45; % 

% sc = scatterm(SP_Lt, SP_Ln,MkSzRng,'MarkerEdgeColor',data_col,'MarkerFaceColor','w');
% sc.Children.MarkerFaceAlpha = 0.5; 
% sc = scatterm(SP_Lt, SP_Ln, MkSzRng,data_col,'Marker','+' ...
%     ,'LineWidth',1);
% sc.Children.MarkerFaceAlpha = 0.5; 

% others in EU via Heiko et al 2016 Behav Ecol protandry
HS_Lt = [54+11/60 53.380977 52.953602 58+5/60 51.06 ];
HS_Ln = [7+55/60 -3.225359 1.073456 6+47/60 2.39 ];

% sc = scatterm(HS_Lt, HS_Ln,MkSzRng,'MarkerEdgeColor',data_col,'MarkerFaceColor','w');
% sc.Children.MarkerFaceAlpha = 0.5; 
% sc = scatterm(HS_Lt, HS_Ln,MkSzRng,data_col,'Marker','x' ...
%     ,'LineWidth',1);
% sc.Children.MarkerFaceAlpha = 0.5; 

% Koksijde <- c(2.39, 51.06) #Belgium
% Helgoland <- c(7+55/60,54+11/60)
% Hilbre <- c(-3.225359,53.380977)
% Norfolk <- c(1.073456, 52.953602)
% Ottenby <- c(16+23.944/60,56+11.838/60)
% Rybachy <- c(20+51/60,55+9/60)
% Ventotene <- c(13.41832, 40.78791)
% Portovene <- c(9.83877, 44.03798)
% Anacapri <- c(14.22262, 40.55649)
% Lista <- c(6+47/60, 58+5/60)

% Then Jamie found one (!) EURING bird with long wing in fall
JM_Lt =  43.6464;
% JM_Ln = -1.4319;
% sc = scatterm(JM_Lt, JM_Ln,MkSzRng,'MarkerEdgeColor',data_col,'MarkerFaceColor','w');
% sc.Children.MarkerFaceAlpha = 0.5; 
% sc = scatterm(JM_Lt, JM_Ln,MkSzRng,data_col,'Marker','x' ...
%     ,'LineWidth',1);
% sc.Children.MarkerFaceAlpha = 0.5; 

% keyboard
% geolocator trajectory Iqaluit Africa Bairlein 2012

[IqUK_Lt,IqUK_Ln] = track('rh',[63.7 56],[-68.5 -6]);
[UKAf_Lt,UKAf_Ln] = track('rh',[56 25],[-6 -12.5]);

% IqAf_Lt = [63.7 59 56 40 25];
% IqAf_Ln = [-68.5 -45 -6 -9 -12.5];
hold
sc1 = plotm(IqUK_Lt, IqUK_Ln,'Color',data_col,'LineStyle','--','LineWidth',1.25);
sc2 = plotm(UKAf_Lt, UKAf_Ln,'Color',data_col,'LineStyle','--','LineWidth',1.25);

% and now for Ottoson Gld ringing
[Gd_SP_Lt,Gd_SP_Ln] = track('rh',[64.33 43.25],[-52 -2.75]);
[SPAf_Lt,SPAf_Ln] = track('rh',[43.25 25],[-2.75 -12.5]);

% IqAf_Lt = [63.7 59 56 40 25];
% IqAf_Ln = [-68.5 -45 -6 -9 -12.5];
sc1 = plotm(Gd_SP_Lt, Gd_SP_Ln,'Color',data_col,'LineStyle','--','LineWidth',1.25);
sc2 = plotm(SPAf_Lt, SPAf_Ln,'Color',data_col,'LineStyle','--','LineWidth',1.25);