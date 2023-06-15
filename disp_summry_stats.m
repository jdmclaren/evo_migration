disp('geometric mean and standard geometric range success (1900:2015)')
geomnS = geomean(successful);
round_gmn_S = round(100*geomnS);
geostdS = geostd(successful);
% disp([num2str(round_gmn_S/100) ' (' ...
%     num2str(round(100*(geomnS/geostdS))/100) ...
%   ' - ' num2str(round(100*(geomnS*geostdS))/100) ')'])
disp([num2str(round_gmn_S/100) ' (' ...
    num2str(round(100*(min(successful)))/100) ...
  ' - ' num2str(round(100*(max(successful)))/100) ')'])

if calibr_comp(1) ~= 2 && calibr_comp(1) ~= 3
    
    mn_alf = mod(circ_mean(succ_heads*pi/180)*180/pi,360);
    std_alf = mod(circ_std(succ_heads*pi/180)*180/pi,360);
    mn_zgs = mod(circ_mean(zugkns*pi/180)*180/pi,360);
    std_zgs = mod(circ_std(zugkns*pi/180)*180/pi,360);
    disp('geometric mean and std headings')
    disp('geometric mean and std zugkn headings')
    if std_alf < 1
        disp([num2str(round(mn_alf)) ' (' ...
            num2str(round(10*(std_alf))/10) ')'])
        
    else
        
            disp([num2str(round(mn_alf)) ' (' ...
        num2str(round((std_alf))) ')'])
        
    end
    
    if std_zgs < 1
        
        disp([num2str(round(mn_zgs)) ' (' ...
            num2str(round(10*(std_zgs))/10) ')'])
        
    else
        
        disp([num2str(round(mn_zgs)) ' (' ...
            num2str(round((std_zgs))) ')'])
        
    end        
    
else
    
    for iz = 1:n_zugs+1
        
        disp(['inclination projection ' num2str(iz)])
        % take abs(coeff) forcalibr_comp(1) == 3 since it has sign(init heads)
        inclns = (calibr_comp(1)==3)*abs(atan(coeff_succ(:,iz))) + ...
            (calibr_comp(1)==4)*acos(coeff_succ(:,iz));
        mn_gam = circ_mean(inclns)*180/pi;
        std_gam = circ_std(inclns)*180/pi;
        if std_gam < 1
            disp([num2str(round(mn_gam)) ' (' ...
                num2str(round(10*(std_gam))/10) ')'])
        else
            disp([num2str(round(mn_gam)) ' (' ...
            num2str(round(std_gam)) ')'])
        end

    end

end