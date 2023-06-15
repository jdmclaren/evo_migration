function [nSubPolys , LonSub_poly, LatSub_poly] = ...
det_subPolys(Lon_poly,Lat_poly)

% determine number of disjoint subpolys 
nSubPolys = sum(isnan(Lon_poly));

for iSubPol = 1:nSubPolys

    LatSub_poly{iSubPol} = Lat_poly((iSubPol-1)*6 + (1:5));
    LonSub_poly{iSubPol} = Lon_poly((iSubPol-1)*6 + (1:5));

end