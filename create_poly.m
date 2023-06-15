function [Lon_poly, Lat_poly] = create_poly(Verts)

nVerts = numel(Verts);

Lon_poly = [];
Lat_poly = [];

for iV = 1:nVerts
    
    Lon_1 = Verts{iV}(1);
    Lon_2 = Verts{iV}(2);
    Lon_poly = [Lon_poly Lon_1 Lon_1 Lon_2 Lon_2 Lon_1 NaN];
    Lat_1 = Verts{iV}(3);
    Lat_2 = Verts{iV}(4);
    Lat_poly = [Lat_poly Lat_1 Lat_2 Lat_2 Lat_1 Lat_1 NaN];
    
end
    