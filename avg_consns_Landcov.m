for kLand = 1:1
    
    y = imread(['Consensus_full_class_' num2str(kLand+9) '.tif'], 'tif');
%     y = imread('Consensus_reduced_class_8.tif', 'tif');
    for jLat = 1:146
        for iLon = 1:360

            yy(jLat,iLon,kLand) = median(y((jLat-1)*120+(1:120), ...
                (iLon-1)*120+(1:120)),'all');

        end
    end
    
end