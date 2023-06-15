function [interval] = randinterval(low, high,N_inds)
% returns a random number within specified interval

% previously used for random start locations from polygon


interval = low + ( (high - low) * rand(N_inds,1) );


% x = datasample(locs,N_inds,'Replace',false,'Weights',w);
% dif = high-low;
% locs = low + 0:dif/(N_inds-1):high;
% 
% w = locs-min(locs)+1;
% keyboard