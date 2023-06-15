function [meanps, medianps, lq_ps, uq_ps] = bin_mean_prbs(ps,bins, ...
    av_dbins,pbins)

for ib = 1:numel(av_dbins)
    
    try
    is_bin =  bins == av_dbins(ib);
    meanps(ib) = mean(ps(is_bin));
    if std(ps(bins == ib)) > 0 && numel(ps(is_bin)) > 3
            qs = quantile(ps(is_bin),pbins);
            lq_ps(ib) = qs(1);
             medianps(ib) = qs(2);
             uq_ps(ib) = qs(3);
    else
        lq_ps(ib) = 0;
        uq_ps(ib) = 0;
        medianps(ib) = median(ps(is_bin));
    end
    catch 
        keyboard
    end
end

% for ib = numel(bins):79
%     
%      meanps(ib) = 0;
%     
% end   