%function [CrossSum, CrossSum_NSample] = compute_LagCrossCorrSum(ts1, ts2)
%%% Do not output CrossSum_NSample: can be a prior told by the length of ts1

function CrossSum = compute_LagCrossCorrSum(ts1, ts2)
    % assert(length(ts1) == length(ts2));
    ts_len = length(ts1);
    
    CrossSum = zeros(ts_len, 1);
    % CrossSum_NSample = zeros(ts_len, 1);    % Equal to [ts_len:-1:1]'
    
    for lag = 0:(ts_len-1)
        NSample_lag = ts_len-lag;
        
        ts1_u = ts1(1:NSample_lag);     %u: unshifted
        ts2_s = ts2(1+lag:ts_len);      %s: shifted
 
        CrossSum(lag+1) = sum(ts1_u.*ts2_s);
        % CrossSum_NSample(lag+1) = NSample_lag;  % = length(ts1_u)
    end
    
    % assert(all(CrossSum_NSample == [ts_len:-1:1]'));
end
