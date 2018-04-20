% Correlates with the initial data only
function CrossSum = compute_LagCrossCorrSum_init(ts1, ts2)
    % assert(length(ts1) == length(ts2));
    ts_len = length(ts2);
    CrossSum = zeros(ts_len, 1);
    % CrossSum_NSample = zeros(ts_len, 1);    % Equal to [ts_len:-1:1]'
    
    CrossSum = ts1(1).*ts2;
    
    % assert(all(CrossSum_NSample == [ts_len:-1:1]'));
end
