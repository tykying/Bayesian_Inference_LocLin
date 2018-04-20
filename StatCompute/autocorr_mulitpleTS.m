% Plot autocorrelations of multiple realisation of a time-series
function acf = autocorr_mulitpleTS(x, numLags)
    [ts_len, NSample] = size(x);
    
    assert(ts_len > numLags);
    
    acf = zeros(numLags+1, 1);
    
    x_mean = mean(x(:));
    
    for lag = 0: numLags
        x_unlagged = x(1:(end-lag), :)-x_mean;
        x_lagged = x((1+lag):end, :)-x_mean;
        
        assert(all(size(x_unlagged) == size(x_lagged)));
        
        auto_sample = (x_unlagged).*(x_lagged);
        
        acf(lag+1) = mean(auto_sample(:));
    end

    acf = acf./acf(1);
end
