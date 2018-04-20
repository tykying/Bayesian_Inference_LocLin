function AutoCorr = Compute_ShiftCorr_R2Data(dWx, dWy, AutoCorr_MaxLength)

% 4 components: (x_unshifted, y_unshifted) x (x_shifted, y_shifted)
% 1st: corrcoef(x_unshifted, x_shifted)
% 2nd: corrcoef(x_unshifted, y_shifted)
% 3rd: corrcoef(y_unshifted, x_shifted)
% 4th: corrcoef(y_unshifted, y_shifted)
% choose a SampInt such that all these correlations are approximately zero
AutoCorr = zeros(AutoCorr_MaxLength+1, 4);  % Hardcoded = 4

for Dn = 0:AutoCorr_MaxLength
    dWx_unshifted = dWx(1:end-Dn);
    dWy_unshifted = dWy(1:end-Dn);

    dWx_shifted  = dWx(1+Dn:end);
    dWy_shifted  = dWy(1+Dn:end);

    assert(length(dWx_unshifted) == length(dWx_shifted));
    
    % Manually compute the corrcoef: no idea how the built-in one works for
    % vector data
    
    % Component-wise consideration: correlation between 4 scalar R.V.
    Cxuxs = cov(dWx_unshifted, dWx_shifted);
    Cxuys = cov(dWx_unshifted, dWy_shifted);
    Cyuxs = cov(dWy_unshifted, dWx_shifted);
    Cyuys = cov(dWy_unshifted, dWy_shifted);
    
    Rxuxs = Cxuxs(1,2)/sqrt((Cxuxs(1,1)*Cxuxs(2,2)));
    Rxuys = Cxuys(1,2)/sqrt((Cxuys(1,1)*Cxuys(2,2)));
    Ryuxs = Cyuxs(1,2)/sqrt((Cyuxs(1,1)*Cyuxs(2,2)));
    Ryuys = Cyuys(1,2)/sqrt((Cyuys(1,1)*Cyuys(2,2)));

    R = [Rxuxs, Rxuys, Ryuxs, Ryuys];
    AutoCorr(Dn+1, :) = R;  % Hardcoded: Number of Columns = 4
    
    % Build-in functions
    % dW_unshifted = [dWx_unshifted, dWy_unshifted];
    % dW_shifted = [dWx_shifted, dWy_shifted];
    % R = corrcoef(dW_unshifted, dW_shifted)   % <-- No idea how it works

    % R = corrcoef(dWx_unshifted, dWx_shifted)
    % assert(R(1,2) == Cxx(1,2)/sqrt((Cxx(1,1)*Cxx(2,2))))  % True
end