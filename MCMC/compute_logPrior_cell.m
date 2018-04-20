function logPrior_cell = compute_logPrior_cell(Prior_Param, theta_cellk)
% % diffusivity components: prop 1/detK^2 [Actually performed badly compared with the Uni one!]
if Prior_Param.PriorType == 1
    % Jeffrey's Prior
    % Fields_cell(3:5) = [K11, K22, K12]
    detK = (theta_cellk(7)*theta_cellk(8))^2;
    
    % Jeffrey's Prior
    logPrior_cell = -0.5*(2+1).*log(detK);

elseif Prior_Param.PriorType == 2
    b = theta_cellk(1);
    nu = theta_cellk(3:4);
    tr = theta_cellk(6);
    sigma = theta_cellk(7:8);
    
    b_BND = Prior_Param.b;
    nu_BND = Prior_Param.nu;
    tr_BND = Prior_Param.tr;
    sigma_BND = Prior_Param.sigma;

    % Check if lies in the range
    b_TF = (b < b_BND(1)) || (b > b_BND(2));
    nu_TF = any(nu < nu_BND(1)) || any(nu > nu_BND(2));
    tr_TF = any(tr < tr_BND(1)) || any(tr > tr_BND(2));
    sigma_TF = any(sigma < sigma_BND(1)) || any(sigma > sigma_BND(2));
    
    % Uniform Prior
    if (b_TF || nu_TF || tr_TF || sigma_TF)
        %logPrior_cell = -1/(eps*10);
        logPrior_cell = -inf;
    else
        logPrior_cell = 0;
    end
end

end
