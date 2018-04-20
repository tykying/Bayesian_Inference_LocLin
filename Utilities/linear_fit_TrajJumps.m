function [theta_KPolar, Fields, GradxFields, GradyFields] = linear_fit_TrajJumps(TrajJumps)
NVars = 2 + 4 + 3;  % Hardcoded
% u, v; gradu(2), gradv(2); sigma1, sigma2, phi

% Assume x0, y0 centred to (0, 0)

Ncells = length(TrajJumps);

Fields = zeros(Ncells, 5);
GradxFields = zeros(Ncells, 5);
GradyFields = zeros(Ncells, 5);
theta_KPolar = zeros(Ncells, NVars);

% PWL Mean Velocity and Diffusivity
parfor cell_k = 1:Ncells  % Parallelisable but seems takes even longer for data traffic
    if (mod(log2(cell_k), 4) == 0)
        disp(['initialise_theta: Fitting Jumps in cell', num2str(cell_k)]);
    end
    
    TrajJumps_cellk = TrajJumps(cell_k);
    
    % Fit a constant plane to the velocity field
    %alpha_CellData = [TrajJumps_cellk.alphax, TrajJumps_cellk.alphay];
    x_CellData = TrajJumps_cellk.x0;  % Shifted wrt CellCen already
    y_CellData = TrajJumps_cellk.y0;
    diffx_CellData = TrajJumps_cellk.diffx;
    diffy_CellData = TrajJumps_cellk.diffy;
    h_CellData = TrajJumps_cellk.h;
    
    u_CellData = diffx_CellData./h_CellData;
    v_CellData = diffy_CellData./h_CellData;
    
    % Fit a piecewise linear plane to the velocity field
    % CellCen = Mesh(cell_k).cell_centre';
    %CellCen = Mesh(cell_k).cell_centre'*0;  % Centred
    CellCen = [0, 0];  % Centred
    
    %   %%% Only when spatial structure is involved
    X_inCell = [x_CellData, y_CellData];
    
    % Decomposing instantaneous velocity by linear fitting
    [u_sf, u_goodness] = fit(X_inCell, u_CellData, 'poly11');
    [v_sf, v_goodness] = fit(X_inCell, v_CellData, 'poly11');
    
    % Variance = (sse of u)/N = 2*kappa/h
    Kappa_uu = 0.5*u_goodness.sse*mean(h_CellData)/length(h_CellData);
    Kappa_vv = 0.5*v_goodness.sse*mean(h_CellData)/length(h_CellData);
    
    %         % For debugging
    %         uv_inCell = [u_CellData, v_CellData];
    %         UV_inCell = [feval(u_sf, X_inCell), feval(v_sf, X_inCell)];
    %         uv_misfit = uv_inCell - UV_inCell;
    %         uv_sse = sum(uv_misfit.^2, 1);  % =[u_goodness.sse, v_goodness.sse]
    %         uv_mse = uv_sse/length(h_CellData);
    %         uv_var = mean(uv_inCell.^2, 1) - mean(uv_inCell, 1).^2;
    %         mean(uv_inCell, 1)
    %         var(uv_inCell, 1)  % < uv_mse; But of same magnitude
    %         [u_goodness.sse, v_goodness.sse] ./length(h_CellData);
    
    % Obtain cell-averaged values
    U_fit = [feval(u_sf, CellCen), feval(v_sf, CellCen)];
    
    Gradu_fit = [u_sf.p10, u_sf.p01];
    Gradv_fit = [v_sf.p10, v_sf.p01];
    
    if Kappa_uu >= Kappa_vv
        Sigma_fit = [sqrt(Kappa_uu), sqrt(Kappa_vv), 0];
    else
        Sigma_fit = [sqrt(Kappa_vv), sqrt(Kappa_uu), pi/2];
    end
    
    Fields_cellk = [U_fit, Sigma_fit];
    GradxFeilds_cellk = [u_sf.p10, v_sf.p10, 0, 0, 0];
    GradyFeilds_cellk = [u_sf.p01, v_sf.p01, 0, 0, 0];
    
    Fields(cell_k, :) = Fields_cellk;
    GradxFields(cell_k, :) = GradxFeilds_cellk;
    GradyFields(cell_k, :) = GradyFeilds_cellk;
    theta_KPolar(cell_k, :) = [U_fit, u_sf.p10, u_sf.p01, v_sf.p10, v_sf.p01, Sigma_fit];
end