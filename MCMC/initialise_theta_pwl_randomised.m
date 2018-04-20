%% Take the shortcut to compute the sample mean; piecewise constant diffusivity field
function [theta_Init, theta_Scale] = initialise_theta_pwl_randomised(MCMC_Param, Mesh, TrajJumps)

% Initialisation
Ncells = length(Mesh);

if (strcmp(MCMC_Param.InfScheme, 'Local') || ...
        strcmp(MCMC_Param.InfScheme, 'EA_PWC') || ...
        strcmp(MCMC_Param.InfScheme, 'EA_PWL') || ...
        strcmp(MCMC_Param.InfScheme, 'LocLinPolar') || ...
        strcmp(MCMC_Param.InfScheme, 'LocIsoK') )
    NVars = 2 + 4 + 3;  % Hardcoded
    % u, v; gradu(2), gradv(2); sigma1, sigma2, phi
    
    theta_Init = zeros(Ncells, NVars);
    theta_Scale = zeros(Ncells, NVars);
    
    % PWL Mean Velocity and Diffusivity
    for cell_k = 1:Ncells  % Parallelisable but seems takes even longer for data traffic
        if (mod(log2(cell_k), 2) == 0)
            disp(['initialise_theta: Fitting Jumps in cell', num2str(cell_k)]);
        end
        
        % Fit a constant plane to the velocity field
        %alpha_CellData = [TrajJumps(cell_k).alphax, TrajJumps(cell_k).alphay];
        x_CellData = TrajJumps(cell_k).x0;  % Shifted wrt CellCen already
        y_CellData = TrajJumps(cell_k).y0;
        diffx_CellData = TrajJumps(cell_k).diffx;
        diffy_CellData = TrajJumps(cell_k).diffy;
        h_CellData = TrajJumps(cell_k).h;
        
        u_CellData = diffx_CellData./h_CellData;
        v_CellData = diffy_CellData./h_CellData;
        
        % Fit a piecewise linear plane to the velocity field
        % CellCen = Mesh(cell_k).cell_centre';
        CellCen = Mesh(cell_k).cell_centre'*0;  % Centred
        
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
        
        % Save to output
        if (strcmp(MCMC_Param.InfScheme, 'Local') || ...
                strcmp(MCMC_Param.InfScheme, 'EA_PWC') || ...
                strcmp(MCMC_Param.InfScheme, 'EA_PWL') )
            % 1) Scale = Magnitude of initial guess
            theta_Scale(cell_k, :) = abs([U_fit, Gradu_fit, Gradv_fit, Sigma_fit]);
            theta_Scale(cell_k, 3) = 0.5*abs(Gradu_fit(1)-Gradv_fit(2));
            theta_Scale(cell_k, 6) = 0.5*abs(Gradu_fit(1)-Gradv_fit(2));
            
            % 2) Setup theta_Init
            theta_Init(cell_k, 1:6) = theta_Scale(cell_k, 1:6) .* randn([1,6]);
            
            % Incompressiblity
            dudx = theta_Init(cell_k, 3);
            dvdy = theta_Init(cell_k, 6);
            mag = 0.5*abs(dudx-dvdy);

            theta_Init(cell_k, 3) = sign(dudx)*mag;
            theta_Init(cell_k, 6) = -sign(dudx)*mag;   
            
            % Flatten the mean velocities, slope
            theta_Init(cell_k, 1:2) = theta_Init(cell_k, 1:2)*1E-2;
            theta_Init(cell_k, 3:6) = theta_Init(cell_k, 3:6)*1E-4;
            
            % S.P.D. kappa
            theta_Init(cell_k, 7:8) = sort(theta_Scale(cell_k, 7:8) .* rand([1,2]), 'descend');  % Guaranteed positive
            %theta_Init(cell_k, 9) = -0.5*pi + rand(1)*pi;
            theta_Init(cell_k, 9) = rand(1)*pi*0.5;

%             % Double check the initial conditions
%             theta_Init_restricted = restrict_theta_cell_U1K0(theta_Init(cell_k, :));
%             assert(all(theta_Init_restricted == theta_Init(cell_k, :)));
            
        elseif (strcmp(MCMC_Param.InfScheme, 'LocLinPolar'))
            % 1) Scale = Magnitude of initial guess
            theta_Scale(cell_k, 1) = norm(U_fit, 2);    % Magnitude of b
            theta_Scale(cell_k, 2) = atan(U_fit(2)/U_fit(1));
            
            % Use the mean of A11 and -A22 to define A_phi2 and A_sigmahS
            A_phi2_fromA11 = atan( 2*Gradu_fit(1)/(Gradu_fit(2)+Gradv_fit(1)));
            A_phi2_fromA22 = atan(-2*Gradv_fit(2)/(Gradu_fit(2)+Gradv_fit(1)));
            %
            % A_phi2 = 0.5*(A_phi2_fromA11+A_phi2_fromA22);
            % A_sigmahS = sign(Gradu_fit(1)) * 0.5*(abs(Gradu_fit(1))+abs(Gradv_fit(2))) /sin(A_phi2);   % Half Sum of A_sigma_1 and A_sigma_2
            
            % Defining A_sigmahD only requires A12 and A21
            A_sigmahD = 0.5*(Gradu_fit(2) - Gradv_fit(1));   % Half difference of A_sigma_1 and A_sigma_2: by matching with 0.5*(A21-A12)
            
            A_sigmahScos2phi = 0.5*(Gradu_fit(2) + Gradv_fit(1));   % Half Sum of A_sigma_1 and A_sigma_2
            A_sigmahSsin2phi = 0.5*(Gradu_fit(1) + Gradv_fit(2));   % Half Sum of A_sigma_1 and A_sigma_2
            
            A_sigmahS = norm(A_sigmahScos2phi, A_sigmahSsin2phi);
            A_2phi = atan( A_sigmahSsin2phi/A_sigmahScos2phi);
            
            theta_Scale(cell_k, 3:5) = [abs(A_sigmahS), abs(A_sigmahD), A_2phi];
            theta_Scale(cell_k, 6) = A_2phi + 0.5*pi;  % Not used
            
            theta_Scale(cell_k, 7:9) = Sigma_fit;
            
            % 2) Setup theta_Init
            theta_Init(cell_k, 1) = abs( theta_Scale(cell_k, 1) .* randn(1) ); % Magnitude; Need to be positive
            theta_Init(cell_k, 2) = -pi + rand(1)*2*pi;     % Angle for velocity
            theta_Init(cell_k, 3:4) = theta_Scale(cell_k, 3:4) .* randn([1,2]); % Two parameters for A
            theta_Init(cell_k, 5) = -0.5*pi + rand(1)*pi;   % Angle for matrix A
            theta_Init(cell_k, 6) = theta_Init(cell_k, 5) + 0.5*pi;
            
            % Flatten the mean velocities, slope
            theta_Init(cell_k, 1) = theta_Init(cell_k, 1)*1E-2;
            theta_Init(cell_k, 3:4) = theta_Init(cell_k, 3:4)*1E-4;
            
            % S.P.D. kappa
            theta_Init(cell_k, 7:8) = sort(theta_Scale(cell_k, 7:8) .* rand([1,2]), 'descend');  % Guaranteed positive
            theta_Init(cell_k, 9) = -0.5*pi + rand(1)*pi;
            
            % Double check the initial conditions
            theta_Init_restricted = restrict_theta_cell_LocLinPolar(theta_Init(cell_k, :));
            assert(all(theta_Init_restricted == theta_Init(cell_k, :)));
            
        elseif (strcmp(MCMC_Param.InfScheme, 'LocIsoK'))
            % 1) Scale = Magnitude of initial guess
            theta_Scale(cell_k, :) = abs([U_fit, Gradu_fit, Gradv_fit, Sigma_fit]);
            theta_Scale(cell_k, 3) = 0.5*abs(Gradu_fit(1)-Gradv_fit(2));
            theta_Scale(cell_k, 6) = 0.5*abs(Gradu_fit(1)-Gradv_fit(2));
            
            % 2) Setup theta_Init
            theta_Init(cell_k, 1:6) = theta_Scale(cell_k, 1:6) .* randn([1,6]);
            
            % Incompressiblity
            dudx = theta_Init(cell_k, 3);
            dvdy = theta_Init(cell_k, 6);
            mag = 0.5*abs(dudx-dvdy);
            
            theta_Init(cell_k, 3) = sign(dudx)*mag;
            theta_Init(cell_k, 6) = -sign(dudx)*mag;
            
            % Flatten the mean velocities, slope
            theta_Init(cell_k, 1:2) = theta_Init(cell_k, 1:2)*1E-2;
            theta_Init(cell_k, 3:6) = theta_Init(cell_k, 3:6)*1E-4;
            
            % S.P.D. kappa
            theta_Init(cell_k, 7:8) = mean(theta_Scale(cell_k, 7:8)) .* rand(1);  % Guaranteed positive
            theta_Init(cell_k, 9) = -0.5*pi + rand(1)*pi;
            
            % Double check the initial conditions
            theta_Init_restricted = restrict_theta_cell_LocIsoK(theta_Init(cell_k, :));
            assert(all(theta_Init_restricted == theta_Init(cell_k, :)));
            
%             % For debugging
%             fprintf('cell %d: Matching between linear fit and incompressibility: \n', cell_k);
%             fprintf('A_phi2_fromA11 =%f, A_phi2_fromA22 =%f. \n',A_phi2_fromA11, A_phi2_fromA22);
%             fprintf('A_phi2 =%f, A_sigmahS =%f, A_sigmahD =%f. \n',A_2phi, A_sigmahS, A_sigmahD);
%             fprintf('du/dx =%f, dv/dy =%f, du/dy =%f, dv/dx =%f. \n',Gradu_fit(1), Gradv_fit(2), Gradu_fit(2), Gradv_fit(1));
        end
        
    end
end

end