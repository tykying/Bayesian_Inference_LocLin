% Set-up Fine Fields from theta_mean
function [FineFields, grid_vis_f] = convert_theta_into_FineFields(theta_mean, Mesh_Struct)
% Make sure no gradient is used in diffusivity field
assert(size(theta_mean, 2) == 9)

RectMesh_Param = Mesh_Struct.RectMesh_Param;

theta_mean_ij = convert_cellk_to_cellij( theta_mean, RectMesh_Param );

% Convert to 2D arrays
U_vis = squeeze(theta_mean_ij(:, :, 1));
V_vis = squeeze(theta_mean_ij(:, :, 2));

dUdx_vis = squeeze(theta_mean_ij(:, :, 3));
dUdy_vis = squeeze(theta_mean_ij(:, :, 4));
dVdx_vis = squeeze(theta_mean_ij(:, :, 5));
dVdy_vis = squeeze(theta_mean_ij(:, :, 6));

K_sigma1_vis = squeeze(theta_mean_ij(:, :, 7));
K_sigma2_vis = squeeze(theta_mean_ij(:, :, 8));
K_phi_vis = squeeze(theta_mean_ij(:, :, 9));

[KXX_vis, KYY_vis, KXY_vis] = Kpolar_to_Kcart_vectorised(K_sigma1_vis, K_sigma2_vis, K_phi_vis);

FineGrid_Resolution = 256+1;
[U_f, X_f, Y_f] = griddata_interpolation_pwl(U_vis, dUdx_vis, dUdy_vis, RectMesh_Param, FineGrid_Resolution);
[V_f, X_f, Y_f] = griddata_interpolation_pwl(V_vis, dVdx_vis, dVdy_vis, RectMesh_Param, FineGrid_Resolution);


% Ad-hoc: assume no gradients in Diffusivity Field
[K_sigma1_f, X_f, Y_f] = griddata_interpolation_pwl(K_sigma1_vis, zeros(size(K_sigma1_vis)), zeros(size(K_sigma1_vis)), RectMesh_Param, FineGrid_Resolution);
[K_sigma2_f, X_f, Y_f] = griddata_interpolation_pwl(K_sigma2_vis, zeros(size(K_sigma2_vis)), zeros(size(K_sigma2_vis)), RectMesh_Param, FineGrid_Resolution);
[K_phi_f, X_f, Y_f] = griddata_interpolation_pwl(K_phi_vis, zeros(size(K_phi_vis)), zeros(size(K_phi_vis)), RectMesh_Param, FineGrid_Resolution);

[KXX_f, X_f, Y_f] = griddata_interpolation_pwl(KXX_vis, zeros(size(KXX_vis)), zeros(size(KXX_vis)), RectMesh_Param, FineGrid_Resolution);
[KYY_f, X_f, Y_f] = griddata_interpolation_pwl(KYY_vis, zeros(size(KYY_vis)), zeros(size(KYY_vis)), RectMesh_Param, FineGrid_Resolution);
[KXY_f, X_f, Y_f] = griddata_interpolation_pwl(KXY_vis, zeros(size(KXY_vis)), zeros(size(KXY_vis)), RectMesh_Param, FineGrid_Resolution);

% Set up grid_vis_f; something is missing
[Mesh_f, RectMesh_Param_f] = setup_RectMeshStruct(size(X_f, 1), size(Y_f, 2), RectMesh_Param.SpDis_ranges_min, RectMesh_Param.SpDis_ranges_max);

% Output
grid_vis_f = setup_grid_vis(Mesh_f, RectMesh_Param_f);

FineFields= struct('u', u, 'v', v, 'K_sigma1', K_sigma1_f, ...
                    'K_sigma2', K_sigma2_f, 'K_phi', K_phi_f, ...
                    'KXX', KXX_f, 'KYY', KYY_f, 'KXY', KXY_f);

end
