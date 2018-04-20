% Set-up Fine Fields from theta_Fields
function [FineFields, grid_vis_f] = convert_theta_Fields_into_FineFields(theta_Fields, Mesh_Struct)
% Make sure no gradient is used in diffusivity field
assert(size(theta_Fields, 2) == 9)

RectMesh_Param = Mesh_Struct.RectMesh_Param;

theta_Fields_ij = convert_cellk_to_cellij( theta_Fields, RectMesh_Param );

% Convert to 2D arrays
U_vis = squeeze(theta_Fields_ij(:, :, 1));
V_vis = squeeze(theta_Fields_ij(:, :, 2));

dUdx_vis = squeeze(theta_Fields_ij(:, :, 3));
dUdy_vis = squeeze(theta_Fields_ij(:, :, 4));
dVdx_vis = squeeze(theta_Fields_ij(:, :, 5));
dVdy_vis = squeeze(theta_Fields_ij(:, :, 6));

KXX_vis = squeeze(theta_Fields_ij(:, :, 7));
KYY_vis = squeeze(theta_Fields_ij(:, :, 8));
KXY_vis = squeeze(theta_Fields_ij(:, :, 9));

[KXX_vis, KYY_vis, KXY_vis] = Kpolar_to_Kcart_vectorised(KXX_vis, KYY_vis, KXY_vis);

FineGrid_Resolution = 60;
%FineGrid_Resolution = 50;
%FineGrid_Resolution = 128+1;   %<-this will create multiple-valued nodes
[U_f, X_f, Y_f] = griddata_interpolation_pwl(U_vis, dUdx_vis, dUdy_vis, RectMesh_Param, FineGrid_Resolution);
[V_f, X_f, Y_f] = griddata_interpolation_pwl(V_vis, dVdx_vis, dVdy_vis, RectMesh_Param, FineGrid_Resolution);

% Remark:
% if some grid points of FineGrid overlap with the coarse grid boundary,
% U_f, V_f will have an extra dimension.
% E.g. take FineGrid_Resolution = 128+1;
    

%% Ad-hoc: for U1K0 assume no gradients in Diffusivity Field
[KXX_f, X_f, Y_f] = griddata_interpolation_pwl(KXX_vis, zeros(size(KXX_vis)), zeros(size(KXX_vis)), RectMesh_Param, FineGrid_Resolution);
[KYY_f, X_f, Y_f] = griddata_interpolation_pwl(KYY_vis, zeros(size(KYY_vis)), zeros(size(KYY_vis)), RectMesh_Param, FineGrid_Resolution);
[KXY_f, X_f, Y_f] = griddata_interpolation_pwl(KXY_vis, zeros(size(KXY_vis)), zeros(size(KXY_vis)), RectMesh_Param, FineGrid_Resolution);

% Set up grid_vis_f; something is missing
[Mesh_f, RectMesh_Param_f] = setup_RectMeshStruct(size(X_f, 1), size(Y_f, 2), RectMesh_Param.SpDis_ranges_min, RectMesh_Param.SpDis_ranges_max);

% Output
grid_vis_f = setup_grid_vis(Mesh_f, RectMesh_Param_f);

FineFields= struct('u', U_f, 'v', V_f, 'KXX', KXX_f, 'KYY', KYY_f, 'KXY', KXY_f);

end
