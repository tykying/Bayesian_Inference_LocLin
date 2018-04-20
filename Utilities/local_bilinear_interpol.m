function [q, Gradq] = local_bilinear_interpol(Q, alpha)
%% Bilinear interpolation of estimated kappa off the grid

% Original for 2D;
% Now written in a manner that is applicable to 1D also 

% Specifications:

% Input:
% cell_dx = [cell_dx; cell_dy]: grid size
% alpha = [alpha_x; alpha_y]: located relative position

% Q(i, j) = Q at (x=x_i, y=y_j)
% GradQ = Q(x=x_i+1/2, y=y_j+1/2) matrix of size Nx_cell x Ny_cell x 2; dQ_x(i, j) = GradQ(x=x_i+1/2, y=y_j+1/2, 1)

% The structured grid function is defined also on the boundaries [0, 1]*[0, 1]

%% Computation

assert(all(size(Q) == [2, 2]))

% c.f. Wiki Bilinear Nonlinear
% x_N = [1-alpha(1), alpha(1)];
% y_N = [1-alpha(2); alpha(2)];
% 
% q = x_N*Q*y_N;

% Just take GradQ as the gradient
a00 = Q(1,1);
a10 = Q(2,1)-Q(1,1);
a01 = Q(1,2)-Q(1,1);
a11 = Q(2,2)+Q(1,1)-(Q(2,1)+Q(1,2));

q =  a00 + a10*alpha(1) + a01*alpha(2) + a11*alpha(1)*alpha(2);
Gradq = [a10 + a11*alpha(2), a01 + a11*alpha(1)];
end
