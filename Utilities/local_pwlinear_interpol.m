function [q, Gradq] = local_pwlinear_interpol(cell_dx, alpha, Q, GradQ)
%% Bilinear interpolation of estimated kappa off the grid

% Original for 2D;
% Now written in a manner that is applicable to 1D also 

% Specifications:

% Input:
% cell_dx = [cell_dx; cell_dy]: grid size
% alpha = [alpha_x; alpha_y]: located relative position

% Q = Q(x=x_i+1/2, y=y_j+1/2): value of Q at cell centre
% GradQ = Q(x=x_i+1/2, y=y_j+1/2) matrix of size Nx_cell x Ny_cell x 2; dQ_x(i, j) = GradQ(x=x_i+1/2, y=y_j+1/2, 1)

% The structured grid function is defined also on the boundaries [0, 1]*[0, 1]

%% Computation

% Obtain q based on the cell averaged and prescribed slope

% Orginal 2D Version: 
% q = Q + GradQ(1)*((alpha_x-0.5)*cell_dx) + GradQ(2)*((alpha_y-0.5)*cell_dy);
% q = Q + GradQ(1)*((alpha(1)-0.5)*cell_dx(1)) + GradQ(2)*((alpha(2)-0.5)*cell_dx(2));

assert(all(size(Q) == [1]))
assert(all(size(GradQ) == [1,2]))
assert(all(size(cell_dx) == [1,2]))
assert(all(size(alpha) == [1,2]))

% Obtain q based on the cell averaged and prescribed slope
q = Q + dot(GradQ, (alpha-0.5).*cell_dx);

% Just take GradQ as the gradient
Gradq = GradQ;

end
