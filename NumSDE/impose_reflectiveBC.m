function [x_new, y_new] = impose_reflectiveBC(x, y, ranges_min, ranges_max)
%% Initialise_particle_positions% Specifications:

% Input:
% x_new, y_new

% Output:
% x0: 2D array, size(1, nparticles); Initial x cooridnates of all particles
% y0: 2D array, size(1, nparticles); Initial y cooridnates of all particles

% Simple reflective BC: also conserve KE
% Symmetrically reflect how much it goes beyond the boundary
%% computation

% Along x-direction
x_new = x;
y_new = y;

% Reflecting along the left boundary
OutB_ind = (x < ranges_min(1));
% refx = ranges_min(1) - x(OutB_ind);
% x_new(OutB_ind) = ranges_min(1) + refx;
x_new(OutB_ind) = 2*ranges_min(1) - x(OutB_ind);

% Reflecting along the right boundary
OutB_ind = (x > ranges_max(1));
% refx = x(OutB_ind) - ranges_max(1);
% x_new(OutB_ind) = ranges_max(1) - refx;
x_new(OutB_ind) = 2*ranges_max(1) - x(OutB_ind);

% Reflecting along the left boundary
OutB_ind = (y < ranges_min(2));
% refy = ranges_min(2) - y(OutB_ind);
% y_new(OutB_ind) = ranges_min(2) + refy;
y_new(OutB_ind) = 2*ranges_min(2) - y(OutB_ind);
 
% Reflecting along the right boundary
OutB_ind = (y > ranges_max(2));
% refy = y(OutB_ind) - ranges_max(2);
% y_new(OutB_ind) = ranges_max(2) - refy;
y_new(OutB_ind) = 2*ranges_max(2) - y(OutB_ind) ;

end