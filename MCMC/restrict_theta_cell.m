function [theta_prop_cell] = restrict_theta_cell(theta_prop_cell)
% Restrict positive parameters
theta_prop_cell(:, [1,7:8]) = abs(theta_prop_cell(:, [1,7:8]));

end
