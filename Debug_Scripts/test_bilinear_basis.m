% Incompressible field: Quadratic streamfunction
grid_x = linspace(0, 1, 10);
grid_y = linspace(0, 1, 10);

[X,Y] = meshgrid(grid_x,grid_y)

a11 = 1; 
a22 = 0;
a12 = 1;
b1 = -1;
b2 = 1;

Psi = 0.5*a11*X.*X + 0.5*a22*Y.*Y + a12*X.*Y + b1*X + b2*Y

contourf(X,Y, Psi)