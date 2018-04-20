function Point_AssoGridInd = assoicate_points_with_GridInd_uniform(x, bin_size, grid_startpt, Nx_cell)
x_Divided = (x-grid_startpt)./bin_size;

Point_AssoGridInd = zeros(2, length(x_Divided));
for i = 1:length(x_Divided)
    if mod(x_Divided(i), 1) == 0
        % Grid Point lying on a boundary
        Point_AssoGridInd(1, i) = (x_Divided(i));
        Point_AssoGridInd(2, i) = (x_Divided(i)+1);
    else
        % Interior Grid Point
        Point_AssoGridInd(:, i) = ceil(x_Divided(i));
    end
end
Point_AssoGridInd(Point_AssoGridInd==0) = 1;
Point_AssoGridInd(Point_AssoGridInd==(Nx_cell+1)) = Nx_cell;
end