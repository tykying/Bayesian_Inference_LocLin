% grid_FD = struct('X', X', 'Y', Y');  % Use (x_i, y_j) for (i,j)
function fX_FD = feval_FDBilinear(f_FD, grid_FD, X)
    assert(size(X, 2) == 2);  % Each row = coordinates X(1, :)=[x_1, y_1]
    fX_FD = zeros(size(X, 1), 1);

    X_mesh = grid_FD.X;
    Y_mesh = grid_FD.Y;    
    
    % Ensure in (x_i, y_j) for (i, j) format
    assert(all(X_mesh(:,1) == X_mesh(:,2))); 
    assert(all(diff(X_mesh(:,1)) == diff(X_mesh(:,2)))); 

    assert(all(Y_mesh(1,:) == Y_mesh(2,:))); 
    assert(all(diff(Y_mesh(1,:)) == diff(Y_mesh(2,:))));
    
    grid_boundary_x = X_mesh(:, 1);
    grid_boundary_y = Y_mesh(1, :);
    
    for item_X = 1:size(X, 1)
        [x_CellInd, x_alpha] = search_sorted_nonuniform(X(item_X, 1), grid_boundary_x);
        [y_CellInd, y_alpha] = search_sorted_nonuniform(X(item_X, 2), grid_boundary_y);
        
        alpha = [x_alpha, y_alpha];
        f_FD_clipped = f_FD(x_CellInd:x_CellInd+1, y_CellInd:y_CellInd+1);
        
        fX_FD(item_X) = local_bilinear_interpol(f_FD_clipped, alpha);
    end
end