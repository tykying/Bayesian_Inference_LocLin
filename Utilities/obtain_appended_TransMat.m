function [TransDen, pi_Stat, RatioRemain, RatioRemainNeigh] = obtain_appended_TransMat(TransMat, Mesh_Struct)
        % Unfold the variables
        RectMesh_Param = Mesh_Struct.RectMesh_Param;
                
        % Sum along row i = total number of Jumps starting from cell_i
        TrajJumps_celli = sum(TransMat, 2);
        
        % Reduce it into Prob Matrix: Pij = Pr(j|i); sum along each row=1
        TransDen = TransMat./TrajJumps_celli;
        
        % Power method <=> Physical method to obtain the stationary distribution (pretty stupid...)
        % pi*Pij = pi
        pi_init = rand(1, size(TransMat, 1));
        pi_old = pi_init;
        for iteration = 1:50000
            pi = pi_old*TransDen;
            
            %disp(max(pi./pi_old))
            pi_old = pi;
        end
        
        pi_Stat = pi/sum(pi);
        
        % Proportion of particles remained in the same cell = diagonal terms of Pij
        RatioRemain = diag(TransDen);
        
        % Proportion of particles Transit to neighbouring cells (5-points
        % stencil)
        Nx_cell = RectMesh_Param.Nx_cell;
        Ny_cell = RectMesh_Param.Ny_cell;
        cellij_to_k = RectMesh_Param.cellij_to_k;
        cellk_to_ij = RectMesh_Param.cellk_to_ij;
        
        RatioRemainNeigh = zeros(size(RatioRemain));
        for cell_k = 1:length(RatioRemainNeigh)
            % Brainless method: first convert back to 2D array and only care +-1
            % along i and j directions
            
            % Note: basically it is just the sum of some terms along the row of the
            % sparse
            TransDen_from_k = TransDen(cell_k, :);
            
            % Note:
            % TransDen_from_k_xy(cell_x, cell_y) store the prob from cell_k to (cell_x, cell_j)
            
            cell_xy = cellk_to_ij(cell_k);
            cell_x = cell_xy(1);
            cell_y = cell_xy(2);
            
            % Ensure not out of boundaries; 9-points stencil
            for cell_dx = [-1, 0, 1]
                for cell_dy = [-1, 0, 1]
                    cell_xP = cell_x + cell_dx;
                    cell_yP = cell_y + cell_dy;
                    
                    if ( (cell_xP <= Nx_cell) &&  (cell_yP <= Ny_cell) ...
                            && (cell_xP >= 1) &&  (cell_yP >= 1) )
                        cell_l = cellij_to_k(cell_xP, cell_yP);
                        RatioRemainNeigh(cell_k) = RatioRemainNeigh(cell_k) + TransDen_from_k(cell_l);
                    end
                end
            end

        end
end
