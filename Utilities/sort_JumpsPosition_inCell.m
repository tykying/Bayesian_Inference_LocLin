function TrajJumps_inCell = sort_JumpsPosition_inCell(x, y, diffx, diffy, GridInd, difft, cell_k)
        GridInd_inCell = (GridInd == cell_k);
                
        %   %%% Only when spatial structure is involved
        x_inCell = x(GridInd_inCell);
        y_inCell = y(GridInd_inCell);
        
        u_inCell = diffx(GridInd_inCell);
        v_inCell = diffy(GridInd_inCell);
        
        difft_inCell = difft(GridInd_inCell);
                
        TrajJumps_inCell = struct('x', x_inCell, 'y', y_inCell, ...
        	'diffx', u_inCell, 'diffy', v_inCell, 'difft', difft_inCell, ...
            'GridInd', GridInd_inCell);
end