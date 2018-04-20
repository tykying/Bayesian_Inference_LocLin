function [x0, y0, nparticles] = initialise_particle_positions(npart_pc_DIR, ncell_DIR, PART_ranges_min, PART_ranges_max)
%% Initialise_particle_positions% Specifications:

% Input:
% npart_pc_DIR: vector of size 2; number of particles per cell along each direction
% ncell_DIR: vector of size 2; number of particles per cell along each direction
% ranges_min, ranges_max: Ranges of the rectangular box in which the particles are placed

% Output:
% x0: 2D array, size(1, nparticles); Initial x cooridnates of all particles
% y0: 2D array, size(1, nparticles); Initial y cooridnates of all particles


%% computation
field_dim = 2;

cellsize_DIR = zeros(1, field_dim);
partsep_DIR = zeros(1, field_dim);

onehalf = 0.5;

% Calculate nparticles
nparticles = 1;
for i = 1: field_dim
    nparticles = nparticles * (npart_pc_DIR(i) * ncell_DIR(i));
    cellsize_DIR(i) = (PART_ranges_max(i) - PART_ranges_min(i))/(ncell_DIR(i));
    partsep_DIR(i) = cellsize_DIR(i)/(npart_pc_DIR(i));
end

x0 = zeros(1, nparticles);
y0 = zeros(1, nparticles);

nparticles_counter = 0;
% For readibility, hardcode for dim=2
% Lopp over all cells
for j_cell = 1: ncell_DIR(2)
    for i_cell = 1: ncell_DIR(1)
        % Within each cell
        for j_part_pc = 1: npart_pc_DIR(2)
            for i_part_pc = 1: npart_pc_DIR(1)
                nparticles_counter = nparticles_counter + 1;
                
                x0(nparticles_counter) = (i_cell-1)*cellsize_DIR(1) + (i_part_pc-onehalf)*partsep_DIR(1)+PART_ranges_min(1);
                y0(nparticles_counter) = (j_cell-1)*cellsize_DIR(2) + (j_part_pc-onehalf)*partsep_DIR(2)+PART_ranges_min(2);
            end
        end
    end
end

% Consistency check
assert(nparticles_counter == nparticles);

end