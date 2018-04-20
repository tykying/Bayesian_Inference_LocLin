function theta_store  = unzip_theta(theta_zipped)
    [Ncells, NVars, NSamples_zipped] = size(theta_zipped);
        
    Nsteps_pc = (NSamples_zipped-1)*NVars;
       
    theta_store = zeros(Ncells, NVars, Nsteps_pc+1);
    
    % step_z = 1: initial conditions
    for step_z = 1:NSamples_zipped
        for cmp = 1:NVars
            % Deal with initial condition: step_z-1 and step+1
            steps_list = convert_step_z_to_step_list(cmp, NVars, (step_z-1));
            steps_list = steps_list + 1;
            
            steps_list = steps_list((steps_list>=1) & (steps_list<=(Nsteps_pc+1)));
            theta_store(:, cmp, steps_list) = theta_zipped(:,cmp, step_z);
        end
    end
end