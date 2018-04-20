function steps_list = convert_step_z_to_step_list(cmp, NVars, step_z)
    assert(cmp <= NVars);
    assert(cmp >= 1);

    steps_init = cmp + (step_z-1)*NVars;
    steps_list = steps_init + [0:(NVars-1)];
end