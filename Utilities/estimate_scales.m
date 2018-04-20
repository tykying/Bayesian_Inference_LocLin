function Est_Scales = estimate_scales(TrajData_SubSampled)
x = TrajData_SubSampled.x_SubSampled;
y = TrajData_SubSampled.y_SubSampled;
ts_list = TrajData_SubSampled.ts_list_SubSampled;


% 0. Crude preliminary estimation of scales of velocity and diffusivity
% u_scale, v_scale: from r.m.s estimate [Seem to overestimate actual scale]
u_scale = mean( sqrt(mean((diff(x).^2),2))./diff(ts_list) );
v_scale = mean( sqrt(mean((diff(y).^2),2))./diff(ts_list) );

% ReS: Rescale the data: minus the row mean, then divide by sqrt(2*h)
% Idea: Normalisation
% For Dx = u dt + N(0, 2kh)
% At each time step, approximate the drift by ensemble mean, i.e. u dt = mean(Dx)
% Then divide by sqrt(2h)
% Finally compute the variance of the rescaled/normalised data
% [Seem to slightly underestimate the actual prescribed scale]
dx_ReS = (diff(x, 1)-mean(diff(x), 2))./(2*sqrt(diff(ts_list)));
dy_ReS = (diff(y, 1)-mean(diff(y), 2))./(2*sqrt(diff(ts_list)));

kappa_scale = var(dx_ReS(:)) + var(dy_ReS(:));

L = ceil(max([abs(x(:)); abs(y(:))]));

Est_Scales = struct('u_scale', u_scale, 'v_scale', v_scale, 'kappa_scale', kappa_scale, 'L', L);
end