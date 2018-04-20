function [sigmas, phi] = sort_sigmas_restrict_phi(sigmas, phi)
% Sort sigma1 and sigma2
if sigmas(1) < sigmas(2)
    phi = phi + 0.5*pi;
    
    sigmas(1:2) = sigmas(2:-1:1);
end

% Restrict to 0 to 2*pi
while (phi >= 2*pi) || (phi < 0)
    phi = phi - sign(phi)*2*pi;
end

% Restrict to 0 to pi
if (phi >= pi)
    phi = phi - pi;
end
end
