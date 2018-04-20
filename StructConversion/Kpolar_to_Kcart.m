function [Kxx, Kyy, Kxy, R]  = Kpolar_to_Kcart(K_sigma1, K_sigma2, K_phi)

R = [[cos(K_phi); sin(K_phi)], [-sin(K_phi); cos(K_phi)]];
K = R*diag([K_sigma1^2, K_sigma2^2])*R';

Kxx = K(1,1);
Kyy = K(2,2);
Kxy = K(1,2);

end