function [Kxx, Kyy, Kxy]  = Kpolar_to_Kcart_vectorised(K_sigma1, K_sigma2, K_phi)

% Tested; Compared with Maple
% Formula is correct
% Assumed K_sigma1 > K_sigma2; K_phi between -pi/2 to pi/2

c = cos(K_phi);
s = sin(K_phi);

Kxx = (c .*K_sigma1).^2 + (s .* K_sigma2).^2;
Kyy = (s .*K_sigma1).^2 + (c .* K_sigma2).^2;
Kxy = c.*s.*(K_sigma1.^2 - K_sigma2.^2);

end