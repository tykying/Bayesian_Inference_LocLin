function [K_sigma1, K_sigma2, K_phi]  = Kcart_to_Kpolar(Kxx, Kyy, Kxy)
K = [Kxx, Kxy; Kxy, Kyy];
[V, D] = eig(K);
eigval = diag(D);
MajorAxeInd = find(eigval == max(eigval), 1);
MinorAxeInd = 1*(MajorAxeInd==2) + 2*(MajorAxeInd==1);

sigma1_CellMean = sqrt(eigval(MajorAxeInd));
sigma2_CellMean = sqrt(eigval(MinorAxeInd));

MajorAxeEigvec = V(:, MajorAxeInd);

phi_CellMean = atan2(MajorAxeEigvec(2), MajorAxeEigvec(1));  % angle for eigvec_1

if abs(phi_CellMean) > pi*0.5
    phi_CellMean = atan2(-MajorAxeEigvec(2), -MajorAxeEigvec(1));
end

K_sigma1 = sigma1_CellMean;
K_sigma2 = sigma2_CellMean;
K_phi = phi_CellMean;
end