function [Fields, GradxFields, GradyFields] = convert_theta_to_Fields(theta)

if size(theta,3) == 0
    Fields = zeros(size(theta,1), 5);   % 5 = U,V, K11,K22,K12
    GradxFields = zeros(size(theta,1), 5);   % 5 = U,V, K11,K22,K12
    GradyFields = zeros(size(theta,1), 5);   % 5 = U,V, K11,K22,K12   
        
    b = [theta(:, 1).* cos(theta(:, 2)), theta(:, 1).* sin(theta(:, 2))];
    nu_1 = theta(:, 3) + theta(:, 4);
    nu_2 = theta(:, 3) - theta(:, 4);

    A11 = -theta(:, 3).*sin(2*theta(:, 5));
    A12 = theta(:, 3).*2.*cos(theta(:, 5)).*cos(theta(:, 5)) - nu_2;
    A21 = theta(:, 3).*2.*cos(theta(:, 5)).*cos(theta(:, 5)) - nu_1;

    % trA_h = 0.5 * trace(A)
    trA_h = theta(:, 6);
    
    [K11, K22, K12] = Kpolar_to_Kcart_vectorised(theta(:,7), theta(:,8), theta(:,9));
    
    Fields(:,:) = [b, [K11, K22, K12]];
    
    GradxFields(:, 1:2) = [trA_h+A11, A21];
    GradyFields(:, 1:2) = [A12, trA_h-A11];
else
    Fields = zeros(size(theta,1), 5, size(theta,3));   % 5 = U,V, K11,K22,K12
    GradxFields = zeros(size(theta,1), 5, size(theta,3));   % 5 = U,V, K11,K22,K12
    GradyFields = zeros(size(theta,1), 5, size(theta,3));   % 5 = U,V, K11,K22,K12
        
    
    b = [theta(:, 1, :).* cos(theta(:, 2, :)), theta(:, 1, :).* sin(theta(:, 2, :))];
    nu_1 = theta(:, 3, :) + theta(:, 4, :);
    nu_2 = theta(:, 3, :) - theta(:, 4, :);
    
    A11 = -theta(:, 3, :).*sin(2*theta(:, 5, :));
    A12 = theta(:, 3, :).*2.*cos(theta(:, 5, :)).*cos(theta(:, 5, :)) - nu_2;
    A21 = theta(:, 3, :).*2.*cos(theta(:, 5, :)).*cos(theta(:, 5, :)) - nu_1;

    trA_h = theta(:, 6, :);
    
    [K11, K22, K12] = Kpolar_to_Kcart_vectorised(theta(:,7,:), theta(:,8,:), theta(:,9,:));

    for k = 1:size(theta,3)
        Fields(:,:,k) = [b(:,:,k), [K11(:,:,k), K22(:,:,k), K12(:,:,k)]];
        
        GradxFields(:, 1:2, k) = [trA_h(:,:,k) + A11(:,:,k), A21(:,:,k)];
        GradyFields(:, 1:2, k) = [A12(:,:,k), trA_h(:,:,k) - A11(:,:,k)];
    end
end



end