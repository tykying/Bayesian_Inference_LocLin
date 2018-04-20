function logLik_cell = compute_logLik_jump_Fields_cell_linearSDE(TrajJumps_cell, Fields_cell, GradxFields_cell, GradyFields_cell)
%Jtheta: contains u(Jump), kappa(Jump). Interpolated by the right scheme of inference.
logLik_jump = zeros(size(TrajJumps_cell.diffx));

% ASSUME UNIFORM TIME STEP HERE
Jh = TrajJumps_cell.h(1);


% A = [Gradu; Gradv]; b = [u0, v0]
%disp('Fields:')
A = [[GradxFields_cell(1), GradyFields_cell(1)]; [GradxFields_cell(2), GradyFields_cell(2)]];
b = [Fields_cell(1); Fields_cell(2)];

kappa = [Fields_cell(3), Fields_cell(5); Fields_cell(5), Fields_cell(4)];

% disp('Fields:')
% A
% b
% kappa
% det(A)

%%
% Check real vs complex eigenvalues of A
detA = det(A);
trA = trace(A);

p = 0.5*trA;
det_qsq_d = trA^2 - 4*detA;
q = 0.5*sqrt(abs(det_qsq_d));
pq = p*q;
ppq = p+q;
pmq = p-q;
% ppq_sq = p^2 + q^2;
% pmq_sq = p^2 - q^2;
ppq_sq = 0.25*(trA^2 + abs(det_qsq_d));
pmq_sq = 0.25*(trA^2 - abs(det_qsq_d));
q_sq = 0.25*abs(det_qsq_d);

if all(A(:)==0)
    % PWC Velocity Fields
    % disp('zero A')
    mut_X0 = b.*Jh;
    
    JSigma11 = 2*kappa(1,1)*Jh;  % Correct this part
    JSigma22 = 2*kappa(2,2)*Jh;
    JSigma12 = 2*kappa(1,2)*Jh;
    
    detSigma = JSigma11*JSigma22 - JSigma12*JSigma12;
    invdetSigma = 1/detSigma;
    
    Misfit_x = TrajJumps_cell.diffx' - mut_X0(1,:);
    Misfit_y = TrajJumps_cell.diffy' - mut_X0(2,:);
else
    % PWL Velocity Fields
    if abs(trA) < 1E-15
        % Case 1: Incompressible
        if abs(detA) > 1E-20
            % Case 1a: Non-vanishing detA
            rtdetA = sqrt(abs(detA));
            
            %invdetA = 1/detA;
            %invrtdetA = 1/rtdetA;
            
            if detA > 0
                c = cos(rtdetA*Jh); s = sin(rtdetA*Jh);
                %disp('positive det(A)')
            else
                c = cosh(rtdetA*Jh); s = sinh(rtdetA*Jh);
                %disp('negative det(A)')
            end
            
            cA_sq_2int = (Jh-s*c/rtdetA)/detA;
            cI_sq_2int = (Jh+s*c/rtdetA);
            cA_cI_2int = s*s/abs(detA);
            
            expAt_mI = (c-1)*eye(2) + (s/rtdetA)*A;
            expAt_int = (s/rtdetA)*eye(2) + ((1-c)/detA)*A;
            expAt_2kappa_expAt_int = cI_sq_2int*kappa + cA_cI_2int*(A*kappa+kappa*A') + cA_sq_2int*A*kappa*A';
        else
        % Case 1b: Incompressible
            expAt_mI = Jh*A;
            expAt_int = Jh*eye(2) + 0.5*Jh*Jh*A;
            expAt_2kappa_expAt_int = 2*Jh*kappa + Jh*Jh*(A*kappa+kappa*A') + (2/3)*Jh^3*A*kappa*A';
        end
    else
        % Case 2: Compressible
        if det_qsq_d > 0
            %disp('real lambda')
            
            % Real eigenvalues of A
            if abs(detA) < eps
                %disp('One eigenvalue = 0');
                
                % One of the eigenvalues is zero
                e2p = exp(Jh*2*p);
                
                cI = 1;
                cA = (e2p - 1)/(2*p);
                
                cI_int = Jh;
                cA_int = (-2*Jh*e2p-1)/(4*p^2);
                
                cI_sq_int = Jh;
                cIcA_int = (-2*Jh*p + e2p - 1)/(4*p^2);
                cA_sq_int = (exp(4*Jh*p) + 4*Jh*p - 4*e2p + 3)/(16*p^3);
            else
                %disp('Distinct Real');
                eppq = exp(Jh*(ppq));
                epmq = exp(Jh*(pmq));
                e2ppq = exp(Jh*2*(ppq));
                e2pmq = exp(Jh*2*(pmq));
                e2p = exp(Jh*2*p);
                
                cI = ( -pmq*eppq + ppq*epmq )/(2*q);
                cA = ( eppq - epmq )/(2*q);
                
                dem_fac = (pmq_sq)*2*q;
                
                cI_int = (-pmq^2*eppq + ppq^2*epmq - 4*pq)/dem_fac;
                cA_int = ( pmq*eppq - ppq*epmq + 2*q)/dem_fac;
                
                cI_sq_int = ( p*(pmq^3*e2ppq + ppq^3*e2pmq) - 2*pmq_sq^2*e2p - 10*pq^2 + 2*q_sq^2 )/( 4*pq* dem_fac);
                cIcA_int = ( -pmq^2*e2ppq - ppq^2*e2pmq + 2*pmq_sq*e2p + 4*q_sq )/( 4*q* dem_fac );
                cA_sq_int = ( p*pmq*e2ppq + p*ppq*e2pmq - 2*pmq_sq*e2p - 2*q_sq )/( 4*pq* dem_fac );
            end
            
        else
            %disp('Complex lambda')
            % Complex eigenvalues of A
            epsq = exp(Jh*p).*sin(Jh*q);
            epcq = exp(Jh*p).*cos(Jh*q);
            
            e2p = exp(2*Jh*p);
            e2ps2q = e2p.*sin(2*Jh*q);
            e2pc2q = e2p.*cos(2*Jh*q);
            
            cI = ( -p*epsq + q*epcq)/q;
            cA = epsq/q;
            
            dem_fac = q*ppq_sq;
            
            cI_int = (-pmq_sq*epsq - 2*pq*(1-epcq))/(dem_fac);
            cA_int = ( p*epsq + q*(1-epcq))/(dem_fac);
            
            cI_sq_int = ( -p^2*(pmq_sq-2*q^2)*e2pc2q - pq*(2*p^2+pmq_sq)*e2ps2q + ppq_sq^2*e2p -5*pq^2 -q_sq^2)/( 4*pq*dem_fac);
            cIcA_int = ( pmq_sq*e2pc2q + 2*pq*e2ps2q - ppq_sq*e2p + 2*q_sq )/( 4*q*dem_fac);
            cA_sq_int = (-p^2*e2pc2q - pq*e2ps2q + ppq_sq*e2p - q_sq)/( 4*pq*dem_fac);
        end
        
        % vec_new = [cI, cA, cI_int, cA_int, cI_sq_int, cIcA_int, cA_sq_int]
        
        expAt_mI = (cI-1)*eye(2) + cA*A;
        expAt_int = cI_int*eye(2) + cA_int*A;
        expAt_2kappa_expAt_int = cI_sq_int*2*kappa + cIcA_int*2*(A*kappa+kappa*A') + cA_sq_int*2*A*kappa*A';
     end
        
     JSigma11 = expAt_2kappa_expAt_int(1,1);
     JSigma22 = expAt_2kappa_expAt_int(2,2);
     JSigma12 = expAt_2kappa_expAt_int(1,2);
     
     detSigma = det(expAt_2kappa_expAt_int);
     invdetSigma = 1/detSigma;
     
     
     X0 = [TrajJumps_cell.x0'; TrajJumps_cell.y0'];
     mut_X0 = expAt_mI*X0 + expAt_int*b;
     
     Misfit_x = TrajJumps_cell.diffx' - mut_X0(1,:);
     Misfit_y = TrajJumps_cell.diffy' - mut_X0(2,:);
end

% -0.5 = Prefactor of -1/2 in Normal distribution
ExpoCost = (-0.5).*(JSigma22.*Misfit_x.^2 - 2*JSigma12.*Misfit_x.*Misfit_y + JSigma11.*Misfit_y.^2) * invdetSigma;

logLik_jump = ExpoCost - 0.5.*log(detSigma);

logLik_cell = sum(logLik_jump);

if isreal(logLik_cell) == 0
    disp('Check Matrix: non-real logLik_cell!');
    A
    kappa
    expAt_2kappa_expAt_int
    detSigma
    detA
    trA
    [cI, cA, cI_int, cA_int, cI_sq_int, cIcA_int, cA_sq_int]
    
    % Set to exp(-inf) = 0
    logLik_cell = -inf;
end


end
