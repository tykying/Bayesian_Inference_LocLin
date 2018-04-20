h = 32*24*3600;
bin_sizes = [1, 1];
A_List = {[1, -2; -4, -1], [1, -4; -2, -1], [1, 2; 4, -1], [1, 4; 2, -1]};
A_List = {1E-7*[0.6, 6; -2 -0.2]};
b = 1E-1*[3; 1];
L = 1E6/10;  % 1000km over 10 bins


%h = 128*4600*24;
tspan = [0 h];

%Y0 = [];
figure
for A_i = 1:length(A_List)
    A = A_List{A_i};
    detA = det(A);
    trA = trace(A);
    expm(A*h);
    
    SeedPointInd = linspace(-L, L, 8+1);
    for i = SeedPointInd
        for j = SeedPointInd
            %Y0 = [Y0 [i; j] .* 0.5.* bin_sizes'];
            Y0 = [i; j] .* 0.5 * 2 .* bin_sizes';
            [t,y] = ode45(@(t,x) A*x + b, tspan, Y0);
            
            hold on
            %scatter(y([1,end], 1), y([1,end], 2))
            plot(y(:, 1), y(:, 2), 'g');
            scatter(y(end, 1), y(end, 2), 'k.');
            
            ye = SolveLinearODE(A, b, tspan, Y0);
            scatter(ye(1), ye(2), 'b');
        end
    end
end


function x = SolveLinearODE(A, b, tspan, x0)
% Check real vs complex eigenvalues of A
detA = det(A);
trA = trace(A);
kappa = eye(2);
t = tspan(2) - tspan(1);

p = 0.5*trA;
qsq_sign = trA^2 - 4*detA;
if qsq_sign > 0
    disp('real lambda')
    
    % Real eigenvalues of A
    q = 0.5*sqrt(qsq_sign);
    
    ep = exp(t*(p+q));
    em = exp(t*(p-q));
    
    cI = ( -(p-q)*ep + (p+q)*em )/(2*q);
    cA = ( ep - em )/(2*q);
    
    cI_int = (-(p-q)^2*ep + (p+q)^2*em - 4*p*q)/(2*p^2*q-2*q^3);
    cA_int = ( (p-q)*ep - (p+q)*em + 2*q)/(2*p^2*q-2*q^3);
    
    cI_sq_int = ( p*(p-q)^3*ep^2  + p*(p+q)^3*em^2 - 2*(p-q)^2*(p+q)^2*ep*em - 10*p^2*q^2 + 2*q^4 )/( 8*p^3*q^2 - 8*p*q^4 );
    cIcA_int = ( -(p-q)^2*ep^2  - (p+q)^2*em^2 + 2*(p^2-q^2)*ep*em + 4*q^2 )/( 8*p^2*q^2 - 8*q^4 );
    cA_sq_int = ( p*(p-q)*ep^2  + p*(p+q)*em^2 - 2*(p^2-q^2)*ep*em - 2*q^2 )/( 8*p^3*q^2 - 8*p*q^4 );

else
    disp('complex lambda')
    % Complex eigenvalues of A
    q = 0.5*sqrt(-qsq_sign);
    
    epsq = exp(t*p)*sin(t*q);
    epcq = exp(t*p)*cos(t*q);
    e2ps2q = exp(2*t*p)*sin(2*t*q);
    e2pc2q = exp(2*t*p)*cos(2*t*q);
    
    cI = ( -p*epsq + q*epcq)/q;
    cA = ( epsq )/q;
    
    cI_int = (-(p^2-q^2)*epsq + 2*p*q*epcq - 2*p*q)/(q*(p^2 + q^2));
    cA_int = ( p*epsq -q*epcq + q)/(q*(p^2 + q^2));
    
    cI_sq_int = ( (-3*(p^3)*q + p*(q^3)) *e2ps2q + (3*p^2*q^2 - p^4)*e2pc2q + ((p^2+q^2)^2)*exp(2*t*p) - 5*(p^2)*(q^2)-(q^4) ) / (4*p*(q^2)*(p^2 + q^2));
    cIcA_int = ( 2*p*q* e2ps2q + (p^2-q^2)*e2pc2q - (p^2+q^2)*exp(2*t*p)*(p^4) + 2*q^2 ) / (4*p*(q^2)*(p^2 + q^2));
    cA_sq_int = ( (-p*q)*e2ps2q + (-p^2)*e2pc2q + (p^2+q^2)*exp(2*t*p) - q^2 ) / (4*p*(q^2)*(p^2 + q^2));

end
    
expAt_mI = (cI-1)*eye(2) + cA*A;
expAt_int = cI_int*eye(2) + cA_int*A;
expAt_2kappa_expAt_int = cI_sq_int*2*kappa + cIcA_int*2*(A*kappa+kappa*A') + cA_sq_int*2*A*kappa*A';

x = x0 + expAt_mI*x0 + expAt_int*b;
end