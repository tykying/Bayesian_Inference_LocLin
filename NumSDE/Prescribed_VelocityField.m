%% Prescription of Velocity Field

% Input required:
% veloc_testcase: specify which velocity field to use
% Pe: Peclet Number, to determine the scale of velocity field

% Output:
% u(x, y): velocity field
function veloc_fldStruct = Prescribed_VelocityField(Profile, param)
Pe = param.Pe;
kappa_scale = param.kappa_scale;
U_scale = param.U_scale;
GradU_scale = param.GradU_scale;


if contains(Profile, 'sincos')
    psi_scale = u_scale;
    psi =  @(x, y) -psi_scale*(1/2*pi)*sin(2*pi*x).*sin(2*pi*y);
    
    u1 = @(x, y) (u_scale)* sin(2*pi*x).*cos(2*pi*y);
    u2 = @(x, y) (v_scale)* cos(2*pi*x).*sin(2*pi*y);
    Gradu1 = @(x, y) (u_scale)* [ (+2*pi)* cos(2*pi*x).*cos(2*pi*y) , (-2*pi)* sin(2*pi*x).*sin(2*pi*y)];
    Gradu2 = @(x, y) (v_scale)* [ (-2*pi)* sin(2*pi*x).*cos(2*pi*y) , (+2*pi)* cos(2*pi*x).*cos(2*pi*y)];
    
    u = @(x, y) [u1(x,y); u2(x,y)];    % Assume x,y is row vector [Python-friendly]
    Gradu = @(x,y) [Gradu1(x,y); Gradu2(x,y)];
    
    Gradu_scale = (u_scale)*(+2*pi);
    Gradv_scale = (v_scale)*(+2*pi);
end

if contains(Profile, 'noisy_sincos')
    U_small_scale = 0.1*U_scale;  % Relative scale (to u_scale) of the small kappa perturbation (high frequency)
    
    %psi_scale = u_scale;
    %psi =  @(x, y) -(psi_scale)*( (1/2*pi)*sin(2*pi*x).*sin(2*pi*y) +  (u_small_scale/(24*pi))*sin(24*pi*x).*sin(24*pi*y) );
    
    u1 = @(x, y) (U_scale(1))*( sin(2*pi*x).*cos(2*pi*y) + U_small_scale(1)*sin(24*pi*x).*cos(24*pi*y) );
    u2 = @(x, y) (U_scale(2))*( cos(2*pi*x).*sin(2*pi*y) + U_small_scale(2)*cos(24*pi*x).*sin(24*pi*y) );
    
    u = @(x, y) [u1(x,y); u2(x,y)];    % Assume x,y is row vector [Python-friendly]
    Gradu = @(x,y) [1, 1; ...
        1, 01];  % Not yet impli
    
end


if contains(Profile, 'vortex')
    
    % psi_scale = u_scale;
    % psi =  @(x, y) -(psi_scale)*( (1/2*pi)*sin(2*pi*x).*sin(2*pi*y);
    
    u1 = @(x, y) (u_scale)*( -y./(x.^2+y.^2) );
    u2 = @(x, y) (v_scale)*(  x./(x.^2+y.^2) );
    Gradu1 = @(x, y) (u_scale)* [ (2.*y.*x)./(x.^2+y.^2).^2 , (y.^2-x.^2)./(x.^2+y.^2).^2 ];
    Gradu2 = @(x, y) (v_scale)* [ (y.^2-x.^2)./(x.^2+y.^2).^2 , (-2*x.*y)./(x.^2+y.^2).^2 ];
    
    u = @(x, y) [u1(x,y); u2(x,y)];    % Assume x,y is row vector [Python-friendly]
    Gradu = @(x,y) [Gradu1(x,y); Gradu2(x,y)];
    
    Gradu_scale = (u_scale)/L;
    Gradv_scale = (v_scale)/L;
    
end


if contains(Profile, 'pwc_shear')
    
    % psi_scale = u_scale;
    % psi =  @(x, y) -(psi_scale)*( (1/2*pi)*sin(2*pi*x).*sin(2*pi*y);
    
    u1 = @(x, y) (u_scale)*(1*(y>0.5) + -2*(y<=0.5) );
    u2 = @(x, y) (v_scale)*(  x*0 );
    Gradu1 = @(x, y) (u_scale)* [ 0 , 0 ];
    Gradu2 = @(x, y) (v_scale)* [ 0 , 0 ];
    
    u = @(x, y) [u1(x,y); u2(x,y)];    % Assume x,y is row vector [Python-friendly]
    Gradu = @(x,y) [Gradu1(x,y); Gradu2(x,y)];
    
    Gradu_scale = (u_scale)/L;
    Gradv_scale = (v_scale)/L;
    
end

if contains(Profile, 'taylor_green')
    
    % psi_scale = u_scale;
    % psi =  @(x, y) -(psi_scale)*( (1/2*pi)*sin(2*pi*x).*sin(2*pi*y);
    
    u1 = @(x, y) (u_scale)*(-cos(y).*sin(x));
    u2 = @(x, y) (v_scale)*( sin(y).*cos(x));
    Gradu1 = @(x, y) (u_scale)* [ -cos(y).*cos(x) , sin(y).*sin(x) ];
    Gradu2 = @(x, y) (v_scale)* [ -sin(y).*sin(x) , cos(y).*cos(x) ];
    
    u = @(x, y) [u1(x,y); u2(x,y)];    % Assume x,y is row vector [Python-friendly]
    Gradu = @(x,y) [Gradu1(x,y); Gradu2(x,y)];
    
    Gradu_scale = (u_scale)/L;
    Gradv_scale = (v_scale)/L;
    
end

if contains(Profile, 'mean_drift')
    
    % psi_scale = u_scale;
    % psi =  @(x, y) -(psi_scale)*( (1/2*pi)*sin(2*pi*x).*sin(2*pi*y);
    
    u1 = @(x, y) (u_scale)*(0.25)*ones(size(x));
    u2 = @(x, y) (v_scale)*(0)*ones(size(x));
    Gradu1 = @(x, y) (u_scale)* [ 0, 0 ];
    Gradu2 = @(x, y) (v_scale)* [ 0, 0 ];
    
    u = @(x, y) [u1(x,y); u2(x,y)];    % Assume x,y is row vector [Python-friendly]
    Gradu = @(x,y) [Gradu1(x,y); Gradu2(x,y)];
    
    Gradu_scale = (u_scale)/L;
    Gradv_scale = (v_scale)/L;
    
end

if contains(Profile, 'shear')
    % psi_scale = u_scale;
    % psi =  @(x, y) -(psi_scale)*( (1/2*pi)*sin(2*pi*x).*sin(2*pi*y);
    dudy = GradU_scale(1);
    
    u1 = @(x, y) (dudy*y);
    u2 = @(x, y) (0.*x);
    Gradu1 = @(x, y) [ 0 , dudy ];
    Gradu2 = @(x, y) [ 0 , 0 ];
    
    u = @(x, y) [u1(x,y); u2(x,y)];    % Assume x,y is row vector [Python-friendly]
    Gradu = @(x,y) [Gradu1(x,y); Gradu2(x,y)];
end

if contains(Profile, 'linear_shear')
    % psi_scale = u_scale;
    % psi =  @(x, y) -(psi_scale)*( (1/2*pi)*sin(2*pi*x).*sin(2*pi*y);
    dudx = 1; dudy = 1;
    dvdx = 3; dvdy = -dudx * GradU_scale(1)/GradU_scale(2);
    bx = 0;
    by = -2;
    
    % Append param struct
    param.A = [GradU_scale(1)*[dudx, dudy]; GradU_scale(2)*[dvdx dvdy]];
    param.b = [U_scale(1)*bx; U_scale(2)*by];
    
    A = param.A;
    b = param.b;
    
    u = @(x, y) A*[x; y] + b;    % Assume x,y is row vector [Python-friendly]
    Gradu = @(x,y) A;
end



if contains(Profile, 'rigid_vortex')
    % psi_scale = u_scale;
    % psi =  @(x, y) -(psi_scale)*( (1/2*pi)*sin(2*pi*x).*sin(2*pi*y);
    dudx = 0; dudy = 1;
    dvdx = -1; dvdy = -dudx * GradU_scale(1)/GradU_scale(2);
    bx = 0;
    by = 0;
    
    % Append param struct
    param.A = [GradU_scale(1)*[dudx, dudy]; GradU_scale(2)*[dvdx dvdy]];
    param.b = [U_scale(1)*bx; U_scale(2)*by];
    
    A = param.A;
    b = param.b;
    
    u = @(x, y) A*[x; y] + b;    % Assume x,y is row vector [Python-friendly]
    Gradu = @(x,y) A;
end


if contains(Profile, 'rigid_vortex_rescaled')
    % psi_scale = u_scale;
    % psi =  @(x, y) -(psi_scale)*( (1/2*pi)*sin(2*pi*x).*sin(2*pi*y);
    dudx = 0; dudy = 1;
    dvdx = -2; dvdy = -dudx * GradU_scale(1)/GradU_scale(2);
    bx = 0;
    by = 0;
    
    % Append param struct
    param.A = [GradU_scale(1)*[dudx, dudy]; GradU_scale(2)*[dvdx dvdy]];
    param.b = [U_scale(1)*bx; U_scale(2)*by];
    
    A = param.A;
    b = param.b;
    
    u = @(x, y) A*[x; y] + b;    % Assume x,y is row vector [Python-friendly]
    Gradu = @(x,y) A;
end


if contains(Profile, 'taylor_green_noisy')
    noise_scale = 0.1;
    % psi_scale = u_scale;
    % psi =  @(x, y) -(psi_scale)*( (1/2*pi)*sin(2*pi*x).*sin(2*pi*y);
    
    Lx = 3;
    Ly = 1;
    Lx_noise = 18;
    Ly_noise = 12;
    
    
    u1 = @(x, y) (U_scale(1))*(-cos(Ly*y).*sin(Lx*x));
    u2 = @(x, y) (U_scale(2))*( sin(Ly*y).*cos(Lx*x));
    
    Gradu1 = @(x, y) (U_scale(1))* [ -Lx*cos(Ly*y).*cos(Lx*x) , Ly*sin(Ly*y).*sin(Lx*x) ];
    Gradu2 = @(x, y) (U_scale(2))* [ -Lx*sin(Ly*y).*sin(Lx*x) , Ly*cos(Ly*y).*cos(Lx*x) ];
    
    
    u1_noise = @(x, y) noise_scale*u1(Lx_noise/Lx*x, Ly_noise/Ly*y);
    u2_noise = @(x, y) noise_scale*u2(Lx_noise/Lx*x, Ly_noise/Ly*y);
    
    % Recall: contrapositive and covariant basis
    Gradu1_noise = @(x, y) noise_scale*Gradu1(Lx_noise/Lx*x, Ly_noise/Ly*y) .* [Lx_noise/Lx, Ly_noise/Ly];
    Gradu2_noise = @(x, y) noise_scale*Gradu2(Lx_noise/Lx*x, Ly_noise/Ly*y) .* [Lx_noise/Lx, Ly_noise/Ly];
    
    % Append param struct
    param.Lx = Lx;
    param.Ly = Ly;
    param.Lx_noise = Lx_noise;
    param.Ly_noise = Ly_noise;
    
    u = @(x, y) [u1(x,y); u2(x,y)] + [u1_noise(x,y); u2_noise(x,y)];    % Assume x,y is row vector [Python-friendly]
    Gradu = @(x,y) [Gradu1(x,y); Gradu2(x,y)] + [Gradu1_noise(x,y); Gradu2_noise(x,y)];
end


if contains(Profile, 'childress_soward')
    % psi_scale = u_scale;
    % psi =  @(x, y) sin(Lxc*x)*sin(Lyc*y) + lambda*cos(Lxs*x)*cos(Lys*y)
    
    lambda = 3/4;
    Lxs = 4/3;
    Lys = 1;
    Lxc = 1/2;
    Lyc = 3/2;
    
    u1 = @(x, y) (U_scale(1))*(-Lys*sin(Lxs*x).*cos(Lys*y) + lambda* Lyc*cos(Lxc*x).*sin(Lyc*y));
    u2 = @(x, y) (U_scale(2))*( Lxs*cos(Lxs*x).*sin(Lys*y) - lambda* Lxc*sin(Lxc*x).*cos(Lyc*y));
    
    ux = @(x, y) Lxc*cos(Lxc*x)*Lyc.*cos(Lyc*y) - lambda* Lxs*sin(Lxs*x).*Lys*sin(Lys*y);
    uy = @(x, y)   sin(Lxc*x)*Lyc^2.*sin(Lyc*y) + lambda* cos(Lxs*x)*Lys^2.*cos(Lys*y);
    vx = @(x, y)  -Lxc^2*sin(Lxc*x).*sin(Lyc*y) - lambda* Lxs^2*cos(Lxs*x).*cos(Lys*y);
    vy = @(x, y) Lxc*cos(Lxc*x)*Lyc.*cos(Lyc*y) + lambda* Lxs*sin(Lxs*x).*Lys*sin(Lys*y);
    
    Gradu1 = @(x, y) (U_scale(1))* [ ux(x,y) , uy(x,y) ];
    Gradu2 = @(x, y) (U_scale(2))* [ vx(x,y) , vy(x,y) ];
    
    
    % Append param struct
    param.Lxs = Lxs;
    param.Lys = Lys;
    param.Lxc = Lxc;
    param.Lyc = Lyc;
    param.lambda = lambda;
    
    u = @(x, y) [u1(x,y); u2(x,y)];    % Assume x,y is row vector [Python-friendly]
    Gradu = @(x,y) [Gradu1(x,y); Gradu2(x,y)];
end

if contains(Profile, 'taylor_green_with_shear')
    % U_scale controls the mean flow
    % GradU_scale controls the shear strength
    % psi_scale controls the vortix strength
    % Lx, Ly controls the lengthscale of the vortices
    
    psi_scale = param.vortix_strength;  % Should be equal to 1 by default
    %psi_scale = 1;
    
    Lx = 1;
    Ly = 1;
    
    Au1 = GradU_scale(1) * [0, 1];
    Au2 = GradU_scale(2) * [0, 0];
    assert(Au1(1) + Au2(2) == 0)  % Incompressibility
    
    b1 = U_scale(1) * [0];
    b2 = U_scale(2) * [0];
    
    psi = @(x, y) psi_scale*sin(Lx*x).*sin(Ly*y) + 0.5*(-Au1(2)*y.^2 + Au2(1)*x.^2) -2*Au1(1)*x.*y - b1*y + b2*x;
    
    
    u1 = @(x, y) (psi_scale*Ly)*(-cos(Ly*y).*sin(Lx*x)) + Au1(1).*x + Au1(2).*y + b1;
    u2 = @(x, y) (psi_scale*Lx)*( sin(Ly*y).*cos(Lx*x)) + Au2(1).*x + Au2(2).*y + b2;
    
    Gradu1 = @(x, y) (psi_scale*Ly)* [ -Lx*cos(Ly*y).*cos(Lx*x) , Ly*sin(Ly*y).*sin(Lx*x) ] + Au1;
    Gradu2 = @(x, y) (psi_scale*Lx)* [ -Lx*sin(Ly*y).*sin(Lx*x) , Ly*cos(Ly*y).*cos(Lx*x) ] + Au2;
    
    % Append param struct
    param.Lx = Lx;
    param.Ly = Ly;
    param.Au1 = Au1;
    param.Au2 = Au2;
    param.b1 = b1;
    param.b2 = b2;
    
    u = @(x, y) [u1(x,y); u2(x,y)];    % Assume x,y is row vector [Python-friendly]
    Gradu = @(x,y) [Gradu1(x,y); Gradu2(x,y)];
    
end



if contains(Profile, 'tg_w_mean')
    % Special case of taylor_green_with_shear
    % Guarantee homogenisation to work
    psi_scale = param.vortix_strength;
    
    Lx = 1;
    Ly = 1;
    
    b1 = U_scale(1);
    b2 = U_scale(2);
    
    psi = @(x, y) psi_scale*sin(Lx*x).*sin(Ly*y) - b1*y + b2*x;
    
    u1 = @(x, y) (psi_scale*Ly)*(-cos(Ly*y).*sin(Lx*x)) + b1;
    u2 = @(x, y) (psi_scale*Lx)*( sin(Ly*y).*cos(Lx*x)) + b2;
    
    Gradu1 = @(x, y) (psi_scale*Ly)* [ -Lx*cos(Ly*y).*cos(Lx*x) , Ly*sin(Ly*y).*sin(Lx*x) ];
    Gradu2 = @(x, y) (psi_scale*Lx)* [ -Lx*sin(Ly*y).*sin(Lx*x) , Ly*cos(Ly*y).*cos(Lx*x) ];
    
    % Append param struct
    param.Lx = Lx;
    param.Ly = Ly;
    param.b1 = b1;
    param.b2 = b2;
    
    u = @(x, y) [u1(x,y); u2(x,y)];    % Assume x,y is row vector [Python-friendly]
    Gradu = @(x,y) [Gradu1(x,y); Gradu2(x,y)];
    
end



if contains(Profile, 'cosine')
    % Special case of taylor_green_with_shear
    % Guarantee homogenisation to work
    %psi_scale = param.vortix_strength;
    
    Lx = 0;
    Ly = pi;
    
    assert(U_scale(2) == 0);
    
    psi_scale = U_scale(1)/Ly;
    
    psi = @(x, y) -psi_scale.*sin(Ly*y);
    
    u1 = @(x, y) (psi_scale*Ly).*cos(Ly*y);
    u2 = @(x, y) psi_scale*0.*y;
    
    Gradu1 = @(x, y) (psi_scale*Ly)* [ 0 , -(psi_scale*Ly*Ly).*sin(Ly*y) ];
    Gradu2 = @(x, y) (psi_scale*Lx)* [ 0 , 0 ];
    
    % Append param struct
    param.Lx = Lx;
    param.Ly = Ly;
    
    u = @(x, y) [u1(x,y); u2(x,y)];    % Assume x,y is row vector [Python-friendly]
    Gradu = @(x,y) [Gradu1(x,y); Gradu2(x,y)];
    
end


if contains(Profile, 'linear_jet_CplxEigen')
    if contains(Profile, '_ISO')
        A_base = [5, 6; -2 -1];
    else
        A_base = [0.6, 6; -2 -0.2];
    end
    A = GradU_scale' .* A_base;
    
    b = param.b_mag .* [cos(param.b_angle); sin(param.b_angle)];
    
    % Append param struct
    param.A = A;
    param.b = b;
    
    u = @(x, y) A*[x; y] + b;    % Assume x,y is row vector [Python-friendly]
    Gradu = @(x,y) A;
end

if contains(Profile, 'linear_jet_RealEigen')
    if contains(Profile, '_ISO')
        A_base = [5, 6; 2 -1];
    else
        A_base = [0.6, 6; 2 -0.2];
    end
    A = GradU_scale' .* A_base;
    
    b = param.b_mag .* [cos(param.b_angle); sin(param.b_angle)];
    
    % Append param struct
    param.A = A;
    param.b = b;
    
    u = @(x, y) A*[x; y] + b;    % Assume x,y is row vector [Python-friendly]
    Gradu = @(x,y) A;
end

if contains(Profile, 'linear_jet_CplxEigen_InC')
    if contains(Profile, '_ISO')
        A_base = [3, 6; -2 -3];
    else
        A_base = [0.8, 6; -2 -0.8];
    end
    A = GradU_scale' .* A_base;
    
    b = param.b_mag .* [cos(param.b_angle); sin(param.b_angle)];
    
    % Append param struct
    param.A = A;
    param.b = b;
    
    u = @(x, y) A*[x; y] + b;    % Assume x,y is row vector [Python-friendly]
    Gradu = @(x,y) A;
end

if contains(Profile, 'linear_jet_RealEigen_InC')
    if contains(Profile, '_ISO')
        A_base = [3, 6; 2 -3];
    else
        A_base = [0.8, 6; 2 -0.8];
    end
    A = GradU_scale' .* A_base;
    
    b = param.b_mag .* [cos(param.b_angle); sin(param.b_angle)];
    
    % Append param struct
    param.A = A;
    param.b = b;
    
    u = @(x, y) A*[x; y] + b;    % Assume x,y is row vector [Python-friendly]
    Gradu = @(x,y) A;
end

% Wrap up for output
veloc_fldStruct = struct('u', u, 'Gradu', Gradu, 'Profile', Profile, 'param', param);


end