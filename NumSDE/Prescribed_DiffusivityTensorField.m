%% Prescription of (Analytic) Diffusivity field

% Remark: Assumed domain = [0, 1]

% Input required:
% kappa_testcase: specify which kappa field to use
% kappa_scale: scale of kappa field

% Output:
% kappa(x, y): scalar diffusivity field
% Gradkappa(x, y): gradient of diffusivity field
% sqrt2kappa: sqrt(2*kappa(x))

% Remark: Assumed domain = [0, 1]

% Accept x,y to be row vectors listing the positions of all particles

%% Computation
function kappa_fldStruct = Prescribed_DiffusivityTensorField(Profile, param)
kappa_scale = param.kappa_scale;


if strcmp(Profile, 'const')   % Constant diffusion   
    BaseTensor = [1, 0; 0, 1];

    kappa  = @(x, y) kappa_scale.* BaseTensor;
    Divkappa_r1 = @(x, y) 0.*x;
    Divkappa_r2 = @(x, y) 0.*y;
    Divkappa = @(x, y) [Divkappa_r1(x,y); Divkappa_r2(x,y)];
    
    Divkappa_scale = kappa_scale;
    param.kappa = kappa_scale.*BaseTensor;
end

if strcmp(Profile, 'const_skew')
    BaseTensor = [3, 1; 1, 2];
    
    kappa  = @(x, y) kappa_scale.* BaseTensor;
    Divkappa_r1 = @(x, y) 0.*x;
    Divkappa_r2 = @(x, y) 0.*y;
    Divkappa = @(x, y) [Divkappa_r1(x,y); Divkappa_r2(x,y)];
    
    Divkappa_scale = kappa_scale;
    param.kappa = kappa_scale.*BaseTensor;
end

if strcmp(Profile, 'planar_sinusoidal')
    Lx1 = 1; Ly1 = 1;
    Lx2 = 1/8; Ly2 = 1/8;
    Lxp = -1; Lyp = 1;

    sigma1 = @(x, y) (1/2)*(sin(Lx1*x)*cos(Ly1*y)+3/2);
    sigma2 = @(x, y) (1/2)*(cos((Lx2+Ly2)*(x+y))+1)*sigma1(x, y);
    phi = @(x, y) (pi/2)*sin(Lxp*x+Lyp*y);
    
    % Formula from Kpolar_to_Kcart_vectorised
    Kxx = @(x, y) (cos(phi(x,y)) .*sigma1(x,y)).^2 + (sin(phi(x,y)) .* sigma2(x,y)).^2;
    Kyy = @(x, y) (sin(phi(x,y)) .*sigma1(x,y)).^2 + (cos(phi(x,y)) .* sigma2(x,y)).^2;
    Kxy = @(x, y)  cos(phi(x,y)) .* sin(phi(x,y)) .* (sigma1(x,y).^2 - sigma2(x,y).^2);

    BaseTensor = @(x, y) [Kxx(x,y), Kxy(x,y); Kxy(x,y), Kyy(x,y)];

    kappa  = @(x, y) kappa_scale.* BaseTensor(x, y);
    Divkappa_r1 = @(x, y) 0.*x;  % Not yet implement
    Divkappa_r2 = @(x, y) 0.*y;  % Not yet implement
    Divkappa = @(x, y) [Divkappa_r1(x,y); Divkappa_r2(x,y)];  % Not yet implement
    
    Divkappa_scale = kappa_scale;
    param.Lx1 = Lx1;
    param.Ly1 = Ly1;
    param.Lx2 = Lx2;
    param.Ly2 = Ly2;
    param.Lxp = Lxp;
    param.Lyp = Lyp;
end


if strcmp(Profile, 'sinusoidal')
    Lx1 = 1; Ly1 = 1;
    Lx2 = 1/8; Ly2 = 1/8;
    Lxp = -1; Lyp = 1;

    sigma1 = @(x, y) (1/2)*(sin(Lx1*x)*cos(Ly1*y)+3/2);
    sigma2 = @(x, y) (1/2)*(cos((Lx2+Ly2)*(x+y))+1)*sigma1(x, y);
    phi = @(x, y) (pi/2)*sin(Lxp*x+Lyp*y);
    
    % Formula from Kpolar_to_Kcart_vectorised
    Kxx = @(x, y) (cos(phi(x,y)) .*sigma1(x,y)).^2 + (sin(phi(x,y)) .* sigma2(x,y)).^2;
    Kyy = @(x, y) (sin(phi(x,y)) .*sigma1(x,y)).^2 + (cos(phi(x,y)) .* sigma2(x,y)).^2;
    Kxy = @(x, y)  cos(phi(x,y)) .* sin(phi(x,y)) .* (sigma1(x,y).^2 - sigma2(x,y).^2);

    BaseTensor = @(x, y) [Kxx(x,y), Kxy(x,y); Kxy(x,y), Kyy(x,y)];

    kappa  = @(x, y) kappa_scale.* BaseTensor(x, y);
    Divkappa_r1 = @(x, y) 0.*x;  % Not yet implement
    Divkappa_r2 = @(x, y) 0.*y;  % Not yet implement
    Divkappa = @(x, y) [Divkappa_r1(x,y); Divkappa_r2(x,y)];  % Not yet implement
    
    Divkappa_scale = kappa_scale;
    param.Lx1 = Lx1;
    param.Ly1 = Ly1;
    param.Lx2 = Lx2;
    param.Ly2 = Ly2;
    param.Lxp = Lxp;
    param.Lyp = Lyp;
end

if strcmp(Profile, 'cg_taylor_green_with_mean_drift')   % from experiment        
    BaseTensor = [0.4346, -0.0051; -0.0051, 0.2947];
    
    kappa  = @(x, y) BaseTensor;
    Divkappa_r1 = @(x, y) 0.*x;
    Divkappa_r2 = @(x, y) 0.*y;
    Divkappa = @(x, y) [Divkappa_r1(x,y); Divkappa_r2(x,y)];
    
    Divkappa_scale = kappa_scale;
end

if strcmp(Profile, 'pwc_shear')    
    BaseTensor1 = [1, 1.5; 1.5, 3];
    BaseTensor2 = [5, -3; -3, 2];

    kappa  = @(x, y) kappa_scale*((x>0.5).* BaseTensor2 + (x<=0.5).* BaseTensor1);
    Divkappa_r1 = @(x, y) 0.*x;
    Divkappa_r2 = @(x, y) 0.*y;
    Divkappa = @(x, y) [Divkappa_r1(x,y); Divkappa_r2(x,y)];
    
    Divkappa_scale = kappa_scale;
end


if strcmp(Profile, 'quadratic_x')    
    BaseTensor1 = [1, 1.5; 1.5, 3];
    
    kappa = @(x, y) kappa_scale*(BaseTensor1*x*x);
    Divkappa_r1 = @(x, y) 2.*x;
    Divkappa_r2 = @(x, y) 3.*x;  %(1.5*2x)
    Divkappa = @(x, y) [Divkappa_r1(x,y); Divkappa_r2(x,y)];
    
    Divkappa_scale = kappa_scale;
end

kappa_fldStruct = struct('kappa', kappa, 'Divkappa', Divkappa,  'Profile', Profile, 'param', param);
end