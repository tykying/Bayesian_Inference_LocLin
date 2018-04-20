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
function kappa_fldStruct = Prescribed_DiffusivityField(kappa_scale, kappa_testcase)
BaseTensor = [1, 0; 0, 1];

if strcmp(kappa_testcase, 'sincos')   % Variable diffusion
    kappa  = @(x, y) kappa_scale*(1 -0.9*sin(1*pi*x).*cos(1*pi*y));
    Dxkappa = @(x, y) kappa_scale*( -0.9*1*pi*cos(1*pi*x).*cos(1*pi*y));
    Dykappa = @(x, y) kappa_scale*( -0.9*1*pi*sin(1*pi*x).*(-sin(1*pi*y)));
    Gradkappa = @(x, y) [Dxkappa(x,y); Dykappa(x,y)]';
    
    Gradkappa_scale = kappa_scale*pi;
end

if strcmp(kappa_testcase, 'noisy_sincos')   % Variable diffusion
    kappa_small_scale = 0.1;  % Relative scale (to kappa) of the small kappa perturbation (high frequency)
    
    kappa  = @(x, y) kappa_scale*(1 -0.9*sin(1*pi*x).*cos(1*pi*y) + kappa_small_scale*cos(2*pi*x).*sin(2*pi*y));
    Dxkappa = @(x, y) kappa_scale*( -0.9*1*pi*cos(1*pi*x).*cos(1*pi*y) + kappa_small_scale*2*pi*(-sin(2*pi*x)).*sin(2*pi*y) );
    Dykappa = @(x, y) kappa_scale*( -0.9*1*pi*sin(1*pi*x).*(-sin(1*pi*y)) + kappa_small_scale*2*pi*cos(2*pi*x).*cos(2*pi*y) );
    Gradkappa = @(x, y) [Dxkappa(x,y); Dykappa(x,y)]';
    
    Gradkappa_scale = kappa_scale*pi;
end


if strcmp(kappa_testcase, 'const')   % Variable diffusion    
    kappa  = @(x, y) kappa_scale;
    Dxkappa = @(x, y) 0;
    Dykappa = @(x, y) 0;
    Gradkappa = @(x, y) [Dxkappa(x,y); Dykappa(x,y)]';
    
    Gradkappa_scale = kappa_scale;
end

if strcmp(kappa_testcase, 'pwc_shear')
    kappa  = @(x, y) kappa_scale*(2*(x>0.5) + 1*(x<=0.5) );
    Dxkappa = @(x, y) 0;
    Dykappa = @(x, y) 0;
    Gradkappa = @(x, y) [Dxkappa(x,y); Dykappa(x,y)]';
    
    Gradkappa_scale = kappa_scale;
    
    BaseTensor = [1, 1.5; 1.5, 3];
end


sqrt2kappa = @(x, y) sqrt(2*kappa(x, y));   % Assume kappa: equal, diagonal matrix


kappa_fldStruct = struct('kappa', kappa, 'Gradkappa', Gradkappa, 'kappa_scale', kappa_scale, 'Gradkappa_scale', Gradkappa_scale, 'BaseTensor', BaseTensor);

end