function PlotEllipse(ax, ellipse_Param)

%ellipse_Param = struct('a', 5, 'b', 3, 'phi', pi/4, 'centre', [1;0])

a = ellipse_Param.a;
b = ellipse_Param.b;
phi = ellipse_Param.phi;
centre = ellipse_Param.centre;
PlotSpec = ellipse_Param.PlotSpec;

assert(a >= b)

e1 = [cos(phi); sin(phi)];
e2 = [-sin(phi); cos(phi)];

t = linspace(-pi, pi, 120);

alpha = cos(t);
beta = sin(t);

X = alpha.*e1.*a + beta.*e2.*b;

u = X(1,:) + centre(1);
v = X(2,:) + centre(2);
plot(ax, u, v, 'Color', PlotSpec.Color, 'LineStyle', PlotSpec.LineStyle);
