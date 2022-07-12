x = linspace(0, 2, 50);
y = linspace(0, .7, 60);
[X,Y] = ndgrid(x,y);
val = X.^2.*sin(pi*X).*sin(2*pi*Y/.7)-Y.^3.*sin(4*pi*X).*sin(5*pi*Y/.7);
xData = reshape(X, numel(X), 1);
yData = reshape(Y, numel(Y), 1);
values = reshape(val, numel(val), 1);
% surf(X,Y,val)
%% analytical equations
dx = @(X,Y) 2*X.*sin(pi*X).*sin((20*pi*Y)/7) + X.^2.*pi.*cos(pi*X).*sin((20*pi*Y)/7) - 4*Y.^3.*pi.*cos(4*pi*X).*sin((50*pi*Y)/7);
dy = @(X,Y) (20*X.^2.*pi.*cos((20*pi*Y)/7).*sin(pi*X))/7 - 3*Y.^2.*sin(4*pi*X).*sin((50*pi*Y)/7) - (50*Y.^3.*pi.*cos((50*pi*Y)/7).*sin(4*pi*X))/7;
dx3 = @(X,Y) 6.*pi.*cos(pi.*X).*sin((20.*pi.*Y)/7) - 6.*X.*pi^2.*sin(pi.*X).*sin((20.*pi.*Y)/7) - X.^2.*pi^3.*cos(pi.*X).*sin((20.*pi.*Y)/7) + 64.*Y.^3.*pi^3.*cos(4.*pi.*X).*sin((50.*pi.*Y)/7);
dy3 = @(X,Y) (125000*Y.^3.*pi.^3.*cos((50*pi*Y)/7).*sin(4*pi*X))/343 - (8000.*X.^2.*pi.^3.*cos((20*pi*Y)/7).*sin(pi*X))/343 - 6.*sin(4*pi*X).*sin((50*pi*Y)/7) + (22500.*Y.^2.*pi.^2.*sin(4*pi*X).*sin((50*pi*Y)/7))/49 - (900.*Y.*pi.*cos((50*pi*Y)/7).*sin(4*pi*X))/7;
[X,Y] = ndgrid(x,y);
dxv = dx3(X,Y);
dyv = dy3(X,Y);
%% new method
cfit = ChebyshevSurfaceFitter([23,24],xData, yData, values);
csurf = cfit.fitSurface;
%% quintic splines
sfit = quinticBSplineSurfaceFitter([xData, yData], values, {'open', 'open'}, [9, 9]);
sfit.fit_spline_surfaces();
sSurf = sfit.output_solved_spline_evaluator();

%%
figure(1)
subplot(1,3,1)
csurf.plot_deriv(x,y, [3,0])
title('Chebyshev - New Method')
subplot(1,3,2)
surf(X,Y,dxv)
title('Analytical')
subplot(1,3,3)
sSurf.plot_derivative(x,y,3,0)
title('Quintic Splines')
%%
figure(2)
subplot(1,3,1)
csurf.plot_deriv(x,y, [0,3])
title('Chebyshev - New Method')
subplot(1,3,2)
surf(X,Y,dyv)
title('Analytical')
subplot(1,3,3)
sSurf.plot_derivative(x,y,0,3)
title('Quintic Splines')
%%
figure(3)
subplot(1,3,1)
csurf.plot(x,y)
title('Chebyshev - New Method')
subplot(1,3,2)
surf(X,Y,val)
title('Analytical')
subplot(1,3,3)
sSurf.plot_spline(x,y)
title('Quintic Splines')