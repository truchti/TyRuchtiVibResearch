%% test for Chebyshev Fit of 1D sinusoidal equation and derivatives
x = linspace(0, 1.2, 86);
f = @(x) sin(2*pi*x/1.2).*x.^2;
df = @(x) 2*x.*sin((5*pi*x)/3) + (5*pi*x.^2.*cos((5*pi*x)/3))/3;
ddf = @(x) 2*sin((5*pi*x)/3) - (25*pi^2*x.^2.*sin((5*pi*x)/3))/9 + (20*pi*x.*cos((5*pi*x)/3))/3;
dddf = @(x) 10*pi*cos((5*pi*x)/3) - (50*pi^2*x.*sin((5*pi*x)/3))/3 - (125*pi^3*x.^2.*cos((5*pi*x)/3))/27;

fx = f(x);
dfx = df(x);
d2fx = ddf(x);
d3fx = dddf(x);
scatter(x,fx)

%%


b = ChebyshevPolyFit(12,x,fx);
b.find_coefficients
% b.plot_data_and_fit
fitC = b.fitPoly;
scatter(x,d3fx)
hold on;
fitC.plot_deriv(x,3)
hold off;