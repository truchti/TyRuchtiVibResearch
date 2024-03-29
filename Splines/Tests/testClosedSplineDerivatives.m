% Test quintic splines on a closed loop

% set up points around a circle
theta = linspace(0, 2*pi, 181)';
r = .43;
x = r*cos(theta);
y = r*sin(theta);
% set up function evalutation at points around circle
cc = 2;
w = cos(cc*theta);
dw = -cc*sin(cc*theta);
d2w = -cc^2*cos(cc*theta);
d3w = cc^3*sin(cc*theta);
% create single spline to match function 
parameters = theta; data = w;
%
a = quinticBSplineFitter(parameters, data, 'closed', 15);
a.fit_spline
b = a.output_solved_spline_evaluator();
ppts = linspace(0,2*pi, 181);

fxs = b.evaluate_at_parameters(ppts);
dfxs = b.evaluate_at_parameters(ppts, 1);
d2fxs = b.evaluate_at_parameters(ppts, 2);
d3fxs = b.evaluate_at_parameters(ppts, 3);
figure(1)
subplot(2,2,1)
plot(theta, w, ':k', ppts, fxs, 'r')
title('disp')
subplot(2,2,2)
plot(theta, dw, ':k', ppts, dfxs, 'r')
title('1st')
subplot(2,2,3)
plot(theta, d2w, ':k', ppts, d2fxs, 'r')
title('2nd')
subplot(2,2,4)
plot(theta, d3w, ':k', ppts, d3fxs, 'r')
title('3rd')

%% effects of adding minor noise
ndata = data + .01*randn(size(data));
a = quinticBSplineFitter(parameters, ndata, 'closed', 15);
a.fit_spline
b = a.output_solved_spline_evaluator();
ppts = linspace(0,2*pi, 181);

fxs = b.evaluate_at_parameters(ppts);
dfxs = b.evaluate_at_parameters(ppts, 1);
d2fxs = b.evaluate_at_parameters(ppts, 2);
d3fxs = b.evaluate_at_parameters(ppts, 3);
figure(2)
subplot(2,2,1)
plot(theta, w, ':k', ppts, fxs, 'r')
title('disp')
subplot(2,2,2)
plot(theta, dw, ':k', ppts, dfxs, 'r')
title('1st')
subplot(2,2,3)
plot(theta, d2w, ':k', ppts, d2fxs, 'r')
title('2nd')
subplot(2,2,4)
plot(theta, d3w, ':k', ppts, d3fxs, 'r')
title('3rd')


