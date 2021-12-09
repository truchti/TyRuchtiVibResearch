
clear variables
close all
%% Set up 
%plate Parameters
width = .6;% width (x)
height = 1.0;% height (y)
xres = .01;
yres = .01;
thk = .002;
%material Parameters
E = 2.00e11;
dens = 7850.0;
poi = .3;
eta = .001;
% force parameters
forceMag = 10;
loc1 = [.15 .85];
loc2 = [.45 .15];
fWidth = .002;
frequency = 50;
modes = 40;
%%
%set up mesh for plate
x = 0:xres:width; 
y = 0:yres:height;
[X,Y] = meshgrid(x,y);
D = E*thk^3/(12*(1-poi^2));
%freq
omega = frequency*2*pi;
% force magnitude
F = forceMag;
F2 = -forceMag;
%location of force 1
x1 = loc1(1);
y1 = loc1(2);
x2 = loc2(1);
y2 = loc2(2);
% x1 = loc1(1) - fWidth/2;% bound of x load
% x2 = loc1(1) + fWidth/2;
% y1 = loc1(2) - fWidth/2; % bound of y load
% y2 = loc1(2) + fWidth/2;
% %bounds of force 2
% x3 = loc2(1) - fWidth/2;
% x4 = loc2(1) + fWidth/2;
% y3 = loc2(2) - fWidth/2;
% y4 = loc2(2) + fWidth/2;
%%
% place holder variables initialized
w1 = zeros(length(y), length(x));
% dw1dt = zeros(length(y), length(x));
w2 = zeros(length(y), length(x));
% dw2dt = zeros(length(y), length(x));
%%
w1temp = 0;
w2temp = 0;
dwdt1temp = 0;
dwdt2temp = 0;
C1 = 4*F/(dens*thk*width*height);
C2 = 4*F2/(dens*thk*width*height);
for i = 1:length(x) %for all x coordinates
    for j = 1:length(y) %for all y coordinates
        for m = 1:modes %summation for m 1 to 40
            for n = 1:modes %summation for n 1 to 40
                %calculate omega mn 
                omn = pi^2*((m/width)^2+(n/height)^2)*sqrt(D/(dens*thk)); 
                %calculate w (Z displacement) and dwdt (z velocity)
                %calculate the constant part 
                aplitude1 =C1*...
                    sin(m*pi*x1/width)*...(cos(m*pi*x1/width)-cos(m*pi*x2/width))...                
                    sin(n*pi*y1/height)/...(cos(n*pi*y1/height)-cos(n*pi*y2/height))...
                    ((omn^2-omega^2)^2+eta^2*omn^4);
                amplitude2 = C2*...
                    sin(m*pi*x2/width)*...(cos(m*pi*x3/width)-cos(m*pi*x4/width))...
                    sin(n*pi*y2/height)/...(cos(n*pi*y3/height)-cos(n*pi*y4/height))...
                    ((omn^2-omega^2)^2+eta^2*omn^4);
                %calc plate shape for each mode
                disp =sin(m*pi*x(i)/width)* sin(n*pi*y(j)/height)*(omn^2-omega^2-1i*eta*omn^2);
                w1temp = disp*aplitude1+ w1temp;
                w2temp = disp*amplitude2 + w2temp;
            end
        end
        w1(j,i) = w1temp;
        w1temp = 0;
        w2(j,i) = w2temp;
        w2temp = 0;
    end
end
dw1dt = w1*1i*omega;
dw2dt = w2*1i*omega;
%%
w = w1+w2;
dwdt = dw1dt+dw2dt;

%% Quintic spline interpolation
%Creates additional points based on quintic splines to smooth surface
%interpolate across one dimension
% this creates splines across the surface by connecting all the y
% coordianes (x for each spline is constant)
for i =1:length(x)
    %creates a quintic spline structure for every row 
    c(i) = spapi(5,Y(:,i),w(:,i)); %#ok<*SAGROW>
    cv(i) =spapi(5,Y(:,i),dwdt(:,i));
    %interpolate additional points from spline structures (increases number
    %of columns)
    %Ztemp(:,i)= fnval(y,c(i));
end
%%
% surf(Xtemp,Ytemp,Ztemp)
% Creates spline objects for the surface where y is constant
for i =1:length(y)
    %create splines for each column of the expanded data set
    d(i) = spapi(5,X(i,:),w(i,:));
    dv(i) = spapi(5,X(i,:),dwdt(i,:));
    %interpolate additional points (increases rows)
    %Ztemp2(i,:) = fnval(x,d(i));
end
%clear X Y Z Xtemp Ytemp Ztemp i k n c 
%% Der function
%dw/dx
[dwdx, splinesdx] = surSplinePD2(d,x,1);
%d2wdx2
[dwddx, splinesddx] = surSplinePD2(splinesdx,x,1);
%d3wdx3
[dwdddx, splinesdddx] = surSplinePD2(splinesddx,x,1);
%dw/dy    
[dwdy, splinesdy] = surSplinePD2(c,y,2);
%d2wdy2
[dwddy,splinesddy] = surSplinePD2(splinesdy,y,2);
%d3wdy3
[dwdddy,splinesdddy] = surSplinePD2(splinesddy,y,2);

%% Cross Derive position
dwdydx = surSplineCrossder(d,X,Y);
%% 3rd order crossDerivs
for i = 1:length(y)
    xsplineddy(i) = spapi(5, X(i,:), dwddy(i,:));
end
for i = 1:length(x)
    ysplineddx(i) = spapi(5, Y(:,i), dwddx(:,i));
end
dwddydxtest = surSplinePD2(xsplineddy, x,1);
dwddxdytest = surSplinePD2(ysplineddx, y,2);
dwddydx = surSplineCrossder(splinesdx,X,Y);
dwddxdy = transpose(surSplineCrossder(splinesdy,Y',X'));
%% spacial derive velocity
[dwdxdt, splinesdxdt] = surSplinePD2(dv,x,1);
[dwdydt, splinesdydt] = surSplinePD2(cv,y,2);
Mx = D * (dwddx + poi*dwddy);
My = D * (dwddy + poi*dwddx);
Qx = D * (dwdddx+dwddydxtest);
Qy = D * (dwdddy+dwddxdytest);
Mxy = D * (1-poi) *dwdydx;


%%
qx = .5*D*real((dwdddx+dwddydx).*conj(dwdt)+ ...
     (dwddx+poi.*(dwddy)).*conj(dwdxdt)+...
     (1-poi).*dwdydx.*conj(dwdydt));
qy = .5*D*real((dwdddy+dwddxdy).*conj(dwdt)+...
    (dwddy+poi.*(dwddx)).*conj(dwdydt)+...
    (1-poi).*dwdydx.*conj(dwdydt));
% %%
 quiver(X(1:end,(1:end)),Y(1:end,(1:end)),qx(1:end,(1:end)),qy(1:end,(1:end)),'LineWidth',.75,'MaxHeadSize',1.5,'Color','k')
% axis equal
% axis([0 width 0 height])
% st = ' Hz';
% f= round(100.*frequency)/100;
% title([num2str(f) st])
% savefig([num2str(f) st '.fig'])
%% resultant plots




%  surf(X,Y,real(Mxy));
%     view(0,90)
%     axis equal
%     colorbar
