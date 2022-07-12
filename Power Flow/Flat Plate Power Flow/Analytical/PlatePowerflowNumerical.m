% Quintic Spline interpolation
clear variables
close all
%% Set up 
%plate dimensions
a = 1;% width (x)
b = .7;% height (y)
x1 = .49;% bound of x load
x2 = .51;
y1 = .34; % bound of y load
y2 = .36;
x = 0:.005:a; 
y = 0:.005:b;
[X,Y] = meshgrid(x,y);
thk = .002;
%properties and constants
dens = 7869.0/(a*b); %mass per area
F = 5; %force into plate
E = 2.040e11; %modulus of material
poi = .3; %poissons ration
eta = .1;%structural damping
D = E*thk^3/(12*(1-poi^2));
im = sqrt(-1);
%%
frequency =90;% in Hz
omega = frequency*2*pi; %frequency in rad/sec
% temporary variables for summations
wtemp = 0;
dwdttemp = 0;
w = zeros(length(x), length(y));
dwdt = zeros(length(x), length(y));
Q = 4*F/(dens*thk*pi^2);
for i = 1:length(x) %for all x coordinates
    for j = 1:length(y) %for all y coordinates
        for m = 1:25 % modes in x direction
            for n = 1:25 %modes in y direction
                %calculate omega mn 
                omn = pi^2*((m/a)^2+(n/b)^2)*sqrt(D/(dens*thk));
                OMN(m,n) = omn;
                denom = ((omn^2-omega^2)^2+eta^2*omn^4);
                DENOM(m,n) = denom;
                %calculate w (Z displacement) and dwdt (z velocity)
                %calculate the constant part 
                R = 1/(m*n)*(cos(m*pi*x1/a)-cos(m*pi*x2/a))*(cos(n*pi*y1/b)-cos(n*pi*y2/b))*sin(m*pi*x(i)/a)*sin(n*pi*y(j)/b);
                %calculate the iterative portion of the positon and
                %velocity
                wtemp =    R*(          omn^2-omega^2-im*eta*omn^2) /denom + wtemp;
                dwdttemp = R*(im*omega*(omn^2-omega^2-im*eta*omn^2))/denom + dwdttemp;
            end
        end
        w(i,j) = Q*wtemp;
        wtemp = 0;
        dwdt(i,j) = dwdttemp;
        dwdttemp = 0;
    end
end
%%
w = w';
dwdt = dwdt';
wr = real(w);
dwrdt = real(dwdt);
surf(X,Y,real(w))
% savefig('Plate_Wxy')
surf(X,Y,real(dwdt))
% savefig('Plate_Vel')
%% Quintic spline interpolation
%Creates additional points based on quintic splines to smooth surface
%interpolate across one dimension
for i =1:length(x)
    %creates a quintic spline structure for every row 
    c(i) = spapi(5,Y(:,i),w(:,i));
    cv(i) =spapi(5,Y(:,i),dwdt(:,i));
    %interpolate additional points from spline structures (increases number
    %of columns)
    Ztemp(:,i)= fnval(y,c(i));
end
%%
% surf(Xtemp,Ytemp,Ztemp)
for i =1:length(y)
    %create splines for each column of the expanded data set
    d(i) = spapi(5,X(i,:),w(i,:));
    dv(i) = spapi(5,X(i,:),dwdt(i,:));
    %interpolate additional points (increases rows)
    Ztemp2(i,:) = fnval(x,d(i));
end
%clear X Y Z Xtemp Ytemp Ztemp i k n c 
%% spatial derivatives
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

%% Spacial cross derivatives
dwdydx = surSplineCrossder(d,X,Y);
dwddydx = surSplineCrossder(splinesdx,X,Y);
dwddxdy = transpose(surSplineCrossder(splinesdy,Y',X'));
%% spacial  time derivatives
[dwdxdt, splinesdxdt] = surSplinePD2(dv,x,1);
[dwdydt, splinesdydt] = surSplinePD2(cv,y,2);
%%
qx = .5*D*real((dwdddx+dwddydx).*conj(dwdt)-(dwddx+poi.*(dwddy)).*conj(dwdxdt)-(1-poi).*dwdydx.*conj(dwdydt));
qy = .5*D*real((dwdddy+dwddydx).*conj(dwdt)-(dwddy+poi.*(dwddx)).*conj(dwdydt)-(1-poi).*dwdydx.*conj(dwdydt));
%%Plot powerflow
quiver(X(1:3:end,(1:3:end)),Y(1:3:end,(1:3:end)),qx(1:3:end,(1:3:end)),qy(1:3:end,(1:3:end)),'LineWidth',.75,'MaxHeadSize',1.5,'Color','k')
axis equal
axis([0 1 0 .7])

Qmag  = (qx.^2+qy.^2).^(1/2);

% surf(X,Y,qx)
% savefig('Plate_qx')
% surf(X,Y,qy)
% savefig('Plate_qy')
% surf(X,Y,Qmag)
% savefig('Plate_Q')








