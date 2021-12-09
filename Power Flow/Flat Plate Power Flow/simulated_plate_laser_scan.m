%force_plate_simulated_laser_output


%% Cylinder Parameters
E = 200e9; %Pa
nu = 0.297;%0.29;  %Poisson's ratio
rho = 7860; %kg/m3
eta = .001;
thickness = 0.01;  %Plate thickness
width = 0.7; % x
height = 1.0; % y
D = E*thickness^3/(12*(1-nu^2));
% Mesh Parameters
xinc = width/20; %x-axis spatial increment
yinc = height/20; %Theta-axis spatial increment in degrees
x = 0:xinc:height;  %x-axis grid spacing
y = 0:yinc:width;%-yinc;  %y-axis grid spacing
Nmodesx = 20;  %Number of x modes in the modal summation
Nmodesy = 20;  %Number of y modes in the modal summation
mx = 1:1:Nmodesx; %x mode counting array
ny = 1:1:Nmodesy; %y mode counting array
% Force Parameters
xStar1 = .25*height; yStar1 = -60*pi/180;
xStar2 = .75*height; yStar2 = -240*pi/180;
xf1c = pi*xStar1/height;
xf2c = pi*xStar2/height;
F1 = 10;  %Magnitude of the first force in N.
F2 = 10;  %Magnitude of the second force in N.
ForceFreq = 150;  %Forcing frequency in Hz
OmegaF1 = 2*pi*ForceFreq; %Forcing frequency in rad/sec
OmegaF2 = OmegaF1;  %Second forcing frequency equals the first
% Time Parameters
%Put the time and phase part in the solution.
tinc = 1/ForceFreq/20;%0.0002; %time increment
tmax = 1/ForceFreq;%1.25e-4;%0.005; %max time. (make sure these are set to give a full cylce based on the forcing freq.
time = 0:tinc:tmax; %time array
dfreq = 1/tmax; %frequency axis spacing
Faxis = dfreq*(0:1:length(time)); %Frequency axis
SpikeN = round(ForceFreq/dfreq); %Determining the line number associated with the forcing frequency
% Viscous Damping
Lam = rho*thickness*eta/OmegaF1; %Damping ########################
% Initialize sizes for speed
W1 = zeros(length(x),length(y),length(time));
W2 = zeros(size(W1));
WtotalF = zeros(length(x),length(y));
% SaveData = zeros(length(x)*length(y),4);
WtotalD = zeros(length(x), length(y), length(time));
%% Natural Frequencies and Phase
Omega = pi^2*((mx'/width).^2+(ny/height).^2)*sqrt(D/(rho*thickness));
Phi1 = atan(2*Lam*(OmegaF1./Omega)./(1-(OmegaF1./Omega).^2));
Phi2 = Phi1-170*pi/180;


%% Equation to find W(x,y)

for j = 1:length(x)
    for k = 1:length(y)
        for t = 1:length(time)
            for M = 1:Nmodesx
                for N = 0:Nmodesy-1
                    
                    WtotalD(j,k,t) = W1(j,k,t) + W2(j,k,t);
                end
            end
        end
    end
end
%% Compute the FFT and determine the spectral real and imaginary terms
for j = 1:length(x)
    for k = 1:length(theta)
        WF = fft(WtotalD(j,k,:));
        WtotalF(j,k) = WF(SpikeN);
    end
end
%% Convert to rectangular coordinates
% xtrans = reshape(repmat(cos(theta), length(x),1), length(theta)*length(x),1);
% ytrans = reshape(repmat(sin(theta), length(x),1), length(theta)*length(x),1);
% xcoord = radius*xtrans;
% ycoord = radius*ytrans;
% thetaCoor = repmat(theta, length(x), 1);
% zcoord = repmat(x,length(theta),1)';
% myX = reshape(xcoord,numel(xcoord), 1);
% myY = reshape(ycoord,numel(ycoord), 1);
% myZ = reshape(zcoord,numel(zcoord), 1);
% radialDisplacement = reshape(WtotalF, numel(WtotalF), 1);
% dispX = radialDisplacement.*xtrans;
% rX = real(dispX);
% iX = imag(dispX);
% dispY = radialDisplacement.*ytrans;
% rY = real(dispY);
% iY = imag(dispY);
% rZ = 0*rY;
% iZ = 0*iY;
%% Save Data
dataFormat = "%d\t%d\t%d\t%d\t%d\t%d\t%d\n";
dataFile = "C:\Users\ME\Documents\Simulated Data\";
mkdir(strcat(dataFile, fileName, "_data"));
old = cd(strcat(dataFile, fileName, "_data"));
for point = 1:length(xcoord)
    num = sprintf("%03d", point);
    savePointFileName = strcat(fileName, "_Simulated_", num, ".txt");
    fid = fopen(savePointFileName,'wt');
    %Build the data file
    fprintf(fid, strcat("Source File Name:\t", fileName, "\n"));
    fprintf(fid, strcat("Point Index:\t", num2str(point), "\t[x = %d\ty = %d\tz = %d]\n"), myX(point), myY(point), myZ(point));
    fprintf(fid, "Component:\tRoot\n");
    fprintf(fid, "Signal: FFT - Vib 3D Displacement - Real & Imag.\n\n");
    fprintf(fid, "Frequency\tReal X\tReal Y\tReal Z\tImaginary X\t Imaginary Y\tImaginary Z\n");
    fprintf(fid, "[ Hz ]\t[ m ]\t[ m ]\t[ m ]\t[ m ]\t[ m ]\t[ m ]\n");
    fprintf(fid, dataFormat, ForceFreq, rX(point), rY(point), rZ(point), iX(point), iY(point), iZ(point));
    fclose(fid);
end
cd(old)
