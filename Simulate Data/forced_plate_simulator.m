classdef forced_plate_simulator < handle
    properties
        material
        force1
        force2
        geometry
        mesh
    end
    properties(Dependent)
        fileName
        naturalFrequencies
    end
    properties (Hidden = true, Dependent = true)
        heights
        widths
        D
        lambda
        time
        SpikeN
    end
    properties (Hidden = true)
        timesteps = 20
        nOmega= [];
        complexDampingPhase = [];
        WtotalD
        WtotalF
        data
        saveFolder
        dWdx
        dWdy
        dWdx2
        dWdxy
        dWdy2
        dWdx3
        dWdx2y
        dWdxy2
        dWdy3
        Oxdot
        Oydot
    end
    methods
        function show_full_analytical_power_flow(obj)
            [Mx, My, Mxy, Qx, Qy] = obj.calculate_moments_and_shears();
            [wdotS, OxdotS, OydotS] = obj.conjugate_velocities();
            qx = 1/2*real(Qx.*wdotS + Mx.*OydotS - Mxy.*OxdotS);
            qy = 1/2*real(Qy.*wdotS - My.*OxdotS + Mxy.*OydotS);
            [W, H] =meshgrid(obj.widths, obj.heights);
            quiver(W,H, qx, qy);
            
        end
        function show_velocity(obj, type)
            [wd, Odx, Ody] = obj.conjugate_velocities();
            switch type
                case 'wdot'
                    value = wd;
                case 'Oxdot'
                    value = Odx;
                case 'Oydot'
                    value = Ody;
            end
            surf(obj.widths, obj.heights, zeros*real(value), real(value))
            view(0,90)
            colorbar
            colormap jet
        end
        function show_resultant(obj, type)
            [Mx, My, Mxy, Qx, Qy] = obj.calculate_moments_and_shears();
            switch type
                case 'Mx'
                    value = Mx;
                case 'My'
                    value = My;
                case 'Mxy' 
                    value = Mxy;
                case 'Qx'
                    value = Qx;
                case 'Qy'
                    value = Qy;
            end
            surf(obj.widths, obj.heights, zeros*real(value), real(value))
            view(0,90)
            colorbar
            colormap jet
            
        end
        function obj = forced_plate_simulator(geometryOrFreq, forces, material, mesh)
            if nargin < 1
                obj.create_default();
            elseif nargin < 2
                obj.create_default_at_frequency(geometryOrFreq);
            else
                obj.geometry = geometryOrFreq;
                obj.force1 = forces(1);
                obj.force2 = forces(2);
                obj.material = material;
                obj.mesh = mesh;
            end
        end
        function animate_displacement_in_time(obj, imagData, loops)
            if nargin < 3 
                loops = 3;
            end
            if nargin < 2 
                imagData = false;
            end
            [H, W] = meshgrid(obj.heights, obj.widths);
            if imagData
                values = imag(obj.WtotalD);
            else
                values = real(obj.WtotalD);
            end            
            
            figure(1)
            s = surf(W,H,values(:,:,1),values(:,:,1));
            s.EdgeAlpha = 0.0;
            s.FaceColor = 'interp';
            colorbar
            clim = max(max(max(values)));
            zlim([-clim, clim]);
            caxis([-clim, clim])
            sz = size(values,3);
            for i = 1:loops*sz
                s.ZData = values(:,:,mod(i,sz)+1);
                s.CData = values(:,:,mod(i,sz)+1);
                pause(0.1);
            end
        end 
        function simulate_and_save(obj)
            obj.simulate_data();
            obj.save_data();
        end
        function simulate_data(obj)
            obj.calculate_mode_frequencies_and_phases();
            obj.calculate_displacement();
            obj.calculate_spectral_terms();
            obj.collectData();
        end
        function save_data(obj)
            if ~isempty(obj.data)
                obj.write_velocity_file();
                obj.write_displacement_file();
            else
                warning("data has not been simulated")
            end
        end
        function disable_force(obj, force_number)
            switch force_number
                case 1
                    obj.force1.enabled = false;
                case 2
                    obj.force2.enabled = false;
            end
        end
        function enable_force(obj, force_number)
            switch force_number
                case 1
                    obj.force1.enabled = true;
                case 2
                    obj.force2.enabled = true;
            end
        end
        function props = export_plate_power_flow_properties(obj)
            props = Plate_Properties(obj.geometry.width, obj.geometry.height, obj.geometry.thickness, obj.material.E, obj.material.poisson);
        end
    end
    methods (Access = private)
        function calculate_displacement(obj)
            if obj.force2.enabled
                obj.calculate_2_force_displacement;
            else
                obj.calculate_1_force_displacement;
            end
        end
        function calculate_1_force_displacement(obj)
            %precalculate values that will be used a lot
            W1 = zeros(length(obj.widths), length(obj.heights), length(obj.time));
            Phi1 = obj.complexDampingPhase + obj.force1.phase;
            timeV = obj.time;
            F1W = obj.force1.omega;
            Widths = obj.widths;
            Height = obj.heights;
            pix1_a = pi*obj.force1.location(1)/obj.geometry.width;
            piy1_b = pi*obj.force1.location(2)/obj.geometry.height;
            Cf1 = 4*obj.force1.magnitude/(pi^2*obj.material.density*obj.geometry.thickness);
            % loops over M modes along width
            for M = 1:obj.mesh.width_modes
                sinX1 = sin(M*pix1_a);
                sinWidths = (sin(M*pi*Widths/obj.geometry.width))';
                Cf1_m = Cf1/M;
                % loops over N modes along height
                for N = 1:obj.mesh.height_modes
                    Cf1_mn = Cf1_m/N;
                    sinY1 = sin(N*piy1_b);
                    sinHeights = sin(N*pi*Height/obj.geometry.height);
                    denom = obj.calculate_denominator(M,N);
                    complexDampCoeff = obj.calculate_complex_damping_coefficient(M,N);
                    allScalars1 = Cf1_mn*sinX1*sinY1/denom;
                    gridEvals = (sinWidths*sinHeights)*complexDampCoeff;
                    amp1 = allScalars1*gridEvals;  % because sinMxj is column and cosNtk is row straight multiplying give us complete matrix that is X by Theta
                    times1(1,1,:) = (cos(timeV*F1W-Phi1(M,N))-1i*sin(timeV*F1W-Phi1(M,N)));% create time vector with phase shift as a 1 by 1 by t vector
                    W1 = W1 + amp1.*times1; %because the first two dims of times are 1 the .* causes the entire matrix amp to be multiplied and by each value of times and gives an M by N by t array
                end
            end
            obj.WtotalD = W1;
        end
        function calculate_2_force_displacement(obj)
            %precalculate values that will be used a lot
            W1 = zeros(length(obj.widths), length(obj.heights), length(obj.time));
            W2 = zeros(length(obj.widths), length(obj.heights), length(obj.time));
            dwx = zeros(length(obj.widths), length(obj.heights), length(obj.time));
            dwy = zeros(length(obj.widths), length(obj.heights), length(obj.time));
            dwx2 = zeros(length(obj.widths), length(obj.heights), length(obj.time));
            dwxy = zeros(length(obj.widths), length(obj.heights), length(obj.time));
            dwy2 = zeros(length(obj.widths), length(obj.heights), length(obj.time));
            dwx3 = zeros(length(obj.widths), length(obj.heights), length(obj.time));
            dwx2y = zeros(length(obj.widths), length(obj.heights), length(obj.time));
            dwxy2 = zeros(length(obj.widths), length(obj.heights), length(obj.time));
            dwy3 = zeros(length(obj.widths), length(obj.heights), length(obj.time));
            Odotx = zeros(length(obj.widths), length(obj.heights), length(obj.time));
            Odoty = zeros(length(obj.widths), length(obj.heights), length(obj.time));
            Phi1 = obj.complexDampingPhase + obj.force1.phase; % combine complex phase lag and phase shift
            Phi2 = obj.complexDampingPhase + obj.force2.phase;
            timeV = obj.time;
            F1W = obj.force1.omega;
            F2W = obj.force2.omega;
            Widths = obj.widths;
            Height = obj.heights;
            pix1_a = pi*obj.force1.location(1)/obj.geometry.width;
            piy1_b = pi*obj.force1.location(2)/obj.geometry.height;
            pix2_a = pi*obj.force2.location(1)/obj.geometry.width;
            piy2_b = pi*obj.force2.location(2)/obj.geometry.height;
            Cf1 = 4*obj.force1.magnitude/(pi^2*obj.material.density*obj.geometry.thickness);
            Cf2 = 4*obj.force2.magnitude/(pi^2*obj.material.density*obj.geometry.thickness);
            
            % loops over M modes along width
            for M = 1:obj.mesh.width_modes
                sinX1 = sin(M*pix1_a);
                sinX2 = sin(M*pix2_a);
                sinWidths = (sin(M*pi*Widths/obj.geometry.width))';
                dWx = (M*pi/obj.geometry.width)*(cos(M*pi*Widths/obj.geometry.width))';
                dWx2 = -(M*pi/obj.geometry.width)^2*(sin(M*pi*Widths/obj.geometry.width))';
                dWx3 = -(M*pi/obj.geometry.width)^3*(cos(M*pi*Widths/obj.geometry.width))';
                Cf1_m = Cf1/M;
                Cf2_m = Cf2/M;
                % loops over N modes along height
                for N = 1:obj.mesh.height_modes
                    Cf1_mn = Cf1_m/N;
                    Cf2_mn = Cf2_m/N;
                    sinY1 = sin(N*piy1_b);
                    sinY2 = sin(N*piy2_b);
                    sinHeights = sin(N*pi*Height/obj.geometry.height);
                    dWy = (N*pi/obj.geometry.height) * cos(N*pi*Height/obj.geometry.height);
                    dWy2 = -(N*pi/obj.geometry.height)^2 * sin(N*pi*Height/obj.geometry.height);
                    dWy3 = -(N*pi/obj.geometry.height)^3 * cos(N*pi*Height/obj.geometry.height);
                    denom = obj.calculate_denominator(M,N);
                    complexDampCoeff = obj.calculate_complex_damping_coefficient(M,N);
                    allScalars1 = Cf1_mn*sinX1*sinY1/denom;
                    allScalars2 = Cf2_mn*sinX2*sinY2/denom; 
                    gridEvals = (sinWidths*sinHeights)*complexDampCoeff;
                    amp1 = allScalars1*gridEvals;  % because sinMxj is column and cosNtk is row straight multiplying give us complete matrix that is X by Theta
                    amp2 = allScalars2*gridEvals;
                    K1 = allScalars1*complexDampCoeff;
                    K2 = allScalars2*complexDampCoeff;
                    dxamp1 = K1*(dWx*sinHeights);
                    dxamp2 = K2*(dWx*sinHeights);
                    dyamp1 = K1 * (sinWidths*dWy);
                    dyamp2 = K2 * (sinWidths*dWy);
                    dx2amp1 = K1 *(dWx2*sinHeights);
                    dx2amp2 = K2 *(dWx2*sinHeights);
                    dxyamp1 = K1 *(dWx*dWy);
                    dxyamp2 = K2 *(dWx*dWy);
                    dy2amp1 = K1 *(sinWidths*dWy2);
                    dy2amp2 = K2* (sinWidths*dWy2) ;
                    dx3amp1 = K1 *(dWx3*sinHeights);
                    dx3amp2 = K2 *(dWx3*sinHeights);
                    dx2yamp1 = K1 *(dWx2*dWy);
                    dx2yamp2 = K2 *(dWx2*dWy);
                    dxy2amp1 = K1 *(dWx*dWy2);
                    dxy2amp2 = K2 *(dWx*dWy2); 
                    dy3amp1 = K1 *(sinWidths*dWy3);
                    dy3amp2 = K2 *(sinWidths*dWy3);
                    times1(1,1,:) = (cos(timeV*F1W-Phi1(M,N))-1i*sin(timeV*F1W-Phi1(M,N)));% create time vector with phase shift as a 1 by 1 by t vector
                    times2(1,1,:) = (cos(timeV*F2W-Phi2(M,N))-1i*sin(timeV*F2W-Phi2(M,N)));
                    dtimes1(1,1,:) = 1i*F1W*(cos(timeV*F1W-Phi1(M,N))-1i*sin(timeV*F1W-Phi1(M,N)));
                    dtimes2(1,1,:) = 1i*F2W*(cos(timeV*F2W-Phi1(M,N))-1i*sin(timeV*F1W-Phi1(M,N)));
                    W1 = W1 + amp1.*times1; %because the first two dims of times are 1 the .* causes the entire matrix amp to be multiplied and by each value of times and gives an M by N by t array
                    W2 = W2 + amp2.*times2;
                    dwx = dwx + dxamp1.*times1 + dxamp2.*times2;
                    dwy = dwy + dyamp1.*times1 + dyamp2.*times2;
                    dwx2 = dwx2 + dx2amp1.*times1 + dx2amp2.*times2;
                    dwxy = dwxy + dxyamp1.*times1 + dxyamp2.*times2;
                    dwy2 = dwy2 + dy2amp1.*times1 + dy2amp2.*times2;
                    dwx3 = dwx3 + dx3amp1.*times1 + dx3amp2.*times2;
                    dwx2y = dwx2y + dx2yamp1.*times1 + dx2yamp2.*times2;
                    dwxy2 = dwxy2 + dxy2amp1.*times1 + dxy2amp2.*times2;
                    dwy3 = dwy3 + dy3amp1.*times1 + dy3amp2.*times2;
                    Odotx = Odotx + dyamp1.*dtimes1 + dyamp2.*dtimes2;
                    Odoty = Odoty + dxamp1.*dtimes1 + dyamp2.*dtimes2;
                end
            end
            obj.WtotalD = W1+W2;
            obj.dWdx = dwx;
            obj.dWdy = dwy;
            obj.dWdx2 = dwx2;
            obj.dWdxy = dwxy;
            obj.dWdy2 = dwy2;
            obj.dWdx3 = dwx3;
            obj.dWdx2y = dwx2y;
            obj.dWdxy2 = dwxy2;
            obj.dWdy3 = dwy3;
            obj.Oxdot = Odotx;
            obj.Oydot = Odoty;
        end
        function calculate_spectral_terms(obj)
            timeDim = 3;% time is the third dimension of the data
            fullFFT = 1/(size(obj.WtotalD,3))*fft(obj.WtotalD, [], timeDim); % take the FFT over time dimension
            dxFFT = 1/(size(obj.WtotalD,3))*fft(obj.dWdx, [], timeDim);
            dyFFT = 1/(size(obj.WtotalD,3))*fft(obj.dWdy, [], timeDim);
            dx2FFT = 1/(size(obj.WtotalD,3))*fft(obj.dWdx2, [], timeDim);
            dxyFFT = 1/(size(obj.WtotalD,3))*fft(obj.dWdxy, [], timeDim);
            dy2FFT = 1/(size(obj.WtotalD,3))*fft(obj.dWdy2, [], timeDim);
            dx3FFT = 1/(size(obj.WtotalD,3))*fft(obj.dWdx3, [], timeDim);
            dx2yFFT = 1/(size(obj.WtotalD,3))*fft(obj.dWdx2y, [], timeDim);
            dxy2FFT = 1/(size(obj.WtotalD,3))*fft(obj.dWdxy2, [], timeDim);
            dy3FFT = 1/(size(obj.WtotalD,3))*fft(obj.dWdy3, [], timeDim);
            OxFFT = 1/(size(obj.WtotalD,3))*fft(obj.Oxdot , [], timeDim);
            OyFFT = 1/(size(obj.WtotalD,3))*fft(obj.Oydot, [], timeDim);
%             obj.WtotalF = zeros(size(fullFFT(:,:,1)));
            %identify max location avoid edges because they are supported
            %and have no displacement
            [~, indx] = max(squeeze(fullFFT(2,2,:)));
            obj.WtotalF = fullFFT(:,:,indx);
            obj.dWdx = dxFFT(:,:,indx);
            obj.dWdy = dyFFT(:,:,indx);
            obj.dWdx2 = dx2FFT(:,:,indx);
            obj.dWdxy = dxyFFT(:,:,indx);
            obj.dWdy2 = dy2FFT(:,:,indx);
            obj.dWdx3 = dx3FFT(:,:,indx);
            obj.dWdx2y = dx2yFFT(:,:,indx);
            obj.dWdxy2 = dxy2FFT(:,:,indx);
            obj.dWdy3 = dy3FFT(:,:,indx);
            obj.Oxdot= OxFFT(:,:,indx);
            obj.Oydot = OyFFT(:,:,indx);
        end
        function write_velocity_file(obj)
            obj.write_file('vel')
        end
        function write_displacement_file(obj)
            obj.write_file('disp')
        end
        function write_file(obj, type)
            switch type
                case 'vel'
                    tstr = 'Velocity.txt';
                case 'disp'
                    tstr = 'Displacement.txt';
            end
            old = cd(obj.saveFolder);
            DataFileName = strcat(obj.fileName, "_Simulated_", tstr);
            fid = fopen(DataFileName,'wt');
            obj.write_header(fid, type, obj.force1.frequency)
            obj.write_data(fid, obj.data, type);
            fclose(fid);
            cd(old);
        end
        function calculate_mode_frequencies_and_phases(obj)
            scale = pi^2*sqrt(obj.D/(obj.material.density*obj.geometry.thickness));
            obj.nOmega = zeros(obj.mesh.height_modes,obj.mesh.width_modes);
            obj.complexDampingPhase = zeros(size(obj.nOmega));
            for m = 1:obj.mesh.height_modes
                for n = 1:obj.mesh.width_modes
                    obj.nOmega(m,n) = ((m/obj.geometry.height)^2+(n/obj.geometry.width)^2)*scale;
                    obj.complexDampingPhase(m,n) = atan(2*obj.lambda*(obj.force1.omega/obj.nOmega(m,n))/(1-(obj.force1.omega/obj.nOmega(m,n))^2));
                end
            end
        end
        function create_default(obj)
            obj.geometry = plateGeometry();
            obj.force1 = sinusoidalForce(3);
            obj.force2 = sinusoidalForce(4);
            obj.material = simMaterial();
            obj.mesh = meshRectDetails();
            obj.timesteps = 20;
            obj.saveFolder = "C:\Users\ME\Desktop\Simulated Data\";
        end
        function create_default_at_frequency(obj, frequency)
            obj.geometry = plateGeometry();
            obj.force1 = sinusoidalForce( [.05, .1], 40, frequency, 0);
            obj.force2 = sinusoidalForce( [.35, .4], 40, frequency, 179*pi/180);
            obj.material = simMaterial();
            obj.mesh = meshRectDetails();
            obj.timesteps = 20;
            obj.saveFolder = "C:\Users\ME\Desktop\Simulated Data\";
        end
        function denom = calculate_denominator(obj,M,N)
            denom =(obj.nOmega(M,N)^2-obj.force1.omega^2)^2+obj.material.eta^2*obj.nOmega(M,N)^4;
        end
        function imagCoef = calculate_complex_damping_coefficient(obj, M,N)
            imagCoef = (1-1i*obj.material.eta)*obj.nOmega(M,N)^2-obj.force1.omega^2;
        end
        function collectData(obj)
            % put into column vectors
            [W, H] = meshgrid(obj.widths, obj.heights);
            tempData.xCoor =reshape(W, numel(W),1);
            tempData.yCoor =reshape(H, numel(H), 1);
            tempData.zCoor = ones(size(tempData.yCoor));
            % no displacement in X or Y
            tempData.rX = zeros(size(tempData.xCoor));
            tempData.iX = zeros(size(tempData.xCoor));
            tempData.rY = zeros(size(tempData.yCoor));
            tempData.iY = zeros(size(tempData.yCoor));
            tempData.rXv = zeros(size(tempData.xCoor));
            tempData.iXv = zeros(size(tempData.xCoor));
            tempData.rYv = zeros(size(tempData.yCoor));
            tempData.iYv = zeros(size(tempData.yCoor));
            %
            fullDisp = reshape(obj.WtotalF, numel(obj.WtotalF), 1);
            fullVel = 1i*obj.force1.omega*fullDisp;
            tempData.rZ = real(fullDisp);
            tempData.iZ = imag(fullDisp);
            tempData.rZv = real(fullVel);
            tempData.iZv = imag(fullVel);
            obj.data = tempData;            
        end
        function [Mx, My, Mxy, Qx, Qy]= calculate_moments_and_shears(obj)
            nu = obj.material.poisson;
            Mx = -obj.D * (obj.dWdx2 + nu*obj.dWdy2);  
            My = -obj.D * (obj.dWdy2 + nu*obj.dWdx2);
            Mxy = -obj.D * (1-nu) * obj.dWdxy;
            Qy = -obj.D * (obj.dWdx3 + obj.dWdxy2);
            Qx = -obj.D * (obj.dWdx2y + obj.dWdy3);
        end
        function [wdotS, OxdotS, OydotS] = conjugate_velocities(obj)
            wdotS =  conj(1i*obj.force1.omega*obj.WtotalF);
            OxdotS = conj(obj.Oxdot);
            OydotS = conj(obj.Oydot); 
        end
    end
    methods % getters
        function value = get.fileName(obj)
            date = datestr(floor(now));
            date = strrep(date,'-','_');
            value = strcat('plate_simulation_2_forces_', date, '_', num2str(obj.force1.frequency), '_Hz');
        end
        function value = get.D(obj)
            value = obj.material.E*obj.geometry.thickness^3/(12*(1-obj.material.poisson));
        end
        function value = get.heights(obj)
            value = linspace(0,obj.geometry.height,obj.mesh.height_divisions+1);
        end
        function value = get.widths(obj)
            value = linspace(0, obj.geometry.width, obj.mesh.width_divisions+1);
        end
        function value = get.SpikeN(obj)
            dfreq = 1/max(obj.time);
            value = round(obj.force1.frequency/dfreq);
        end
        function value = get.time(obj)
            t =  linspace(0, 1/obj.force1.frequency, obj.timesteps+1);
            value = t(1:end-1);
        end
        function value = get.lambda(obj)
            value = obj.material.density*obj.geometry.thickness*obj.material.eta/obj.force1.omega;
        end
        function value = get.naturalFrequencies(obj)
            value = obj.nOmega/(2*pi);
        end
    end
    methods (Static)
        function write_header(fid, type, frequency)
            switch type
                case 'vel'
                    Type = 'Velocity';
                    unit = 'm/s';
                case 'disp'
                    Type = 'Displacement';
                    unit = 'm';
            end
            fprintf(fid, strcat("Source File Name:\t", "Simulated Data\n"));
            fprintf(fid, strcat("Signal: FFT - Vib 3D ", Type, " - Real & Imag.\n\n"));
            fprintf(fid, strcat("Band No.:\t 1\nFrequency:\t", num2str(frequency), " Hz\n\n\n"));
            fprintf(fid, "Interpolated:\tYes\nFiltered:\tYes\n\n");
            fprintf(fid, strcat("Index\tX\tY\tZ\tReal X [", unit,"]\tReal Y [", unit,"]\tReal Z [", unit,"]\tImaginary X [", unit,"]\t Imaginary Y [", unit,"]\tImaginary Z [", unit,"]\n"));
        end
        function write_data(fid, data, type)
            dataFormat = strcat('%d',repmat('\t%g', 1, 9), '\n');
            X = data.xCoor;  Y = data.yCoor;   Z = data.zCoor;
            switch type
                case 'vel'
                    rX = data.rXv;  rY = data.rYv;   rZ = data.rZv;
                    iX = data.iXv;  iY = data.iYv;   iZ = data.iZv;
                case 'disp'
                    rX = data.rX;  rY = data.rY;   rZ = data.rZ;
                    iX = data.iX;  iY = data.iY;   iZ = data.iZ;
            end
            for point = 1:length(X)
                fprintf(fid, dataFormat, point, X(point), Y(point), Z(point), rX(point), rY(point), rZ(point), iX(point), iY(point), iZ(point));
            end
        end
    end
end
