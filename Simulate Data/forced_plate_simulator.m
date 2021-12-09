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
        timesteps
        nOmega= [];
        complexDampingPhase = [];
        WtotalD
        WtotalF
        data
        saveFolder
    end
    methods
        function in_testing(obj)
            obj.calculate_mode_frequencies_and_phases
            obj.calculate_2_force_displacement
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
            obj.calculate_2_force_displacement();
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
    end
    methods (Access = private)
        function calculate_2_force_displacement(obj)
            %precalculate values that will be used a lot
            W1 = zeros(length(obj.widths), length(obj.heights), length(obj.time));
            W2 = zeros(length(obj.widths), length(obj.heights), length(obj.time));
            Phi1 = obj.complexDampingPhase + obj.force1.phase;
            Phi2 = obj.complexDampingPhase + obj.force2.phase;
            timeV = obj.time;
            F1W = obj.force1.omega;
            F2W = obj.force2.omega;
            Widths = obj.widths;
            Height = obj.heights;
            pix1_a = pi*obj.force1.location(1)/obj.geometry.height;
            piy1_b = pi*obj.force1.location(2)/obj.geometry.width;
            pix2_a = pi*obj.force2.location(1)/obj.geometry.height;
            piy2_b = pi*obj.force2.location(2)/obj.geometry.width;
            Cf1 = 4*obj.force1.magnitude/(pi^2*obj.material.density*obj.geometry.height);
            Cf2 = 4*obj.force2.magnitude/(pi^2*obj.material.density*obj.geometry.height);
            % loops over M modes along width
            for M = 1:obj.mesh.width_modes
                sinX1 = sin(M*pix1_a);
                sinX2 = sin(M*pix2_a);
                sinWidths = (sin(M*pi*Widths/obj.geometry.width))';
                Cf1_m = Cf1/M;
                Cf2_m = Cf2/M;
                % loops over N modes along height
                for N = 1:obj.mesh.height_modes
                    Cf1_mn = Cf1_m/N;
                    Cf2_mn = Cf2_m/N;
                    sinY1 = sin(N*piy1_b);
                    sinY2 = sin(N*piy2_b);
                    sinHeights = sin(N*pi*Height/obj.geometry.height);
                    denom = obj.calculate_denominator(M,N);
                    complexDampCoeff = obj.calculate_complex_damping_coefficient(M,N);
                    allScalars1 = Cf1_mn*sinX1*sinY1/denom;
                    allScalars2 = Cf2_mn*sinX2*sinY2/denom; 
                    gridEvals = (sinWidths*sinHeights)*complexDampCoeff;
                    amp1 = allScalars1*gridEvals;  % because sinMxj is column and cosNtk is row straight multiplying give us complete matrix that is X by Theta
                    amp2 = allScalars2*gridEvals;
                    times1(1,1,:) = (cos(timeV*F1W-Phi1(M,N))-1i*sin(timeV*F1W-Phi1(M,N)));% create time vector with phase shift as a 1 by 1 by t vector
                    times2(1,1,:) = (cos(timeV*F2W-Phi2(M,N))-1i*sin(timeV*F2W-Phi2(M,N)));
                    W1 = W1 + amp1.*times1; %because the first two dims of times are 1 the .* causes the entire matrix amp to be multiplied and by each value of times and gives an M by N by t array
                    W2 = W2 + amp2.*times2;
                end
            end
            obj.WtotalD = W1+W2;
        end
        function calculate_spectral_terms(obj)
            timeDim = 3;% time is the third dimension of the data
            fullFFT = 1/(size(obj.WtotalD,3))*fft(obj.WtotalD, [], timeDim); % take the FFT over time dimension
            obj.WtotalF = zeros(size(fullFFT(:,:,1)));
            %identify max location avoid edges because they are supported
            %and have no displacement
            [~, indx] = max(squeeze(fullFFT(2,2,:)));
            obj.WtotalF = fullFFT(:,:,indx);
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
            tempData.rX = zeros(size(tempData.yCoor));
            tempData.iX = zeros(size(tempData.yCoor));
            tempData.rY = zeros(size(tempData.yCoor));
            tempData.iY = zeros(size(tempData.yCoor));
            tempData.rXv = zeros(size(tempData.yCoor));
            tempData.iXv = zeros(size(tempData.yCoor));
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
            dataFormat = strcat('%d',repmat('\t%0.8f', 1, 9), '\n');
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
