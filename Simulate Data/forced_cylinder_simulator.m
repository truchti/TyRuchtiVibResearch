classdef forced_cylinder_simulator< handle
    % This class simulates the displacement field of a cylinder subjected
    % to 2 point loads both driving at the same frequency but they can be
    % phased independently. 
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
        xyzForce1
        xyzForce2
        trzForce1
        trzForce2
    end
    properties (Hidden = true, Dependent = true)
        xs
        thetas
        xf1c
        xf2c
        C1
        C2
        lambda
        time
        SpikeN
        sortedNFs
    end
    properties (Hidden = true)
        timesteps =20;
        nOmega= [];
        complexDampingPhaseLag = [];
        WtotalD
        WtotalF
        data
        saveFolder
    end
    methods
        function obj = forced_cylinder_simulator(geometry, forces, material, mesh)
            %% creates the object that simulates the displacement of a cylinder to 2 forces at a gived frequency
            if nargin < 1
                obj.create_default();
            elseif nargin < 2
                obj.create_default_at_frequency(geometry);
            else
                obj.geometry = geometry;
                obj.force1 = forces(1);
                obj.force2 = forces(2);
                obj.material = material;
                obj.mesh = mesh;
            end
        end
        function simulate_and_save(obj)
            obj.simulate_data();
            obj.save_data();
        end
        function simulate_data(obj)
            obj.calculate_mode_frequencies();
            obj.calculate_complex_damping_phase_lag();
            obj.calculate_displacement_through_space_and_time();
            obj.calculate_spectral_terms();
            obj.calculate_rectangular_data();
        end
        function save_data(obj)
            % this function creates both the velocity and displacement
            % files that mimic the SLDV output from the simulated displacements 
            if ~isempty(obj.data)
                obj.write_velocity_file();
                obj.write_displacement_file();
            else
                warning("data has not been simulated")
            end
        end
        function disable_force(obj, force_number)
            % turn off a force from being applied
            switch force_number
                case 1
                    obj.force1.enabled = false;
                case 2
                    obj.force2.enabled = false;
            end
        end
        function enable_force(obj, force_number)
            % turn on a force so that it affects displacement
            switch force_number
                case 1
                    obj.force1.enabled = true;
                case 2
                    obj.force2.enabled = true;
            end
        end
        function show_radial_displacement(obj)
            subplot(1,2,1)
            surf(obj.thetas, obj.xs, zeros(length(obj.xs), length(obj.thetas)),real(obj.WtotalF))
            a = [obj.force1.location; obj.force2.location];
            hold on;
            scatter(wrapTo2Pi(a(:,1)), a(:,2), 'r')
            colorbar;
            view(0,90)
            hold off;
            subplot(1,2,2)
            surf(obj.thetas, obj.xs, zeros(length(obj.xs), length(obj.thetas)),imag(obj.WtotalF))
            hold on;
            scatter(wrapTo2Pi(a(:,1)), a(:,2), 'r')
            colorbar;
            view(0,90)
            hold off;
         end
        function animate_displacement_in_time(obj, imagData, loops)
            if nargin < 3 
                loops = 3;
            end
            if nargin < 2 
                imagData = false;
            end
            [T, H] = meshgrid([obj.thetas 0], obj.xs);
            x = obj.geometry.radius *cos(T);
            y = obj.geometry.radius * sin(T);
            z = H;
            if imagData
                values = imag(obj.WtotalD);
            else
                values = real(obj.WtotalD);
%                 ivalues = imag(obj.WtotalD);
            end
            loopedValues = zeros(size(values) + [0 1 0]);
%             iloopedValues = zeros(size(ivalues) + [0 1 0]);
            for i = 1:size(values,3)
                loopedValues(:,:,i) = [values(:,:,i)  values(:,1,i)];
%                 iloopedValues(:,:,i) = [ivalues(:,:,i) ivalues(:,1,i)];
            end
            figure(1)
            s = surf(x,y,z,loopedValues(:,:,1));
            axis equal
            s.EdgeAlpha = 0.0;
            s.FaceColor = 'interp';
            colorbar
            clim = max(max(max(values)));
            caxis([-clim, clim])
            sz = size(values,3);
            for i = 1:loops*sz
                s.CData = loopedValues(:,:,mod(i,sz)+1);
%                 s2.CData = iloopedValues(:,:,mod(i, sz)+1);
                pause(0.1);
            end
        end             
        function show_cyl_radial_disp(obj, imagine)
            if nargin < 2 
                imagine = false;
            end
            rad = obj.geometry.radius;
            [T, H] = meshgrid([obj.thetas 0], obj.xs);
            x = rad*cos(T);
            y = rad*sin(T);
            z = H;
            f = [obj.force1.location; obj.force2.location];
            xForces = rad*cos(f(:,1));
            yForces = rad*sin(f(:,1));
            zForces = f(:,2);
            if imagine
                value = imag(obj.WtotalF(:,:,1));
            else
                value = real(obj.WtotalF(:,:,1));
            end
            value = [value value(:,1)];
            s = surf(x,y,z,value);
            s.EdgeAlpha = 0.0;
            s.FaceColor = 'interp';
            colorbar
            hold on;
            scatter3(xForces, yForces, zForces, 'dr')
            hold off;
        end
    end
    methods (Access = private)
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
        function calculate_mode_frequencies(obj, useSoedel)
            % This function calculates the natural frequencies of the modes
            % of the cylinder it also calculates the phase associated
            % with the each mode because of the damping
            if nargin <2 || ~useSoedel
                obj.rao_nat_freq(obj.mesh.longitude_modes, obj.mesh.theta_modes);
            else
                obj.soedel_nat_freq(obj.mesh.longitude_modes, obj.mesh.theta_modes);
            end            
        end
        function calculate_complex_damping_phase_lag(obj)
            zpart = obj.lambda/(2*obj.material.density*obj.geometry.thickness); % damping Coeff
            fomg1 = obj.force1.omega; % force Freq
            for m = 1:obj.mesh.longitude_modes
                for n = 0:obj.mesh.theta_modes
                    zeta = zpart/obj.nOmega(m,n+1); % soedel eq 8.3.3 (pg. 212)
                    omegRatio = fomg1/obj.nOmega(m,n+1);
                    obj.complexDampingPhaseLag(m,n+1) = atan2(2*zeta*(omegRatio),(1-(omegRatio)^2)); % sodel eq 8.5.6 pg 215
                end
            end
        end
        function omgs = rao_nat_freq(obj, mModes, nModes)   
            E = obj.material.E;
            nu = obj.material.poisson;
            rho = obj.material.density;
            h = obj.geometry.thickness;
            R = obj.geometry.radius;
            L = obj.geometry.height;
            
            mu = h^2/(12*R^2);
            a1 = (1-nu)/2;
            a2 = (1+nu)/2;
            O = (1-nu^2)*R^2*rho/E;
            %preallocate omega^2 array
            w_squared = zeros(1, mModes, nModes+1);
            for m = 1:mModes
                y = m*pi*R/L;
                for n = 0:nModes
                    b1 = -y^2*(1+a1+2*n^2*mu) - n^2*(1+a1) - y^4*mu - n^4*mu - 1;
                    b2 = y^6*mu*(1+a1) + y^4*(a1+3*n^2*mu+3*a1*n^2*mu) + ...
                        y^2*(1 + n^2 + a1^2*n^2 - a2^2*n^2 - nu^2 + a1 + 3*a1*n^4*mu + 3*n^4*mu) + ...
                        n^6*mu*(1+a1) + n^4*a1 + n^2*a1;
                    
                    b3 = -y^8*a1*mu - y^6*n^2*mu*(1+2*a1-a2^2+a1^2) + ...
                        -y^4*(a1 + 2*a1*n^4*mu + 2*a1^2*n^4*mu + 2*n^4*mu - 2*a2^2*n^4*mu - a2*nu^2) + ...
                        -y^2*(n^6*mu*(1 + a1^2 + 2*a1 - a2^2) + n^2*(-a2^2 + a1^2 + 2*a2*nu - nu^2)) - n^8*a1*mu;
                    % mode associated with transverse deflection is lowest of roots
                    w_squared(1,m,n+1) = min(roots([O^3, b1*O^2, b2*O, b3]));
                end
            end
            if nargout ==0
                obj.nOmega = squeeze(sqrt(w_squared));
            else
                omgs = squeeze(sqrt(w_squared));
            end
        end
        function omgs = soedel_nat_freq(obj,mModes, nModes)
            E = obj.material.E;
            mu = obj.material.poisson;
            rho = obj.material.density;
            h = obj.geometry.thickness;
            a = obj.geometry.radius;
            L = obj.geometry.height;
            % stiffnesses
            K = E*h/(1-mu^2); % soedel eq 2.5.10; Membrane stiffness
            D = E*h^3/(12*(1-mu^2)); % soedel eq 2.5.16 Bending Stiffness            
            w_squared = zeros(3,mModes, nModes+1);
            for m = 1:mModes
                for n = 0:nModes
                    % eqs 5.5.67-72 k matrix
                    k11 = K*((m*pi/L)^2 + (1-mu)/2*(n/a)^2);
                    k12 = K*(1+mu)/2*(m*pi/L)*n/a;
                    k13 = mu*K/a*(m*pi/L);
                    k22 = (K+ D/a^2)*(((1-mu)/2)*(m*pi/L)^2+ (n/a)^2);
                    k23 = -K*n/a^2-((D*n)/a^2)*((m*pi/L)^2+(n/a)^2);
                    k33 = D*((m*pi/L)^2 + (n/a)^2)^2+K/a^2;
                    
                    ph = (rho*h);
                    % coefficients of determinate
                    a1 = -1/(ph)*(k11+k22+k33); % eq 5.5.74
                    a2 = 1/(ph)^2*(k11*k33+k22*k33+k11*k22-k23^2-k12^2-k13^2); % eq 5.5.75
                    a3 = 1/(ph)^3*(k11*k23^2+k22*k13^2+k33*k12^2 + 2*k12*k23*k13-k11*k22*k33); % eq 5.5.76
                    
                    % solution to polynomial w^6 + a1*w^4 + a2*w^2 + a3
                    %verified alpha
                    alpha = acos((27*a3+2*a1^3-9*a1*a2)/(2*(a1^2-3*a2)^(3/2)));
                    %verified w_squared
                    w_squared(1,m,n+1) = -2/3 * (a1^2-3*a2)^(1/2)*cos(alpha/3)-a1/3;
                    w_squared(2,m,n+1) = -2/3 * (a1^2-3*a2)^(1/2)*cos((alpha+2*pi)/3)-a1/3;
                    w_squared(3,m,n+1) = -2/3 * (a1^2-3*a2)^(1/2)*cos((alpha+4*pi)/3)-a1/3;
                    %test if w^2 is actually roots as calcualted
                    %         w2test = roots([1, a1, a2, a3]);
                end
            end
            w = w_squared.^(1/2);
            if nargout == 0
                obj.nOmega = squeeze(w(1,:,:));
            else
                omgs = squeeze(w(1,:,:));
            end
        end
        function calculate_single_force_displacement(obj, force_num)
            % if one of the forces is disabled only calculate the
            % displacement from that force
            W = zeros(length(obj.xs), length(obj.thetas), length(obj.time));
            timeV = obj.time;
            Thetas = obj.thetas;
            Xs = obj.xs;
            switch force_num
                case 1
                    Phi = obj.complexDampingPhaseLag + obj.force1.phase;
                    FW = obj.force1.omega;
                    tstar = obj.force1.tStar;
                    C = obj.C1;
                    xfc = obj.xf2c;
                case 2
                    Phi = obj.complexDampingPhaseLag + obj.force2.phase;
                    FW = obj.force2.omega;
                    tstar = obj.force2.tStar;
                    C = obj.C2;
                    xfc = obj.xf2c;
            end
            for M = 1:obj.mesh.longitude_modes
                CsinM = C*sin(M*xfc);
                for N = 0:obj.mesh.theta_modes-1
                    denom = obj.calculate_denominator(force_num,M,N);
                    sinMxj = sin(M*pi*Xs/obj.geometry.height);
                    cosNtk = cos(N*(Thetas-tstar));
                    amplitude = (sinMxj' * cosNtk) *CsinM/denom;  
                    timeValues(1,1,:) = (cos(timeV*FW-Phi(M,N+1))+1i*sin(timeV*FW-Phi(M,N+1)));
                    W = W + amplitude.*timeValues;
                end
            end
            obj.WtotalD = W;
        end
        function calculate_displacement_through_space_and_time(obj)
            if obj.force1.enabled
                if obj.force2.enabled
                    obj.calculate_2_force_displacement();
                else
                    obj.calculate_single_force_displacement(1);
                end
            else
                obj.calculate_single_force_displacement(2);
            end                
        end
        function calculate_2_force_displacement(obj)
            % this calculates the displacement function in terms of height
            % theta and time
            % precalculate values that will be used a lot
            W1 = zeros(length(obj.xs), length(obj.thetas), length(obj.time));
            W2 = zeros(length(obj.xs), length(obj.thetas), length(obj.time));
            Phi1 = obj.complexDampingPhaseLag + obj.force1.phase;
            Phi2 = obj.complexDampingPhaseLag + obj.force2.phase;
            timeV = obj.time;
            F1W = obj.force1.omega;
            F2W = obj.force2.omega;
            t1star = obj.force1.tStar;
            t2star = obj.force2.tStar;
            Thetas = obj.thetas;
            Xs = obj.xs;
            c1 = obj.C1;
            c2 = obj.C2;
            xfc1 = obj.xf1c;
            xfc2 = obj.xf2c;
            % loops over M modes in longitudinal Direction
            for M = 1:obj.mesh.longitude_modes
                CsinM1 = c1*sin(M*xfc1);
                CsinM2 = c2*sin(M*xfc2);
                sinMxj = (sin(M*pi*Xs/obj.geometry.height))';
                % loops over N modes in theta Direction
                for N = 0:obj.mesh.theta_modes-1
                    denom1 = obj.calculate_denominator(1, M,N);
                    denom2 = obj.calculate_denominator(2, M,N);
                    cosNtk1 = cos(N*(Thetas-t1star));
                    cosNtk2 = cos(N*(Thetas-t2star));
                     % because sinMxj is column and cosNtk is row straight multiplying give us complete matrix that is X by Theta
                    amp1 = sinMxj * (cosNtk1 * CsinM1/denom1); 
                    amp2 = sinMxj * (cosNtk2 * CsinM2/denom2);
                    % create time vector with phase shift as a 1 by 1 by t vector
                    times1(1,1,:) = (cos(timeV*F1W-Phi1(M,N+1))-1i*sin(timeV*F1W-Phi1(M,N+1)));
                    times2(1,1,:) = (cos(timeV*F2W-Phi2(M,N+1))-1i*sin(timeV*F2W-Phi2(M,N+1)));
                    %because the first two dimensionss of the multiplication are 1 the .* causes the entire matrix, amp1, to be multiplied by each value in array, times, and results in an M by N by t array
                    W1 = W1 + amp1.*times1; 
                    W2 = W2 + amp2.*times2;
                end
            end
            obj.WtotalD = W1+W2;
        end
        function calculate_spectral_terms(obj)
            % this fuction calculates the FFT of the diplacement field to
            % get the complex values for displacement
            
            timeDim = 3;% time is the third dimension of the data
            fullFFT = 1/(size(obj.WtotalD,3))*fft(obj.WtotalD, [], timeDim); % take the FFT over time dimension
%             obj.WtotalF = zeros(size(fullFFT(:,:,1)));
            %identify max location
            [~, indx] = max(squeeze(fullFFT(2,2,:))); %% first row  is simply supported so it will always be zero
            obj.WtotalF = fullFFT(:,:,indx);
            %check accuracy of FFT
%             plot(obj.time, real(squeeze(obj.WtotalD(m,n,:))), 'k', obj.time, abs(obj.WtotalF(m,n))*cos(obj.time*obj.force1.omega +angle(obj.WtotalF(m,n))),'--r')
        end
        function denom = calculate_denominator(obj,forceNum,M,N)
            if forceNum == 1
                omegaRatiosqrd = (obj.force1.omega/obj.nOmega(M,N+1))^2;
            elseif forceNum == 2
                omegaRatiosqrd = (obj.force2.omega/obj.nOmega(M,N+1))^2;
            else
                error("Invalid Force Number");
            end
            zeta = obj.lambda/(2*obj.material.density*obj.geometry*thickness*obj.nOmega(M,N+1));
            denom = obj.nOmega(M,N+1)^2.*((1 - omegaRatiosqrd)^2 + (4*zeta^2*omegaRatiosqrd))^0.5;
            if ~N
                denom = denom*2;
            end
        end
        function create_default(obj)
            obj.geometry = cylinderGeometry();
            obj.force1 = sinusoidalForce(1);
            obj.force2 = sinusoidalForce(2);
            obj.material = simMaterial();
            obj.mesh = meshDetails();
            obj.timesteps = 20;
            obj.saveFolder = "C:\Users\ME\Desktop\Simulated Data\";
        end
        function create_default_at_frequency(obj, frequency)
            obj.geometry = cylinderGeometry();
            obj.force1 = sinusoidalForce( [.1, -60*pi/180], 500, frequency, 0);
            obj.force2 = sinusoidalForce( [.4, -240*pi/180], 500, frequency, 179*pi/180);
            obj.material = simMaterial();
            obj.mesh = meshDetails();
            obj.timesteps = 20;
            obj.saveFolder = "C:\Users\ME\Desktop\Simulated Data\";
        end
        function calculate_rectangular_data(obj)
            % converts the displacement from the cylindrical coordinate
            % system it is calculated in to the cartesian coordinate system
            % which the SLDV measures in
            xtrans = reshape(repmat(cos(obj.thetas), length(obj.xs),1), length(obj.thetas)*length(obj.xs),1);
            ytrans = reshape(repmat(sin(obj.thetas), length(obj.xs),1), length(obj.thetas)*length(obj.xs),1);
            xcoord = obj.geometry.radius*xtrans;
            ycoord = obj.geometry.radius*ytrans;
            zcoord = repmat(obj.xs,length(obj.thetas),1)';
            tdata.myX = reshape(xcoord,numel(xcoord), 1);
            tdata.myY = reshape(ycoord,numel(ycoord), 1);
            tdata.myZ = reshape(zcoord,numel(zcoord), 1);
            radialDisplacement = reshape(obj.WtotalF, numel(obj.WtotalF), 1);
            dispX = radialDisplacement.*xtrans;
            scaling = 1; 
            %scaling = 1000; % disp in mm
            tdata.rX = real(dispX)*scaling;
            tdata.iX = imag(dispX)*scaling;
            dispY = radialDisplacement.*ytrans;
            tdata.rY = real(dispY)*scaling;
            tdata.iY = imag(dispY)*scaling;
            velX = 1i*obj.force1.omega*dispX;
            tdata.rXv = real(velX)*scaling;
            tdata.iXv = imag(velX)*scaling;
            velY = 1i*obj.force1.omega*dispY;
            tdata.rYv = real(velY)*scaling;
            tdata.iYv = imag(velY)*scaling;
            tdata.rZ = 0*tdata.rY;
            tdata.iZ = 0*tdata.iY;
            tdata.rZv = 0*tdata.rZ;
            tdata.iZv = 0*tdata.rZ;
            obj.data = tdata;
        end
    end
    methods % getters
        function value = get.fileName(obj)
            date = datestr(floor(now));
            date = strrep(date,'-','_');
            value = strcat('cylinder_simulation_2_forces_', date, '_', num2str(obj.force1.frequency), '_Hz');
        end
        function value = get.xs(obj)
            value = linspace(0,obj.geometry.height,obj.mesh.longitude_divisions+1);
        end
        function value = get.thetas(obj)
            temp = linspace(0, 2*pi, obj.mesh.theta_divisions+1);
            value = temp(1:end-1);
        end
        function value = get.SpikeN(obj)
            dfreq = 1/max(obj.time);
            value = round(obj.force1.frequency/dfreq);
        end
        function value = get.xf1c(obj)
            value = pi*obj.force1.xStar/obj.geometry.height;
        end
        function value = get.xf2c(obj)
            value = pi*obj.force2.xStar/obj.geometry.height;
        end
        function value = get.C1(obj)
            if obj.force1.enabled
                M = obj.force1.magnitude;
            else
                M = 0;
            end
            rho = obj.material.density;
            thick = obj.geometry.thickness;
            r = obj.geometry.radius;
            l = obj.geometry.height;
            value = 2*M/(rho*thick*r*l*pi);
        end
        function value = get.C2(obj)
            if obj.force2.enabled
                M = obj.force2.magnitude;
            else 
                M = 0;
            end
            rho = obj.material.density;
            thick = obj.geometry.thickness;
            r = obj.geometry.radius;
            l = obj.geometry.height;
            value = 2*M/(rho*thick*r*l*pi);
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
        function value = get.xyzForce1(obj)
            F = obj.force1.location;
            x = obj.geometry.radius*cos(F(1));
            y = obj.geometry.radius*sin(F(1));
            z = F(2);
            value = [x, y, z];
        end
        function value = get.xyzForce2(obj)
            F = obj.force2.location;
            x = obj.geometry.radius*cos(F(1));
            y = obj.geometry.radius*sin(F(1));
            z = F(2);
            value = [x, y, z];
        end
        function value = get.trzForce1(obj)
            F = obj.force1.location;
            F(1) = wrapTo2Pi(F(1));
            value = [F(1), obj.geometry.radius, F(2)];
        end
        function value = get.trzForce2(obj)
            F = obj.force2.location;
            F(1) = wrapTo2Pi(F(1));
            value = [F(1), obj.geometry.radius, F(2)];
        end
        function value = get.sortedNFs(obj)
            value = sort(reshape(obj.naturalFrequencies,[numel(obj.nOmega),1]));
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
            X = data.myX;  Y = data.myY;   Z = data.myZ;
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