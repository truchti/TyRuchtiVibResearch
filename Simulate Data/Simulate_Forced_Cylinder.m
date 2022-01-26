classdef Simulate_Forced_Cylinder< handle
    % This class simulates the displacement field of a cylinder subjected
    % to 2 point loads both driving at the same frequency but they can be
    % phased independently. 
    properties
        material
        forces = sinusoidalForce.empty()
        geometry
        mesh
        saveFolder = 'C:\Users\ME\Desktop\RandomData\';
        timesteps =20;
        nOmega= [];
        complexDampingPhaseLag
    end
    properties(Dependent)
        fileName
        naturalFrequencies
    end
    properties (Hidden = true, Dependent = true)
        xs
        thetas
        lambda
        time
        sortedNFs
    end
    properties (Hidden = true)
        data
        WtotalD = [];
        WtotalF
        K
        D
        PFCylProps
        useSoedel = false;
    end
    methods
        function obj = Simulate_Forced_Cylinder(geometry, forces, material, mesh)
            %% creates the object that simulates the displacement of a cylinder in response to applied forces
            if nargin < 1
                obj.create_default();
            elseif nargin < 2
                obj.create_default_at_frequency(geometry);
            else
                obj.geometry = geometry;
                obj.forces = forces;
                obj.material = material;
                obj.mesh = mesh;
            end
        end
        function simulate_and_save(obj)
            obj.simulate_data();
            obj.save_data();
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
            s = surf(x,y,z,loopedValues(:,:,1));
            axis equal
            s.EdgeAlpha = 0.0;
            s.FaceColor = 'interp';
            colorbar
            clim = max(max(max(values)));
            scl = 1;
            caxis([-scl*clim, scl*clim])
            sz = size(values,3);
            view(10, 25)
            for i = 1:loops*sz
                s.CData = loopedValues(:,:,mod(i,sz)+1);
%                 s2.CData = iloopedValues(:,:,mod(i, sz)+1);
                pause(0.1);
            end
        end             
        function simulate_data(obj)
            obj.calculate_mode_frequencies();
            obj.calculate_complex_damping_phase_lag();
            obj.calculate_spatial_temporal_displacement();
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
        function PFCP = export_power_flow_cylinder_properties(obj)
            if isempty(obj.PFCylProps)
                obj.create_PFCylProps()
            end
            PFCP = obj.PFCylProps;
        end
    end
    methods (Access = private)
        %% calculations based of Sodel Eqs
        function calculate_spatial_temporal_displacement(obj)
            %% calculates the displacement of a cylinder due to a harmonic force
            % variables to hold properties that are constant in loops
            ph = obj.material.density*obj.geometry.thickness;
            ts = obj.time;  L = obj.geometry.height;    Thetas = obj.thetas;    Xs = obj.xs;
            obj.WtotalD = [];
            for frc = obj.forces()
                %for each force applyied to the cylinder 
                W = zeros(length(Xs), length(Thetas), length(ts));
                omg_F = frc.omega; % get frequency of force (should be same for all forces)
                for m = 1:obj.mesh.longitude_modes
                    % components that only depend on m-number of mode
                    sinTermVector = sin(m*pi*frc.xStar/L) * sin(m*pi*Xs/L);
                    for n = 0:obj.mesh.theta_modes-1
                        % mode contribution factor
                        Nmn = obj.calculate_Nk(m,n);
                        % complex Damping phase lag 
                        phimn = obj.complexDampingPhaseLag(m,n+1);
                        % f(w) see sodel chapter 8 pg. 238 8.14.19
                        fwmn = obj.force_at_omega(frc, m, n);
                        % constants for multiplying 
                        consts = frc.magnitude/(ph*Nmn*fwmn);
                        % terms associated with n mode number
                        cosTermVector = cos(n*(Thetas-frc.tStar));
                        % calculate amplitude
                        amplitude = (consts*sinTermVector')*cosTermVector;
                        % calculate time vector
                        timeVector(1,1,:) = complex(cos(omg_F*ts-phimn-frc.phase), sin(omg_F*ts-phimn-frc.phase));
                        % compute displacement over grid for all time
                        % stepsand add to total displacement
                        W = W + amplitude.*timeVector;
                    end
                end
                if isempty(obj.WtotalD)
                    obj.WtotalD = W;
                else
                    obj.WtotalD = obj.WtotalD + W;
                end
            end
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
        function calculate_mode_frequencies(obj)
            % This function calculates the natural frequencies of the modes
            % of the cylinder it also calculates the phase associated
            % with the each mode because of the damping
            if ~obj.useSoedel
                obj.rao_nat_freq(obj.mesh.longitude_modes, obj.mesh.theta_modes);
            else
                obj.soedel_nat_freq(obj.mesh.longitude_modes, obj.mesh.theta_modes);
            end            
        end
        function calculate_complex_damping_phase_lag(obj)
            %calcuated based off of eq. sodel 8.5.7 and 8.8.28 
            % zeta_k = lambda/(2*rho*h*omega_k) equation 8.3.3
            zpart = obj.lambda/(2*obj.material.density*obj.geometry.thickness); % damping Coeff constant part
            fomg1 = obj.forces(1).omega; % force Freq
            for m = 1:obj.mesh.longitude_modes
                for n = 0:obj.mesh.theta_modes
                    zeta = zpart/obj.nOmega(m,n+1); % soedel eq 8.3.3 (pg. 212)
                    omegRatio = fomg1/obj.nOmega(m,n+1);
                    obj.complexDampingPhaseLag(m,n+1) = atan((2*zeta*(omegRatio))/(1-(omegRatio)^2)); % sodel eq 8.5.6 pg 215
                end
            end
        end
        function [k11, k12, k13, k22, k23, k33] = calculate_matrix_ks(obj, m, n)
            a = obj.geometry.radius;
            L = obj.geometry.height;
            KK = obj.K;
            DD = obj.D;
            mu = obj.material.poisson;
            noa = n/a;
            mpL = m*pi/L;
            noa2 = noa^2;
            mpL2 = mpL^2;
            %eqs 5.5.67 - 5.5.72
            k11 = KK * (mpL2 + (1-mu)/2 * noa2);
            k12 = KK * (1+mu)/2 * mpL * noa;
            k13 = mu*KK/a * (mpL);
            k22 = (KK + DD/a^2) * ((1-mu)/2*mpL2 + noa2);
            k23 = -(KK/a)*noa - (DD/a*noa)*(mpL2 + noa2);
            k33 = DD*(mpL2 + noa2)^2 + KK/a^2;
        end   
        function [AoCi, BoCi] = calculate_mode_ratio(obj, m, n)
            p = obj.material.density;
            h = obj.geometry.thickness;
            phw2 = p*h*(obj.nOmega(m,n+1))^2;
            [k11, k12, k13, k22, k23, ~] = obj.calculate_matrix_ks(m,n);
            denom = ((phw2-k11)*(phw2-k22)-k12^2);
            %equation5.5.85
            AoCi = -(k13*(phw2 - k22)- k12*k23)/denom;
            %equation 5.5.86
            BoCi = -(k23*(phw2 - k11)- k12*k13)/denom;
        end
        function fw = force_at_omega(obj, frc, m, n)
            wmn = obj.nOmega(m,n+1);
            w = frc.omega;
            zeta = obj.lambda/(2*obj.material.density*obj.geometry.thickness*obj.nOmega(m,n+1));
            fw = wmn^2 * sqrt((1-(w/wmn)^2)^2 + 4*zeta^2*(w/wmn)^2);
        end
        function Nk = calculate_Nk(obj, m, n)
            [AoCmn, BoCmn] = obj.calculate_mode_ratio(m,n);
            L = obj.geometry.height;
            a = obj.geometry.radius;
            % equation 8.14.30
            if n == 0
                Nk = (AoCmn^2 +1) * L*a*pi;
            else
                Nk = (AoCmn^2 + BoCmn^2 + 1) * L*a*pi/2;
            end            
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
            velX = 1i*obj.forces(1).omega*dispX;
            tdata.rXv = real(velX)*scaling;
            tdata.iXv = imag(velX)*scaling;
            velY = 1i*obj.forces(1).omega*dispY;
            tdata.rYv = real(velY)*scaling;
            tdata.iYv = imag(velY)*scaling;
            tdata.rZ = 0*tdata.rY;
            tdata.iZ = 0*tdata.iY;
            tdata.rZv = 0*tdata.rZ;
            tdata.iZv = 0*tdata.rZ;
            obj.data = tdata;
        end
        %% data file functions
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
            obj.write_header(fid, type, obj.forces(1).frequency)
            obj.write_data(fid, obj.data, type);
            fclose(fid);
            cd(old);
        end
        function create_PFCylProps(obj)
            r = obj.geometry.radius;
            t = obj.geometry.thickness;
            h = obj.geometry.height;
            E = obj.material.E;
            p = obj.material.poisson;
            fLs =[];
            for i = 1:length(obj.forces)
                fLs = [fLs; obj.forces(i).location]; %#ok<AGROW>
            end
            obj.PFCylProps = PFCylProperties(r, t, h, E, p, fLs);
        end
        %% natural freq functions
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
            Kk = E*h/(1-mu^2); % soedel eq 2.5.10; Membrane stiffness
            Dd = E*h^3/(12*(1-mu^2)); % soedel eq 2.5.16 Bending Stiffness            
            w_squared = zeros(3,mModes, nModes+1);
            for m = 1:mModes
                for n = 0:nModes
                    % eqs 5.5.67-72 k matrix
                    k11 = Kk*((m*pi/L)^2 + (1-mu)/2*(n/a)^2);
                    k12 = Kk*(1+mu)/2*(m*pi/L)*n/a;
                    k13 = mu*Kk/a*(m*pi/L);
                    k22 = (Kk+ Dd/a^2)*(((1-mu)/2)*(m*pi/L)^2+ (n/a)^2);
                    k23 = -Kk*n/a^2-((Dd*n)/a^2)*((m*pi/L)^2+(n/a)^2);
                    k33 = Dd*((m*pi/L)^2 + (n/a)^2)^2+Kk/a^2;
                    
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
        %% setup functions
        function create_default(obj)
            obj.geometry = cylinderGeometry();
            obj.forces(1) = sinusoidalForce(1);
            obj.forces(2) = sinusoidalForce(2);
            obj.material = simMaterial();
            obj.mesh = meshDetails();
            obj.timesteps = 20;
            obj.saveFolder = "C:\Users\ME\Desktop\Simulated Data\";
        end
        function create_default_at_frequency(obj, frequency)
            obj.geometry = cylinderGeometry();
            obj.forces(1) = sinusoidalForce( [.025, 0*pi/180], 10, frequency, 0);
            obj.forces(2) = sinusoidalForce( [.135, 180*pi/180], 10, frequency, 179*pi/180);
            obj.material = simMaterial();
            obj.mesh = meshDetails();
            obj.timesteps = 20;
            obj.saveFolder = "C:\Users\ME\Desktop\Simulated Data\";
        end
    end
    methods % getters
        function value = get.fileName(obj)
            date = datestr(floor(now));
            date = strrep(date,'-','_');
            value = strcat('cylinder_simulation_2_forces_', date, '_', num2str(obj.forces(1).frequency), '_Hz');
        end
        function value = get.xs(obj)
            value = linspace(0,obj.geometry.height,obj.mesh.longitude_divisions+1);
        end
        function value = get.thetas(obj)
            temp = linspace(0, 2*pi, obj.mesh.theta_divisions+1);
            value = temp(1:end-1);
        end
        function value = get.time(obj)
            t =  linspace(0, 1/obj.forces(1).frequency, obj.timesteps+1);
            value = t(1:end-1);
        end
        function value = get.lambda(obj)
            value = obj.material.density*obj.geometry.thickness*obj.material.eta*obj.forces(1).omega;
        end
        function value = get.naturalFrequencies(obj)
            value = obj.nOmega/(2*pi);
        end
        function value = get.sortedNFs(obj)
            value = sort(reshape(obj.naturalFrequencies,[numel(obj.nOmega),1]));
            value = value(1:20);
        end
        function value = get.K(obj)
            h = obj.geometry.thickness;
            E = obj.material.E;
            mu = obj.material.poisson;
            value = E*h/(1-mu^2);
        end
        function value = get.D(obj)
            h = obj.geometry.thickness;
            E = obj.material.E;
            mu = obj.material.poisson;
            value = E*h^3/(12*(1-mu^2));
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