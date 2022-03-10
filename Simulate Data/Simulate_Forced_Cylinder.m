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
        WtotalV
        WtotalFV
        dwdx
        dwd0
        d2wdx2
        d2wdxd0
        d2wd02
        d3wdx3
        d3wdx2d0
        d3wdxd02
        d3wd03
        dwdxT
        dwd0T
        d2wdx2T
        d2wdxd0T
        d2wd02T
        d3wdx3T
        d3wdx2d0T
        d3wdxd02T
        d3wd03T
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
        function plot_flat(obj, type, isReal)
            if nargin< 3 
                isReal = true;
            end
            switch type
                case { 'Nt', 'Nl', 'Ntl', 'Mt', 'Ml', 'Mtl', 'Qt', 'Ql'}
                    obj.plot_flat_resultants(type, isReal);
                otherwise
                    obj.plot_flat_disp_or_derive(type, isReal);
            end            
        end
        function plot_flat_disp_or_derive(obj, deriv, isReal)
            if nargin < 3
                isReal = true;
            end
            value = obj.get_value_or_derivative(deriv);
            if ~isReal
                Values = imag(value);
            else
                Values = real(value);
            end
            [T, H] = meshgrid(obj.thetas, obj.xs);
            surf(T*obj.geometry.radius, H, -.0001*zeros(size(H)), Values, 'EdgeAlpha', 0);
            view(0,90)
            colorbar
            colormap jet
        end
        function plot_flat_resultants(obj, resultant, isReal)
            if nargin < 3
                isReal = true;
            end
            value = obj.calculate_resultant(resultant);
            if ~isReal
                Values = imag(value);
            else
                Values = real(value);
            end
            [T, H] = meshgrid(obj.thetas, obj.xs);
            surf(T*obj.geometry.radius, H, -.0001*zeros(size(H)), Values, 'EdgeAlpha', 0);
            view(0,90)
            colorbar
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
            % that is the radial dispalcement is a function of (long, theta,
            % t)
            % here the first set of indices represents longituinal mesh
            % point, 2nd theta point and 3rd time point
            % variables to hold properties that are constant in loops
            ph = obj.material.density*obj.geometry.thickness;
            L = obj.geometry.height; 
            Xs = obj.xs; % used to evaluate function at all longitudinal points
            Thetas = obj.thetas;  % evaluate at all theta points
            ts = obj.time;       % evaluate at all time points
            obj.WtotalD = []; % ensure current displacement cleared
            for frc = obj.forces()
                %for each force applyied to the cylinder
                %% pre-allocate all to zeros for each force
                W = zeros(length(Xs), length(Thetas), length(ts));
                dx = zeros(length(Xs), length(Thetas), length(ts));
                d0 = zeros(length(Xs), length(Thetas), length(ts));
                dx2 = zeros(length(Xs), length(Thetas), length(ts));
                dxd0 = zeros(length(Xs), length(Thetas), length(ts));
                d02 = zeros(length(Xs), length(Thetas), length(ts));
                dx3 = zeros(length(Xs), length(Thetas), length(ts));
                dx2d0 = zeros(length(Xs), length(Thetas), length(ts));
                dxd02 = zeros(length(Xs), length(Thetas), length(ts));
                d03 = zeros(length(Xs), length(Thetas), length(ts));
                %%
                omg_F = frc.omega; % get frequency of force (should be same for all forces)
                % for every longitudinal mode
                for m = 1:obj.mesh.longitude_modes
                    % terms that only depend on m-number of mode
                    M = m*pi/L;
                    S = sin(m*pi*frc.xStar/L);
                    %longitudinal part of function and derivatives
                    longitudinalFunction =  S * sin(M*Xs);
                    LD1 = S *  M   * cos(M*Xs);
                    LD2 = S * -M^2 * sin(M*Xs);
                    LD3 = S * -M^3 * cos(M*Xs);
                    % for every theta mode
                    for n = 0:obj.mesh.theta_modes-1
                        % mode contribution factor
                        Nmn = obj.calculate_Nk(m,n);
                        % complex Damping phase lag
                        phimn = obj.complexDampingPhaseLag(m,n+1);
                        % f(w) see sodel chapter 8 pg. 238 8.14.19
                        fwmn = obj.force_at_omega(frc, m, n);
                        % constants for multiplying
                        consts = frc.magnitude/(ph*Nmn*fwmn);
                        % theta portion of function and derivatives
                        cosTermVector = cos(n*(Thetas-frc.tStar));
                        cD1 =  -n  * sin(n*(Thetas-frc.tStar));
                        cD2 = -n^2 * cos(n*(Thetas-frc.tStar));
                        cD3 =  n^3 * sin(n*(Thetas-frc.tStar));
                        % calculate amplitude and derivative amplitudes
                        amplitude = (consts*longitudinalFunction')*cosTermVector;
                        Adx = consts*LD1'*cosTermVector;
                        Ad0 = consts*longitudinalFunction'*cD1;
                        Adx2 = consts*LD2'*cosTermVector;
                        Adxd0 = consts*LD1'*cD1;
                        Ad02 = consts*longitudinalFunction'*cD2;
                        Adx3 = consts*LD3'*cosTermVector;
                        Adx2d0 = consts*LD2'*cD1;
                        Adxd02 = consts*LD1'*cD2;
                        Ad03 =  consts*longitudinalFunction'*cD3;
                        % calculate time vector % phase is the combination
                        % of the phase shift and the complex phase lag
                        % for all time steps
                        timeVector(1,1,:) = complex(cos(omg_F*ts-phimn-frc.phase), sin(omg_F*ts-phimn-frc.phase));
                        % compute displacement over grid for all time
                        % steps and add to total displacement
                        W = W + amplitude.*timeVector;
                        dx = dx + Adx.*timeVector;
                        d0 = d0 + Ad0.*timeVector;
                        dx2 = dx2 + Adx2.*timeVector;
                        dxd0 = dxd0 + Adxd0.*timeVector;
                        d02 = d02 + Ad02.*timeVector;
                        dx3 = dx3 + Adx3.*timeVector;
                        dx2d0 = dx2d0 + Adx2d0.*timeVector;
                        dxd02 = dxd02 + Adxd02.*timeVector;
                        d03 = d03 + Ad03.*timeVector;
                    end
                end
                % add modal contribution to total deflection
                if isempty(obj.WtotalD)%if first time through there is nothing to add to
                    obj.WtotalD = W;
                    obj.dwdx = dx;
                    obj.dwd0 = d0;
                    obj.d2wdx2 = dx2;
                    obj.d2wdxd0 = dxd0;
                    obj.d2wd02 = d02;
                    obj.d3wdx3 = dx3;
                    obj.d3wdx2d0 = dx2d0;
                    obj.d3wdxd02 = dxd02;
                    obj.d3wd03 = d03;
                
                else % add to existing displacement
                    obj.WtotalD = obj.WtotalD + W;
                    obj.dwdx = obj.dwdx + dx;
                    obj.dwd0 = obj.dwd0 + d0;
                    obj.d2wdx2 = obj.d2wdx2 + dx2;
                    obj.d2wdxd0 = obj.d2wdxd0 + dxd0;
                    obj.d2wd02 = obj.d2wd02 + d02;
                    obj.d3wdx3 = obj.d3wdx3 + dx3;
                    obj.d3wdx2d0 = obj.d3wdx2d0 + dx2d0;
                    obj.d3wdxd02 = obj.d3wdxd02 + dxd02;
                    obj.d3wd03 = obj.d3wd03 + d03;
                end
            end
            % compute velocity from displacement by multiplying by
            % imaginary angular velocity
            obj.WtotalV = obj.WtotalD * 1i*obj.forces(1).omega;
        end
        function calculate_spectral_terms(obj)
            % this fuction calculates the FFT of the diplacement field to
            % get the complex values for displacement
            timeDim = 3;% time is the third dimension of the data
            % take the fft of each of the data sets 
            fullFFT = 1/(size(obj.WtotalD,3))*fft(obj.WtotalD, [], timeDim); % take the FFT over time dimension
            fullFFTVel = 1/(size(obj.WtotalV,3))*fft(obj.WtotalV, [], timeDim);
            dx = 1/(size(obj.WtotalD,3))*fft(obj.dwdx, [], timeDim);
            d0 = 1/(size(obj.WtotalD,3))*fft(obj.dwd0, [], timeDim);
            dx2 = 1/(size(obj.WtotalD,3))*fft(obj.d2wdx2, [], timeDim);
            dxd0 = 1/(size(obj.WtotalD,3))*fft(obj.d2wdxd0, [], timeDim);
            d02 = 1/(size(obj.WtotalD,3))*fft(obj.d2wd02, [], timeDim);
            dx3 = 1/(size(obj.WtotalD,3))*fft(obj.d3wdx3, [], timeDim);
            dx2d0 = 1/(size(obj.WtotalD,3))*fft(obj.d3wdx2d0, [], timeDim);
            dxd02 = 1/(size(obj.WtotalD,3))*fft(obj.d3wdxd02, [], timeDim);
            d03 = 1/(size(obj.WtotalD,3))*fft(obj.d3wd03, [], timeDim);
            %identify FFT peak index
            [~, indx] = max(squeeze(fullFFT(2,2,:))); %% first row is simply supported so it will always be zero so check inside the mesh
            % get the FFT at the peak value for each
            obj.WtotalF = fullFFT(:,:,indx);
            obj.WtotalFV = fullFFTVel(:,:,indx);
            obj.dwdxT = dx(:,:,indx);
            obj.dwd0T = d0(:,:,indx);
            obj.d2wdx2T = dx2(:,:,indx);
            obj.d2wdxd0T = dxd0(:,:,indx);
            obj.d2wd02T = d02(:,:,indx);
            obj.d3wdx3T = dx3(:,:,indx);
            obj.d3wdx2d0T = dx2d0(:,:,indx);
            obj.d3wdxd02T = dxd02(:,:,indx);
            obj.d3wd03T = d03(:,:,indx);
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
            %Soedel eqs 5.5.67 - 5.5.72
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
            % Soedel equation 5.5.85
            AoCi = -(k13*(phw2 - k22)- k12*k23)/denom;
            % Soedel equation 5.5.86
            BoCi = -(k23*(phw2 - k11)- k12*k13)/denom;
        end
        function fw = force_at_omega(obj, frc, m, n)
            % this equation calculates the equivalent force for calcuating
            % the contribution of a harmonic load at each mode shape
            wmn = obj.nOmega(m,n+1); % mode shape frequency
            w = frc.omega; % forcing frequency
            %damping ratio
            zeta = obj.lambda/(2*obj.material.density*obj.geometry.thickness*obj.nOmega(m,n+1));
            % equivalent force 
            fw = wmn^2 * sqrt((1-(w/wmn)^2)^2 + 4*zeta^2*(w/wmn)^2);
        end
        function Nk = calculate_Nk(obj, m, n)
            % calculate the contribution factor from the Mode ratio
            [AoCmn, BoCmn] = obj.calculate_mode_ratio(m,n);
            L = obj.geometry.height;
            a = obj.geometry.radius;
            % Soedel equation 8.14.30
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
            %             radialV = 1i*obj.forces(1).omega * radialDisplacement;
            radialVelocity = reshape(obj.WtotalFV, numel(obj.WtotalFV), 1);
            dispX = radialDisplacement.*xtrans;
            scaling = 1;
            tdata.rX = real(dispX)*scaling;
            tdata.iX = imag(dispX)*scaling;
            dispY = radialDisplacement.*ytrans;
            tdata.rY = real(dispY)*scaling;
            tdata.iY = imag(dispY)*scaling;
            velX = radialVelocity.*xtrans;
            tdata.rXv = real(velX)*scaling;
            tdata.iXv = imag(velX)*scaling;
            velY = radialVelocity.*ytrans;
            tdata.rYv = real(velY)*scaling;
            tdata.iYv = imag(velY)*scaling;
            tdata.rZ = 0*tdata.rY;
            tdata.iZ = 0*tdata.iY;
            tdata.rZv = 0*tdata.rZ;
            tdata.iZv = 0*tdata.rZ;
            obj.data = tdata;
        end
        function value = get_value_or_derivative(obj, deriv)
            % get the values for the dispalcement, velocity, or up to 3rd order spatial
            % derivatives of the displacement
            switch deriv
                case {'dx'}
                    value = obj.dwdxT;
                case {'d0'}
                    value = obj.dwd0T;
                case {'dx2'}
                    value = obj.d2wdx2T;
                case {'dxd0'}
                    value = obj.d2wdxd0T;
                case {'d02'}
                    value = obj.d2wd02T;
                case {'dx3'}
                    value = obj.d3wdx3T;
                case {'dx2d0'}
                    value = obj.d3wdx2d0T;
                case {'dxd02'}
                    value = obj.d3wdxd02T;
                case {'d03'}
                    value = obj.d3wd03T;
                case {'rdot'}
                    value = obj.WtotalFV;
                otherwise
                    value = obj.WtotalF;
            end
        end
        function value = calculate_resultant(obj, resultant)
            % calculate a specified resultant using the fact that there is
            % only radial displacement
            r = obj.geometry.radius;
            nu = obj.material.poisson;
            switch resultant
                case {'Nl'}
                    value = obj.D * (nu/r*(obj.get_value_or_derivative('w'))) - obj.K/r*obj.get_value_or_derivative('dx2');
                case {'Nt'}
                    value = obj.D * (obj.get_value_or_derivative('w')/r) + obj.K/r^3 *(obj.get_value_or_derivative('w') + obj.get_value_or_derivative('d02'));
                case {'Nlt'}
                    value = obj.K*(1-nu)/2*(-obj.get_value_or_derivative('dxd0'))/r^2;
                case {'Ntl'}
                    value = obj.K*(1-nu)/2*(-obj.get_value_or_derivative('dxd0'))/r^2;
                case {'Mt'}
                    value = obj.K *(obj.get_value_or_derivative('w')/r^2 +obj.get_value_or_derivative('d02')/r^2 + nu*obj.get_value_or_derivative('dx2'));
                case {'Ml'}
                    value = obj.K*(nu/r^2*obj.get_value_or_derivative('d02') + obj.get_value_or_derivative('dx2'));
                case {'Mlt'}
                    value = obj.K*(1-nu)/r*(obj.get_value_or_derivative('dxd0'));
                case {'Mtl'}
                    value = obj.K*(1-nu)/r*(obj.get_value_or_derivative('dxd0'));
                case {'Ql'}
                    dMtl_dth_a =  obj.K*(1-nu)*(obj.get_value_or_derivative('dxd02')/r^2);
                    dMl_dz =  obj.K*(nu/r^2*obj.get_value_or_derivative('dxd02') + obj.get_value_or_derivative('dx3'));
                    value = dMl_dz + dMtl_dth_a;
                case {'Qt'}
                    dMt_dth_a = obj.K *(obj.get_value_or_derivative('d0')/r^3 +obj.get_value_or_derivative('d03')/r^3 + nu/r*obj.get_value_or_derivative('dx2d0'));
                    dMlt_dz = obj.K*(1-nu)/r*(obj.get_value_or_derivative('dx2d0'));
                    value = dMt_dth_a + dMlt_dz;
                otherwise
                    obj.get_value_or_derivative
                    value = [];
            end
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
            obj.forces(1) = sinusoidalForce( [.025, 90*pi/180], -10, frequency, 0);
            obj.forces(2) = sinusoidalForce( [.135, 270*pi/180], -10, frequency, 179*pi/180);
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
            dataFormat = strcat('%d',repmat('\t%0.8g', 1, 9), '\n');
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