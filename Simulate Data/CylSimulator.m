classdef CylSimulator< handle
    % coordinates should be (R, Theta, Z)
    % This class simulates the displacement field of a cylinder subjected
    % to 2 point loads both driving at the same frequency but they can be
    % phased independently.
    properties
        % for all properties t represents theta and x represents
        % longitudinal
        material
        forces
        geometry
        mesh
        saveFolder = 'C:\Users\tysum\Documents\SimulatedData';
        phaseCorrection = true
    end
    properties(Dependent)
        fileName
        thetas
        longs
        naturalFrequencies
    end
    properties (Hidden = true)
        timesteps =50;
        nOmega= [];
        complexDampingPhaseLag = [];
        useSoedel = false;
        data % structure to hold rectangular coordinate transforms
        %% deriv vars
        w
        wdot
        dwdL, dwd0
        d2wdL2, d2wdL0, d2wd02
        d3wdL3, d3wdL20, d3wdL02, d3wd03
        d2wdLt, d2wd0t
        %% test vars
        denomArray
        ampArray
    end
    properties (Access = private)
          %% per force vars
          wforce
          dw
          dL, d0
          dL2, dL0, d02
          dL3, dL20, dL02, d03
          d0t, dLt
          %% names of properties
          names = {'w', 'wdot', 'dwdL', 'dwd0', 'd2wdL2', 'd2wdL0', 'd2wd02', 'd3wdL3', 'd3wdL20', 'd3wdL02', 'd3wd03', 'd2wdLt', 'd2wd0t'};
          perForceNames = {'wforce', 'dw', 'dL', 'd0', 'dL2', 'dL0', 'd02', 'dL3', 'dL20', 'dL02', 'd03', 'd0t', 'dLt'};
          %% animation vars
          loopedValues
          currentLoopStep
          surfX, surfY, surfZ
          animationTimer
          aniPlotHandle
          surfHandle
    end
    %% dimensional derivatives
    properties(Hidden = true, Dependent = true)
        % length coors
        dwdx, dwdy
        d2wdx2, d2wdxy, d2wdy2
        d3wdx3, d3wdx2y, d3wdxy2, d3wdy3
        d2wdxt, d2wdyt
    end
    methods
        function value = get.dwdx(obj)
            value = obj.dwdL;
        end
        function value = get.dwdy(obj)
            r = obj.geometry.radius;
            value = 1/r* obj.dwd0;
        end
        function value = get.d2wdx2(obj)
            value = obj.d2wdL2;
        end
        function value = get.d2wdxy(obj)
            r = obj.geometry.radius;
            value = 1/r* obj.d2wdL0;
        end
        function value = get.d2wdy2(obj)
            r = obj.geometry.radius;
            value = 1/r^2* obj.d2wd02;
        end
        function value = get.d3wdx3(obj)
            value = obj.d3wdL3;
        end
        function value = get.d3wdx2y(obj)
            r = obj.geometry.radius;
            value = 1/r* obj.d3wdL20;
        end
        function value = get.d3wdxy2(obj)
            r = obj.geometry.radius;
            value = 1/r^2* obj.d3wdL02;
        end
        function value = get.d3wdy3(obj)
            r = obj.geometry.radius;
            value = 1/r^3* obj.d3wd03;
        end
        function value = get.d2wdxt(obj)
            value = obj.d2wdLt;
        end
        function value = get.d2wdyt(obj)
            r = obj.geometry.radius;
            value = 1/r * obj.d2wd0t;
        end
    end
    %% general methods
    methods
        function obj = CylSimulator(geometry, forces, material, mesh)
            %% creates the object that simulates the displacement of a cylinder to 2 forces at a gived frequency
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
        %% run simulations
        function simulate_and_save(obj)
            obj.simulate_data();
            obj.save_data();
        end
        function simulate_with_specified_mode_frequencies(obj)
            obj.calculate_complex_damping_phase_lag();
            obj.calculate_displacements();
            obj.convert_to_rectangular_coordinates();
        end
        function simulate_data(obj)
            obj.calculate_mode_frequencies();
            obj.calculate_complex_damping_phase_lag();
            obj.calculate_displacements();
            obj.convert_to_rectangular_coordinates();
        end
        function [fRao, fSoe] = compare_both_frequency_calcs(obj)
            omgRao = obj.rao_nat_freq(obj.mesh.longitude_modes, obj.mesh.theta_modes);
            omgSoe = obj.soedel_nat_freq(obj.mesh.longitude_modes, obj.mesh.theta_modes);
            fRao = omgRao/(2*pi);
            fSoe = omgSoe/(2*pi);
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
        %% force adjusting
        function disable_force(obj, force_number)
            % turn off a force from being applied
            if force_number <= length(obj.forces)
                obj.forces(force_number).enabled = false;
            end
        end
        function enable_force(obj, force_number)
            % turn on a force so that it affects displacement
            if force_number <= length(obj.forces)
                obj.forces(force_number).enabled = true;
            end
        end
        %% plotting
        function phandle = plot_cyl_component(obj, type, isReal)
            if nargin< 3
                isReal = true;
            end
            [Th, Lng, value] = obj.get_plot_values(type, isReal);
            Th = [Th, Th(:,1)];
            Lng = [Lng, Lng(:,1)];
            value = [value, value(:,1)];
            X = cos(Th)*obj.geometry.radius;
            Y = sin(Th)*obj.geometry.radius;
            Z = Lng;
            phandle = surf(X,Y,Z,value, 'EdgeAlpha', .25);
%             axis equal
%             colormap jet
            colorbar
        end
        function plot_RI_cyl_component(obj, type)
            figure('units','normalized','outerposition',[0.05 0.05 .9 .9])
            s1 = subplot(1,2,1);
            obj.plot_cyl_component(type, true);
            shading interp
            s1.Title.String = [s1.Title.String ' Real'];
            s2 = subplot(1,2,2);
            obj.plot_cyl_component(type, false);
            shading interp
            s2.Title.String = [s2.Title.String ' Imag'];
            obj.link_subplots(s1, s2);
        end
        function pHandle = plot_component(obj, type, isReal, truScale)
            if nargin < 3 
                isReal = true;
            end
            if nargin< 4
                truScale = true;
            end
            [Th, Lng, values] = obj.get_plot_values(type, isReal);
            if truScale
                Th = Th*obj.geometry.radius;
            end
            pHandle = surf(Th, Lng, values, 'EdgeAlpha', 0);
            Z = max(max(values));
%             W = min(min(values));
            obj.plot_forces(Z);
            xlabel('Tangential'); ylabel('Height'); title(type);
            colormap jet; colorbar
            axis equal; view(0, 90)
         end
        function vMax = plot_RI_component(obj, type, truScale)
            if nargin < 3
                truScale = true;
            end
%             figure('units','normalized','outerposition',[0.05 0.05 .5 .35])
            s1 = subplot(2,1,1);
            obj.plot_component(type, true, truScale);
            s1.Title.String = [s1.Title.String ' Real'];
            s2 = subplot(2,1,2);
            obj.plot_component(type, false, truScale);
            s2.Title.String = [s2.Title.String ' Imag'];
            obj.link_subplots(s1, s2);
            vMax = s2.CLim;
        end
        function plot_power_flow(obj, type, truScale)  
            if nargin < 2
                type = '';
            end
            if nargin < 3
                truScale = true;
            end
            [X,Y,U,V] = obj.get_quiver_values(type);
            if truScale
                X = X*obj.geometry.radius;
            end
            obj.plot_shadow_displacement(X,Y)
            hold on;
            quiver(X,Y,U,V);
            obj.plot_pf_forces()
            hold off;
        end
        function plot_dim_power_flow(obj, type)
            [Cir, Len] = meshgrid(obj.geometry.radius*obj.thetas, obj.longs);
            if nargin < 2
                [qLen,qCir] = obj.calculate_dimensional_power_flow;
            else
                [qLen,qCir] = obj.calculate_power_flow_part(type);
            end  
            quiver(Cir, Len, qCir, qLen);
            obj.plot_pf_forces()
        end
        function plot_shadow_displacement(obj,x,y)
            z = abs(obj.w);
            surf(x,y,zeros(size(z)),z, 'FaceAlpha', .1, 'EdgeAlpha', 0)
        end
        function plot_power_contours(obj)
            [X,Y,U,V] = obj.get_quiver_values('');
            X = X*obj.geometry.radius;
            Z = sqrt(real(U).^2 + real(V).^2);
            contour(X,Y,Z)            
        end
        function maxis = calculate_maximums(obj, pt) 
            mTypes = {'w', 'wdot', 'ML', 'M0', 'ML0', 'NL', 'N0', 'NL0', 'QL', 'Q0', 'B0', 'BL'};
            for i = 1:length(mTypes)
                if i == 2 || i > 10
                    isReal = false;
                else
                    isReal = true;
                end
                [x, y, values] = obj.get_plot_values(mTypes{i}, isReal);
                if nargin >1 
                    
                end
                maxis.(mTypes{i}) = max(values, [], 'all');
            end                
        end
        %% animation
        function animate(obj)
            obj.cyl_disp_animation_loop()
            if isempty(obj.animationTimer)
                obj.setupTimer;
            end
            start(obj.animationTimer)
        end
        function animate_flat(obj, duration)
            obj.flat_disp_animation_loop()
            if nargin>1 
                obj.setupTimer(duration);
            elseif isempty(obj.animationTimer) 
                obj.setupTimer();
            end
            start(obj.animationTimer)
        end
        function value = sortedNFs(obj, num)
            if nargin < 2
                num = Inf;
            end
            value = sort(reshape(obj.naturalFrequencies,[numel(obj.nOmega),1]));
            if length(value) > num
                value = value(1:num);
            end
        end
    end
    methods (Access = private)
        %% Natural frequency 
        function calculate_mode_frequencies(obj, useSoedel)
            % This function calculates the natural frequencies of the modes
            % of the cylinder it also calculates the phase lag associated
            % with the each mode because of the material damping
            if nargin < 2; useSoedel = true; end
            if  ~useSoedel
                obj.rao_nat_freq(obj.mesh.longitude_modes, obj.mesh.theta_modes);
            else
                obj.soedel_nat_freq(obj.mesh.longitude_modes, obj.mesh.theta_modes);
            end
        end
        function calculate_complex_damping_phase_lag(obj)
            eta = obj.material.eta;
%             zpart = obj.lambda/(2*obj.material.density*obj.geometry.thickness); % damping Coeff
            fomg1 = obj.forces(1).omega; % force Freq
            %% find closest mode freq to forcing freq
            diffMatrix = abs( obj.nOmega - fomg1);
            minDiff = min(min(diffMatrix));
            indx = find(minDiff == diffMatrix);
            % if force is really close (<0.1% off) to natural frequency shift phase so
            % that that mode has fully real displacement
            if obj.phaseCorrection && (minDiff/obj.nOmega(indx) < .001) 
                phaseAdjust = atan2(eta,(1-(fomg1/obj.nOmega(indx))^2));
            else
                phaseAdjust = 0;
            end
            for m = 1:obj.mesh.longitude_modes
                for n = 0:obj.mesh.theta_modes-1
%                     zeta = zpart/obj.nOmega(m,n+1); % soedel eq 8.3.3 (pg. 212)
%                     zeta = eta*obj.nOmega(m,n+1)/(2*fomg1); % soedel eq 14.3.18
                    % using relationship on lambda and zeta the numerator
                    % 2(zeta)(w/wk) equals eta
                    if obj.mesh.longitude_modes == 1
                        omegRatio = fomg1/obj.nOmega(n+1);
                    else
                        omegRatio = fomg1/obj.nOmega(m,n+1); 
                    end
                    obj.complexDampingPhaseLag(m,n+1) = atan2(eta,(1-(omegRatio)^2))-phaseAdjust; % sodel eq 8.5.6 pg 215
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
            mu = obj.material.poisson;
            rho = obj.material.density;
            h = obj.geometry.thickness;
            a = obj.geometry.radius;
            L = obj.geometry.height;
            % stiffnesses
            K = obj.calc_D; %Soedel eq 2.5.10 for K is the same as Flugge for D Flugge 5.8a
            D = obj.calc_K; % soedel eq 2.5.16 same as Flugge for K see flugge 5.8b
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
            omega3 = w_squared.^(1/2);
            if nargout == 0
                obj.nOmega = squeeze(omega3(1,:,:));
            else
                omgs = squeeze(omega3(1,:,:));
            end
        end
        %% calculate_displacements
        function calculate_displacements(obj)
            % if one of the forces is disabled only calculate the
            % displacement from that force
            % see Soedel Pg 225 equation number 8.8.29
%             timeV = obj.calc_time_vector();
            Ts = obj.thetas;
            Ls = obj.longs;
            obj.clean_all_simulated_variables(length(Ls), length(Ts));
%             obj.w = zeros(length(Ls), length(Ts)); % clear all displacement info
            for i = 1:length(obj.forces) % for each force on the cylinder
                if ~obj.forces(i).enabled % if force is disabled skip it
                    continue;
                end
                frc = obj.forces(i);
                obj.clean_per_force_variables(length(Ls), length(Ts));
                tstar = frc.tStar; % theta location of Force 
                % need to shift so that the main mode frequency is at
                % relative 0 phase shift
                Phi = frc.phase + obj.complexDampingPhaseLag;
                phaseMulti = cos(Phi) + 1i*sin(Phi);
                lStar = frc.xStar; %  longitudinal location of force
                L = obj.geometry.height;
                % longitude by theta mode
%                 obj.denomArray = zeros(obj.mesh.longitude_modes, obj.mesh.theta_modes);
                obj.ampArray = zeros(obj.mesh.longitude_modes, obj.mesh.theta_modes);
                for M = 1:obj.mesh.longitude_modes
                    s = M*pi/L; % constant inside of sin terms
                    CM = sin(s*lStar); % m mode force position multiplier
                    % m mode grid evaluations
                    gL = sin(s*Ls); % m mode grid evaluations
                    dgL = s * cos(s*Ls); % derivatives
                    d2gL = s^2 * -sin(s*Ls);
                    d3gL = s^3 * -cos(s*Ls);
                    for N = 0:obj.mesh.theta_modes-1
                        %% n mode grid evaluations
                        if obj.nOmega(M,N+1) == 0 %% if I have zeroed out this mode
                            continue;
                        end
                        fO = cos(N*(Ts-tstar)); 
                        dfO = N * -sin(N*(Ts-tstar));
                        d2fO = N^2 * -cos(N*(Ts-tstar));
                        d3fO = N^3 * sin(N*(Ts-tstar));
                        %% mode constants
                        denom = obj.calculate_denominator(frc,M,N);
                        obj.denomArray(M,N+1) = denom;
                        if numel(phaseMulti) > 1
                            Consts = CM/denom*phaseMulti(M,N+1);
                        else
                            Consts = CM/denom*phaseMulti;
                        end
                        obj.ampArray(M,N+1) = Consts;
                        %% compute various derivs
                        obj.wforce = obj.wforce + (Consts * gL')*fO;
                        obj.dw = obj.dw + 1i*frc.omega * (Consts*gL')*fO;
                        obj.dL = obj.dL + (Consts * dgL')*fO;
                        obj.d0 = obj.d0 + (Consts * gL')*dfO;
                        obj.dL2 = obj.dL2 + (Consts * d2gL')*fO;
                        obj.dL0 = obj.dL0 + (Consts * dgL')*dfO;
                        obj.d02 = obj.d02 + (Consts * gL')*d2fO;
                        obj.dL3 = obj.dL3 + (Consts * d3gL')*fO;
                        obj.dL20 = obj.dL20 + (Consts * d2gL')*dfO;
                        obj.dL02 = obj.dL02 + (Consts * dgL')*d2fO;
                        obj.d03 = obj.d03 + (Consts * gL')*d3fO;
                        obj.d0t = obj.d0t + 1i*frc.omega * (Consts*gL')*dfO;
                        obj.dLt = obj.dLt + 1i*frc.omega * (Consts*dgL')*fO;
                    end
                end
                % force magnitude constant
                C = obj.calculate_C_from_force(frc);
                % add current force contribution to displacement
                obj.add_scale_force_contribution_to_total(C);
            end
            obj.clean_per_force_variables(0,0);
        end            
        %% need to be verified still
        function values = calculate_resultant(obj, type)
            switch type
                case 'M0' 
                    [values, ~,~,~] = obj.calculate_moment_resultants();
                case 'ML' 
                    [~, values, ~,~] = obj.calculate_moment_resultants();
                case 'ML0' 
                    [~,~,values, ~] = obj.calculate_moment_resultants();
                case 'M0L' 
                     [ ~,~,~, values]= obj.calculate_moment_resultants();
                case 'N0'
                     [ values, ~,~,~ ]= obj.calculate_normal_resultants();
                case 'NL' 
                    [ ~,values,~,~] = obj.calculate_normal_resultants();
                case 'N0L' 
                    [ ~,~,values,~] = obj.calculate_normal_resultants();
                case 'NL0' 
                    [~,~,~,values] = obj.calculate_normal_resultants();
                case 'Q0' 
                    [values, ~] = obj.calculate_shear_resultants();
                case 'QL'
                    [~, values] = obj.calculate_shear_resultants();
            end
        end
        function [N0, NL, N0L, NL0] = calculate_normal_resultants(obj)
            D = obj.calc_D();
            K = obj.calc_K();
            v = obj.material.poisson;
            a = obj.geometry.radius;
            %since this is analytical data the u and v displacements are 0
            %so their derivatives are also zero
            dudL = 0; dvd0 = 0; dvdL= 0; dud0 = 0;
            NL = D*(dudL+ v/a *(obj.w + dvd0)) - K/a*obj.d2wdL2;
            N0 = D*(1/a*dvd0 + v*dudL+obj.w/a) + K/a^3 * ( obj.w + obj.d2wd02);
            N0L = D*(1-v)/2*(dvdL + 1/a*dud0)+ K*(1-v)/2 * (obj.d2wdL0/a^2 + dud0/a^3);
            NL0 = D*(1-v)/2*(dvdL + 1/a*dud0)+ K*(1-v)/2 * (-obj.d2wdL0/a^2 + dvdL/a^2);
        end
        function [M0, ML, M0L, ML0] = calculate_moment_resultants(obj)
            K = obj.calc_K();
            v = obj.material.poisson;
            a = obj.geometry.radius;
            %since this is analytical data the u and v displacements are 0
            %so their derivatives are also zero
            dudL = 0; dvd0 = 0; dvdL= 0; dud0 = 0;
            ML = K*(v/a^2*obj.d2wd02 + obj.d2wdL2 - 1/a*dudL - v/a^2*dvd0);
            M0 = K*(obj.w/a^2 + obj.d2wd02/a^2 + v*obj.d2wdL2);
            M0L = K*(1-v)*(1/a*obj.d2wdL0 - 1/a*dvdL);
            ML0 = K*(1-v)*(1/a*obj.d2wdL0 - 1/(2*a)*dvdL + 1/(2*a^2)*dud0);
        end
        function [Q0, QL] = calculate_shear_resultants(obj)
            K = obj.calc_K();
            v = obj.material.poisson;
            a = obj.geometry.radius;
            dudL2 = 0; dvdL0 = 0; dud02 = 0; dvdL2 = 0;
            QL = K* (v/a*obj.d3wdL02 + obj.d3wdL3- 1/a*dudL2 - v/a^2 * dvdL0 + ...
                    ((1-v)*(obj.d3wdL02/a^2 - 1/(2*a^2) * dvdL0 + 1/(2*a^3)*dud02)));
            Q0 = K* (1/a^3*obj.dwd0 + 1/a^3*obj.d3wd03 + v/a* obj.d3wdL02 + (1-v)/a*(obj.d3wdL20- dvdL2));
        end
        function [q0, qL] = calculate_power_flow(obj)
            [q0norm, qLnorm] = obj.calculate_power_flow_normal;
            [q0bend, qLbend] = obj.calculate_power_flow_bending;
            [q0twst, qLtwst] = obj.calculate_power_flow_twisting;
            [q0shear, qLshear] = obj.calculate_power_flow_shear;
            q0 = q0norm - q0bend - q0twst + q0shear;
            qL = qLnorm - qLbend - qLtwst + qLshear;
        end
        %% compute power flow using ext flex and curv
        function [qX, qY] = calculate_dimensional_power_flow(obj)
            [qXe, qYe] = obj.calculate_power_flow_extensional();
            [qXf, qYf] = obj.calculate_power_flow_flextural();
            [qXc, qYc] = obj.calculate_power_flow_curvature();
            qX = qXe + qXf + qXc;
            qY = qYe + qYf + qYc;
        end
        function [qX, qY] = calculate_power_flow_extensional(obj)
            D = obj.calc_D();
            dudx = zeros(size(obj.w)); 
            dvdy = dudx; dudy = dudx; dvdx = dudx;
            udot = dudx; vdot = dudx;
            nu = obj.material.poisson;
            qX = -1/2*D*real((dudx+ nu*dvdy).*conj(udot) - (1-nu)/2*(dudy+dvdx).*conj(vdot));
            qY = -1/2*D*real((nu*dudx+ dvdy).*conj(vdot) - (1-nu)/2*(dudy+dvdx).*conj(udot));
        end
        function [qX, qY] = calculate_power_flow_flextural(obj)
            K = obj.calc_K;
            nu = obj.material.poisson;
            [Bdotx, Bdoty] = obj.calculate_beta_dots;
            qX = -K/2*real((obj.d3wdx3 + obj.d3wdxy2).*conj(obj.wdot) - (obj.d2wdx2 + nu*obj.d2wdy2) .* conj(Bdoty) - (1-nu)*obj.d2wdxy.*conj(Bdotx));
            qY = -K/2*real((obj.d3wdx2y + obj.d3wdy3).*conj(obj.wdot) - (nu*obj.d2wdx2 + obj.d2wdy2) .* conj(Bdotx) - (1-nu)*obj.d2wdxy.*conj(Bdoty));
        end
        function [qX, qY] = calculate_power_flow_curvature(obj)
            udot = zeros(size(obj.w)); vdot = udot;
            dudxy = udot; dudy = udot; dvdx = udot; dvdxy = udot;
            dudy2 = udot; dvdx2 = udot;
            ax = Inf;
            ay = obj.geometry.radius;
            D = obj.calc_D; nu = obj.material.poisson; K = obj.calc_K();
            [Bdotx, Bdoty] = obj.calculate_beta_dots();
            qX = -D/ax*nu*obj.w.*udot ...
                + -K*(nu*(obj.w/ax^2 + obj.d2wdy2).*conj(udot) ...
              - (1-nu)/2*(obj.d2wdxy - 1/ax*dvdx).*conj(vdot) ...
                      -(1-nu)/2*(dudy2-dvdxy).*conj(obj.wdot) ... 
                             - nu/ax*obj.dwdx.*conj(obj.wdot) ...
                                   + nu/ax*obj.w.*conj(Bdoty) ...
                          + (1-nu)/2*(dudy-dvdx).*conj(Bdotx));
                                            
            qY = -D/ay*obj.w.*vdot...
                + -K*((1/ay*dudy + obj.d2wdxy).*conj(udot) ...
                                 + obj.d2wdx2.*conj(vdot) ...
                  -(1-nu)/2*(dudxy-dvdx2).*conj(obj.wdot) ...
                          - 1/ay*obj.dwdy.*conj(obj.wdot) ...
                      + (1-nu)/2*(dudy-dvdx).*conj(Bdoty) ...
                              + (1/ay)*obj.w.*conj(Bdotx));
                         
            qX = -1/2*real(qX);
            qY = -1/2*real(qY);
        end
        %% compute power flow using normals moments and shears
        function [q0p, qLp] = calculate_power_flow_part(obj, type)
            switch type
                case {'normal', 'norm', 'N', 'n'}
                    [q0p, qLp] = obj.calculate_power_flow_normal();
                case {'bending', 'bend', 'B', 'b'}
                    [q0p, qLp] = obj.calculate_power_flow_bending();
                case{'twisting', 'twist', 'twst', 'T', 't'}
                    [q0p, qLp] = obj.calculate_power_flow_twisting();
                case{'shear', 'shr', 'S', 's'}
                    [q0p, qLp] = obj.calculate_power_flow_shear();
                case{'e', 'ext', 'extentional', 'E'}
                    [q0p, qLp] = obj.calculate_power_flow_extensional();
                case {'f', 'flex', 'flexural', 'F'}
                    [q0p, qLp] = obj.calculate_power_flow_flextural();
                case {'c', 'curve', 'curvature', 'C'}
                    [q0p, qLp] = obj.calculate_power_flow_curvature();
            end
        end
        function [q0p, qLp] = calculate_power_flow_normal(obj)
            [N0, NL, N0L, NL0] = obj.calculate_normal_resultants();
            [uSd, vSd, ~] = obj.calculate_conj_velocities();
            q0p = -1/2*real(N0.*uSd + N0L.*vSd);
            qLp = -1/2*real(NL.*vSd + NL0.*uSd);
        end
        function [q0p, qLp] = calculate_power_flow_bending(obj)
            [M0, ML, ~, ~] = obj.calculate_moment_resultants();
            [Beta0dotS, BetaLdotS] = obj.calculate_conj_ang_velocities();
            q0p = -1/2*real(M0.*BetaLdotS);
            qLp = -1/2*real(ML.*Beta0dotS);           
        end
        function [q0p, qLp] = calculate_power_flow_twisting(obj)
            [~, ~, M0L, ML0] = obj.calculate_moment_resultants();
            [Beta0dotS, BetaLdotS] = obj.calculate_conj_ang_velocities();
            q0p = -1/2*real(M0L .* Beta0dotS);
            qLp = -1/2*real(ML0 .* BetaLdotS);
        end
        function [q0p, qLp] = calculate_power_flow_shear(obj)
            [Q0, QL] = obj.calculate_shear_resultants();
            [~, ~, wSd] = obj.calculate_conj_velocities();
            q0p = -1/2*real(Q0.*wSd);
            qLp = -1/2*real(QL.*wSd);
        end
        %% calculate velocities
        function [ud, vd, wd] = calculate_velocities(obj)
            wd = obj.wdot;
            ud = zeros(size(obj.wdot));
            vd = zeros(size(obj.wdot));
        end
        function [uSd, vSd, wSd] = calculate_conj_velocities(obj)
            [ud, vd, wd] = obj.calculate_velocities();
            wSd = conj(wd);
            uSd = conj(ud);
            vSd = conj(vd);
        end
        function [BdotO, BdotL] = calculate_ang_velocities(obj)
            BdotO = obj.d2wdLt;
            BdotL = obj.d2wd0t/obj.geometry.radius;
        end
        function [BdotX, BdotY] = calculate_beta_dots(obj)
            BdotX = obj.d2wdyt;
            BdotY = obj.d2wdxt;
        end
        function val = get_ang_vel(obj, type)
            if contains(type, 'L')
                [~, val] = obj.calculate_ang_velocities();
            else
                [val, ~] = obj.calculate_ang_velocities();
            end
        end
        function [beta0dotS, betaLdotS] = calculate_conj_ang_velocities(obj)
            [Bdot0, BdotL] = obj.calculate_ang_velocities();
            beta0dotS = conj(Bdot0);
            betaLdotS = conj(BdotL);
        end
        %% plotting helpers
        function plot_pf_forces(obj, axisHandle)
            if nargin < 2
                axisHandle = gca;
            end
            for i = 1:length(obj.forces)
                if obj.forces(i).enabled
                    Y = obj.forces(i).xStar;
                    X = obj.forces(i).tStar* obj.geometry.radius;
                    hold on;
                    scatter(axisHandle, X,Y, 'x', 'MarkerEdgeColor', [0 0 0], 'SizeData', 100)
                    hold off;
                end
            end
        end
        function plot_forces(obj, Z, axisHandle)
            if nargin < 3
                axisHandle = gca;
            end
            for i = 1:length(obj.forces)
                if obj.forces(i).enabled
                    Y = obj.forces(i).xStar;
                    X = obj.forces(i).tStar* obj.geometry.radius;
                    hold on;
                    scatter3(axisHandle, X,Y,Z,'MarkerEdgeColor', [0 0 0] , 'MarkerFaceColor', [.1,.1,.1], 'SizeData', 50)
                    hold off;
                end
            end
        end
        function [Th, Lng, values] = get_plot_values(obj, type, isReal)
            [Th, Lng] = meshgrid(obj.thetas, obj.longs);
            switch type
                case {'w', 'W'}
                    values = obj.w;
                case {'wDot', 'WDot', 'Wdot', 'wdot', 'dwdt'}
                    values = obj.wdot;
                case {'dwdL', 'dwd0', 'd2wdL2', 'd2wdL0', 'd2wd02', 'd3wdL3', 'd3wdL20', 'd3wdL02', 'd3wd03', 'd2wdLt', 'd2wd0t'}
                    values = obj.(type);
                case {'M0', 'ML', 'ML0', 'M0L', 'N0', 'NL', 'N0L', 'NL0', 'Q0', 'QL'}
                    values = obj.calculate_resultant(type);
                case {'B0', 'BL'}
                    values = obj.get_ang_vel(type);
                otherwise
                    warning('Requested info not availble: Check type entered');
                    fprintf('Possible Types are: w, wDot, dwdL, dwd0, d2wdL2, d2wdL0, d2wd02, d3wdL3, d3wdL20, d3wdL02, d3wd03, d2wdLt, d2wd0t or resultants L/0')
                    values = nan(size(obj.w));
            end
            if isReal
                values = real(values);
            else
                values = imag(values);
            end
        end
        function [X,Y,U,V] = get_quiver_values(obj, type)
            [X,Y] = meshgrid(obj.thetas, obj.longs);
            if nargin < 2 || isempty(type)
                [U,V] = obj.calculate_power_flow();
            else
                [U, V] = obj.calculate_power_flow_part(type);
            end                
        end 
        function clear_animation_timer(obj)
            if ~isempty(obj.animationTimer)
                stop(obj.animationTimer);
                delete(obj.animationTimer);
                obj.animationTimer = [];
            end
        end
        function setupTimer(obj, duration)
            dt = .033; % about 30 Hz
            if nargin >1
                steps = floor(duration/dt);
            else
                steps = obj.timesteps * 5;
            end
            obj.animationTimer = timer('TimerFcn',@(~,~) obj.plot_next_animation_frame,...
                                    'ExecutionMode','fixedRate',...
                                    'Period',dt,...
                                    'TasksToExecute',steps);
        end
        function plot_next_animation_frame(obj)
            if size(obj.surfZ,3) > 1
                Z = obj.surfZ(:,:,1);
            else
                Z = obj.surfZ;
            end
            if isempty(obj.aniPlotHandle) || ~isgraphics(obj.aniPlotHandle)
                obj.aniPlotHandle = figure();
                obj.aniPlotHandle.DeleteFcn = @(~,~) obj.clear_animation_timer();
                obj.surfHandle = surf(obj.surfX, obj.surfY, Z, obj.loopedValues(:,:,obj.currentLoopStep));
                axis equal
                obj.surfHandle.EdgeAlpha = 0.0;
                obj.surfHandle.FaceColor = 'interp';
                colorbar
                clim = max(max(max(obj.loopedValues)));
                caxis([-clim, clim]) 
            else
                obj.surfHandle.ZData = Z;
                obj.surfHandle.CData = obj.loopedValues(:,:,obj.currentLoopStep);
            end
            sz = size(obj.loopedValues,3);
            obj.currentLoopStep = mod(obj.currentLoopStep+1, sz)+1; 
        end
        function cyl_disp_animation_loop(obj, useReal)
            if isgraphics(obj.aniPlotHandle)
                delete(obj.aniPlotHandle);
            end
            obj.currentLoopStep = 1;
            if nargin < 2 
                useReal = true;
            end
            [T, H] = meshgrid([obj.thetas 0], obj.longs);
            obj.surfX = obj.geometry.radius *cos(T);
            obj.surfY = obj.geometry.radius * sin(T);
            obj.surfZ = H;
            t = obj.calc_time_vector();
            ft(1,1,:) = exp(1i*obj.forces(1).omega*t);
            values = obj.w .*ft;
            if useReal
                values = real(values);
            else
                values = imag(values);
            end
            obj.loopedValues = zeros(size(values) + [0 1 0]);
            for i = 1:size(values,3)
                obj.loopedValues(:,:,i) = [values(:,:,i)  values(:,1,i)];
            end            
        end
        function flat_disp_animation_loop(obj, useReal)
            if isgraphics(obj.aniPlotHandle)
                delete(obj.aniPlotHandle);
            end
            obj.currentLoopStep = 1;
            if nargin < 2 
                useReal = true;
            end
            [T, H] = meshgrid(obj.thetas, obj.longs);
            obj.surfX = obj.geometry.radius * T;
            obj.surfY = H;
            t = obj.calc_time_vector();
            ft(1,1,:) = exp(1i*obj.forces(1).omega*t);
            values = obj.w .*ft;
            if useReal
                values = real(values);
            else
                values = imag(values);
            end
            obj.surfZ = values;
            obj.loopedValues = values;      
        end
        %% helpers
        function clean_all_simulated_variables(obj, Lsize, Tsize)
            varNames = obj.names;
            for v = 1:length(varNames)
                obj.(varNames{v}) = zeros(Lsize, Tsize);
            end
        end
        function clean_per_force_variables(obj, Lsize, Tsize)
            varNames = obj.perForceNames;
            for v = 1:length(varNames)
                obj.(varNames{v}) = zeros(Lsize, Tsize);
            end
        end
        function add_scale_force_contribution_to_total(obj, C)
            fNames = obj.perForceNames;
            tNames = obj.names;
            for i = 1:length(tNames)
                obj.(tNames{i}) = obj.(tNames{i}) + C*obj.(fNames{i});
            end
        end
        function denom = calculate_denominator(obj, frc, M, N)
            if obj.mesh.longitude_modes == 1
                wmn = obj.nOmega(N+1);
            else
                wmn = obj.nOmega(M, N+1); % N is naturally zero based so the one shifts to the correct index
            end
            wratsqr = (frc.omega/wmn)^2; 
%             if obj.dampingOn
%                 if obj.constDamp
%                     zeta = obj.material.eta;
%                 else
                    zeta = obj.material.eta*wmn/(2*frc.omega); %% see soedel equ 8.3.3 and 14.2.21
%                 end
%             else
%                 zeta = 0.0;
%             end
            denom = wmn^2.*((1 - wratsqr)^2 + (4*zeta^2*wratsqr))^0.5;
            if N == 0
                denom = denom*2;
            end
        end
        function C = calculate_C_from_force(obj, force)
            Magnitude = force.magnitude;
            rho = obj.material.density;
            thick = obj.geometry.thickness;
            r = obj.geometry.radius;
            l = obj.geometry.height;
            C = 2*Magnitude/(rho*thick*r*l*pi);
        end
        function lStar = normalize_force_long_position(obj, frc)
            lStar = pi*frc.xStar/obj.geometry.height;        
        end
        function timeVector = calc_time_vector(obj)
            t =  linspace(0, 1/obj.forces(1).frequency, obj.timesteps+1);
            timeVector = t(1:end-1); % get rid of last to prevent aliasing
        end
        
        function value = lambda(obj)
            value = obj.material.density*obj.geometry.thickness*obj.material.eta*obj.forces(1).omega;
        end
        function D = calc_D(obj)
            E = obj.material.E;
            mu = obj.material.poisson;
            h = obj.geometry.thickness;
            D =  E*h/(1-mu^2); % soedel eq 2.5.10
        end
        function K = calc_K(obj)
            E = obj.material.E;
            mu = obj.material.poisson;
            h = obj.geometry.thickness;
            K = E*h^3/(12*(1-mu^2)); % soedel eq 2.5.16
        end
        %% default makers
        function create_default(obj)
            obj.geometry = cylinderGeometry();
            obj.forces = [sinusoidalForce(1), sinusoidalForce(2)];
            obj.material = simMaterial();
            obj.mesh = meshDetails();
            obj.timesteps = 100;
%             obj.saveFolder = "C:\Users\ME\Desktop\Simulated Data\";
        end
        function create_default_at_frequency(obj, frequency)
            obj.geometry = cylinderGeometry();
            obj.forces = [sinusoidalForce( [.1, -60*pi/180], 500, frequency, 0), sinusoidalForce( [.4, -240*pi/180], 500, frequency, 179*pi/180)];
            obj.material = simMaterial();
            obj.mesh = meshDetails();
            obj.timesteps = 20;
%             obj.saveFolder = "C:\Users\ME\Desktop\Simulated Data\";
        end
        %% coordinate conversion
        function convert_to_rectangular_coordinates(obj)
            % converts the displacement from the cylindrical coordinate
            % system it is calculated in to the cartesian coordinate system
            % which the SLDV measures in
            xtrans = reshape(repmat(cos(obj.thetas), length(obj.longs),1), length(obj.thetas)*length(obj.longs),1);
            ytrans = reshape(repmat(sin(obj.thetas), length(obj.longs),1), length(obj.thetas)*length(obj.longs),1);
            xcoord = obj.geometry.radius*xtrans;
            ycoord = obj.geometry.radius*ytrans;
            zcoord = repmat(obj.thetas,length(obj.thetas),1)';
            tdata.myX = reshape(xcoord,numel(xcoord), 1);
            tdata.myY = reshape(ycoord,numel(ycoord), 1);
            tdata.myZ = reshape(zcoord,numel(zcoord), 1);
            radialDisplacement = reshape(obj.w, numel(obj.w), 1);
            dispX = radialDisplacement.*xtrans;
            tdata.rX = real(dispX);
            tdata.iX = imag(dispX);
            dispY = radialDisplacement.*ytrans;
            tdata.rY = real(dispY);
            tdata.iY = imag(dispY);
            velX = 1i*obj.forces(1).omega*dispX;
            tdata.rXv = real(velX);
            tdata.iXv = imag(velX);
            velY = 1i*obj.forces(1).omega*dispY;
            tdata.rYv = real(velY);
            tdata.iYv = imag(velY);
            tdata.rZ = 0*tdata.rY;
            tdata.iZ = 0*tdata.iY;
            tdata.rZv = 0*tdata.rZ;
            tdata.iZv = 0*tdata.rZ;
            obj.data = tdata;
        end
        %% file writing functions
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
    end
    methods % getters
        function value = get.fileName(obj)
            date = datestr(floor(now));
            date = strrep(date,'-','_');
            value = strcat('cylinder_simulation_2_forces_', date, '_', num2str(obj.forces(1).frequency), '_Hz');
        end
        function value = get.longs(obj)
            value = linspace(0,obj.geometry.height,obj.mesh.longitude_divisions+1);
        end
        function value = get.thetas(obj)
            temp = linspace(0, 2*pi, obj.mesh.theta_divisions+1);
            value = temp(1:end-1);
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
        function link_subplots(s1, s2)
            Z1 = s1.ZLim; Z2 = s2.ZLim;
            Z = [min([Z1(1), Z2(1)]) max([Z1(2), Z2(2)])];
            s1.ZLim = Z;
            Link = linkprop([s1,s2],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim','ZLim'});
            setappdata(gcf, 'StoreTheLink', Link);
        end
    end
end