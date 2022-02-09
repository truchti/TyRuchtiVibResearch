classdef CylindricalPowerFlowResultants
    properties
        % object that can evaluate splines and their derivatives for the
        % components of the displacment fields
        splines
        % hold information on the physical cylinder
        cylProps
        % the parameter points where the reaction forces are to be
        % calculated
        xi % the parameter that represents distance in the tangential direction
        eta % the parameter that represents distance in the longitudinal direction
    end
    properties (Dependent)
        
    end
    properties (Hidden = true)
        D % stiffness coeff
        K % rigidity coeff
        a % radius        
        mu % poissons
        % longitudinal displacement derivatives
        long
        longVel
        dlong_dz
        dlong_dth
        %tangential displacement derivatives
        tang
        tangVel
        dtang_dz
        dtang_dth
        % radial displacement and 1st 2nd and 3rd order derivatives
        rad
        radVel
        drad_dth
        drad_dz
        d2rad_dz2
        d2rad_dth2
        d2rad_dthdz
        d3rad_dth2dz
        d3rad_dz3
        d3rad_dth3
        d3rad_dthdz2
        d2tang_dz2 
        d2tang_dthdz
        d2long_dz2
        d2long_dth2
    end
    methods
        function obj = CylindricalPowerFlowResultants(splines, cylProps, xi, eta)
            % check performed by Cylindrical_Power_Flow
            obj.splines = splines;
            obj.cylProps = cylProps;
            if nargin > 2
                obj.xi = xi;
                obj.eta = eta;
            end
            obj.a = cylProps.radius;
            obj.mu = cylProps.poissons;
        end
        function [Nt, Nl, Ntl, Nlt, Mt, Ml, Mtl, Mlt, Qt, Ql, obj] = calculate_resultants(obj,xi,eta)
            if nargin > 1
                obj = obj.set_xi_and_eta(xi,eta);
            end
            obj = obj.calculate_D_coeff();
            obj = obj.calculate_K_coeff();
            obj = obj.calculate_necessary_derivatives();
            [Nt, Nl, Ntl, Nlt] = obj.calculate_normal_resultants();
            [Mt, Ml, Mtl, Mlt] = obj.calculate_moment_resultants();
            [Qt, Ql] = obj.calculate_shear_resultantss();
        end
        function plot_derivatives(obj, derivative)
            switch derivative
                case 'dlong_dz'
                    values = obj.dlong_dz;
                case 'dlong_dth'
                    values = obj.dlong_dth;
                case 'dtang_dz'
                    values = obj.dtang_dz;
                case 'dtang_dth'
                    values = obj.dtang_dth;
                case 'drad_dz'
                    values = obj.drad_dz;
                case 'drad_dth'
                    values = obj.drad_dth;
                case 'd2rad_dz2'
                    values = obj.d2rad_dz2;
                case 'd2rad_dth2'
                    values = obj.d2rad_dth2;
                case 'd2rad_dthdz'
                    values = obj.d2rad_dthdz;
                case 'rad'
                    values = obj.rad;
                case 'd3rad_dth2dz'
                    values = obj.d3rad_dth2dz;
                case 'd3rad_dz3'
                    values = obj.d3rad_dz3;
                case 'd3rad_dth3'
                    values = obj.d3rad_dth3;
                case 'd3rad_dthdz2'
                    values = obj.d3rad_dthdz2;
            end
            rValues = real(values);
            iValues = imag(values);
            [Eta, Xi] = meshgrid(obj.eta, obj.xi);
            subplot(1,2,1)
            surf(Xi, Eta, rValues)
            title('Real Values')
            subplot(1,2,2)
            surf(Xi, Eta, iValues)
            title('Imaginary Values')
        end
        function plot_displacement(obj, displacement)
            switch displacement
                case { 'r'; 'rad'; 'radial'}
                    values = obj.rad;
                case {'t'; 'theta'}
                    values = obj.tang;
                case {'l'; 'long'; 'longitude'; 'z'}
                    values = obj.long;
            end
            rValues = real(values);
            iValues = imag(values);
            [Eta, Xi] = meshgrid(obj.eta, wrapTo2Pi(obj.xi));
            X = cos(Xi)*obj.cylProps.radius;
            Y = sin(Xi)*obj.cylProps.radius;
            
%             figure
%             subplot(1,2,1)
            surf(obj.cylProps.radius*Xi, Eta, zeros(size(Eta)),  iValues, 'EdgeColor', 'none')
            view(0,90);
            colormap jet
%             surf(X, Y, Eta, rValues)
            title('Real Values')
%             axis equal
            colorbar
%             subplot(1,2,2)
%             surf(X, Y, Eta, iValues)
%             title('Imaginary Values')
%             axis equal
%             colorbar
        end
        function plot_velocity(obj, direc)
            switch direc
                case { 'rdot'; 'raddot'; 'radialdot'}
                    values = obj.radVel;
                case {'tdot'; 'thetadot'}
                    values = obj.tangVel;
                case {'ldot'; 'longdot'; 'longitudedot'; 'zdot'}
                    values = obj.longVel;
            end
            rValues = real(values);
            iValues = imag(values);
            [Eta, Xi] = meshgrid(obj.eta, wrapTo2Pi(obj.xi));
            X = cos(Xi)*obj.cylProps.radius;
            Y = sin(Xi)*obj.cylProps.radius;
            
%             figure
%             subplot(1,2,1)
            surf(obj.cylProps.radius*Xi, Eta, zeros(size(Eta)),  rValues, 'EdgeColor', 'none')
            view(0,90);
%             surf(X, Y, Eta, rValues, 'EdgeColor', 'none')
            title('Real Values')
            colormap jet
%             axis equal
            colorbar
%             subplot(1,2,2)
%             surf(X, Y, Eta, iValues)
%             title('Imaginary Values')
%             axis equal
%             colorbar
        end
        function plot_flat_displacement(obj, displacement)
            switch displacement
                case { 'r'; 'rad'; 'radial'}
                    values = obj.rad;
                case {'t'; 'theta'}
                    values = obj.tang;
                case {'l'; 'long'; 'longitude'; 'z'}
                    values = obj.long;
            end
            rValues = real(values);
            iValues = imag(values);
            [Eta, Xi] = meshgrid(obj.eta, obj.xi);
            figure(1)
            surf(Xi, Eta, zeros(size(Eta)), rValues)
            title('Real Values')
            xlabel('Theta')
            ylabel('Height')
            view(0,90);
            colorbar
            figure(2)
            surf(Xi, Eta, zeros(size(Eta)), iValues)
            title('Imaginary Values')
            xlabel('Theta')
            ylabel('Height')
            view(0,90);
            colorbar
        end
        function plot_resultant(obj, resultant)
            switch resultant
                case 'Nt'
                    [values, ~, ~, ~] = obj.calculate_normal_resultants();
                case 'Nl'
                    [~, values, ~, ~] = obj.calculate_normal_resultants();
                case 'Ntl'
                    [~, ~, values, ~] = obj.calculate_normal_resultants();
                case 'Nlt'
                    [~, ~, ~, values] = obj.calculate_normal_resultants();
                case 'Mt'
                    [values, ~, ~, ~] = obj.calculate_moment_resultants();
                case 'Ml'
                    [~, values, ~, ~] = obj.calculate_moment_resultants();
                case 'Mtl'
                    [~, ~, values, ~] = obj.calculate_moment_resultants();
                case 'Mlt'
                    [~, ~, ~, values] = obj.calculate_moment_resultants();
                case 'Qt'
                    [values, ~] = obj.calculate_shear_resultantss();
                case 'Ql'
                    [~, values] = obj.calculate_shear_resultantss();
            end
            rValues = real(values);
            iValues = imag(values);
            [Eta, Xi] = meshgrid(obj.eta, obj.cylProps.radius*wrapTo2Pi(obj.xi));
%             subplot(1,2,1)
            surf(obj.cylProps.radius*Xi, Eta, zeros(size(Eta)),  rValues, 'EdgeColor', 'none')
            title('Real Values')
            view(0,90);
%             colorbar
%             subplot(1,2,2)
%             surf(Xi, Eta, zeros(size(Eta)), iValues)
%             title('Imaginary Values')
%             view(0,90);
            colorbar
            colormap('jet')
        end
    end
    methods (Hidden = true)
        function obj = set_xi_and_eta(obj,xi,eta)
            obj.xi = xi;
            obj.eta = eta;
        end
        function obj = calculate_D_coeff(obj)
            % D = E*h/(1-v^2)
            obj.D = obj.cylProps.E*obj.cylProps.thickness/(1-obj.mu^2);
        end
        function obj = calculate_K_coeff(obj)
            % K = Eh^3/(12*(1-v^2))
            obj.K = obj.cylProps.E*obj.cylProps.thickness^3/(12*(1-obj.mu^2));
        end
        function obj = calculate_necessary_derivatives(obj)
            obj.long = obj.get_complex_deriv('long');
            obj.dlong_dz = obj.get_complex_deriv('dlong_dz');
            obj.dlong_dth = obj.get_complex_deriv('dlong_dth');
            obj.d2long_dz2 = obj.get_complex_deriv('d2long_dz2');
            obj.d2long_dth2 = obj.get_complex_deriv('d2long_dth2');
            
            obj.tang = obj.get_complex_deriv('tang');
            obj.dtang_dz = obj.get_complex_deriv('dtang_dz');
            obj.dtang_dth = obj.get_complex_deriv('dtang_dth');
            obj.d2tang_dz2 = obj.get_complex_deriv('d2tang_dz2');
            obj.d2tang_dthdz = obj.get_complex_deriv('d2tang_dthdz');
            
            obj.rad = obj.get_complex_deriv('rad');
%             obj.drad_dz = obj.get_complex_deriv('drad_dz');
            obj.drad_dth = obj.get_complex_deriv('drad_dth,');
            
            obj.d2rad_dz2 =obj.get_complex_deriv('d2rad_dz2');
            obj.d2rad_dth2 =  obj.get_complex_deriv('d2rad _dth2');
            obj.d2rad_dthdz =obj.get_complex_deriv('d2rad_dthdz');
            
            obj.d3rad_dz3 =obj.get_complex_deriv('d3rad_dz3');
            obj.d3rad_dthdz2 = obj.get_complex_deriv('d3rad_dthdz2');
            obj.d3rad_dth2dz = obj.get_complex_deriv('d3rad_dth2dz');
            obj.d3rad_dth3 = obj.get_complex_deriv('d3rad_dth3');
            obj = obj.calculate_velocities();
        end
        function obj = calculate_velocities(obj)
            lV1 = obj.splines.evaluate_at_parameters(obj.xi, obj.eta, 'vel', false, 3, [0 0]);
            lV2 = obj.splines.evaluate_at_parameters(obj.xi, obj.eta, 'vel', true, 3, [0 0]);
            obj.longVel = complex(lV1, lV2);
            tV1 = obj.splines.evaluate_at_parameters(obj.xi, obj.eta, 'vel', false, 1, [0 0]);
            tV2= obj.splines.evaluate_at_parameters(obj.xi, obj.eta, 'vel', true, 1, [0 0]);
            obj.tangVel = complex(tV1, tV2);
            rV1 = obj.splines.evaluate_at_parameters(obj.xi, obj.eta, 'vel', false, 2, [0 0]);
            rV2 = obj.splines.evaluate_at_parameters(obj.xi, obj.eta, 'vel', true, 2, [0 0]);
            obj.radVel = complex(rV1, rV2);
        end
        function complexValues = get_complex_deriv(obj, deriv)
            [measure, dimension, derivative] = obj.get_params_from_deriv(deriv);
            realValues = obj.splines.evaluate_at_parameters(obj.xi, obj.eta, measure, false, dimension, derivative);
            imagValues = obj.splines.evaluate_at_parameters(obj.xi, obj.eta, measure, true, dimension, derivative);
            complexValues = complex(realValues, imagValues);
        end
        function [Nt, Nl, Ntl, Nlt] = calculate_normal_resultants(obj)
            %% see lab notebook pgs 12-14 where u is longitudinal v is theta and w is radial displacements
            nu = obj.cylProps.poissons;
            r = obj.a;
            Nl = obj.D * (obj.dlong_dz + nu/r*(obj.dtang_dth +obj.rad)) - obj.K/r*obj.d2rad_dth2;
            Nt = obj.D * (1/r*obj.dtang_dth + nu*obj.dlong_dz + obj.rad/r) + obj.K/r^3 *(obj.rad + obj.d2rad_dth2);
            Ntl = obj.D*(1-nu)/2*(obj.dtang_dz+1/r*obj.dlong_dth) + obj.K*(1-nu)/2 *(obj.d2rad_dthdz/r^2 + obj.dlong_dth/r^3);
            Nlt = obj.D*(1-nu)/2*(obj.dtang_dz+1/r*obj.dlong_dth) + obj.K*(1-nu)/2*(obj.dtang_dz - obj.d2rad_dthdz)/r^2;
        end
        function [Mt, Ml, Mtl, Mlt] = calculate_moment_resultants(obj)
            %% see lab notebook pgs 12-14 where u is longitudinal v is theta and w is radial displacements
            r = obj.a;
            nu = obj.cylProps.poissons;
            Mt= obj.K *(obj.rad/r^2 +obj.d2rad_dth2/r^2 + nu*obj.d2rad_dz2);
            Ml = obj.K*(nu/r^2*obj.d2rad_dth2 + obj.d2rad_dz2 -1/r * obj.dlong_dz - nu/r^2*obj.dtang_dth);
            Mlt = obj.K*(1-nu)/r*(obj.d2rad_dthdz-obj.dtang_dz);
            Mtl = obj.K*(1-nu)/2*(2*obj.d2rad_dthdz/r+ obj.dtang_dz/r + obj.dlong_dth/r^2);
        end
        function [Qtang, Qlong] = calculate_shear_resultantss(obj)
            %% see lab notebook pgs 12-14 where u is longitudinal v is theta and w is radial displacements
            r = obj.a;
            nu = obj.cylProps.poissons;
            dMt_dth_a = obj.K *(obj.drad_dth/r^3 +obj.d3rad_dth3/r^3 + nu/r*obj.d3rad_dthdz2);
            dMtl_dth_a =  obj.K*(1-nu)/2*(2*obj.d3rad_dth2dz/r^2+ obj.d2tang_dthdz/r^2 + obj.d2long_dth2/r^3);
            dMl_dz =  obj.K*(nu/r^2*obj.d3rad_dth2dz + obj.d3rad_dz3 -1/r * obj.d2long_dz2 - nu/r^2*obj.d2tang_dthdz);
            dMlt_dz = obj.K*(1-nu)/r*(obj.d3rad_dthdz2-obj.d2tang_dz2);
           
            Qlong = dMl_dz + dMtl_dth_a;
            Qtang = dMt_dth_a + dMlt_dz;
        end
    end
    methods (Static = true, Hidden = true)
        function [measure, dimension, derivative] = get_params_from_deriv(deriv)
            measure = 'disp';
            parts = strsplit(deriv, '_');
            if contains(parts{1}, 'rad')
                dimension = 2;
            elseif contains(parts{1}, 'long')
                dimension = 3;
            elseif contains(parts{1}, 'tang')
                dimension = 1;
            else
                error('Wrong Dimension for resultant calculations')
            end
            if length(parts) == 2
                inxz = strfind(parts{2}, 'z');
                if isempty(inxz)
                    dzs = 0;
                elseif length(parts{2}) > inxz  && ~isnan(str2double(parts{2}(inxz+1)))
                    dzs = str2double(parts{2}(inxz+1));
                else
                    dzs = 1;
                end
                inxt = strfind(parts{2}, 'th');
                if isempty(inxt)
                    dths = 0;
                elseif length(parts{2}) > inxt + 1 && ~isnan(str2double(parts{2}(inxt+2)))
                    dths = str2double(parts{2}(inxt+2));
                else
                    dths = 1;
                end
                derivative = [dths, dzs];
            else
                derivative = [0,0];
            end
        end
    end
end