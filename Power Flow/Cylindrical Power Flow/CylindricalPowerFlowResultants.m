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
        dlong_dz
        dlong_dth
        %tangential displacement derivatives
        dtang_dz
        dtang_dth
        % radial displacement and 1st 2nd and 3rd order derivatives
        rad
        drad_dth
        drad_dz
        d2rad_dz2
        d2rad_dth2
        d2rad_dthdz
        d3rad_dth2dz
        d3rad_dz3
        d3rad_dth3
        d3rad_dthdz2
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
            [Nt, Nl, Ntl, Nlt] = obj.calculate_normals();
            [Mt, Ml, Mtl, Mlt] = obj.calculate_moments();
            [Qt, Ql] = obj.calculate_shears();
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
        function plot_resultant(obj, resultant)
            switch resultant
                case 'Nt'
                    [values, ~, ~, ~] = obj.calculate_normals();
                case 'Nl'
                    [~, values, ~, ~] = obj.calculate_normals();
                case 'Ntl'
                    [~, ~, values, ~] = obj.calculate_normals();
                case 'Nlt'
                    [~, ~, ~, values] = obj.calculate_normals();
                case 'Mt'
                    [values, ~, ~, ~] = obj.calculate_moments();
                case 'Ml'
                    [~, values, ~, ~] = obj.calculate_moments();
                case 'Mtl'
                    [~, ~, values, ~] = obj.calculate_moments();
                case 'Mlt'
                    [~, ~, ~, values] = obj.calculate_moments();
                case 'Qt'
                    [values, ~] = obj.calculate_shears();
                case 'Ql'
                    [~, values] = obj.calculate_shears();
            end
            rValues = real(values);
            iValues = imag(values);
            [Eta, Xi] = meshgrid(obj.eta, obj.xi);
            figure(1)
            surf(Xi, Eta, rValues)
            title('Real Values')
            figure(2)
            surf(Xi, Eta, iValues)
            title('Imaginary Values')
        end
    end
    methods (Hidden = true)
        function obj = set_xi_and_eta(obj,xi,eta)
            obj.xi = xi;
            obj.eta = eta;
        end
        function obj = calculate_D_coeff(obj)
            obj.D = obj.cylProps.E*obj.cylProps.thickness/(1-obj.mu^2);
        end
        function obj = calculate_K_coeff(obj)
            obj.K = obj.cylProps.E*obj.cylProps.thickness^3/(12*(1-obj.mu^2));
        end
        function obj = calculate_necessary_derivatives(obj)
            obj.dlong_dz = obj.a*obj.get_complex_deriv('dlong_dz');
            obj.dlong_dth = obj.get_complex_deriv('dlong_dth');
            
            obj.dtang_dz = obj.a * obj.get_complex_deriv('dtang_dz');
            obj.dtang_dth = obj.get_complex_deriv('dtang_dth');
            
            obj.rad = obj.get_complex_deriv('rad');
            obj.drad_dz = obj.a * obj.get_complex_deriv('drad_dz');
            obj.drad_dth = obj.get_complex_deriv('drad_dth,');
            
            obj.d2rad_dz2 = obj.a^2 *obj.get_complex_deriv('d2rad_dz2');
            obj.d2rad_dth2 =  obj.get_complex_deriv('d2rad _dth2');
            obj.d2rad_dthdz = obj.a * obj.get_complex_deriv('d2rad_dthdz');
            
            obj.d3rad_dz3 = obj.a^3*obj.get_complex_deriv('d3rad_dz3');
            obj.d3rad_dthdz2 = obj.a^2*obj.get_complex_deriv('d3rad_dthdz2');
            obj.d3rad_dth2dz = obj.a*obj.get_complex_deriv('d3rad_dth2dz');
            obj.d3rad_dth3 = obj.get_complex_deriv('d3rad_dth3');
        end
        function complexValues = get_complex_deriv(obj, deriv)
            [measure, dimension, derivative] = obj.get_params_from_deriv(deriv);
            realValues = obj.splines.evaluate_at_parameters(obj.xi, obj.eta, measure, false, dimension, derivative);
            imagValues = obj.splines.evaluate_at_parameters(obj.xi, obj.eta, measure, true, dimension, derivative);
            complexValues = complex(realValues, imagValues);
        end
        function [Nt, Nl, Ntl, Nlt] = calculate_normals(obj)
            Nt = obj.D/obj.a*(obj.dtang_dth+obj.rad+obj.mu*obj.dlong_dz) + obj.K/obj.a^3*(obj.rad + obj.d2rad_dth2);
            Nl= obj.D/obj.a*(obj.dlong_dz+obj.mu*(obj.dtang_dth+obj.rad)) - obj.K/obj.a^3*(obj.d2rad_dz2);
            temp = (obj.dlong_dth+obj.dtang_dz); 
            Ntl = obj.D*(1-obj.mu)/(2*obj.a) * temp + obj.K/obj.a^3*(1-obj.mu)/2*(obj.dlong_dth+obj.d2rad_dthdz);
            Nlt = obj.D*(1-obj.mu)/(2*obj.a) * temp + obj.K/obj.a^3*(1-obj.mu)/2*(obj.dtang_dz-obj.d2rad_dthdz);
        end
        function [Mt, Ml, Mtl, Mlt] = calculate_moments(obj)
            Kovera2 = obj.K/(obj.a^2);
            Mt= Kovera2 * (obj.rad + obj.d2rad_dth2 + obj.mu*obj.d2rad_dz2);
            Ml= Kovera2 * (obj.d2rad_dz2 + obj.mu*obj.d2rad_dth2- obj.dlong_dz - obj.mu*obj.dtang_dth);
            Mlt = Kovera2 * (1-obj.mu) * (obj.d2rad_dthdz - obj.dtang_dz);
            Mtl = Kovera2 * (1-obj.mu) * (obj.d2rad_dthdz + 1/2*obj.dlong_dth-1/2*obj.dtang_dz);
        end
        function [Qtang, Qlong] = calculate_shears(obj)
            Mtdot = obj.K/obj.a^2*(obj.d3rad_dz3 + obj.mu*obj.d3rad_dth2dz);
            Mlprm = obj.K/obj.a^2*(obj.d3rad_dth3 + obj.mu*obj.d3rad_dthdz2);
            Mltprm = obj.K*(1-obj.mu)/obj.a^2*obj.d3rad_dth2dz;
            Mtldot = obj.K*(1-obj.mu)/obj.a^2*obj.d3rad_dthdz2;
            Qtang = (Mtdot+Mltprm)/obj.a;
            Qlong = (Mlprm+Mtldot)/obj.a;
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
                elseif length(parts{2}) > inxz +1 && isnumeric(parts{2}(inxz+1))
                    dzs = parts{2}(inxz+1);
                else
                    dzs = 1;
                end
                inxt = strfind(parts{2}, 'th');
                if isempty(inxt)
                    dths = 0;
                elseif length(parts{2}) > inxt + 2 && isnumeric(parts{2}(inxt+2))
                    dths = parts{2}(inxt+2);
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