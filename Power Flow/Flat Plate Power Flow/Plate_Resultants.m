classdef Plate_Resultants
    properties
        % object that can evaluate splines and their derivatives for the
        % components of the displacment fields
        splines
        % hold information on the physical cylinder
        plateProps
        % the parameter points where the reaction forces are to be
        % calculated
        xi % the parameter that represents distance in the tangential direction
        eta % the parameter that represents distance in the longitudinal direction
    end
    properties (Hidden = true)
        D % stiffness coeff
        mu % poissons
        % normal displacement and velocity
        w
        wDot
        % angluar velocities
        dw2_dydt
        dw2_dxdt
        % 1st order spatial derivs
        dw_dx
        dw_dy
        % 2nd order spatial derivs
        d2w_dx2
        d2w_dxdy
        d2w_dy2
        % 3rd order spatial derivs
        d3w_dx3
        d3w_dx2dy
        d3w_dxdy2
        d3w_dy3
    end
    methods
        function obj = Plate_Resultants(splines, plateProps, xi, eta)
            % check performed by Cylindrical_Power_Flow
            obj.splines = splines;
            obj.plateProps = plateProps;
            if nargin > 2
                obj.xi = xi;
                obj.eta = eta;
            end
            obj.mu = plateProps.poissons;
        end
        function [Mx, My, Mxy, Qx, Qy, obj] = calculate_resultants(obj,xi,eta)
            if nargin > 1
                obj = obj.set_xi_and_eta(xi,eta);
            end
            obj = obj.calculate_D_coeff();
            obj = obj.calculate_necessary_derivatives();
            [Mx, My, Mxy] = obj.calculate_moment_resultants();
            [Qx, Qy] = obj.calculate_shear_resultantss();
        end
        function [wDot, OxDot, OyDot] = get_velocities(obj, xi, eta)
            if nargin > 1
                obj = obj.set_xi_and_eta(xi,eta);
            end
            wDot = obj.wDot;
            OxDot = obj.dw2_dydt;
            OyDot = obj.dw2_dxdt;
        end
        function [wDotS, OxDotS, OyDotS] = get_conjugate_velocities(obj, xi,eta)
            if nargin > 1 
                obj = obj.set_xi_and_eta(xi,eta);
            end
            wDotS = conj(obj.wDot);
            OxDotS = conj(obj.dw2_dydt); % Oxdot = d2w_dydt
            OyDotS = conj(-obj.dw2_dxdt); % Oydot = - d2w_dxdt
        end
        function [X,Y,Z] = plot_derivatives(obj, derivative, isReal, is3D)
            if nargin< 3
                isReal = true;
            end
            if nargin < 4
                is3D = false;
            end
            switch derivative
                case 'dw_dy'
                    values = obj.dw_dy;
                case 'dw_dx'
                    values = obj.dw_dx;
                case 'd2w_dy2'
                    values = obj.d2w_dy2;
                case 'd2w_dx2'
                    values = obj.d2w_dx2;
                case 'd2w_dxdy'
                    values = obj.d2w_dxdy;
                case 'w'
                    values = obj.w;
                case 'd3w_dx2dy'
                    values = obj.d3w_dx2dy;
                case 'd3w_dy3'
                    values = obj.d3w_dy3;
                case 'd3w_dx3'
                    values = obj.d3w_dx3;
                case 'd3w_dxdy2'
                    values = obj.d3w_dxdy2;
            end
            if isReal
                plotValues = real(values);
                tri = '_Real';
            else
                plotValues = imag(values);
                tri = '_Imag';
            end
            [X, Y, Z] = obj.plot_values(plotValues, is3D);
            title([derivative, tri])
        end
        function [X,Y,Z] = plot_displacement(obj, isReal, is3D, trim)
            disp = obj.w;
            if nargin < 4
                trim = false;
            end
            if nargin < 3 
                is3D = false;
            end
            if isReal
                plotValues = real(disp);
                titleString = 'w Real Values';
            else
                plotValues = imag(disp);
                titleString = 'w Imag Values';
            end
            if trim
                plotValues= plotValues(2:end-1, 2:end-1);
            end
            if nargout < 1
                obj.plot_values(plotValues, is3D, trim);
            else
                [X,Y,Z] = obj.plot_values(plotValues, is3D, trim);
            end
            title(titleString)
        end
        function [X,Y,Z] = plot_velocity(obj,isReal, is3D, trim)
            vel = obj.wDot;
            if nargin < 4 
                trim = false;
            end
            if nargin < 3
                is3D = false;
            end
            if isReal
                plotValues = real(vel);
                titleString = 'wdot Real Values';
            else
                plotValues = imag(vel);
                titleString = 'wdot Imag Values';
            end
            if trim 
                plotValues = plotValues(2:end-1, 2:end-1);
            end
            [X,Y,Z] = obj.plot_values(plotValues, is3D, trim);
            title(titleString)
        end
        function [X,Y,Z] = plot_values(obj, plotValues, is3D, trim)
            [Eta, Xi] = meshgrid(obj.eta, obj.xi);
            if trim
               Eta = Eta(2:end-1, 2:end-1);
               Xi = Xi(2:end-1, 2:end-1);
            end
            if ~is3D
                surf(Xi, Eta, zeros(size(Eta)), plotValues, 'EdgeColor', 'none')
                view(0,90);
            else
                surf(Eta, Xi, plotValues, plotValues) % flipping Xi and Eta allow for plots to look like Dr.Blotters Dissertation
                xlabel('y')
                ylabel('x')
            end
            X = Eta; Y = Xi; Z = plotValues;
            colormap jet
            colorbar
        end
        function [X,Y,Z] = plot_resultant(obj, resultant, isReal, is3D, trim)
            if nargin < 5
                trim = false;
            end
            if nargin<3
                isReal = true;
            end
            if nargin < 4
                is3D = false;
            end
            switch resultant
                case 'Mx'
                    [values, ~, ~] = obj.calculate_moment_resultants();
                case 'My'
                    [~, values, ~] = obj.calculate_moment_resultants();
                case {'Mxy', 'Myx'}
                    [~, ~, values] = obj.calculate_moment_resultants();
                case 'Qx'
                    [values, ~] = obj.calculate_shear_resultantss();
                case 'Qy'
                    [~, values] = obj.calculate_shear_resultantss();
                otherwise
                    warning('Incorrect input');
            end
            if isReal
                Values = real(values);
                titleString = [resultant, ' Real Values'];
            else
                Values = imag(values);
                titleString = [ resultant, ' Imag Values'];
            end
            if trim
                Values = Values(2:end-1, 2:end-1);
            end
            [X,Y,Z] = obj.plot_values(Values, is3D, trim);
            title(titleString);
        end
    end
    methods (Hidden = true)
        function obj = set_xi_and_eta(obj,xi,eta)
            obj.xi = xi;
            obj.eta = eta;
        end
        function obj = calculate_D_coeff(obj)
            % D = -Eh^3/(12*(1-v^2)) Blotter Disseration 3.77
            obj.D = -obj.plateProps.E*obj.plateProps.thickness^3/(12*(1-obj.mu^2));
        end
        function obj = calculate_necessary_derivatives(obj)
%             [XI, ETA] = ndgrid(obj.xi, obj.eta);
            obj.w = obj.get_complex_deriv('w');
%             surf(XI, ETA, real(obj.w))
            obj.dw_dy = obj.get_complex_deriv('dw_dy');
            obj.dw_dx = obj.get_complex_deriv('dw_dx');
            
            obj.d2w_dy2 =obj.get_complex_deriv('d2w_dy2');
%             surf(XI, ETA, real(obj.d2w_dy2))
            obj.d2w_dx2 =  obj.get_complex_deriv('d2w_dx2');
            obj.d2w_dxdy =obj.get_complex_deriv('d2w_dxdy');
            
            obj.d3w_dy3 =obj.get_complex_deriv('d3w_dy3');
            obj.d3w_dxdy2 = obj.get_complex_deriv('d3w_dxdy2');
            obj.d3w_dx2dy = obj.get_complex_deriv('d3w_dx2dy');
            obj.d3w_dx3 = obj.get_complex_deriv('d3w_dx3');
            obj = obj.calculate_velocity();
            obj = obj.calculate_angular_velocities();
        end
        function obj = calculate_angular_velocities(obj)
            isImag = true;
            dw2_dydtReal = obj.splines.evaluate_at_parameters(obj.xi, obj.eta, 'Velocity', ~isImag, 3, [0 1]);
            dw2_dydtImag = obj.splines.evaluate_at_parameters(obj.xi, obj.eta, 'Velocity', isImag, 3, [0 1]);
            obj.dw2_dydt = complex(dw2_dydtReal, dw2_dydtImag); % dw2_dxdt = d2w/dydt
            dw2_dxdtReal = obj.splines.evaluate_at_parameters(obj.xi, obj.eta, 'Velocity', ~isImag, 3, [1 0]);
            dw2_dxdtImag = obj.splines.evaluate_at_parameters(obj.xi, obj.eta, 'Velocity', isImag, 3, [1 0]);
            obj.dw2_dxdt = complex(dw2_dxdtReal, dw2_dxdtImag); % Oydot  = -d2w/dxdt
        end
        function obj = calculate_velocity(obj) %% check this
            isImag = true;
            wDotReal = obj.splines.evaluate_at_parameters(obj.xi, obj.eta, 'Velocity', ~isImag, 3, [0 0]);
            wDotImag = obj.splines.evaluate_at_parameters(obj.xi, obj.eta, 'Velocity', isImag, 3, [0 0]);
            obj.wDot = complex(wDotReal, wDotImag);
        end
        function complexValues = get_complex_deriv(obj, deriv)
            [measure, dimension, derivative] = obj.get_params_from_deriv(deriv);
            realValues = obj.splines.evaluate_at_parameters(obj.xi, obj.eta, measure, false, dimension, derivative);
            imagValues = obj.splines.evaluate_at_parameters(obj.xi, obj.eta, measure, true, dimension, derivative);
            complexValues = complex(realValues, imagValues);
        end
        function [Mx, My, Mxy] = calculate_moment_resultants(obj)
            %% seeblotter's dissertation 3.74, 76
            nu = obj.plateProps.poissons;
            Mx = obj.D * (obj.d2w_dx2 + nu*obj.d2w_dy2);  
            My = obj.D * (obj.d2w_dy2 + nu*obj.d2w_dx2);
            Mxy = obj.D * (1-nu) * obj.d2w_dxdy;
        end
        function [Qx, Qy] = calculate_shear_resultantss(obj)
            %% see blotter's dissertation 3.75
            Qx = obj.D * (obj.d3w_dx3 + obj.d3w_dxdy2);
            Qy = obj.D * (obj.d3w_dx2dy + obj.d3w_dy3);
        end
    end
    methods (Static = true, Hidden = true)
        function [measure, dimension, derivative] = get_params_from_deriv(deriv)
            measure = 'disp';
            parts = strsplit(deriv, '_');
            dimension = 3;
            if length(parts) == 2
                inxz = strfind(parts{2}, 'y');
                if isempty(inxz)
                    dys = 0;
                elseif length(parts{2}) > inxz  && ~isnan(str2double(parts{2}(inxz+1)))
                    dys = str2double(parts{2}(inxz+1));
                else
                    dys = 1;
                end
                inxt = strfind(parts{2}, 'x');
                if isempty(inxt)
                    dxs = 0;
                elseif length(parts{2}) > inxt && ~isnan(str2double(parts{2}(inxt+1)))
                    dxs = str2double(parts{2}(inxt+1));
                else
                    dxs = 1;
                end
                derivative = [dxs, dys];
            else
                derivative = [0,0];
            end
        end
    end
end