classdef Plate_Power_Flow <handle
    % This class can take in a SLDV_parser, a SLDV_data or no input in the
    % no input case it creates a default SLDV_FastScan parser. If a parser
    % is input then the data is extracted using the parser. If the raw data
    % is processed separately the processed data can be passed instead of a
    % parser
    properties
        laserData
        numControlPoints = [8 8]
        numberOfEvaluationPnts = [50 50];
        repeatedNodes = [6,6];
        plateProperties
        qx
        qy
    end
    properties (Hidden = true)
        % evaluation parameters
        xi
        eta
        % frequency that the plate was forced at
        freq
        % object that allows for evaluating the various surface splines at
        % a give parameter pair
        fieldEvaluator
        % object that is used to solve for the surface splines that best
        % fit the data
        splineFitter
        % object that is used to plot the splines in cartesian coordinates
        % as function of the parameters
        splinePlotter
        % object that is used to calculate the resultant forces based on
        % love shell theory
        PowerFlowResultants
        % provides plotting functionallity for the power flow data
        PFPlotter
        dirty = true;
        iqx
        iqy
        qxbend
        qybend
        qxtwist
        qytwist
        qxshear
        qyshear
    end
    methods
        function obj = Plate_Power_Flow(SLDV_data, properties)
            if isa(SLDV_data, 'SLDV_data')
                obj.laserData = SLDV_data;
            else
                error('Not Yet Implemented')
            end
            if nargin>1 
                obj.plateProperties = properties;
            end
        end
        function [qx, qy] = calculate_power_flow(obj)
            % This function uses the properties and data that it is
            % initialized with to compute the power flow in the x and y directions 
            % it then creates a plotting object with the results so they can be readily visualized
            if obj.dirty
                obj.process_input_data();
            end
            obj.perform_power_flow_computations();
            if nargout == 1
                qx = {obj.qx, obj.qy};
            end
            if nargout == 2
                qx = obj.qx;
                qy = obj.qy;
            end
        end
        %% plotting functions
        function plot_component(obj, component, isReal, trim)
            if nargin < 3
                isReal = true;
                trim = false;
            elseif nargin < 4
                trim = false;
            end
            switch component
                case {'w'}
                    obj.PowerFlowResultants.plot_displacement(isReal);
                case {'Nx'; 'Ny'; 'Nxy'; 'Mx'; 'My'; 'Mxy'; 'Qx'; 'Qy'}
                    obj.PowerFlowResultants.plot_resultant(component, isReal, trim);
                case {'wdot'}
                    obj.PowerFlowResultants.plot_velocity(isReal, trim);
                case {'Oxdot'}
                    obj.plot_ang_velocity(isReal, 'x', trim);
                case {'Oydot'}
                    obj.plot_ang_velocity(isReal, 'y', trim);
                otherwise
                    if contains(component, 'd')
                        obj.PowerFlowResultants.plot_derivatives(component);
                    else
                        warning('Input is not valid')
                    end
            end 
            title(component)
        end
        function [X, Y, Z] = plot_component_3D(obj, component, isReal, trim)
            if nargin < 4
                trim = false;
            end
            if nargin < 3
                isReal = true;
            end
            is3D = true;
            switch component
                case {'w'}
                    [X,Y,Z] = obj.PowerFlowResultants.plot_displacement(isReal, is3D, trim);
                case {'Nx'; 'Ny'; 'Nxy'; 'Mx'; 'My'; 'Mxy'; 'Qx'; 'Qy'}
                    [X,Y,Z] = obj.PowerFlowResultants.plot_resultant(component, isReal, is3D, trim);
                case {'wdot'}
                    [X,Y,Z] = obj.PowerFlowResultants.plot_velocity(isReal, is3D, trim);
                case {'Oxdot'}
                    [X,Y,Z] = obj.plot_ang_velocity(isReal, 'x', trim);
                case {'Oydot'}
                    [X,Y,Z] = obj.plot_ang_velocity(isReal, 'y', trim);
                otherwise
                    if contains(component, 'd')
                        [X,Y,Z] = obj.PowerFlowResultants.plot_derivatives(component);
                    else
                        warning('Input not valid')
                    end
            end 
        end
        function plot_power_flow(obj)
            [Y, X] = meshgrid(obj.eta, obj.xi);
            quiver(X, Y, obj.qx, obj.qy, 'color', [1 0 .2]);
        end
        function plot_power_flow_rings(obj, centerX, centerY, levels, tol)
            [Y, X] = meshgrid(obj.eta, obj.xi);
            Y = Y-centerY;
            X = X-centerX;
            cDist = (Y.^2 + X.^2).^(1/2);
            validIndxs = zeros(size(Y));
            for k = 1:length(levels)
                validIndxs = validIndxs | abs(cDist-levels(k))< tol;
            end
            quiver(X(validIndxs), Y(validIndxs), obj.qx(validIndxs), obj.qy(validIndxs))
        end
        function plot_power_flow_part(obj, part)
            switch part
                case {'bending', 'b', 'bend'}
                    qxp = obj.qxbend;
                    qyp = obj.qybend;
                case {'twisting', 't', 'twist'}
                    qxp = obj.qxtwist;
                    qyp = obj.qytwist;
                case {'shear', 's'}
                    qxp = obj.qxshear;
                    qyp = obj.qyshear;
            end
            [Y, X] = meshgrid(obj.eta, obj.xi);
            quiver(X, Y, qxp, qyp, 'color', [1 0 .2]);
            title(part)
        end
        function plot_displacement_grayscale(obj, unitTheta)
            Xi = linspace(0, 2*pi, 101);
            Xi = Xi(1:end-1);
            Eta = linspace(0, obj.cylinderProperties.height, 80);
            RRadDisp = obj.fieldEvaluator.evaluate_spline_number_at_parameters(Xi, Eta, 2, [0 0]);
            IRadDisp = obj.fieldEvaluator.evaluate_spline_number_at_parameters(Xi, Eta, 5, [0,0]);
            RIRadDisp = RRadDisp+1i*IRadDisp;
%             phases = angle(RIRadDisp(:,2:end-1));
            amp =abs(RIRadDisp);
%             RadDisp = RRadDisp;
            [Z, Theta] = meshgrid(Eta,Xi);
            if unitTheta 
                Theta = obj.cylinderProperties.radius*Theta;
            end
            colormap('gray')
%             s = surf(Theta, Z, zeros(size(amp)), RadDisp);
            s = surf(Theta, Z, zeros(size(amp)), amp);
            s.EdgeAlpha = 0;
            s.FaceAlpha = .2;
            s.FaceColor = 'interp';
        end  
        %% set spline parameters
        function set_control_points(obj, ptsInxi,ptsIneta)
            obj.numControlPoints = [ptsInxi,ptsIneta];
            obj.dirty = true;
        end
        function set_number_of_patches(obj, xpatch, ypatch)
            obj.set_control_points(xpatch+5, ypatch+5);
        end
        function set_number_evaluation_points(obj, xi, eta)
            obj.numberOfEvaluationPnts = [xi,eta];
            obj.dirty = true;
        end
        function set_number_of_repeated_knots(obj, first, last)
            if nargin < 3
                last = first;
            end
            obj.repeatedNodes = [first, last];
        end
    end
    methods (Hidden = true)
        function obj = process_input_data(obj)
            % This function gets the data from the SLDVdata object the parser 
            % fits spline surface to the data. After these operations are complete it
            % creates object to evaluate the derivatives of the splines and
            % calculates the force resultants
            obj.fit_spline_surfaces();
            obj.set_up_surface_and_derivative_evaluators();
            obj.create_surface_plotter();
            obj.create_resultant_calculator();
            obj.dirty = false;
        end
        function plot_ang_velocity(obj, isReal, ang, trim)
            [~, Ox, Oy] = obj.calculate_velocities();
            switch ang
                case {'x'}
                    values = Ox;
                    str1 = '\theta_x dot';
                case {'y'}
                    values = Oy;
                    str1 = '\theta_y dot';
            end
            if isReal
                values = real(values);
                str2 = ' Real';
            else
                values = imag(values);
                str2 = ' Imag';
            end
            if trim 
                values = values(2:end-1, 2:end-1);
            end
            obj.PowerFlowResultants.plot_values(values, true, trim);
            title(strcat(str1, str2));
            xlabel('x')
            ylabel('y')
        end
        function perform_power_flow_computations(obj, xi, eta)
            % this function calculates power flow by evaluating the
            % force resultants and the complex conjugates of the velocites
            if nargin < 3
                if isempty(obj.xi) || isempty(obj.eta)
                    obj.calculate_xi_and_eta();
                end
                xi = obj.xi;
                eta = obj.eta;
            end
            % this gets the resultants from the resultant calculator object
            [Mx, My, Mxy, Qx, Qy, obj.PowerFlowResultants] = obj.PowerFlowResultants.calculate_resultants(xi, eta);
            
            % this calculates the conjugates of the various velocity components
            [wdot, Odotx, Odoty] = obj.calculate_conjugate_velocities(xi,eta);
            
            obj.qxshear = 1/2 * real(Qx.*wdot);
            obj.qxbend = -1/2 * real(Mx.*Odoty) ;
            obj.qxtwist = -1/2 * real(Mxy.*Odotx);
            obj.qyshear = 1/2 * real(Qy.*wdot);
            obj.qybend = -1/2 * real(My.*Odotx);
            obj.qytwist = -1/2 * real(Mxy.*Odoty);
            fullqx = obj.qxshear + obj.qxbend + obj.qxtwist;
            fullqy = obj.qyshear + obj.qybend + obj.qytwist;
            % save both the real and imaginary components of the power flow
            obj.qx = 1/2*real(fullqx);
            obj.qy = 1/2*real(fullqy);
            obj.iqx = 1/2*imag(fullqx);
            obj.iqy = 1/2*imag(fullqy);
        end
        function [wdot, Odotx, Odoty] = calculate_velocities(obj,eta, xi)
            if nargin < 3
                [wdot, Odotx, Odoty] = obj.PowerFlowResultants.get_velocities();
            else
                [wdot, Odotx, Odoty] = obj.PowerFlowResultants.get_velocities(xi,eta);
            end
        end
        function [wdotS, OdotxS, OdotyS] = calculate_conjugate_velocities(obj, xi,eta)
            if nargin < 3
                [wdotS, OdotxS, OdotyS] = obj.PowerFlowResultants.get_conjugate_velocities();
            else
                [wdotS, OdotxS, OdotyS] = obj.PowerFlowResultants.get_conjugate_velocities(xi,eta);
            end
        end
        function calculate_xi_and_eta(obj)
            % this calculates equally spaced evaluation parameter vectors
            % used in calculating and visualizing power flow
            obj.xi = linspace(0, obj.fieldEvaluator.maxXi, obj.numberOfEvaluationPnts(1));
            obj.eta = linspace(0, obj.fieldEvaluator.maxEta, obj.numberOfEvaluationPnts(2));
        end
        function fit_spline_surfaces(obj)
            [params, data, types] = obj.laserData.get_flat_spline_fitting_data();
            for i = 1:2
                if strcmp(types{i}, 'closed')
                    params(:,i) = wrapTo2Pi(params(:,i));
                end
            end
            %create fitter object
            obj.splineFitter = quinticBSplineSurfaceFitter(params, data, types, obj.numControlPoints, obj.repeatedNodes);
            obj.splineFitter.fit_spline_surfaces();
        end
        % create other objects
        function create_resultant_calculator(obj)
            obj.PowerFlowResultants = Plate_Resultants(obj.fieldEvaluator, obj.plateProperties, obj.xi, obj.eta);
        end
        function create_surface_plotter(obj)
            obj.splinePlotter = quinticSplineSurfacePlotter(obj.fieldEvaluator, obj.xi, obj.eta);
        end
        function set_up_surface_and_derivative_evaluators(obj)
            obj.fieldEvaluator = obj.splineFitter.output_solved_spline_evaluator();
            obj.calculate_xi_and_eta();
        end
    end
    methods (Static = true)
        function value = derivativeNumberfromText(txt)
            switch txt
                case {'dx'}
                    value = [0 1];
                case {'d0'}
                    value = [1 0];
                case {'dx2'}
                    value = [0 2];
                case {'dxd0'}
                    value = [1 1];
                case {'d02'}
                    value = [2 0];
                case {'dx3'}
                    value = [0 3];
                case {'dx2d0'}
                    value = [1 2];
                case {'dxd02'}
                    value = [2 1];
                case {'d03'}
                    value = [3 0];
                otherwise
                    value = [0 0];
            end
        end
    end
end