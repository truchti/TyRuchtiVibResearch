classdef Cylindrical_Power_Flow <handle
    % This class can take in a SLDV_parser, a SLDV_data or no input in the
    % no input case it creates a default SLDV_FastScan parser. If a parser
    % is input then the data is extracted using the parser. If the raw data
    % is processed separately the processed data can be passed instead of a
    % parser
    properties
        originalData
        cylindricalData
        numControlPoints = [12 10]
        numberOfEvaluationPnts = [30 20];
        cylinderProperties
        qtheta
        qlong
    end
    properties %(Hidden = true)
        % evaluation properties
        xi
        eta
        % frequency that the cylinder was forced at
        freq
        % object that allows for evaluating the various surface splines at
        % a give parameter pair
        fieldEvaluator
        % object that is used to solve for the surface splines that best
        % fit the data
        splineFitter
        % object that parameterizes the rectangular coordinates into
        % cylindrical coordinates that is it changes the surface from being
        % f(x,y,z) to one that is represented as f(theta, z)
        rectToCylConverter
        % object that is used to plot the splines in cartesian coordinates
        % as function of the parameters
        splinePlotter
        % allows for ploting of the splines in the cylindrical shape
        cylSplinePlotter
        % object that is used to calculate the resultant forces based on
        % love shell theory
        PowerFlowResultants
        % provides plotting functionallity for the power flow data
        PFPlotter
        
        dirty = true;
        iqtheta
        iqlong
    end
    methods
        function obj = Cylindrical_Power_Flow(SLDV_data, cylinderProps)
            if isa(SLDV_data, 'SLDV_data')
                obj.originalData = SLDV_data;
            else
                error('Not Yet Implemented')
            end
            if isa(cylinderProps, 'PFCylProperties')
                obj.cylinderProperties = cylinderProps;
            else
                warning('Cylinder Properties is an incorrect type or missing.\n')
            end
        end
        function [qt, ql] = calculate_power_flow(obj)
            % This function uses the properties and data that it is
            % initialized with to compute the power flow in a cylinder
            % it does this by 
            % 1. converting the input data into the necessary coordinate
            % system
            % 2. computing the tangential and longitudinal components of
            % power flow
            % it then creates a plotting object with the results so they
            % can be readily visualized
            if obj.dirty
                obj.process_input_data();
            end
            obj.compute_cylindrical_power_flow();
            if nargout == 1
                qt = {obj.qtheta, obj.qlong};
            end
            if nargout == 2
                qt = obj.qtheta;
                ql = obj.qlong;
            end
            obj.create_power_flow_plotter();
        end
        function plot(obj, unitTheta, forceLocs)
            if nargin< 3
                forceLocs = [];
            end
            if nargin < 2
                unitTheta = true;
            end
            if ~isempty(obj.PFPlotter)
                 obj.plot_displacement_grayscale(logical(unitTheta));
                 hold on
                obj.PFPlotter.plot_flat(unitTheta*obj.cylinderProperties.radius, forceLocs);
                view(0,90); % make it so that the plane is perpendicular to camera
                hold off
                xlabel("Distance around Circumfrence")
                ylabel("Height")
                frame_h = gcf();
                frame_h.WindowState = 'maximized';
            end
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
            s.FaceAlpha = .9;
            s.FaceColor = 'interp';
        end  
        function plot3d(obj)
            if ~isempty(obj.PFPlotter)
                obj.PFPlotter.plot_as_cylinder(obj.cylinderProperties.radius, obj.cylinderProperties.height);
            end
        end
        function plot_both(obj)
            subplot(1,2,1)
            obj.plot();
            subplot(1,2,2)
            obj.plot3d();
        end
        function probe_disp(obj)
%             type = 'disp';
%             isImaginary = false;
%             dimension = 2;
%             derivative = [0 0];
%             values = obj.fieldEvaluator.evaluate_at_parameters(obj.xi, obj.eta, type, isImaginary, dimension, derivative);
            obj.PFPlotter.plot_flat
        end
        function probe_derivatives(obj, order)
            type = 'disp';
            isImaginary = false;
            dimension = 2;
            derivative = [0 0];
            figure('units','normalized','outerposition',[0 .55 .5 .45])
            obj.splinePlotter.plot_as_cylinder(type, isImaginary, dimension, derivative)
            title('r')
            switch order
                case 1                    
                    derivative = [1 0];
                    figure('units','normalized','outerposition',[.5 .55 .5 .45])
                    obj.splinePlotter.plot_as_cylinder(type, isImaginary, dimension, derivative)
                    title('dr/d\xi')              
                    
                    derivative = [0 1];
                    figure('units','normalized','outerposition',[0 .1 .5 .45])
                    obj.splinePlotter.plot_as_cylinder(type, isImaginary, dimension, derivative)
                    title('dr/d\eta')
                case 2
                    derivative = [2 0];
                    figure('units','normalized','outerposition',[.5 .55 .5 .45])
                    obj.splinePlotter.plot_as_cylinder(type, isImaginary, dimension, derivative)
                    title('d^2r/d\xi^2')
                    derivative = [1 1];
                    figure('units','normalized','outerposition',[.5 .1 .5 .45])
                    obj.splinePlotter.plot_as_cylinder(type, isImaginary, dimension, derivative)
                    title('d^2r/d\xid\eta')
                    derivative = [0 2];
                    figure('units','normalized','outerposition',[0 .1 .5 .45])
                    obj.splinePlotter.plot_as_cylinder(type, isImaginary, dimension, derivative)
                    title('d^2r/d\eta^2')
                case 3
                    derivative = [3 0];
                    figure('units','normalized','outerposition',[0 .55 .5 .45])
                    obj.splinePlotter.plot_as_cylinder(type, isImaginary, dimension, derivative)
                    title('d^3r/d\xi^3')
                    derivative = [2 1];
                    figure('units','normalized','outerposition',[.5 .55 .5 .45])
                    obj.splinePlotter.plot_as_cylinder(type, isImaginary, dimension, derivative)
                    title('d^3r/d\xi^2d\eta')
                    derivative = [1 2];
                    figure('units','normalized','outerposition',[.5 .1 .5 .45])
                    obj.splinePlotter.plot_as_cylinder(type, isImaginary, dimension, derivative)
                    title('d^3r/d\xid\eta^2')
                    derivative = [0 3];
                    figure('units','normalized','outerposition',[0 .1 .5 .45])
                    obj.splinePlotter.plot_as_cylinder(type, isImaginary, dimension, derivative)
                    title('d^3r/d\eta^3')
            end
        end
        
        function set_control_points(obj, xi,eta)
            obj.numControlPoints = [xi,eta];
            obj.dirty = true;
        end
        function set_evaluation_points(obj, xi, eta)
            obj.numberOfEvaluationPnts = [xi,eta];
            obj.dirty = true;
        end
    end
    methods (Hidden = true)
        function obj = process_input_data(obj)
            % This function gets the data from the SLDVdata object the parser creates
            % and converts it to cylindrical coordinates, fits spline
            % surface to the data. After these operations are complete it
            % creates object to evaluate the derivatives of the splines and
            % calculate the force resultants
            obj.convert_to_cylindrical_coordinates();
            obj.fit_spline_surfaces();
            obj.set_up_surface_and_derivative_evaluators();
            obj.create_surface_plotter();
            obj.create_resultant_calculator();
            obj.dirty = false;
        end
        function compute_cylindrical_power_flow(obj, xi, eta)
            % this function calculates power flow by evaluating the
            % force resultants and the complex conjugates of the velocites
            if nargin < 3
                if isempty(obj.xi) || isempty(obj.eta)
                    obj.calculate_xi_and_eta();
                end
                xi = obj.xi;
                eta = obj.eta;
            end
            % this gets the 10 resultants from the resultant calculator object
            [Nt, Nl, Ntl, Nlt, Mt, Ml, Mtl, Mlt, Qt, Ql, obj.PowerFlowResultants] = obj.PowerFlowResultants.calculate_resultants(xi, eta);
            [Eta, Xi] = meshgrid( obj.eta, obj.xi);
            surf(Xi, Eta, real(Qt))
            % this calculates the conjugates of the various velocity
            % components
            [vtheta, vrad, vlong, omegatheta, omegalong] = obj.calculate_conjugate_velocities(xi,eta);
            fullqtheta = -Qt.*vrad + Nt.*vtheta + Ntl.*vlong +Mtl.*omegalong + Mt.*omegatheta;
            fullqlong = -Ql.*vrad + Nlt.*vtheta + Nl.*vlong + Ml.*omegalong + Mlt.*omegatheta;
            % save both the real and imaginary components of the power flow
            obj.qtheta = 1/2*real(fullqtheta);
            obj.qlong = 1/2*real(fullqlong);
            obj.iqtheta =1/2*imag(fullqtheta);
            obj.iqlong = 1/2*imag(fullqlong);
        end
        function  convert_to_cylindrical_coordinates(obj)
                obj.rectToCylConverter = CylCoorParameterizer(obj.originalData);
                obj.rectToCylConverter.calculate_cylindrical_data();
                obj.cylindricalData = obj.rectToCylConverter.export_cylindrical_data;
        end
        function [vt, vr, vz, omegat, omegaz] = calculate_conjugate_velocities(obj, xi,eta)
            if nargin < 2
                xi = obj.xi;
                eta = obj.eta;
            end
            %tang velocity is 1st component of velocity no derivative
            vt = conj(obj.fieldEvaluator.evaluate_complex_value_at_parameters(xi, eta, 'vel', 1, [0 0]));
            % radial velocity is 2nd componentn of velocity no derivative
            vr = conj(obj.fieldEvaluator.evaluate_complex_value_at_parameters(xi, eta, 'vel', 2, [0 0]));
            % long velocity is 3rd component of velocity no derivative
            vz = conj(obj.fieldEvaluator.evaluate_complex_value_at_parameters(xi, eta, 'vel', 3, [0 0]));
            % angluar velocity abt theta is the spatialderivative wrt long of the theta velocity
            omegat = conj(obj.fieldEvaluator.evaluate_complex_value_at_parameters(xi, eta, 'vel', 1, [0 1]));
            % angular velocity abt long i sthe spatialderiavtive wrt theta of the long velocity
            omegaz = conj(obj.fieldEvaluator.evaluate_complex_value_at_parameters(xi, eta, 'vel', 3, [1 0]));
        end
        function calculate_xi_and_eta(obj)
            % this calculates equally spaced evaluation parameter vectors
            % used in calculating and visualizing power flow
            obj.xi = linspace(0, obj.fieldEvaluator.maxXi, obj.numberOfEvaluationPnts(1));
            obj.eta = linspace(0, obj.fieldEvaluator.maxEta, obj.numberOfEvaluationPnts(2));
        end
        function  fit_spline_surfaces(obj)
            [params, data, types] = obj.cylindricalData.get_cylindrical_spline_fitting_data();
            for i = 1:2
                if strcmp(types{i}, 'closed')
                    params(:,i) = wrapTo2Pi(params(:,i));
                end
            end
            %create fitter object
            obj.splineFitter = quinticBSplineSurfaceFitter(params, data, types, obj.numControlPoints);
            obj.splineFitter.fit_spline_surfaces();
%             obj.splineFitter.check_by_viewing_spline_and_data();
        end
        % create other objects
        function create_resultant_calculator(obj)
            obj.PowerFlowResultants = CylindricalPowerFlowResultants(obj.fieldEvaluator, obj.cylinderProperties, obj.xi, obj.eta);
        end
        function create_power_flow_plotter(obj)
            obj.PFPlotter = PowerFlowDataPlotter(obj.xi, obj.eta, obj.qtheta, obj.qlong);
        end
        function create_surface_plotter(obj)
            obj.splinePlotter = quinticSplineCylSurfacePlotter(obj.fieldEvaluator, obj.xi, obj.eta, obj.cylinderProperties.radius);
        end
        function  set_up_surface_and_derivative_evaluators(obj)
            obj.fieldEvaluator = obj.splineFitter.output_solved_spline_evaluator();
            obj.calculate_xi_and_eta();
        end
    end
end