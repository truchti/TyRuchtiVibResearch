classdef PowerFlowCylindrical <handle
    % This class can take in a SLDV_parser, a SLDV_data or no input in the
    % no input case it creates a default SLDV_FastScan parser. If a parser
    % is input then the data is extracted using the parser. If the raw data
    % is processed separately the processed data can be passed instead of a
    % parser
    %cylindrical coordinates are in the format of (Radial, Theta, Z)
    % it will then fit the cylinder data and create a
    % chebyshev-polynomial/quinticspline tensor surface representation of
    % the surface displacement and velocity fields. These fields can then
    % be differentiated to and used to compute the power. The flow can then
    % be visualized using the animate or plot functionality
    properties
        RawData SLDV_data
        CylProperties PFCylProperties
        CylData SLDV_data_cyl_coors
        ChebOrder 
        QuinticPatches 
    end
    properties %(Hidden = true)
        Surfaces = QuintChebSplineSurface.empty    
        xi 
        eta
        qtheta
        qlong
        numEvalPts = [100, 100];
    end
    properties (Hidden = true)
        processed= false; 
        surfacesFit = false;
    end
    methods
        function obj = PowerFlowCylindrical(SLDV_data, cylinderProps, quinticPatches, chebOrder)
            obj.RawData = SLDV_data;
            obj.CylProperties = cylinderProps;
            
            if nargin < 4 % default Cheb polynomial order
                chebOrder = 30; 
            end
            if nargin < 3 %default number of quintic patches
                quinticPatches = 13; 
            end
            obj.QuinticPatches = quinticPatches;
            obj.ChebOrder = chebOrder;
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
            if ~obj.processed
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
        end
        %% plotting functions
        function plot_component(obj, component, isReal)
            if nargin < 3
                isReal = true;
            end
            switch component
                case {'t'; 'theta'; 'r'; 'rad'; 'l'; 'long', 'tdot'; 'thetadot'; 'rdot'; 'raddot'; 'ldot'; 'longdot'}
                    obj.calc_surface_values(component, isReal); % if the component is just a surface 
                case {'dwd0', 'dwdz', 'd2wd02', 'd2wd0z', 'd2wdz2', 'd3wd03', 'd3wd02z', 'd3wd0z2', 'd3wdz3'}
                    obj.calc_surface(component, isReal); % if component is deriviative of surface
                case {'Nt'; 'Nl'; 'Ntl'; 'Nlt'; 'Mt'; 'Ml'; 'Mtl'; 'Mlt'; 'Qt'; 'Ql'}
                    obj.PowerFlowResultants.plot_resultant(component, isReal)
                case {}
                    obj.PowerFlowResultants.plot_velocity(component, isReal);
                otherwise
                    warning('Input not valid')
            end            
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
            else
                warning("Power flow has not been calculated")
            end
        end
        function plot_part_of_power_flow(obj, part)
            if ~isempty(obj.PFPlotter)
            obj.plot_displacement_grayscale(1);
            hold on
            [Z, Theta] = meshgrid(obj.eta, obj.xi);
            ThetaDist =  obj.cylinderProperties.radius* Theta;
            switch part 
                case {1, 'shear', 'Q', 'Shear'}
                    U = real(obj.shearPowerTang);
                    V = real(obj.shearPowerLong);
                case {2, 'normal', 'Normal', 'N'}
                    U = real(obj.NormalPowerTang);
                    V = real(obj.NormalPowerLong);
                case {3, 'moment', 'Moment', 'M'}
                    U = real(obj.MomentPowerTang);
                    V = real(obj.MomentPowerLong);
            end
            quiver(ThetaDist, Z, U,V)
            view(0,90); % make it so that the plane is perpendicular to camera
            hold off
            xlabel("Distance around Circumfrence")
            ylabel("Height")
            frame_h = gcf();
            frame_h.WindowState = 'maximized';
            else
                warning("Power flow has not been calculated")
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
            s.FaceAlpha = .2;
            s.FaceColor = 'interp';
        end  
%         function plot3d(obj)
%             
%         end
%         function plot_both(obj)
%             subplot(1,2,1)
%             obj.plot();
%             subplot(1,2,2)
%             obj.plot3d();
%         end
    end
    methods (Hidden = true)
        function obj = process_input_data(obj)
            % This function gets the data from the SLDVdata object the parser creates
            % and converts it to cylindrical coordinates, fits spline
            % surface to the data. After these operations are complete it
            % creates object to evaluate the derivatives of the splines and
            % calculate the force resultants
            if ~obj.processed
                obj.convert_to_cylindrical_coordinates();
                obj.processed = true;
            end
            if ~obj.surfacesFit
                obj.fit_surfaces();
            end
        end
        function  fit_surfaces(obj)
            [params, data] = obj.CylData.get_cylindrical_spline_fitting_data();
            Ths = wrapTo2Pi(params(:,1));
            Zs = params(:,2);
            %create fitter object locally
            fitter = QuintChebSplineSurfaceFitter(Ths, Zs, data, obj.QuinticPatches, obj.ChebOrder);
            fitter.fit_spline_surfaces();
            obj.Surfaces = fitter.output_surfaces();
            obj.surfacesFit = true;
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
%             [Eta, Xi] = meshgrid( obj.eta, obj.xi)
            % this calculates the conjugates of the various velocity
            % components
            [vtheta, vrad, vlong, omegatheta, omegalong] = obj.calculate_conjugate_velocities(xi,eta);
            obj.shearPowerLong = -Ql.*vrad;
            obj.NormalPowerLong = Nlt.*vtheta + Nl.*vlong;
            obj.MomentPowerLong = Ml.*omegatheta - Mlt.*omegalong;
            obj.shearPowerTang = -Qt.*vrad ;
            obj.NormalPowerTang = Nt.*vtheta + Ntl.*vlong;
            obj.MomentPowerTang = Mtl.*omegatheta - Mt.*omegalong;
            fullqtheta = -Qt.*vrad + Nt.*vtheta + Ntl.*vlong +Mtl.*omegatheta - Mt.*omegalong;
            fullqlong = -Ql.*vrad + Nlt.*vtheta + Nl.*vlong + Ml.*omegatheta - Mlt.*omegalong;
            % save both the real and imaginary components of the power flow
            obj.qtheta = 1/2*real(fullqtheta);
            obj.qlong = 1/2*real(fullqlong);
            obj.iqtheta =1/2*imag(fullqtheta);
            obj.iqlong = 1/2*imag(fullqlong);
        end
        function  convert_to_cylindrical_coordinates(obj)
            % converter created locally since it does not need to persist
            converter =  CylCoorParameterizer(obj.RawData);
            converter.calculate_cylindrical_data();
            obj.CylData = converter.export_cylindrical_data;
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
            omegaz = 1/obj.cylinderProperties.radius * conj(obj.fieldEvaluator.evaluate_complex_value_at_parameters(xi, eta, 'vel', 3, [1 0]));
        end
        function calculate_xi_and_eta(obj)
            % this calculates equally spaced evaluation parameter vectors
            % used in calculating and visualizing power flow
            obj.xi = linspace(0, obj.fieldEvaluator.maxXi, obj.numEvalPts(1));
            obj.eta = linspace(0, obj.fieldEvaluator.maxEta, obj.numEvalPts(2));
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
                case {'rdot'}
                    value = [0 0];
                otherwise
                    value = [0 0];
            end
        end
    end
end