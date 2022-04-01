classdef Planar_Power_Flow <handle
    % This class can take in a SLDV_data object  
    properties
        laserData
        plateProperties
        numberOfEvaluationPnts = [50 55];
        orders = [8,8]
        qx
        qy
    end
    properties (Hidden = true)
        xPtsCurrent
        yPtsCurrent
        surfaces
        dirty = true;
        qxshear
        qxbend
        qxtwist
        qyshear
        qybend
        qytwist
    end
    properties (Dependent)
        D  
    end
    methods
        function obj = Planar_Power_Flow(SLDV_data, properties)
            if isa(SLDV_data, 'SLDV_data')
                obj.laserData = SLDV_data;
            else
                error('Not Yet Implemented')
            end
            if nargin>1 
                obj.plateProperties = properties;
            end
        end
        function [qx, qy] = calculate_power_flow(obj, xpts, ypts)
            % This function uses the properties and data that it is
            % initialized with to compute the power flow in the x and y directions 
            % it then creates a plotting object with the results so they can be readily visualized
            obj.process_input_data();
            if nargin < 2
                [xpts, ypts] = obj.default_xy_points();
            end
            obj.perform_power_flow_computations(xpts, ypts);
            if nargout == 1
                qx = {obj.qx, obj.qy};
            end
            if nargout == 2
                qx = obj.qx;
                qy = obj.qy;
            end
        end
        function best = optimize_orders(obj,analytical,xpts, ypts)
            h = figure;
            a = subplot(2,2,1);
            [X,Y,dy3true] = analytical.show_3D('dwdy3');
            b = subplot(2,2,3);
            [~, ~, dx3true] = analytical.show_3D('dwdx3');
            bestErr = 1e3;
            c = subplot(2,2,2);
            d = subplot(2,2,4);
            Link = linkprop([a,b,c,d],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim', 'CLim'});
            setappdata(gcf, 'StoreTheLink', Link);
            for ord = 6:2:36                
                obj.set_order_of_fit(ord)
                subplot(2,2,2)
                dy3 = obj.plot_component('dwdy3', true, xpts, ypts, false);
                subplot(2,2,4)
                dx3 = obj.plot_component('dwdx3', true, xpts, ypts, false);
                resX = dx3true(2:end-1)-dx3(2:end-1);
                resY = dy3true(2:end-1)-dy3(2:end-1);
                errsumX = sum(abs(resX),'all');
                errsumY = sum(abs(resY), 'all');
                if errsumX +errsumY < bestErr
                    bestErr = errsumX + errsumY;
                    best = [ord, ord];
                end
            end
            obj.set_order_of_fit(best);
            subplot(2,2,2)
            obj.plot_component('dwdy3', true, xpts, ypts, true);
            subplot(2,2,4)
            obj.plot_component('dwdx3', true, xpts, ypts, true);
        end
        function evaluate_orders(obj, orders)
            qt = true;
            while(qt)
                obj.set_order_of_fit(orders);
                obj.fit_surfaces
                h = figure('units','normalized','outerposition',[0 0 1 1]);
                subplot(1,3,1)
                obj.plot_component('w', true)
                subplot(1,3,2)
                obj.plot_component('dwdx3', true)
                subplot(1,3,3)
                obj.plot_component('dwdy3', true)
                
                orders = input("Input vector of new orders to try");
                if isempty(orders) || ~isnumeric(orders)
                    qt = false;
                end
                if (isobject(h) && isgraphics(h))
                    close(h)
                end
            end
        end
        %% plotting functions
        function plot(obj, trim)
            if nargin <2
                trim = 0;
            end
            if isempty(obj.qx)
                obj.calculate_power_flow(); 
            end
            [X,Y] = ndgrid(obj.xPtsCurrent, obj.yPtsCurrent);
            U = obj.qx;  V = obj.qy;
            if trim 
                X = X(trim:end+1-trim, trim:end+1-trim);
                Y = Y(trim:end+1-trim, trim:end+1-trim);
                U = U(trim:end+1-trim, trim:end+1-trim);
                V = V(trim:end+1-trim, trim:end+1-trim);
            end
            quiver(X,Y,U,V)    
        end        
        function val =  plot_component(obj, component, isReal, xpts, ypts, trim)
            if nargin < 6 
                trim = false;
            end
            if nargin < 4
                [xpts, ypts] = obj.default_xy_points();
            end
            if nargin < 3
                isReal = true;
            end
            switch component
                case {'w','wdot'}
                    val = obj.plot_disp_or_vel(component, isReal, xpts, ypts, trim);
                case {'Qx', 'Qy', 'Mx', 'My', 'Mxy'}
                    val = obj.plot_resultant(component, isReal, xpts, ypts, trim);
                otherwise
                    val = obj.plot_derivative(component, isReal, xpts, ypts, trim);
            end
        end
        function val = plot_resultant(obj, component, isReal, xpts, ypts, trim)
            switch component
                case {'Qx'}
                    [val, ~] = obj.calculate_Qx_and_Qy();
                case {'Qy'}
                    [~, val] = obj.calculate_Qx_and_Qy();
                case {'Mx'}
                    [val, ~,~] = obj.calculate_Mx_My_and_Mxy();
                case {'My'}
                    [~,val, ~] = obj.calculate_Mx_My_and_Mxy();
                case {'Mxy'}
                    [~, ~, val] = obj.calculate_Mx_My_and_Mxy();
            end
            if isReal
                val = real(val);
            else
                val = imag(val);
            end
            plot_surface(val, xpts, ypts, trim)
        end
        function val = plot_disp_or_vel(obj, component, isReal, xpts ,ypts, trim)
            switch component
                case {'w'}
                    val = obj.surfaces.evaluate_complex_displacement_derivative_at_points(xpts, ypts);
                case {'wdot'}
                    [val, ~, ~] = obj.calculate_velocities(xpts, ypts);
                case {'Oxdot'}
                    [~, val, ~] = obj.calculate_velocities(xpts, ypts);
                case {'Oydot'}
                    [~, ~, val] = obj.calculate_velocities(xpts, ypts);
            end
            if isReal 
                val = real(val);
            else
                val = imag(val);
            end
            plot_surface(val, xpts, ypts, trim)
        end
        function val = plot_derivative(obj, component, isReal, xpts, ypts, trim)
            [type, derivs] = convert_component_to_type_and_derivs(component, isReal);
            val = obj.surfaces.evaluate_surface_derivative_at_points(type, xpts, ypts, derivs);
            plot_surface(val, xpts, ypts, trim);
        end
        %% set parameters
        function set_order_of_fit(obj, xOrder, yOrder)
            if nargin < 3
                yOrder = xOrder;
            end
            if length(xOrder) >1
                yOrder = [];
            end
            obj.orders = [xOrder, yOrder];
            obj.fit_surfaces();
        end
        function set_number_evaluation_points(obj, xi, eta)
            obj.numberOfEvaluationPnts = [xi,eta];
            obj.dirty = true;
        end
        %% getters
        function val = get.D(obj)
            val = obj.plateProperties.E*obj.plateProperties.thickness^3/(12*(1-obj.plateProperties.poissons^2));
        end
    end
    methods (Hidden = true)
        function obj = process_input_data(obj)
            % This function gets the data from the SLDVdata object the parser 
            % fits spline surface to the data. After these operations are complete it
            % creates object to evaluate the derivatives of the splines and
            % calculates the force resultants
            obj.fit_surfaces();
            obj.dirty = false;
        end
        %% resultant methods
        function [Qx, Qy] = calculate_Qx_and_Qy(obj, xpts, ypts)
            if nargin < 2
                [xpts, ypts] = obj.default_xy_points();
            end
            [dWdx3, dWdx2y, dWdxy2, dWdy3] = obj.get_all_3rd_derivatives(xpts, ypts);
            Qx = obj.D * (dWdx3 + dWdxy2);
            Qy = obj.D * (dWdx2y + dWdy3);
        end
        function [Mx, My, Mxy] = calculate_Mx_My_and_Mxy(obj, xpts, ypts)
            if nargin < 2
                [xpts, ypts] = obj.default_xy_points();
            end
            [dWdx2, dWdxy, dWdy2] = obj.get_all_2nd_derivatives(xpts, ypts);
            nu = obj.plateProperties.poissons;
            Mx = obj.D * (dWdx2 + nu * dWdy2);
            My = obj.D * (dWdy2 + nu * dWdx2);
            Mxy = obj.D * (1-nu) * dWdxy;
        end
        function [dx2, dxy, dy2] = get_all_2nd_derivatives(obj, xpts, ypts)
            dx2 = obj.surfaces.evaluate_complex_displacement_derivative_at_points(xpts, ypts, [2,0]);
            dxy = obj.surfaces.evaluate_complex_displacement_derivative_at_points(xpts, ypts, [1,1]);
            dy2 = obj.surfaces.evaluate_complex_displacement_derivative_at_points(xpts, ypts, [0,2]);
        end
        function [dx3, dx2y, dxy2, dy3] = get_all_3rd_derivatives(obj, xpts, ypts)
            dx3  = obj.surfaces.evaluate_complex_displacement_derivative_at_points(xpts, ypts, [3,0]);
            dxy2 = obj.surfaces.evaluate_complex_displacement_derivative_at_points(xpts, ypts, [1,2]);
            dx2y = obj.surfaces.evaluate_complex_displacement_derivative_at_points(xpts, ypts, [2,1]);
            dy3  = obj.surfaces.evaluate_complex_displacement_derivative_at_points(xpts, ypts, [0, 3]);
        end
        %% velocity methods
        function [wdot, Odotx, Odoty] = calculate_velocities(obj,xpts, ypts)
            wdot = obj.surfaces.evaluate_complex_velocity_derivative_at_points(xpts, ypts, 0);
            Odotx = obj.surfaces.evaluate_complex_velocity_derivative_at_points(xpts, ypts, [0, 1]);
            Odoty = obj.surfaces.evaluate_complex_velocity_derivative_at_points(xpts, ypts, [1, 0]);
        end
        function [wdotS, OdotxS, OdotyS] = calculate_conjugate_velocities(obj, xpts, ypts)
            [wdot, Odotx, Odoty] = obj.calculate_velocities(xpts, ypts);
            wdotS = conj(wdot);
            OdotxS = conj(Odotx);
            OdotyS = conj(Odoty);
        end
        %% compute power flow
        function perform_power_flow_computations(obj, xpts, ypts)
            % this function calculates power flow by evaluating the
            % force resultants and the complex conjugates of the velocites
            [sx, sy] = obj.calculate_shear_power_flow(xpts, ypts);
            [bx, by] = obj.calculate_bending_power_flow(xpts, ypts);
            [tx, ty] = obj.calculate_twisting_power_flow(xpts, ypts);
            obj.qx = sx + bx + tx;
            obj.qy = sy + by + ty;
            % save points that I currently am using for power flow
            obj.xPtsCurrent = xpts;
            obj.yPtsCurrent = ypts;
        end
        function [qxs, qys] = calculate_shear_power_flow(obj,xpts, ypts)
            [Qx, Qy] = obj.calculate_Qx_and_Qy;
            [wdotStar, ~, ~] = obj.calculate_conjugate_velocities(xpts,ypts);
            qxs = 0.5 * real(Qx .* wdotStar);
            qys = 0.5 * real(Qy .* wdotStar);
        end
        function [qxb, qyb] = calculate_bending_power_flow(obj, xpts, ypts)
            [~,Oxdot, Oydot] = obj.calculate_conjugate_velocities(xpts, ypts);
            [Mx, My, ~] = obj.calculate_Mx_My_and_Mxy(xpts, ypts);
            qxb = -0.5 * real(Mx .* Oydot);
            qyb = -0.5 * real(My .* Oxdot);
        end
        function [qxt, qyt] = calculate_twisting_power_flow(obj, xpts, ypts)
            [~,Oxdot, Oydot] = obj.calculate_conjugate_velocities(xpts, ypts);
            [~, ~, Mxy] = obj.calculate_Mx_My_and_Mxy();
            qxt = -0.5 * real(Mxy .* Oydot);
            qyt = -0.5 * real(Mxy .* Oxdot);
        end
        function [xpts, ypts] = default_xy_points(obj)
            if ~isempty(obj.xPtsCurrent)
                xpts = obj.xPtsCurrent;
                ypts = obj.yPtsCurrent;
            else
                xRng = obj.surfaces.wR.Xrange;
                yRng = obj.surfaces.wR.Yrange;
                xpts = linspace(xRng(1), xRng(2), 99);
                ypts = linspace(yRng(1), yRng(2), 101);
            end
        end
        %% fitting functions
        function fit_surfaces(obj)
            [data_x, data_y] = obj.get_x_and_y_parameters_for_fitting();
            fitter = ChebyshevSurfaceFitter(obj.orders, data_x, data_y, []);
            obj.surfaces = PlanarPowerFlowSurfaceStructure();
            for i = 1:4 % for all non-zero displacements and velocities
                vals = obj.get_measured_values_to_fit(i); % get data for specific part
                surf = fitter.fit_surface_to_new_values(vals); % fit surface
                obj.surfaces.set_surface(i,surf);% store surface
            end
        end
        function [xs, ys] = get_x_and_y_parameters_for_fitting(obj)
            xs = obj.laserData.coordinates(:,1);
            ys = obj.laserData.coordinates(:,2);
        end
        function vals = get_measured_values_to_fit(obj, type)
            switch type
                case 1
                    var = 'realDisp';
                case 2
                    var = 'imagDisp';
                case 3
                    var = 'realVel';
                case 4 
                    var = 'imagVel';
            end
            vals = obj.laserData.(var)(:,3);
        end
    end
end
function [type, derivs] = convert_component_to_type_and_derivs(component, isReal)
    if contains(component, 'dwdot', 'IgnoreCase', true)
        t1 = 'wdot';
    elseif contains(component, 'dw', 'IgnoreCase', true)
        t1 = 'w';
    else
        error('invalid component string')
    end
    derivs = derivs_from_string(component);
    if isReal
        t2 = 'R';
    else
        t2 = 'I';
    end
    type = [t1 t2];
end
function derivs = derivs_from_string(component)
    indx = find(ismember(component,'x'),1);
    indy = find(ismember(component, 'y'),1);
    if isempty(indx)
        xdiv = 0;
    elseif isstrprop(component(indx+1),'digit')
        xdiv = str2double(component(indx+1));
    else
        xdiv = 1;
    end
    if isempty(indy)
        ydiv = 0;
    elseif length(component) > indy && isstrprop(component(indy+1),'digit')
        ydiv = str2double(component(indy+1));
    else
        ydiv = 1;
    end
    
    derivs = [xdiv ydiv];
end
function plot_surface(values, xpts, ypts, trim)
    [Y, X] = ndgrid(ypts,xpts);
    tl = 2;
    if trim
        surf(X(tl:end-tl+1, tl:end-tl+1), Y (tl:end-tl+1, tl:end-tl+1), values(tl:end-tl+1, tl:end-tl+1))
    else
        surf(X,Y, values');
    end
end