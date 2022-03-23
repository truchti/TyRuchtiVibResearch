classdef PowerFlowResidual < handle
   properties
        simulator
        PFcalc
   end
    methods
        function obj = PowerFlowResidual(simulator, powerFlowCalculator)
            obj.simulator = simulator;
            obj.PFcalc = powerFlowCalculator;
            obj.power_flow_same_points_as_simulation()
        end
        function plot_residual(obj, type, isReal)
            if nargin< 3
                isReal = true;
            end
            trim = false;
            [X,Y,Z] = obj.get_power_flow_value(type, isReal, trim);
            [~, ~, Zt] = obj.get_simulation_value(type, isReal);
            surf(X,Y,Z, Z-Zt)
            colorbar
            title(type)
        end
        function side_by_side_with_residual(obj, type, isReal, normalize)
             compareHandle = figure('units','normalized','outerposition',[0 0 1 1]);
            if nargin< 3
                isReal = true;
            end
            if nargin < 4
                normalize = false;
            end
            trim = false;
            [X,Y,Z] = obj.get_power_flow_value(type, isReal, trim);
            [~, ~, Zt] = obj.get_simulation_value(type, isReal);
            hq = subplot(1,2,1);
            surf(X,Y,Z)
            colormap jet
            colorbar
            title([type, ' Quintic Spline'])
            ha = subplot(1,2,2);
            surf(X,Y,Zt)
            colormap jet
            colorbar
            title([type, ' Analytical'])
            hq.CLim = ha.CLim;
            Link = linkprop([ha,hq],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
            setappdata(gcf, 'StoreTheLink', Link);
            residualHandle = figure()
            residual = Zt-Z;
            titleString = 'Residual';
            if normalize
                residual = residual/(max(max(Zt)));
                titleString = [type, ' ', titleString, ' Normalized'];
            end
            surf(X,Y, zeros(size(X)), residual)
            title(titleString)
            colormap jet
            colorbar
            view(0,90)
            
        end
        function [X, Y, Z] = get_power_flow_value(obj, type, isReal, trim)
            h = figure;
            [X, Y, Z] = obj.PFcalc.plot_component_3D(type, isReal, trim);
            close(h)
        end
        function [X, Y, Z] = get_simulation_value(obj, type, isReal)
            h = figure;
            [X, Y, Z] = obj.simulator.show_3D(type, isReal);
            close(h)
        end
        function power_flow_same_points_as_simulation(obj)
            wPts = obj.simulator.mesh.width_divisions+1;
            hPts = obj.simulator.mesh.height_divisions+1;
            if ~(all(obj.PFcalc.numberOfEvaluationPnts == [wPts,hPts]))
                obj.PFcalc.set_number_evaluation_points(wPts,hPts)
                obj.PFcalc.calculate_power_flow()
            end
        end
    end
end