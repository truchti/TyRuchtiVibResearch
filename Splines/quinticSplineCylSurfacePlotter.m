classdef quinticSplineCylSurfacePlotter < quinticSplineSurfacePlotter
    properties
        radius
    end
    methods
        function obj = quinticSplineCylSurfacePlotter(evaluator, xiPoints, etaPoints, radius)
                obj@quinticSplineSurfacePlotter(evaluator, xiPoints, etaPoints);
                obj.radius = radius;
        end
        function plot_as_cylinder(obj, type, isImaginary, dimension, derivative)
            if nargin < 5
                type = 'disp';
                isImaginary = 'false';
                dimension = 2;
                derivative = 0;
            end
            obj.calculate_plotting_points
            values = obj.evaluate_spline_at_plotting_points(type, isImaginary, dimension, derivative);
            [Thetas, Heights] = meshgrid(obj.plottingPoints_xi, obj.plottingPoints_eta);
            xs = obj.radius*cos(Thetas);
            ys = obj.radius*sin(Thetas);
            zs = Heights;
            s = surf(xs, ys, zs, values');
            s.EdgeAlpha =0;
            colorbar
            xlabel('X')
            ylabel('Y')
            zlabel('Z')
        end
    end
    
end