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
            obj.calculate_plotting_points
            values = obj.evaluate_spline_at_plotting_points(type, isImaginary, dimension, derivative);
            [Thetas, Heights] = meshgrid(obj.plottingPoints_xi, obj.plottingPoints_eta);
            xs = obj.radius*cos(Thetas);
            ys = obj.radius*sin(Thetas);
            zs = Heights;
            surf(xs, ys, zs, values')
            colorbar
            xlabel('X')
            ylabel('Y')
            zlabel('Z')
        end
    end
    
end