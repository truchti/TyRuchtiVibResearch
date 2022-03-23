classdef quinticSplineSurfacePlotter < handle
    properties
        splineEvaluator
        numberInterpPoints = [100 101]
    end
    properties (Dependent = true)
        numberOfSurfaces
    end
    properties (Hidden = true)
        plottingPoints_xi
        plottingPoints_eta
        max_xi
        max_eta
    end
    methods
        function obj = quinticSplineSurfacePlotter(evaluator, xiPoints, etaPoints)
            obj.splineEvaluator = evaluator;
            if nargin > 1
                obj.plottingPoints_xi = xiPoints;
                obj.plottingPoints_eta = etaPoints;
                obj.numberInterpPoints = [length(xiPoints), length(etaPoints)];
            end
            obj.calculate_plotting_points;
        end
        function value = get.numberOfSurfaces(obj)
            value = size(obj.splineEvaluator.controlPoints,3);
        end
        function plot_spline(obj, type, isImaginary, dimension, derivative)
            if nargin < 5
                derivative = '';
            end
            if nargin < 4 
                dimension = 2;
            end
            if nargin < 3
                isImaginary = false;
            end
            if nargin < 2 
                type = 'disp';
            end
            if isempty(obj.plottingPoints_xi) || isempty(obj.plottingPoints_eta)
                obj.calculate_plotting_points();
            end
            splineValues = obj.evaluate_spline_at_plotting_points(type, isImaginary, dimension, derivative); 
            surf(obj.plottingPoints_xi, obj.plottingPoints_eta, splineValues')
            xlabel('Width');
            ylabel('Height');
            zlabel('Value');
        end
    end
    methods(Hidden = true)
        function calculate_plotting_points(obj)
            if isempty(obj.max_xi)
                obj.max_xi = max(obj.splineEvaluator.knotVectors{1});
            end
            if isempty(obj.max_eta)
                obj.max_eta = max(obj.splineEvaluator.knotVectors{2});
            end
            obj.plottingPoints_xi = linspace(0, obj.max_xi, obj.numberInterpPoints(1));
            obj.plottingPoints_eta = linspace(0, obj.max_eta, obj.numberInterpPoints(2));
        end
        function values = evaluate_spline_at_plotting_points(obj, type, isImaginary, dimension, derivative)
            if nargin < 5 
                derivative = '';
            end
            if nargin < 4
                dimension = 2;
            end
            if nargin < 3
                isImaginary = false;
            end
            if nargin < 2
                type = 'disp';
            end
            values = obj.splineEvaluator.evaluate_at_parameters(obj.plottingPoints_xi, obj.plottingPoints_eta, type, isImaginary, dimension, derivative);
        end
    end
end