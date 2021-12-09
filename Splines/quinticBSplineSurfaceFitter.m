classdef quinticBSplineSurfaceFitter < matlab.mixin.Copyable
    % This class takes in xi and eta parameter vectors and a 2D array of data
    % values at each of these parameter pairs. It contains methods to
    % fit a spline surface to the data and to export a custom spline object
    properties
        xiParams
        etaParams
        data
        types
        controlNetDimensions
        xiKnot
        etaKnot
    end
    properties (Hidden = true)
        splineEvaluator
        solved = false
    end
    methods
        function obj = quinticBSplineSurfaceFitter(params, data, types, numberControlPoints, knotVectors)
            obj.xiParams = params(:,1);
            obj.etaParams = params(:,2);
            obj.data = data;
            obj.types = types;
            obj.controlNetDimensions = numberControlPoints;
            if nargin > 4
                obj.xiKnot = knotVectors{1};
                obj.etaKnot = knotVectors{2};
            else
                obj.calculate_knot_vectors;
            end
        end
        function [xiB, etaB] = fit_spline_surfaces(obj)
            %% Set up evaluator with all zero control Points
            obj.splineEvaluator = quinticSplineSurfaceEvaluator({obj.xiKnot, obj.etaKnot}, obj.controlNetDimensions, obj.types);
            %evaluate basis functions
            [xiB, etaB] = obj.splineEvaluator.evaluate_basis_functions(obj.xiParams, obj.etaParams);
            for i = size(obj.data,2):-1:1
                    obj.splineEvaluator.controlPoints(:,:,i) = numel(obj.data(:,i))*(xiB\diag(obj.data(:,i))/(etaB'));
            end
            obj.solved = true;
        end  
        function Evaluator = output_solved_spline_evaluator(obj)
            if ~obj.solved
                obj.fit_spline_surfaces;
            end
            Evaluator = copy(obj.splineEvaluator);
        end
        function check_by_viewing_spline_and_data(obj)
            xi = linspace(0, max(obj.xiKnot), 100);
            eta = linspace(0, max(obj.etaKnot), 50);
            [ETA, XI] = meshgrid(eta, xi);
            for n = 1:12
                figure(n)
                values = obj.splineEvaluator.evaluate_spline_number_at_parameters( xi, eta, n);
                scatter3(obj.xiParams, obj.etaParams, obj.data(:,n))
                hold on;
                surf(XI,ETA, values)
                hold off;
            end
        end
    end
    methods (Hidden = true)
        function calculate_knot_vectors(obj)
            obj.xiKnot = obj.calculate_knot_vector_using_type(obj.controlNetDimensions(1), obj.types{1}, max(obj.xiParams));
            obj.etaKnot = obj.calculate_knot_vector_using_type(obj.controlNetDimensions(2), obj.types{2}, max(obj.etaParams));
        end
        function knotVect = calculate_knot_vector_using_type(obj, numPts, type, maxScale)
            if strcmp(type, 'open')
                knotVect = obj.open_loop_knot_vector(numPts, maxScale);
            elseif strcmp(type, 'closed')
                knotVect = obj.closed_loop_knot_vector(numPts,maxScale);
            else
                warning("Invalid Type");
                knotVect = Nan;
            end
        end
    end
    methods (Static = true, Hidden = true)
        function knotVector = closed_loop_knot_vector(npts, units)
            knotVector = (0:6+npts-1)-5;
            if nargin > 1 && units
                knotVector = knotVector./(knotVector(end))* units;
            end
        end
        function knotVector = open_loop_knot_vector(numPts, unit)
            if nargin < 2
                unit = 0;
            end
            npc = numPts+ 6;
            np2 = numPts+2;
            knotVector = zeros(1,npc);
            for i = 2:npc
                if i>6 && i<np2
                    knotVector(i) = knotVector(i-1) + 1;
                else
                    knotVector(i) = knotVector(i-1);
                end
            end
            if unit
                knotVector = knotVector./max(knotVector)*unit;
            end
        end
    end
end