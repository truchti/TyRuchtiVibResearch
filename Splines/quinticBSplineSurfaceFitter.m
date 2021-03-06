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
        boundaryNodeRepeats
    end
    properties (Hidden = true)
        splineEvaluator
        solved = false
    end
    methods
        function obj = quinticBSplineSurfaceFitter(params, data, types, numberControlPoints, numberOfRepeatedNodes)
            obj.xiParams = params(:,1);
            obj.etaParams = params(:,2);
            obj.data = data;
            obj.types = types;
            obj.controlNetDimensions = numberControlPoints;
            if nargin < 5
                obj.boundaryNodeRepeats = [6,6];
            else
                if length(numberOfRepeatedNodes) == 2
                    obj.boundaryNodeRepeats = numberOfRepeatedNodes;
                else
                    obj.boundaryNodeRepeats = repmat(numberOfRepeatedNodes, [1,2]);
                end
            end
            obj.calculate_knot_vectors();
        end
        function [xiB, etaB] = fit_spline_surfaces(obj)
            %% Set up evaluator with all zero control Points
            obj.splineEvaluator = quinticSplineSurfaceEvaluator({obj.xiKnot, obj.etaKnot}, obj.controlNetDimensions, obj.types);
            %evaluate basis functions
            [xiB, etaB] = obj.splineEvaluator.evaluate_basis_functions(obj.xiParams, obj.etaParams);
            for i = size(obj.data,2):-1:1
                tempData = obj.data(:,i);
                if nnz(obj.data(:,i)) < 1 %% if there is no data in out of plane dimension because it is simulated plate don't calculate
                    obj.splineEvaluator.controlPoints(:,:,i) = zeros(obj.controlNetDimensions);
                else
                    if size(obj.data,1) > 9999
                        tempData = sparse(diag(tempData)); %% if would be too large in memory make a sparse matrix since it is diagonal
                    else
                        tempData = diag(tempData);
                    end
                    obj.splineEvaluator.controlPoints(:,:,i) = numel(obj.data(:,i))*(xiB\tempData/(etaB'));
                end
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
            for n = 1:size(obj.data,2)
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
            obj.xiKnot = obj.calculate_knot_vector_using_type(obj.controlNetDimensions(1), obj.types{1}, min(obj.xiParams), max(obj.xiParams));
            obj.etaKnot = obj.calculate_knot_vector_using_type(obj.controlNetDimensions(2), obj.types{2}, min(obj.etaParams), max(obj.etaParams));
        end
        function knotVect = calculate_knot_vector_using_type(obj, numPts, type, minScale, maxScale)
            if strcmp(type, 'open')
                knotVect = obj.open_loop_knot_vector(numPts, minScale, maxScale, obj.boundaryNodeRepeats(1), obj.boundaryNodeRepeats(2));
            elseif strcmp(type, 'closed')
                knotVect = obj.closed_loop_knot_vector(numPts, minScale, maxScale);
            else
                warning("Invalid Type");
                knotVect = Nan;
            end
        end
    end
    methods (Static = true, Hidden = true)
        function knotVector = closed_loop_knot_vector(npts, minVal, maxVal)
%             assert(minVal < 1);
            knotVector = (0:6+npts-1)-5;
            if nargin > 1 && maxVal
                knotVector = knotVector./(knotVector(end))* maxVal;
            end
        end
        function knotVector = open_loop_knot_vector(numPts, minVal, maxVal, numRepeatsStart, numRepeatsEnd)
            npc = numPts+ 6;
            np2 = numPts+2;
            knotVector = zeros(1,npc); % total number of knots is number of control points plus the order in this case 6
            for i = 2:npc
                if i>numRepeatsStart && i<np2 + (6-numRepeatsEnd)
                    knotVector(i) = knotVector(i-1) + 1;
                else
                    knotVector(i) = knotVector(i-1);
                end
            end
            if minVal < 0
                knotVector = (knotVector - max(knotVector))./max(knotVector)*abs(minVal);
            else
                knotVector = knotVector./max(knotVector)*maxVal;
            end
        end
        function knotVector = custom_boundary_open_knot_vectors(numPts, minVal, maxVal, startRepeat, endRepeat)
            npc = numPts + 6;
            np2 = numPts + 2;
            knotVector = zeros(1,npc);
            for i = 2:npc
                if i > startRepeat && i < np2-endRepeat
                    knotVector(i) = knotVector(i-1) + 1;
                else
                    knotVector(i) = knotVector(i-1);
                end
            end
            if minVal < 0
                knotVector = (knotVector - max(knotVector))./max(knotVector)*abs(minVal);
            else
                knotVector = knotVector./max(knotVector)*maxVal;
            end
        end
    end
end