classdef quinticBSplineFitter < matlab.mixin.Copyable
    % This class takes in xi and eta parameter vectors and a 2D array of data
    % values at each of these parameter pairs. It contains methods to
    % fit a spline surface to the data and to export a custom spline object
    properties
        xiParams
        data
        splineType
        numControlPts
        xiKnot
        boundaryNodeRepeats
    end
    properties (Hidden = true)
        splineEvaluator
        solved = false
    end
    methods
        function obj = quinticBSplineFitter(parameters, data, type, numberControlPoints, numberOfRepeatedNodes)
            if isrow(parameters)
                parameters = parameters';
            end
            if size(data,1) == 1
                data = data';
            end
            obj.xiParams = parameters;
            obj.data = data;    obj.splineType = type;
            obj.numControlPts = numberControlPoints;
            if nargin < 5
                obj.boundaryNodeRepeats = 6;
            else
                obj.boundaryNodeRepeats = repmat(numberOfRepeatedNodes, [1,2]);
            end
            obj.calculate_knot_vectors();
        end
        function xiB = fit_spline(obj)
            %% Set up evaluator with all zero control Points
            CPs = zeros(obj.numControlPts, 1); % load dummy control points in initially since they will be solved for
            obj.splineEvaluator = quinticSplineEvaluator(obj.xiKnot, CPs, obj.splineType);
            %evaluate basis functions
            xiB = obj.splineEvaluator.evaluate_basis_functions(obj.xiParams);
            for i = size(obj.data,2):-1:1
                tempData = obj.data(:,i);
                if nnz(obj.data(:,i)) < 1 %% if there is no data in out of plane dimension because it is simulated plate don't calculate
                    obj.splineEvaluator.controlPoints = zeros(obj.controlNetDimensions);
                else
                    obj.splineEvaluator.controlPoints =(xiB\tempData);
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
            obj.xiKnot = obj.calculate_knot_vector_using_type(obj.numControlPts, obj.splineType, min(obj.xiParams), max(obj.xiParams));
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
        function knotVector = closed_loop_knot_vector(npts, ~, maxVal)
            a = 5;
            knotVector = (0:6+npts-1)-a;% five knots before zero 
            if nargin > 1 && maxVal
                knotVector = knotVector./(knotVector(end-a))* maxVal;
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
%         function knotVector = custom_boundary_open_knot_vectors(numPts, minVal, maxVal, startRepeat, endRepeat)
%             npc = numPts + 6;
%             np2 = numPts + 2;
%             knotVector = zeros(1,npc);
%             for i = 2:npc
%                 if i > startRepeat && i < np2-endRepeat
%                     knotVector(i) = knotVector(i-1) + 1;
%                 else
%                     knotVector(i) = knotVector(i-1);
%                 end
%             end
%             if minVal < 0
%                 knotVector = (knotVector - max(knotVector))./max(knotVector)*abs(minVal);
%             else
%                 knotVector = knotVector./max(knotVector)*maxVal;
%             end
%         end
    end
end