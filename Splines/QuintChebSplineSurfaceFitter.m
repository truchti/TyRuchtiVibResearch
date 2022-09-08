classdef QuintChebSplineSurfaceFitter < matlab.mixin.Copyable
    properties
        Ths
        Zs
        Vals
        ChebOrder
        QuinticPatches
        
    end
    properties (Hidden = true)
       
        dummySurf
        exportSurfs
        solved = false
        chebRange
        quinticKnots
    end
    methods
        function obj = QuintChebSplineSurfaceFitter(Ths, Zs, values, numQuintPatches, ChebOrder)
            obj.Ths = Ths;
            obj.Zs = Zs;
            obj.Vals = values;
            obj.QuinticPatches = numQuintPatches;
            obj.ChebOrder = ChebOrder;
            obj.calculate_knot_vector();
            obj.calculate_cheb_range();
        end
        function [Bxi, Teta] = fit_spline_surfaces(obj)
            %% Set up evaluator with all zero control Points
            CN = zeros(obj.QuinticPatches+5, obj.ChebOrder+1);
            obj.dummySurf = QuintChebSplineSurface(CN, obj.quinticKnots, obj.chebRange);
            %evaluate basis functions
            [~, Bxi, Teta] = obj.dummySurf.evaluate(obj.Ths, obj.Zs);
            obj.exportSurfs = QuintChebSplineSurface.empty(0, size(obj.Vals,2));
            for i = size(obj.Vals,2):-1:1
                tempData = obj.Vals(:,i);
                if nnz(obj.Vals(:,i)) < 1 %% if there is no data in out of plane dimension because it is simulated plate don't calculate
                    obj.splineEvaluator.controlPoints(:,:,i) = zeros(obj.controlNetDimensions);
                else
                    if size(obj.Vals,1) > 9999
                        tempData = sparse(diag(tempData)); %% if would be too large in memory so make a sparse matrix since it is diagonal
                    else
                        tempData = diag(tempData);
                    end
                    tempControlNet = numel(obj.Vals(:,i))*(Bxi\tempData/(Teta'));
                    obj.exportSurfs(i) = QuintChebSplineSurface(tempControlNet, obj.quinticKnots, obj.chebRange);
                end
            end
            obj.solved = true;
        end  
        function Surfs = output_surfaces(obj)
            if ~obj.solved
                obj.fit_spline_surfaces;
            end
            for i = length(obj.exportSurfs):-1:1
                Surfs(i) = obj.exportSurfs(i).hard_copy();
            end
        end
%         function check_by_viewing_spline_and_data(obj)
%             xi = linspace(0, max(obj.xiKnot), 100);
%             eta = linspace(0, max(obj.etaKnot), 50);
%             [ETA, XI] = meshgrid(eta, xi);
%             for n = 1:size(obj.Vals,2)
%                 figure(n)
%                 values = obj.splineEvaluator.evaluate_spline_number_at_parameters( xi, eta, n);
%                 scatter3(obj.Ths, obj.Zs, obj.Vals(:,n))
%                 hold on;
%                 surf(XI,ETA, values)
%                 hold off;
%             end
%         end
    end
    methods (Hidden = true)
        function calculate_knot_vector(obj)
            npts = obj.QuinticPatches + 5; % +order(6) -1
            obj.quinticKnots = obj.closed_loop_knot_vector(npts, 2*pi);
        end
        function calculate_cheb_range(obj)
            obj.chebRange = [min(obj.Zs), max(obj.Zs)];
        end
    end
    methods (Static = true, Hidden = true)
        function knotVector = closed_loop_knot_vector(npts, maxVal)

            knotVector = (0:6+npts-1)-5; % start before zero so that 0:2pi maps full circle
            if nargin > 1 && maxVal
                knotVector = knotVector./(knotVector(end))* maxVal;
            end
        end
    end
end