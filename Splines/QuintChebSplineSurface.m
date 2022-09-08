classdef QuintChebSplineSurface <handle
    % this creates a surface where the closed loop is fit using quintic spline
    % patches and the open section is fit using a chebyshev polynomial
    % the first parameter (xi) will always refer to the closed spline and the
    % second parameter (eta) will refer to the open direction fit with the
    % Chebyshev
    properties
        controlNet
        closedKnots
        quinticPatches
        heightRange % used to map the input to range of cheb function (-1:1)
    end
    properties (Dependent = true)
        ChebOrder
        closedCPs 
    end
    properties (Hidden = true)
        quintBasisEvaluator quinticSplineBasisEvaluator
        chebBasisEvaluator  ChebyshevBasisEvaluator
        mapOpenInputToCheb
    end
    methods
        function obj = QuintChebSplineSurface(controlNet, closedKnotVector, openRange)
            if nargin > 0
                obj.controlNet = controlNet;
                obj.closedKnots = closedKnotVector;
                obj.heightRange = openRange;
                obj.create_mapping_function;
                obj.quinticPatches = length(closedKnotVector)-5; 
                obj.quintBasisEvaluator = quinticSplineBasisEvaluator(obj.closedKnots,obj.closedCPs, 'closed');
                obj.chebBasisEvaluator = ChebyshevBasisEvaluator(size(obj.controlNet,2)-1);
            end
        end
        function [fx, Bxi, Teta] = evaluate(obj, Ths, Zs)
            Bxi = obj.quintBasisEvaluator.evaluate_at_parameter(Ths);
            eta = obj.mapOpenInputToCheb(Zs);
            Teta = obj.chebBasisEvaluator.evaluate_at_parameter(eta);
            fx = Bxi*obj.controlNet*Teta';
        end
        function  newSurf = plus(surf1,surf2)
            newSurf = QuintChebSplineSurface.empty();
            if ~ (surf1.ChebOrder == surf2.ChebOrder) || ~(surf1.quinticPatches == surf2.quinticPatches)
                warning("Cannot add surfaces with different orders or patch numbers");
                return;
            end
            if ~(all(approxEqual(surf1.closedKnots, surf2.closedKnots)))
                warning("Cannot add surfaces with different not vectors")
                return;
            end
            if ~(all(approxEqual(surf1.heightRange, surf2.heightRange)))
                warning("Cannot add surfaces with different height ranges")
                return;
            end
            newSurf = QuintChebSplineSurface(surf1.controlNet + surf2.controlNet, surf1.closedKnots, surf1.heightRange);
        end
        function newSurf = uminus(obj)
            newSurf = QuintChebSplineSurface(-obj.controlNet, obj.closedKnots, obj.heightRange);
        end
        function newSurf = minus(surf1, surf2)
            newSurf = surf1 + (-surf2);
        end
        function newSurf = times(surf, scalar)
            if ~isnumeric(scalar) && isnumeric(surf)
                tempscalar = surf;
                surf = scalar;
                scalar = tempscalar;
            end
            newSurf = QuintChebSplineSurface(scalar*surf.controlNet, surf.closedKnots, surf.heightRange);
        end
        function newSurf = mtimes(surf, scalar)
            newSurf = times(surf, scalar);
        end
        function deriv_eval(obj, Ths, Zs)
        end
        %         function plot(obj, Ths, Zs)
        %         end
        function value = get.ChebOrder(obj)
            value = size(obj.controlNet,2)-1;
        end
        function value = get.closedCPs(obj)
            value = size(obj.controlNet,1);
        end
        function create_mapping_function(obj)
            az = min(obj.heightRange);
            bz = max(obj.heightRange);
            obj.mapOpenInputToCheb =  @(z) ((z-az)/(bz-az)-.5)*2;
        end
        function  hc = hard_copy(obj)
            hc = QuintChebSplineSurface(obj.controlNet, obj.closedKnots, obj.heightRange);
        end
    end
end

