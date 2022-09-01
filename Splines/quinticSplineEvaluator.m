classdef quinticSplineEvaluator < matlab.mixin.Copyable
    % this class takes in the knot Vector, order, and control points
    % it has methods for evaluating the surface and its deriviatives for a
    % set of parameters
    properties
        knotVector
        controlPoints
        splineType      
    end
    properties (Hidden = true)
        basisEvaluator
        bxi; dbdxi; dbddxi; dbdddxi
    end
    properties(Dependent)
        maxXi
    end
    methods
        function obj = quinticSplineEvaluator(knotVector, CPs, type)
            obj.knotVector = knotVector;
            obj.controlPoints = CPs;
            obj.splineType = type;
            obj.create_basis_evaluator();
        end
        
        function values = evaluate_at_parameters(obj, xi,  derivative)
            if nargin < 3
                derivative = 0;
            end
            if nargin < 2
                error("Must provide parameters to evaluate at")
            end
            values = obj.evaluate_spline_number_at_parameters(xi, [], derivative);
        end
        function values = evaluate_spline_number_at_parameters(obj, xi, ~, derivative)
            if nargin < 4
                derivative = '';
            end
            if isempty(derivative) || (length(derivative) == 1  && derivative == 0 ) || nnz(derivative) == 0
                bXi = obj.evaluate_basis_functions(xi);
                values = bXi*obj.controlPoints;
            else
                dXi  = obj.evaluate_basis_function_derivatives(xi, derivative);
                scale = 1/mean(diff(obj.knotVector));
                values = (scale^derivative) * dXi * obj.controlPoints;
            end            
        end
        
        function plot_spline_number(obj, xi, eta, number, derivs)
            v = obj.evaluate_spline_number_at_parameters(xi, eta, number, derivs);
%             Zs = max(max(v));
%             Ws = min(min(v));
            [X,Y] = ndgrid(xi,eta);
            surf(X,Y,v, 'EdgeAlpha', 0);
            view(0,90); colormap jet; colorbar;
        end
        function plot_spline(obj, xi)
            v = obj.evaluate_at_parameters(xi);
            plot(xi,v);
        end
        
    end
    methods %% getter
        function value = get.maxXi(obj)
            value = max(obj.knotVector);
        end
    end
    methods (Hidden = true)
        function  xiBasis = evaluate_basis_functions(obj,xi)
            xiBasis =  obj.basisEvaluator.evaluate_at_parameter(xi);
        end
        function xiBasisDeriv = evaluate_basis_function_derivatives(obj, xi, deriv)
            [~, dxi, ddxi, dddxi] = obj.basisEvaluator.evaluate_basis_and_derivatives(xi);
            switch deriv
                case 1
                    xiBasisDeriv = dxi;
                case 2 
                    xiBasisDeriv = ddxi;
                case 3 
                    xiBasisDeriv = dddxi;
            end
        end        
    end
    methods (Access = private)
        function create_basis_evaluator(obj)
            obj.basisEvaluator = quinticSplineBasisEvaluator(obj.knotVector, obj.controlPoints, obj.splineType);
        end
        function validate_or_create_deriv_evaluations(obj, xi, eta)
            if isempty(obj.derivativeEvaluator) || isempty(obj.derivativeEvaluator.controlPoints)|| ~obj.derivativeEvaluator.check_xi_and_eta(xi,eta)
                obj.create_derivative_evaluations(xi,eta)
            end
        end
        function create_derivative_evaluations(obj, xi)
            [obj.bxi, obj.dbdxi, obj.dbddxi, obj.dbdddxi] = obj.basisEvaluator.evaluate_basis_and_derivatives(xi);
%             obj.derivativeEvaluator = splineSurfaceDerivativesStructure(obj.controlPoints, xi, bxi, dbdxi, dbddxi, dbdddxi, eta, beta, dbdeta, dbddeta, dbdddeta);
        end
    end
end
