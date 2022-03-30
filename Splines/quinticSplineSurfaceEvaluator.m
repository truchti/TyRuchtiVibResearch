classdef quinticSplineSurfaceEvaluator < matlab.mixin.Copyable
    % this class takes in the knot Vector, order, and control points
    % it has methods for evaluating the surface and its deriviatives for a
    % set of parameters
    properties
        knotVectors
        controlPoints
        types        
    end
    properties (Hidden = true)
        xiEvaluator
        etaEvaluator
        derivativeEvaluator
    end
    properties(Dependent)
        maxXi
        maxEta
    end
    methods
        function obj = quinticSplineSurfaceEvaluator(knotVectors, CPs, types)
            obj.knotVectors = knotVectors;
            if numel(CPs) ==2
                obj.controlPoints = zeros(CPs);
            else
                obj.controlPoints = CPs;
            end
            obj.types = types;
            obj.create_evaluators();
        end
        
        function values = evaluate_at_parameters(obj, xi, eta, measure, isImaginary, dimension, derivative)
            if nargin < 7
                derivative = 0;
            end
            if nargin < 6
                dimension = 1:3;
            end
            if nargin < 5
                isImaginary = false;
            end
            if nargin < 4
                measure = 'Velocity';
            end
            if nargin < 3
                error("Must provide parameters to evaluate at")
            end
            if size(obj.controlPoints,3)>1
                splineNum = obj.get_number_based_on_types(measure, isImaginary, dimension);
            else
                splineNum = 1;
            end
            values = obj.evaluate_spline_number_at_parameters(xi,eta, splineNum, derivative);
        end
        
        function complexValue = evaluate_complex_value_at_parameters(obj,xi,eta,measure, dimension,derivative)
            if nargin < 6
                derivative = '';
            end
            realPart = obj.evaluate_at_parameters(xi,eta, measure, false,dimension, derivative);
            imagPart = obj.evaluate_at_parameters(xi,eta,measure,true,dimension,derivative);
            complexValue = complex(realPart, imagPart);
        end
        
        function values = evaluate_spline_number_at_parameters(obj, xi, eta, splineNum, derivative)
            if nargin < 5
                derivative = '';
            end
            selectControlPoints = obj.controlPoints(:,:,splineNum);
            if isempty(derivative) || (length(derivative) == 1  && derivative == 0 ) || nnz(derivative) == 0
                [bXi, bEta] = obj.evaluate_basis_functions(xi, eta);
                values = bXi*selectControlPoints*bEta';
            else
                obj.validate_or_create_deriv_evaluator(xi, eta);
                if isnumeric(derivative)
                    dxi = derivative(1);
                    if length(derivative)< 2
                        warning("Since only one order was specified for the derivative assuming deta is 0")
                        deta = 0;
                    else
                        deta = derivative(2);
                    end
                    values = obj.derivativeEvaluator.evaluate_derivative_given_xi_and_eta_orders(dxi, deta, splineNum);
                else
                    values = obj.derivativeEvaluator.evaluate_surface_derivative(derivative, splineNum);
                end
            end            
        end
        function surfaceValues = evaluate_surface_at_parameters(obj, xi, eta)
            [xiValues, etaValues] = obj.evaluate_basis_functions(xi, eta);
            surfaceValues = zeros(length(xi), length(eta), size(obj.controlPoints,3));
            for i = 1:size(obj.controlPoints,3)
                surfaceValues(:,:,i) = xiValues*obj.controlPoints(:,:,i)*etaValues';
            end
        end
        function plot_spline(obj, xi, eta)
            v = obj.evaluate_surface_at_parameters(xi,eta);
            [X,Y] = ndgrid(xi,eta);
            surf(X,Y,v)
        end
        function plot_derivative(obj, xi, eta, dxi, deta)
            v = obj.evaluate_surface_derivative_at_parameter(xi,eta,dxi,deta);
            [X,Y] = ndgrid(xi,eta);
            surf(X,Y,v)
        end
        function values = evaluate_surface_derivative_at_parameter(obj, xi, eta, dxi, deta)
            obj.validate_or_create_deriv_evaluator(xi, eta);
            values = obj.derivativeEvaluator.evaluate_derivative_given_xi_and_eta_orders(dxi,deta);            
        end
    end
    methods 
        function value = get.maxXi(obj)
            value = max(obj.knotVectors{1});
        end
        function value = get.maxEta(obj)
            value = max(obj.knotVectors{2});
        end
    end
    methods (Hidden = true)
        function  [xiBasis, etaBasis] = evaluate_basis_functions(obj,xi,eta)
            xiBasis =  obj.xiEvaluator.evaluate_at_parameter(xi);
            etaBasis = obj.etaEvaluator.evaluate_at_parameter(eta);
        end
    end
    methods (Access = private)
        function numb = get_number_based_on_types(obj, measure, isImaginary, dimension)
            if strcmp(measure, 'disp')
                M = 0;
            else
                M = 6;
            end
            if ~isImaginary
                RI = 0;
            else
                RI = 3;
            end
            switch dimension
                case {'x', 'theta', 1}
                    D = 1;
                case {'y', 'rad', 2}
                    D = 2;
                case {'z', 'height', 3}
                    D = 3;
            end
            numb = M+RI+D;
            if numb > size(obj.controlPoints,3)
                error("Wrong Parameters");
            end
        end
        function create_evaluators(obj)
            obj.xiEvaluator = quinticSplineBasisEvaluator(obj.knotVectors{1}, size(obj.controlPoints, 1), obj.types{1});
            obj.etaEvaluator = quinticSplineBasisEvaluator(obj.knotVectors{2}, size(obj.controlPoints, 2), obj.types{2});
        end
        function validate_or_create_deriv_evaluator(obj, xi, eta)
            if isempty(obj.derivativeEvaluator) || isempty(obj.derivativeEvaluator.controlPoints)|| ~obj.derivativeEvaluator.check_xi_and_eta(xi,eta)
                    obj.create_derivative_evaluator(xi,eta)
            end
        end
        function create_derivative_evaluator(obj, xi, eta)
            [bxi, dbdxi, dbddxi, dbdddxi] = obj.xiEvaluator.evaluate_basis_and_derivatives(xi);
            [beta, dbdeta, dbddeta, dbdddeta] = obj.etaEvaluator.evaluate_basis_and_derivatives(eta);
            obj.derivativeEvaluator = splineSurfaceDerivativesStructure(obj.controlPoints, xi, bxi, dbdxi, dbddxi, dbdddxi, eta, beta, dbdeta, dbddeta, dbdddeta);
        end
    end
end
