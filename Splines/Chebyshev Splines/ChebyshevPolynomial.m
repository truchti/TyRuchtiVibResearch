classdef ChebyshevPolynomial 
    properties
        coefficients
        inputRange
    end
    properties (Dependent = true)
        order
    end
    properties (Hidden = true)
        basisEvaluator
        mapInput
    end
    methods
        function obj = ChebyshevPolynomial(coeffs, range)
            obj.coefficients = coeffs;
            obj.inputRange = range;
            obj = obj.create_mapping_function();
            obj.basisEvaluator = ChebyshevBasisEvaluator(obj.order);
        end
        function [fx, Tn] = eval(obj, pts)
            xis = obj.mapInput(pts);
            Tn = obj.basisEvaluator.evaluate_at_parameter(xis);
            fx = Tn*obj.coefficients;
        end
        function dfdx = deriv_eval(obj, pts, deriv)
            xis = obj.mapInput(pts);
            dTn = obj.basisEvaluator.evaluate_basis_and_derivatives(xis, deriv);
            dfdxi = dTn*obj.coefficients;
            dxidx = mean(diff(xis)./diff(pts)); % assumes linear mapping;
            dfdx = dfdxi*dxidx^(deriv);
        end
        function plot(obj, pts)
            fx = obj.eval(pts);
            plot(pts, fx)
        end
        function plot_deriv(obj,pts, deriv)
            dfx = obj.deriv_eval(pts, deriv);
            plot(pts,dfx)
        end
        function value = get.order(obj)
            value =length(obj.coefficients)-1;
        end
        function obj = create_mapping_function(obj)
            a = min(obj.inputRange);
            b = max(obj.inputRange);
            obj.mapInput = @(x) ((x-a)/(b-a)-.5)*2;
        end
    end
end
