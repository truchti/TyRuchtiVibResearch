classdef splineDerivativesStructure
    properties
        xi
        controlPoints
    end
    properties(Access = private)
        bxi
        dxi
        ddxi
        dddxi
    end
    methods
        function obj = splineSurfaceDerivativesStructure(CPs, xi, bxi, dxi, ddxi, dddxi)
            obj.controlPoints = CPs;
            obj.xi = xi;
            obj.bxi = bxi;
            obj.dxi = dxi;
            obj.ddxi = ddxi;
            obj.dddxi = dddxi;
        end
        function value = evaluate_derivative_order(obj, order, splineNumber)
            if nargin < 4
                splineNumber = 1;
            end
            value = obj.evaluate_derivative(order, splineNumber);            
        end
        function values = evaluate_derivative(obj, derivative)
            if nargin < 3
                splineNumber = 1;
            end
            switch derivative
                case 0
                    values = obj.dxi*obj.controlPoints(:,:,splineNumber)*obj.beta';
                case 1
                    values = obj.bxi*obj.controlPoints(:,:,splineNumber)*obj.deta';
                case 2
                    values = obj.ddxi*obj.controlPoints(:,:,splineNumber)*obj.beta';
                case 3
                    values = obj.dxi*obj.controlPoints(:,:,splineNumber)*obj.deta';
                otherwise
                    warning("Derivative Order not computable");
            end
        end
        function tf = check_xi_and_eta(obj,xi,eta)
            if isequal(obj.xi, xi) && isequal(obj.eta,eta)
                tf = true;
            else
                tf = false;
            end
        end
    end
end
