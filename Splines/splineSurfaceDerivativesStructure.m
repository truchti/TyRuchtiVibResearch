classdef splineSurfaceDerivativesStructure
    properties
        xi
        eta
        controlPoints
    end
    properties(Access = private)
        bxi
        dxi
        ddxi
        dddxi
        beta
        deta
        ddeta
        dddeta
    end
    methods
        function obj = splineSurfaceDerivativesStructure(CPs, xi, bxi, dxi, ddxi, dddxi, eta, beta, deta, ddeta, dddeta)
            obj.controlPoints = CPs;
            obj.xi = xi;
            obj.bxi = bxi;
            obj.dxi = dxi;
            obj.ddxi = ddxi;
            obj.dddxi = dddxi;
            obj.eta = eta;
            obj.beta = beta;
            obj.deta = deta;
            obj.ddeta = ddeta;
            obj.dddeta = dddeta;
        end
        function value = evaluate_derivative_given_xi_and_eta_orders(obj, dxi, deta, splineNumber)
            if nargin < 4
                splineNumber = 1;
            end
            derivStr = obj.convert_derivative_orders_to_string(dxi, deta);
            value = obj.evaluate_surface_derivative(derivStr, splineNumber);            
        end
        function values = evaluate_surface_derivative(obj, derivative, splineNumber)
            if nargin < 3
                splineNumber = 1;
            end
            assert(splineNumber <= size(obj.controlPoints,3))
            switch derivative
                case 'dbdxi'
                    values = obj.dxi*obj.controlPoints(:,:,splineNumber)*obj.beta';
                case 'dbdeta'
                    values = obj.bxi*obj.controlPoints(:,:,splineNumber)*obj.deta';
                case 'ddbddxi'
                    values = obj.ddxi*obj.controlPoints(:,:,splineNumber)*obj.beta';
                case 'ddbdxideta'
                    values = obj.dxi*obj.controlPoints(:,:,splineNumber)*obj.deta';
                case 'ddbddeta'
                    values = obj.bxi*obj.controlPoints(:,:,splineNumber)*obj.ddeta';
                case 'dddbdddxi'
                    values = obj.dddxi*obj.controlPoints(:,:,splineNumber)*obj.beta';
                case 'dddbddxideta'
                    values = obj.ddxi*obj.controlPoints(:,:,splineNumber)*obj.deta';
                case 'dddbdxiddeta'
                    values = obj.dxi*obj.controlPoints(:,:,splineNumber)*obj.ddeta';
                case 'dddbdddeta'
                    values = obj.bxi*obj.controlPoints(:,:,splineNumber)*obj.dddeta';
                otherwise
                    warning("Improper derivative name");
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
    methods (Static = true)
        function derivString = convert_derivative_orders_to_string(dxi,deta)
            p1 = strcat(repmat('d', 1, dxi+deta),'b');
            if dxi >0
                p2 = strcat(repmat('d', 1, dxi), 'xi');
            else
                p2 = '';
            end
            if deta > 0
                p3 = strcat(repmat('d', 1, deta), 'eta');
            else
                p3 = '';
            end
            derivString = strcat(p1, p2, p3);
        end
    end
end
