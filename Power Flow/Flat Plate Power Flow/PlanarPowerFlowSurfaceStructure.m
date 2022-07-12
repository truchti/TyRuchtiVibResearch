classdef PlanarPowerFlowSurfaceStructure < handle
    properties
        wR
        wI
        wdotR
        wdotI
    end
    methods % creation methods
        function obj = PlanarPowerFlowSurfaceStructure()
        end
        function set_surface(obj, type, surf)
            switch type
                case {1, 'wR'}
                    obj.set_real_displacement_surface(surf)
                case {2, 'wI'}
                    obj.set_imag_displacement_surface(surf)
                case {3, 'wdotR'}
                    obj.set_real_velocity_surface(surf)
                case {4, 'wdotI'}
                    obj.set_imag_velocity_surface(surf)
            end
        end
        function set_real_displacement_surface(obj, surf)
            obj.wR = surf;
        end
        function set_imag_displacement_surface(obj, surf)
            obj.wI = surf;
        end
        function set_real_velocity_surface(obj, surf)
            obj.wdotR = surf;
        end
        function set_imag_velocity_surface(obj, surf)
            obj.wdotI = surf;
        end
    end
    methods % evaluation methods
        function vals = evaluate_surface_at_points(obj, type, xpts, ypts)
            vals = obj.(type).eval(xpts, ypts);
        end
        function vals = evaluate_complex_displacement_derivative_at_points(obj, xpts, ypts, derivs)
            if nargin < 4 || nnz(derivs)< 1
                valsR = obj.evaluate_surface_at_points('wR', xpts,ypts);
                valsI = obj.evaluate_surface_at_points('wI', xpts,ypts);
                vals = complex(valsR, valsI);
            else
                valsR = obj.evaluate_surface_derivative_at_points('wR', xpts, ypts, derivs);
                valsI = obj.evaluate_surface_derivative_at_points('wI', xpts, ypts, derivs);
                vals = complex(valsR, valsI);
            end
        end
        function vals = evaluate_complex_velocity_derivative_at_points(obj, xpts, ypts, derivs)
            if nargin < 4 || nnz(derivs)< 1
                valsR = obj.evaluate_surface_at_points('wdotR', xpts,ypts);
                valsI = obj.evaluate_surface_at_points('wdotI', xpts,ypts);
                vals = complex(valsR, valsI);
            else
                valsR = obj.evaluate_surface_derivative_at_points('wdotR', xpts, ypts, derivs);
                valsI = obj.evaluate_surface_derivative_at_points('wdotI', xpts, ypts, derivs);
                vals = complex(valsR, valsI);
            end
        end
        function vals = evaluate_surface_derivative_at_points(obj, type, xpts, ypts, derivs)
            vals = obj.(type).deriv_eval(xpts, ypts, derivs);
        end
        function vals = evaluate_displacement_at_points(obj, xpts, ypts, isReal)
            if nargin< 4
                isReal = true;
            end
            if isReal
                type = 'wR';
            else 
                type = 'wI';
            end
            vals = obj.evaluate_surface_at_points(type, xpts,ypts);
            
        end
        function vals = evaluate_velocity_at_points(obj, xpts, ypts, isReal)
            if nargin< 4
                isReal = true;
            end
            if isReal
                type = 'wdotR';
            else 
                type = 'wdotI';
            end
            vals = obj.evaluate_surface_at_points(type, xpts,ypts);
        end
    end
end