classdef CylinderPowerFlowSurfStruct < handle
    properties
        wR
        wI
        wdotR
        wdotI
        uR
        uI
        udotR
        udotI
        vR
        vI
        vdotR
        vdotI
    end
    methods % creation methods
        function obj = CylinderPowerFlowSurfStruct()
        end
        function set_surface(obj, type, surf)
            if isprop(obj, type)
                obj.(type) = surf;
            else
                warning("No appropriate property")
            end
        end
    end
    methods % evaluation methods
        function vals = evaluate_surface_at_points(obj, type, thPts, lPts)
            if isprop(obj, type)
                vals = obj.(type).eval(thPts, lPts);
            else
                warning("Invalid Input")
            end
        end
        function vals = evaluate_complex_derivative_at_points(obj, thPts, lPts, dir, derivs)
            if nargin < 5 
                derivs = 0;
            end
            valsR = obj.evaluate_surface_derivative_at_points([dir 'R'], thPts, lPts, dir, derivs);
            valsI = obj.evaluate_surface_derivative_at_points([dir 'I'], thPts, lPts, dir, derivs);
            vals = complex(valsR, valsI);
        end
        function vals = evaluate_surface_derivative_at_points(obj, type, thPts, lPts, derivs)
            if isprop(obj, type)
                vals = obj.(type).deriv_eval(thPts, lPts, derivs);
            else
                warning("Invalid Input")
            end
        end
    end
end