classdef cylinderGeometry
    properties
        radius
        height
        thickness
    end
    properties (Dependent = true, Hidden = true)
        mu
    end
    methods
        function obj = cylinderGeometry(rad, height, thick)
            if nargin < 1
                obj = obj.create_default();
            else
                obj.radius = rad;
                obj.height = height;
                obj.thickness = thick;
            end
        end
        function value = get.mu(obj)
            value = obj.thickness^2/(12*obj.radius^2);
        end
        function obj = create_default(obj)
            obj.radius = .02;
            obj.thickness = 2e-3;
            obj.height = .160;
        end
    end
end
