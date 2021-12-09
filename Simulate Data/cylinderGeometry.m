classdef cylinderGeometry
    properties
        radius
        height
        thickness
    end
    properties (Dependent)
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
            obj.radius = .0635;
            obj.thickness = .00163;
            obj.height = .502;
        end
    end
end
