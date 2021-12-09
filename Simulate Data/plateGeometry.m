classdef plateGeometry
    properties
        width
        height
        thickness
    end
    properties (Dependent)
        mu
    end
    methods
        function obj = plateGeometry(width, height, thick)
            if nargin < 1
                obj = obj.create_default();
            else
                obj.width = width;
                obj.height = height;
                obj.thickness = thick;
            end
        end
        function value = get.mu(obj)
            value = obj.thickness^2;
        end
        function obj = create_default(obj)
            obj.width = .400;
            obj.thickness = .00163;
            obj.height = .502;
        end
    end
end
