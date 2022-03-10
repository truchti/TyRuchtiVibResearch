classdef sinusoidalForceCart
    properties
        location
        magnitude
        frequency
        phase
        enabled = true;
    end
    properties (Dependent = true, Hidden = true)
        omega
        xStar
        yStar
    end
    methods
        function obj = sinusoidalForceCart(location, magnitude, frequency, phase)
            if nargin < 2
                switch location
                    case 1
                        obj = obj.create_default_rect_1();
                    case 2
                        obj = obj.create_default_rect_2();
                end
            else
                if nargin < 4
                    phase = 0;
                end                
                obj.location = location;
                obj.magnitude = magnitude;
                obj.frequency = frequency;
                obj.phase = phase;
            end
        end
        function value = get.omega(obj)
            value = obj.frequency*2*pi;
        end
        function value = get.xStar(obj)
            value = obj.location(1);
        end
        function value = get.yStar(obj)
            value = obj.location(2);
        end
        function obj = create_default_rect_1(obj)
            obj.location  = [.05, .1];
            obj.magnitude = 400;
            obj.frequency = 60;
            obj.phase = 0;
        end
        function obj = create_default_rect_2(obj)
            obj.location  = [.35, .4];
            obj.magnitude = 400;
            obj.frequency = 60;
            obj.phase = 178;
        end
    end
end
