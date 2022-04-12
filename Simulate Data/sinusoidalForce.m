classdef sinusoidalForce < handle
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
        tStar
    end
    methods
        function obj = sinusoidalForce(location, magnitude, frequency, phase)
            if nargin < 2
                switch location
                    case 1
                        obj = obj.create_default_1();
                    case 2
                        obj = obj.create_default_2();
                    case 3
                        obj = obj.create_default_rect_1();
                    case 4
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
        function value = get.tStar(obj)
            value = obj.location(2);
        end
        function obj = create_default_1(obj)
            obj.location  = [.025, 0*pi/180];
            obj.magnitude = 1;
            obj.frequency = 3050;
            obj.phase = 0;
        end
        function obj = create_default_2(obj)
            obj.location  = [.135, 180*pi/180];
            obj.magnitude = 1;
            obj.frequency = 3050;
            obj.phase = 170;
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
        function set_angle_loc(obj, angle)
            obj.location(2) = angle;
        end
        function set_height_loc(obj, height)
            obj.location(1) = height;
        end
    end
end
