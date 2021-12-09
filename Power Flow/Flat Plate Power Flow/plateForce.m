classdef plateForce < handle
    
    properties
        location
        magnitude
        width
        frequency
       
    end
    properties (Dependent = true)
        x1
        x2
        y1
        y2
        omega
    end
    methods
        function obj = plateForce(centerLocation, magnitude, size, frequency)
            obj.location = centerLocation;
            obj.magnitude = magnitude;
            obj.width = size;
            obj.frequency = frequency;
        end
        function value = get.x1(obj)
            value = obj.location(1) - obj.width/2;
        end
        function value = get.x2(obj)
            value = obj.location(1) + obj.width/2;
        end
        function value = get.y1(obj)
            value = obj.location(2) - obj.width/2;
        end
        function value = get.y2(obj)
            value = obj.location(2) + obj.width/2;
        end
        function value = get.omega(obj)
            value = obj.frequency*2*pi;
        end
    end
end
    