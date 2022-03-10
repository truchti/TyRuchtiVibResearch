classdef sinusoidalDistForceCart
    properties
        xbounds
        ybounds
        magnitude
        frequency
        phase
        enabled = true;
    end
    properties (Dependent = true, Hidden = true)
%         forceCenter
        omega
    end
    methods
        function obj = sinusoidalDistForceCart(xbounds, ybounds, magnitude, frequency, phase)     
                obj.xbounds = xbounds;
                obj.ybounds = ybounds;
                obj.magnitude = magnitude;
                obj.frequency = frequency;
                obj.phase = phase;
        end
        function value = get.omega(obj)
            value = obj.frequency*2*pi;
        end

    end
end
