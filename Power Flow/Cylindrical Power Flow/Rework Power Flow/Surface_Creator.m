classdef Surface_Creator
    %This class is a general surface class for any spline surface
    %generation and dictate common output methods.  Derived classes 
    
    properties
        CoorPts
        Values
        
    end
    
    methods
        function obj = Surface_Creator(cps, values)
            %UNTITLED11 Construct an instance of this class
            %   Detailed explanation goes here
            obj.CoorPts = cps;
            obj.Values = values;
        end
        
        function outputArg = fit_surface(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
%         function
%         end
    end
end

