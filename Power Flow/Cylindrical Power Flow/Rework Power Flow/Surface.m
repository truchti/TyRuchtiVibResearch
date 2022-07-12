classdef Surface
    %Surface is a prototype class that defines the basic capability of any
    %surface obj. That is a surface must be define in such a way that
    %requesting the value of the surface at parameters is possible. There
    %is a way to visualize the surface over a parameter set and a new
    %surface can be created which is the m by nth spatial derivative of the
    %surface. 
    
    properties
        xiRange {mustBeNumeric}
        etaRange {mustBeNumeric}
    end
    properties(Hidden = true)
        xiKnot
        etaKnot
        controlNet
    end
    
    methods
        function obj = Surface()
        end
        
        function outputArg = plot_surface(obj,inputArg)
            outputArg = obj.Property1 + inputArg;
        end
        function generate_derivative_surface(obj, xiOrder, etaOrder)
        end
        function evaluate_surface(obj)
        end
    end
    methods (Abstract)
        foo(obj)
        
    end
end

