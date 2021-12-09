classdef powerFlowPlate < handle
    properties
        x
        y
        h
        numberOfModes
        meshWidth
        E = 2.040e11; %modulus of material
        poi = .3; %poissons ration
        eta = .1;%structural damping
        areaDen = 7869.0/(a*b); %mass per area
        Force
    end
    properties (Hidden = true)
        D
        Q
        natFreqs 
        xMs
        yMs
        denom 
        cosXmodes 
        cosYmodes 
    end
    properties (Dependent = true, Hidden = true)
        M 
        N
    end
    methods
        function obj = powerFlowPlate(length, width, thickness, modes, meshWidth, force)
            obj.x = width;
            obj.y = length;
            obj.h = thickness;
            obj.numberOfModes = modes;
            obj.meshWidth = meshWidth;
            obj.Force = force;
            obj.calculate_constants();
        end
    end
    methods (Hidden = true)
        function calculate_constants(obj)
            obj.calculate_D()
            obj.calculate_density()
            obj.calculate_Q()
            obj.calculate_natural_frequencies()
            obj.calculate_mode_cosines()
            obj.calculate_denominators()
        end
        function calculate_D(obj)
            obj.D = obj.E*obj.h^3/(12*(1-obj.poi^2));
        end
        function calculate_density(obj)
            obj.areaDen = 7869.0/(obj.x*obj.y);
        end
        function calculate_Q(obj)
            obj.Q = 4*obj.Force.magnitude/(obj.areaDen*obj.h*pi^2);
        end
        function calculate_natural_frequencies(obj)
            obj.natFreqs = pi^2*((obj.xMs'/obj.x).^2+(obj.yMs/obj.y).^2)*sqrt(obj.D/(obj.areaDen*obj.h));
        end
        function calculate_denominators(obj)
             obj.denom = ((obj.natFreqs.^2-obj.Force.frequency^2).^2+obj.eta^2*obj.natFreqs.^4);
        end
        function calculate_mode_cosines(obj)
            obj.cosXmodes = cos(obj.M*pi*x1/a)-cos(obj.M*pi*x2/a);
            obj.cosYmodes = cos(obj.N*pi*y1/b)-cos(obj.N*pi*y2/b);
        end
        function w = calculate_plate_displacement(obj)
        end
        function dwdt = calculate_plate_velocity(obj)
        end
    end
end
