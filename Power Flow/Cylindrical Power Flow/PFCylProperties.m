classdef PFCylProperties
    properties
        radius
        thickness
        height
        E
        poissons
        forceLocations
    end
    methods
        function obj = PFCylProperties(radius, thickness, height, E, poissons, fLs)
            obj.radius = radius;
            obj.thickness = thickness;
            obj.height = height;
            obj.E = E;
            obj.poissons = poissons;
            if nargin > 5
                obj.forceLocations = fLs;
            end
        end
    end
end