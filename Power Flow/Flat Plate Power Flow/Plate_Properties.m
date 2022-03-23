classdef Plate_Properties
    properties
        width
        thickness
        height
        E
        poissons
        forceLocations
    end
    methods
        function obj = Plate_Properties(width, height, thickness, E, poissons, fLs)
            obj.width = width;
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