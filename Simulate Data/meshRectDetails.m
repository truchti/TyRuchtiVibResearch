classdef meshRectDetails
    properties
        height_modes
        width_modes
        height_divisions
        width_divisions
    end
    methods
        function obj = meshRectDetails(modesL, modesT, divsL, divsT)
            if nargin < 1
                obj = obj.create_default();
            else
                obj.height_modes = modesL;
                obj.width_modes = modesT;
                obj.height_divisions = divsL;
                obj.width_divisions= divsT;
            end
        end
        function obj = create_default(obj)
            obj.height_modes = 20;
            obj.width_modes = 20;
            obj.height_divisions = 40;
            obj.width_divisions = 50;
        end
    end
end