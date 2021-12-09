classdef meshDetails
    properties
        longitude_modes
        theta_modes
        longitude_divisions
        theta_divisions
    end
    methods
        function obj = meshDetails(modesL, modesT, divsL, divsT)
            if nargin < 1
                obj = obj.create_default();
            else
                obj.longitude_modes = modesL;
                obj.theta_modes = modesT;
                obj.longitude_divisions = divsL;
                obj.theta_divisions= divsT;
            end
        end
        function obj = create_default(obj)
            obj.longitude_modes = 25;
            obj.theta_modes = 25;
            obj.longitude_divisions = 30;
            obj.theta_divisions = 50;
        end
    end
end