classdef SLDV_data_cyl_coors< SLDV_data
    properties
        radius = 1;
    end
    methods
        function flat_plot_displacement(obj, dim ,imaginary)
            if nargin < 3
                imaginary = false;
            end
            if nargin < 2 
                dim = 2;
            end
            obj.flat_plot('disp', dim, imaginary);
        end
        function plot_cyl_coors(obj)
            obj.coordinates;
            obj.radius;
        end
        function flat_plot_velocity(obj, imaginary)
            if nargin < 2 
                imaginary = false;
            end
            obj.flat_plot('vel', 2, imaginary);
        end
        function [params, data, types] = get_cylindrical_spline_fitting_data(obj)
            params = obj.coordinates(:,1:2:3);
            params(:,2) = params(:,2)-min(params(:,2)); % ensure that the min height is zero
            data = obj.get_array_of_complex_displacement_and_velocity_values();
            types = {'closed', 'open'};
        end
        function obj = eliminate_points_off_of_cylinder(obj)
            CCP = CylCoorParameterizer(obj);
            CCP.calculate_cylindrical_data();
            cylData = CCP.export_cylindrical_data();
            radii = cylData.coordinates(:,2);
            m = mean(radii);    d = std(radii); rng  = [m-1.5*d m+1.5*d];
            valid = radii > rng(1) & radii < rng(2);
            obj.validPoints = logical(valid);
        end
    end
    methods(Hidden = true)
        function flat_plot(obj, type, dimension, imaginary)
            if strcmp(type, 'disp')
                [T,R,H] = obj.separate_displacements(imaginary);
            elseif strcmp(type, 'vel')
                [T,R,H] = obj.separate_velocities(imaginary);
            else
                error('Wrong Type')
            end
            switch dimension
                case 1
                    plotData = T;
                case 2
                    plotData = R;
                case 3
                    plotData = H;
            end
            [T,~,H] = obj.separate_coordinates;
            scatter3(wrapTo2Pi(T),H,plotData)
        end
    end
end