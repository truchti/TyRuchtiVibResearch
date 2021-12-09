classdef SLDV3D_FFT_Data_Cylinder < handle
    properties
        original_data_source
        coordinates_cyl
        displacements_cyl
        velocities_cyl
        frequencies
    end
    properties (Hidden = true)
        cyl_axis
        cyl_radius
        cyl_zero
    end
    properties (Hidden = true, Access = private)
        coors_rect
        disp_rect
        vel_rect
    end
    methods
        function obj = SLDV3D_FFT_Data_Cylinder(fft_data)
            obj.coors_rect = fft_data.coordinates;
            obj.disp_rect = fft_data.displacements;
            obj.vel_rect = fft_data.velocities;
            obj.frequencies = fft_data.frequencies;
            %             obj.convert_to_cylindrical_coors();
        end
        function convert_to_cylindrical_coors(obj)
            [th, r, z] = cart2pol(obj.coors_rect(:,1), obj.coors_rect(:,2), obj.coors_rect(:,3));
            obj.coordinates_cyl = [th, r, z];
        end
        function fit_cylinder_to_geometry_pts(obj)
            [obj.cyl_axis, obj.cyl_radius, obj.cyl_zero] = obj.fit_cylinder_to_points(obj.coors_rect);
            obj.rotate_points();
            obj.translate_coordinates()
            obj.plot_cylinder();
        end
        function rotate_points(obj)
            R = obj.rotation_matrix(obj.cyl_axis);
            obj.coors_rect = (R*obj.coors_rect')';
            obj.disp_rect = (R*obj.disp_rect')';
            obj.vel_rect = (R*obj.vel_rect')';
        end
        function translate_coordinates(obj)
            xshift = mean(obj.coors_rect(:,1));
            yshift = mean(obj.coors_rect(:,2));
            zshift = min(obj.coors_rect(:,3));
            obj.coors_rect = obj.coors_rect - repmat([xshift, yshift, zshift], size(obj.coors_rect,1),1);
        end
        function plot_cylinder(obj)
            scatter3(obj.coors_rect(:,1), obj.coors_rect(:,2), obj.coors_rect(:,3))
        end
    end
    methods (Static = true)
        function R = rotation_matrix(axis)
            a = axis./ norm(axis);
            b= [0; 0; 1];
            v = cross(a,b);
            c = dot(a,b);
            vx = [0 -v(3) v(2);
                v(3) 0 -v(1);
                -v(2) v(1) 0];
            R = eye(3)+ vx + vx^2*(1/(1+c));
        end
        function [axisNormal, cylinderRadius, cylinderCenter] = fit_cylinder_to_points(points)
            % axis(x y z) rad center(x,y,z)
            
            x0 = calculate_guess(points);
            [ub, lb] = calculate_bounds(points);
            
            % ------------Linear constraints------------
            A = [];
            b = [];
            Aeq = [];
            beq = [];
            % ------------Pass Parameters into fmincon------------
            objectivFuntion = @(x) obj(x, points);
            constraintFunction = @(x) con(x, points);
            
            % ------------Call fmincon------------
            options = optimoptions(@fmincon, 'StepTolerance', 1e-10, 'MaxIterations', 1000, 'MaxFunctionEvaluations', 10000);
            [xopt, ~, ~, ~] = fmincon(objectivFuntion, x0, A, b, Aeq, beq, lb, ub, constraintFunction, options);
            axisNormal = xopt(1:3);
            cylinderRadius = xopt(4);
            cylinderCenter = xopt(5:7);
            % ------------Objective and Non-linear Constraints------------
            function [f, c, ceq] = objcon(x, points)
                %design variables
                normal = x(1:3)';
                radius = x(4);
                center = x(5:7)';
                %other analysis variables
                dist = zeros(length(points),1);
                for p = 1:length(points)
                    dist(p) = norm((center-points(p,:))-(dot((center-points(p,:)), normal)*normal));
                end
                error = sum((dist-radius).^2);
                %objective function
                f = error;%what are you minimizing
                %inequality constraints in a column vector c where c<=0
                c = 0;
                %equality contraints (ceq)
                ceq = norm(normal) - 1;
            end
            % ------------Separate obj/con (do not change)------------
            function [f] = obj(x, points)
                [f, ~, ~] = objcon(x, points);
            end
            function [c, ceq] = con(x, points)
                [~, c, ceq] = objcon(x, points);
            end
            function x0 = calculate_guess(points)
                x0= zeros(7,1);
                x0(3) = 1; % assume axis is in z direction
                x0(4) = (max(points(:,1))-min(points(:,1)))/2; %radius is half of distance between max x and min x
                x0(5:7) = mean(points);
            end
            function [ub, lb] = calculate_bounds(points)
                ub = ones(7,1);        lb = -ones(7,1);
                ub(4) = max(max(points)-min(points))/2;
                lb(4) = 1e-3;
                ub(5:7) = repmat(max(max(points)),1,3);
                lb(5:7) = repmat(min(min(points)),1,3);
            end
        end
    end
end

