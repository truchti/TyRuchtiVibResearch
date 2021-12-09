classdef CylCoorParameterizer < handle
    properties
        dataCart
        dataCylindrical
    end
    properties (Hidden = true)
        upperBounds = [1.1,   1.1,  1.1, 0.5,    0.05,  0.05,  0.1];
        lowerBounds =  [-1.1, -1.1, -1.1, 0.01, -0.05, -0.05, -0.1];
        startValues = [0,   0,  1,  0.076,    0.0  0.0  0.0];
        cylAxis
        cylRadius
        cylCenter
        rotatedCart
        forceLocations
    end
    methods
        function obj = CylCoorParameterizer(SLDVdata, forceLocations)
            obj.dataCart = SLDVdata;
            if nargin>1
                obj.forceLocations = forceLocations;
            end
        end
        function calculate_cylindrical_data(obj)
            points = obj.dataCart.coordinates;
            [obj.cylAxis, obj.cylRadius, obj.cylCenter] = obj.fit_cylinder_to_points(points, obj.startValues, obj.lowerBounds, obj.upperBounds);
%             [obj.rotatedCart, obj.forceLocations] = obj.dataCart.translate_and_rotate_data(obj.cylCenter, obj.cylAxis, obj.forceLocations);
            [obj.rotatedCart] = obj.dataCart.translate_and_rotate_data(obj.cylCenter, obj.cylAxis);
            obj.convert_to_cylindrical_coordinates();       
%             obj.plot_to_verify();
        end
        function cylData = export_cylindrical_data(obj)
            if ~isempty(obj.dataCylindrical)
                cylData = obj.dataCylindrical;
            else
                warning('Cylindrical Data must be calculated before exporting');
            end
        end
    end
    methods (Hidden = true)
        function convert_to_cylindrical_coordinates(obj)
            points = obj.rotatedCart.coordinates;
            %put in 1 3 2 since y is the cylinder axis
            [theta, rad, h] = cart2pol(points(:,1), points(:,3), points(:,2)); 
            coors = [theta, rad, h];
            [reD, imD, reV, imV] = obj.convert_measurements_using_theta(theta);
            obj.dataCylindrical = SLDV_data_cyl_coors(coors, reD, imD, reV, imV, obj.rotatedCart.frequency);
        end
        function plot_to_verify(obj)
            a = obj.dataCylindrical;
            figure(1)
            a.flat_plot_displacement(1);
            title('theta Disp')
            figure(2)
            a.flat_plot_displacement(2)
            title('radial Disp')
            figure(3)
            a.flat_plot_displacement(3)
            title('Z disp')
            
            obj.plot_to_compare_disp_magnitudes()
        end
        function plot_to_compare_disp_magnitudes(obj)
            cartDisp = obj.dataCart.realDisp;
            normCart = sqrt(sum(cartDisp.^2,2));
            cylDisp = obj.dataCylindrical.realDisp;
            normCyl = sqrt(sum(cylDisp.^2,2));
            err = normCart- normCyl;
        end
        function [reD, imD, reV, imV] =convert_measurements_using_theta(obj, theta)
            reD = obj.rotatedCart.realDisp;
            imD = obj.rotatedCart.imagDisp;
            reV = obj.rotatedCart.realVel;
            imV = obj.rotatedCart.imagVel;
            for i = 1:length(theta)
                R = obj.calculate_theta_transform(theta(i));
                reD(i,:) = reD(i,:)*R;
                imD(i,:) = imD(i,:)*R;
                reV(i,:) = reV(i,:)*R;
                imV(i,:) = imV(i,:)*R;
            end  
        end
    end
    methods (Hidden = true, Static = true)
        function [cylAxis, cylRadius, cylCenter] =fit_cylinder_to_points(points, x0, lb, ub)
            if nargin <2
                x0 = [0,   0,  1,  0.075,    0.0  0.0  0.0];
            end
            if nargin <3
                ub = [1,   1,  1,  10,    0.1,  0.1,  0.1];
            end
            if nargin < 4
                lb = [-1, -1, -1, 0.01, -0.1, -0.1, -0.1];
            end
            % ------------Linear constraints------------
            A = [];
            b = [];
            Aeq = [];
            beq = [];
            % ------------Pass Parameters into fmincon------------
            objectivFuntion = @(x) objective(x, points);
            constraintFunction = @(x) constraint(x, points);
            % ------------Call fmincon------------
            options = optimoptions(@fmincon, 'StepTolerance', 1e-8, 'MaxIterations', 500, 'MaxFunctionEvaluations', 5000, 'Display', 'off');
            [xopt, ~, ~, ~] = fmincon(objectivFuntion, x0, A, b, Aeq, beq, lb, ub, constraintFunction, options);
            cylAxis = xopt(1:3);
            cylRadius = xopt(4);
            cylCenter = xopt(5:7);
            % ------------Objective and Non-linear Constraints------------
            function [f, c, ceq] = objcon(x, points)
                %design variables
                normal = x(1:3);
                radius = x(4);
                center = x(5:7);
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
            function [f] = objective(x, points)
                [f, ~, ~] = objcon(x, points);
            end
            function [c, ceq] = constraint(x, points)
                [~, c, ceq] = objcon(x, points);
            end
        end
        function R = calculate_theta_transform(theta)
            R = [sin(theta) cos(theta) 0;
                 0 0 1;
                 -cos(theta) sin(theta) 0];
        end
    end
end