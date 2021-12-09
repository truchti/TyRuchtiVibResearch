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
    function [f] = obj(x, points)
        [f, ~, ~] = objcon(x, points);
    end
    function [c, ceq] = con(x, points)
        [~, c, ceq] = objcon(x, points);
    end
    function x0 = calculate_guess(points)
        x0= zeros(7,1);
        x0(3) = 1; % assume axis is in z direction
        x0(4) = (max(points(:,1))-min(points(:,1)))/x; %radius is half of distance between max x and min x
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