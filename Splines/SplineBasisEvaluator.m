classdef SplineBasisEvaluator < handle
    %% This class evaluates the quintic Spline Basis functions for a given knot Vector and parameter
    properties
        knotVector
        numberOfControlPoints
        type
    end
    properties (Access = private)
        N = [ 1/120, -1/24,  1/12, -1/12,  1/24, -1/120;
                0,  1/24,  -1/6,   1/4,  -1/6,   1/24;
                0,  1/12,  -1/6,     0,   1/6,  -1/12;
                0,  1/12,   1/6,  -1/2,   1/6,   1/12;
                0,  1/24,  5/12,     0, -5/12,  -1/24;
                0, 1/120, 13/60, 11/20, 13/60,  1/120];
        D = [0     0     0     0     0     0;
             5     0     0     0     0     0;
             0     4     0     0     0     0;
             0     0     3     0     0     0;
             0     0     0     2     0     0;
             0     0     0     0     1     0];
        order
    end
    methods
        function obj = SplineBasisEvaluator(order, knotVector, numberOfControlPoints, type)
            obj.knotVector = knotVector;
            obj.numberOfControlPoints = numberOfControlPoints;
            obj.type = type;
            obj.order = order;
        end
        function values = evaluate_at_parameter(obj, parameter)
            if strcmp(obj.type, 'open')
                values = obj.open_loop_basis_eval(parameter);
            elseif strcmp(obj.type, 'closed')
                values = obj.closed_loop_basis_eval(parameter);
            end
        end
        function [values, d1, d2, d3] = evaluate_basis_and_derivatives(obj, parameter)
            if strcmp(obj.type, 'open')
                [values, d1, d2, d3] = obj.open_loop_basis_derivatives(parameter);
            elseif strcmp(obj.type, 'closed')
                [values, d1, d2, d3] = obj.closed_loop_basis_derivatives(parameter);
            end
        end
    end
    methods (Hidden = true)
        function values = open_loop_basis_eval(obj, parameter)
            % calculate the values of all basis functions at a given parameter value
            % requires the order of basis functions parameter value, number of points and the knot vector
            npc = obj.numberOfControlPoints+obj.order;
            assert(length(obj.knotVector) <= npc)
            for pt = 1:length(parameter)
                % this evaluates the basis
                temp = zeros(1, npc-1);
                for i = 1:npc-1
                    if parameter(pt) >= obj.knotVector(i) && parameter(pt) < obj.knotVector(i+1)
                        temp(i) = 1;
                    else
                        temp(i) = 0;
                    end
                end
                for k = 2:obj.order
                    for i = 1:npc-k
                        if temp(i) ~=0 && (obj.knotVector(i+k-1)-obj.knotVector(i)) ~= 0
                            d = ((parameter(pt)-obj.knotVector(i))*temp(i))/(obj.knotVector(i+k-1)-obj.knotVector(i));
                        else
                            d = 0;
                        end
                        if temp(i+1) ~= 0 && (obj.knotVector(i+k)-obj.knotVector(i+1)) ~= 0
                            e =((obj.knotVector(i+k)-parameter(pt))*temp(i+1))/(obj.knotVector(i+k)-obj.knotVector(i+1));
                        else
                            e = 0;
                        end
                        temp(i) = d+e;
                    end
                end
                if parameter(pt) >= obj.knotVector(npc)
                    temp(obj.numberOfControlPoints) = 1;
                end
                if ~nnz(temp)
                    temp(1) = 1;
                end
                if ~exist('values','var')
                    values = zeros(length(parameter), obj.numberOfControlPoints);
                    values(1,:) = temp(1:obj.numberOfControlPoints);
                end
                values(pt,:) = temp(1:obj.numberOfControlPoints);
            end
        end
        function values = closed_loop_basis_eval(obj, parameter)
            values = zeros(length(parameter), obj.numberOfControlPoints);
            for param_index = 1:length(parameter)
                %find knots that surround parameter value ( if parameter equals knot value
                %take that knot and the next knot to be the surrounding knots
                leftknotIndex = find(parameter(param_index) > obj.knotVector, 1, 'last');
                if isempty(leftknotIndex)
                    leftknotIndex = 1;
                end
                rightknotIndex = leftknotIndex +1;
                u = (parameter(param_index)-obj.knotVector(leftknotIndex))/(obj.knotVector(rightknotIndex)-obj.knotVector(leftknotIndex));
                % calculate local basis values
                localBu = obj.quintic_periodic_basis_matrix(u);
                % pad localBu to full length
                Bu = padarray(localBu,[0 obj.numberOfControlPoints-length(localBu)], 0,'post');
                % shift to correct location
                shift =leftknotIndex - obj.order;
                values(param_index, :) = circshift(Bu, -shift);
                % put basis values into correct part of full basis values vector
            end
        end
        function [values, d1, d2, d3] = open_loop_basis_derivatives(obj, parameter)
            npts = obj.numberOfControlPoints;
            values = zeros(length(parameter), npts);
            d1 = zeros(length(parameter), npts);
            d2 = zeros(length(parameter), npts);
            d3 = zeros(length(parameter), npts);
            npc = obj.order + obj.numberOfControlPoints;
            for pt = 1:length(parameter)
                temp = zeros(npc,1);
                temp1 = zeros(npc,1);
                temp2 = zeros(npc,1);
                temp3 = zeros(npc,1);
                for i = 1:npc-1
                    if parameter(pt) >= obj.knotVector(i) && parameter(pt)<obj.knotVector(i+1)
                        temp(i) = 1;
                    else
                        temp(i) = 0;
                    end
                end
                if parameter(pt) == obj.knotVector(npts+1)
                    temp(npts) = 1;
                    temp(npts+1) = 0;
                end
                for k = 2:obj.order
                    for i = 1:npc-obj.order
                        if temp(i) ~= 0
                            b1 = ((parameter(pt)-obj.knotVector(i))*temp(i))/(obj.knotVector(i+k-1)-obj.knotVector(i));
                            f1 = temp(i)/(obj.knotVector(i+k-1)-obj.knotVector(i));
                        else
                            b1 = 0;
                            f1 = 0;
                        end
                        if temp(i+1) ~= 0
                            b2 = ((obj.knotVector(i+k)-parameter(pt))*temp(i+1))/(obj.knotVector(i+k)-obj.knotVector(i+1));
                            f2 = -temp(i+1)/(obj.knotVector(i+k)-obj.knotVector(i+1));
                        else
                            b2 = 0;
                            f2 = 0;
                        end
                        if temp1(i) ~= 0
                            f3 = (parameter(pt)-obj.knotVector(i))*temp1(i)/(obj.knotVector(i+k-1)-obj.knotVector(i));
                            s1 = 2*temp1(i)/(obj.knotVector(i+k-1)-obj.knotVector(i));
                        else
                            f3 = 0;
                            s1 = 0;
                        end
                        if temp1(i+1) ~= 0
                            f4 = (obj.knotVector(i+k)-parameter(pt))*temp1(i+1)/(obj.knotVector(i+k)-obj.knotVector(i+1));
                            s2 = -2*temp1(i+1)/(obj.knotVector(i+k)-obj.knotVector(i+1));
                        else
                            f4 = 0;
                            s2 = 0;
                        end
                        if temp2(i) ~= 0
                            s3 = (parameter(pt)-obj.knotVector(i))*temp2(i)/(obj.knotVector(i+k-1)-obj.knotVector(i));
                            t1 = 2*temp2(i)/(obj.knotVector(i+k-1)-obj.knotVector(i));
                        else
                            s3 = 0;
                            t1 = 0;
                        end
                        if temp2(i+1) ~= 0
                            s4 = (obj.knotVector(i+k)-parameter(pt))*temp2(i+1)/(obj.knotVector(i+k)-obj.knotVector(i+1));
                            t2 = -2*temp2(i+1)/(obj.knotVector(i+k)-obj.knotVector(i+1));
                        else
                            s4 = 0;
                            t2 = 0;
                        end
                        if temp3(i) ~= 0
                            t3 = (parameter(pt)-obj.knotVector(i))*temp3(i)/(obj.knotVector(i+k-1)-obj.knotVector(i));
                        else
                            t3 = 0;
                        end
                        if temp3(i+1) ~= 0
                            t4 = (obj.knotVector(i+k)-parameter(pt))*temp3(i+1)/(obj.knotVector(i+k)-obj.knotVector(i+1));
                        else
                            t4 = 0;
                        end
                        temp(i) = b1+b2;
                        temp1(i) = f1+f2+f3+f4; %first derivatives
                        temp2(i) = s1+s2+s3+s4; %second derivatives
                        temp3(i) = t1+t2+t3+t4; %third derivatives
                    end
                end
                values(pt,:) = temp(1:npts);
                d1(pt,:) = temp1(1:npts);
                d2(pt,:) = temp2(1:npts);
                d3(pt,:) = temp3(1:npts);
            end            
        end
        function [values, d1, d2, d3] = closed_loop_basis_derivatives(obj, parameter)
            values = zeros(length(parameter), obj.numberOfControlPoints);
            d1 = zeros(length(parameter), obj.numberOfControlPoints);
            d2 = zeros(length(parameter), obj.numberOfControlPoints);
            d3 = zeros(length(parameter), obj.numberOfControlPoints);
            for param_index = 1:length(parameter)
                %find knots that surround parameter value ( if parameter equals knot value
                %take that knot and the next knot to be the surrounding knots
                leftknotIndex = find(parameter(param_index) > obj.knotVector, 1, 'last');
                if isempty(leftknotIndex)
                    leftknotIndex = 1;
                end
                rightknotIndex = leftknotIndex +1;
                u = (parameter(param_index)-obj.knotVector(leftknotIndex))/(obj.knotVector(rightknotIndex)-obj.knotVector(leftknotIndex));
                % calculate local basis values
                [localBu, ld1, ld2, ld3] = obj.quintic_periodic_basis_matrix_and_derivatives(u);
                % pad localBu to full length
                Bu = padarray(localBu,[0 obj.numberOfControlPoints-length(localBu)], 0,'post');
                deriv1 = padarray(ld1, [0 obj.numberOfControlPoints-length(localBu)], 0,'post');
                deriv2 = padarray(ld2, [0 obj.numberOfControlPoints-length(localBu)], 0,'post');
                deriv3 = padarray(ld3, [0 obj.numberOfControlPoints-length(localBu)], 0,'post');
                % shift to correct location
                shift =leftknotIndex - obj.order;
                % put basis values into correct part of full basis values vector
                values(param_index, :) = circshift(Bu, -shift);
                d1(param_index, :) = circshift(deriv1, -shift);
                d2(param_index, :) = circshift(deriv2, -shift);
                d3(param_index, :) = circshift(deriv3, -shift); 
            end
        end
        function [F] = quintic_periodic_basis_matrix(obj, t_star)
            T = [t_star.^5, t_star.^4, t_star.^3, t_star.^2, t_star, t_star.^0];
            F = T*obj.N;
        end
        function [F, d1, d2, d3] = quintic_periodic_basis_matrix_and_derivatives(obj, t_star)
            T = [t_star.^5, t_star.^4, t_star.^3, t_star.^2, t_star, t_star.^0];
            F = T*obj.N;
            d1 = (T*obj.D)*obj.N;
            d2 = (T*obj.D^2)*obj.N;
            d3 = (T*obj.D^3)*obj.N;
        end
    end
end
