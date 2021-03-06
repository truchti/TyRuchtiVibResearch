classdef ANSYS_181Shell_Mesh < handle
    properties
        elements
        nodes
        isCylinder
    end
    properties (Hidden = true)
        extrapolatedToNodes = false;
        elementResultList = {'N11', 'N22', 'N12', 'M11', 'M22', 'M12', 'Q1', 'Q2'};
        nodeResultList = {'u', 'v', 'w', 'udot', 'vdot', 'wdot'};
        Fitter
        DotSpline
    end
    methods
        function obj = ANSYS_181Shell_Mesh(nodes, elements, eleResultTypes, nodeResultTypes)
            obj.nodes = nodes;
            obj.elements = elements;
            if nargin >2
                obj.elementResultList = eleResultTypes;
            end
            if nargin > 3
                obj.nodeResultList = nodeResultTypes;
            end
            obj.extrapolate_all_resultants_to_nodes();
        end
        function plot_element_result(obj, type)
            for ele = 1:length(obj.elements)
                resultantValue = obj.elements(ele).get_resultant(type);
                [Xs, Ys, Zs] = obj.get_nodal_coordinates_from_element(ele);
                patch(Xs, Ys, Zs, real(resultantValue))
            end
            colormap(gca, jet)
            view(0,35)
            axis equal
            colorbar
        end
        function plot_nodal_result(obj, type, scale)
            if nargin < 3
                scale = 1;
            end
            for e = 1:length(obj.elements)
                [Xs, Ys, ~] = obj.get_nodal_coordinates_from_element(e);
                nodeValues = obj.get_nodal_values_for_an_element(type, e);
                patch(Xs, Ys, real(nodeValues), real(scale*nodeValues));
            end
            colormap(gca, jet)
            view(0,90)
            axis equal
            colorbar
            title(type)
        end
        function plot_imag_nodal_result(obj, type, scale)
            if nargin < 3
                scale = 1;
            end
            for e = 1:length(obj.elements)
                [Xs, Ys, ~] = obj.get_nodal_coordinates_from_element(e);
                nodeValues = obj.get_nodal_values_for_an_element(type, e);
                patch(Xs, Ys, imag(nodeValues), imag(scale*nodeValues));
            end
            colormap(gca, jet)
            view(0,90)
            axis equal
            colorbar
            title(['Imaginary ' type])
        end
        function calculate_and_plot_plate_power_flow(obj, density)
            if nargin < 2
                density = 1;
            end
            [xs, ys, ~] = obj.get_all_node_coordinates();
            [Q1, Q2, M11, M22, M12, ~, ~, ~] = obj.get_all_resultants();
            wdot = obj.get_all_node_value_by_type('wdot');
            thetaDotX = obj.get_all_node_value_by_type('thetaDotX');
            thetaDotY = obj.get_all_node_value_by_type('thetaDotY');
            % calcluate power components
            qx = 1/2* real(Q1.*conj(wdot) + M11.*conj(thetaDotY) + M12.*conj(thetaDotX));
            qy = 1/2* real(Q2.*conj(wdot) + M22.*conj(thetaDotX) + M12.*conj(thetaDotY));
            %plot power flow vectors
            quiver(xs(1:density:end)', ys(1:density:end)', qx(1:density:end), qy(1:density:end),'LineWidth',.75,'MaxHeadSize',1.5,'Color','k')
            axis equal
            title('Power flow from ANSYS Results')
        end
    end
    methods(Hidden = true)
        function [Xs, Ys, Zs] = get_nodal_coordinates_from_element(obj, elementNumber)
            ele = obj.elements(elementNumber);
            allCoors = [];
            for j = 1:length(ele.nodeNums)
                allCoors = [allCoors; obj.nodes(ele.nodeNums(j)).coors];
            end
            if nargout == 1
                Xs = allCoors;
            else
                Xs = allCoors(:,1);
                Ys = allCoors(:,2);
                Zs = allCoors(:,3);
            end
        end
        function values = get_nodal_values_for_an_element(obj, type, element)
            ele = obj.elements(element);
            values = [];
            if any(strcmp(obj.elementResultList,type))
                for i = 1:length(ele.nodeNums)
                    values = [values; obj.nodes(ele.nodeNums(i)).get_resultant(type)]; %#ok<*AGROW>
                end
            else
                for i = 1:length(ele.nodeNums)
                    values = [values; obj.nodes(ele.nodeNums(i)).get_disp_or_vel(type)];
                end
            end
        end
        function extrapolate_all_resultants_to_nodes(obj)
            for r = 1:length(obj.elementResultList)
                obj.extrapolate_single_resultant_to_nodes(obj.elementResultList{r});
            end
            obj.extrapolatedToNodes = true;
        end
        function extrapolate_single_resultant_to_nodes(obj, type)
            for n = 1:length(obj.nodes)
                AEs = obj.nodes(n).associatedElements;
                temp = 0;
                for el = 1:length(AEs)
                    temp = temp+ obj.elements(AEs(el)).get_resultant(type);
                end
                obj.nodes(n).set_extrapolated_resultant(type, temp/length(AEs) )
            end
        end
        function [Xs, Ys, Zs] = get_all_node_coordinates(obj)
            Coors = zeros(length(obj.nodes),3);
            for n = 1:length(obj.nodes)
                Coors(n,:) = obj.nodes(n).coors;
            end
            Xs = Coors(:,1);
            Ys = Coors(:,2);
            Zs = Coors(:,3);
        end
        function [Thetas, Radii, Zs] = get_all_node_coordinates_in_cylindrical_coors(obj)
            [Xs, Ys, Zs] = obj.get_all_node_coordinates();
            [Thetas, Radii, Zs] = cart2pol(Xs, Ys, Zs);
        end
        function values = get_all_node_value_by_type(obj, type)
            values = zeros(size(obj.nodes));
            for n = 1:length(obj.nodes)
                values(n) = obj.nodes(n).get_value_by_type(type);
            end
        end
        function [Q1, Q2, M11, M22, M12, N11, N22, N12] = get_all_resultants(obj)
            Q1 = obj.get_all_node_value_by_type('Q1');
            Q2 = obj.get_all_node_value_by_type('Q2');
            M11 = obj.get_all_node_value_by_type('M11');
            M22 = obj.get_all_node_value_by_type('M22');
            M12 = obj.get_all_node_value_by_type('M12');
            N11 = obj.get_all_node_value_by_type('N11');
            N22 = obj.get_all_node_value_by_type('N22');
            N12 = obj.get_all_node_value_by_type('N12');
        end
        function calculate_theta_dots(obj)          
            [Xs, Ys, ~] =obj.get_all_node_coordinates();
            values = obj.get_all_node_value_by_type('wdot');
            data = [real(values)' imag(values)'];
            % create spline surface and interpolate
            obj.Fitter = quinticBSplineSurfaceFitter([Xs,Ys], data, {"open","open"}, [15,15]);
            obj.Fitter.fit_spline_surfaces();
            obj.DotSpline = obj.Fitter.output_solved_spline_evaluator();
            % take deriviatives wrt x and y
            for n= 1:length(obj.nodes)
                thetaDotXReal = obj.DotSpline.evaluate_spline_number_at_parameters(Xs(n), Ys(n), 1, [0,1]);
                thetaDotXImag = obj.DotSpline.evaluate_spline_number_at_parameters(Xs(n), Ys(n), 2, [0,1]);
                thetaDotX = complex(thetaDotXReal, thetaDotXImag);
                thetaDotYReal = obj.DotSpline.evaluate_spline_number_at_parameters(Xs(n), Ys(n), 1, [1,0]);
                thetaDotYImag = obj.DotSpline.evaluate_spline_number_at_parameters(Xs(n), Ys(n), 2, [1,0]);
                thetaDotY = complex(thetaDotYReal, thetaDotYImag);
                %save values
                obj.nodes(n).add_disp_or_vel('thetaDotX', thetaDotX)
                obj.nodes(n).add_disp_or_vel('thetaDotY', thetaDotY)
            end
        end
    end
end