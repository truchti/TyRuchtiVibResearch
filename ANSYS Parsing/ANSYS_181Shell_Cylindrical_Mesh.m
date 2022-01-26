classdef ANSYS_181Shell_Cylindrical_Mesh < handle
    properties
        elements
        nodes
    end
    properties (Hidden = true)
        extrapolatedToNodes = false;
        elementResultList = {'N11', 'N22', 'N12', 'M11', 'M22', 'M12', 'Q1', 'Q2'};
        % u is displacement in the rad direct
        % v is displacement in theta direction 
        % w is displacement in the long direction
        nodeResultList = {'u', 'v', 'w', 'udot', 'vdot', 'wdot'};
        Fitter
        DotSpline
        qt
        ql
        qtimag
        qlimag
        LMxx
        LMxy
        LNxx
        LNxy
        LQ
        TMxx
        TMxy
        TNxx
        TNxy
        TQ
    end
    methods
        function obj = ANSYS_181Shell_Cylindrical_Mesh(nodes, elements, nodeResultTypes)
            obj.nodes = nodes;
            obj.elements = elements;
            if nargin > 2
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
                [Xs, Ys, Zs] = obj.get_nodal_coordinates_from_element(e);
                nodeValues = obj.get_nodal_values_for_an_element(type, e);
                patch(Xs, Ys, Zs, real(scale*nodeValues));
            end
            colormap(gca, jet)
            view(90, 0)
            axis equal
            colorbar
        end
        function plot_flat_nodal_result(obj,type)
            for e = 1:length(obj.elements)
                [Xs, Ys, Zs] = obj.get_nodal_coordinates_from_element(e);
                thetas = atan2(Ys,Xs);
                Thetas = wrapTo2Pi(thetas);
                if max(Thetas) - min(Thetas) > 5
                    Thetas(Thetas== 0) = 2*pi;
                end
                rads = (Xs.^2 + Ys.^2).^(1/2);
                rad = mean(rads);
                nodeValues = obj.get_nodal_values_for_an_element(type, e);
                patch(Thetas*rad, Zs, -.005*ones(size(rads)), real(nodeValues), 'EdgeColor', 'none');
            end
            colormap(gca, gray)
            alpha(.9)
            axis equal
%             colorbar
        end
        function plot_geometry(obj)
            coors = zeros(length(obj.nodes),3);
            for n = 1:length(obj.nodes)
                coors(n, :) = obj.nodes(n).coors;
            end
            scatter3(coors(:,1), coors(:,2), coors(:,3))
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
        end
        function calculate_cyl_power_flow(obj)
            %resultants
            [Qt, Ql, Mt, Ml, Mtl, Nt, Nl, Ntl] = obj.get_all_resultants();
            [~, radii, ~] = obj.get_all_node_coordinates_in_cylindrical_coors();
            rad = mean(mean(radii));
            Nlt = Ntl;
            Mlt = Mtl;
            %velocities
            wdot = obj.get_all_node_value_by_type('rdot');
            thetadot = obj.get_all_node_value_by_type('thdot');
            vdot = rad*thetadot;
            vlong = obj.get_all_node_value_by_type('ldot');
            omegaTh = obj.get_all_node_value_by_type('betaDotTh');
            omegaL = obj.get_all_node_value_by_type('betaDotLong');
            %power components
            fullqtheta =  - Nt.*vdot - Ntl.*vlong -Qt.*wdot +Mtl.*omegaTh - Mt.*omegaL;
            fullqlong = -Ql.*wdot - Nlt.*thetadot - Nl.*vlong + Ml.*omegaTh - Mlt.*omegaL;
            obj.TQ =-Qt.*wdot;
            obj.LQ =  -Ql.*wdot;
            obj.TMxx = - Mt.*omegaL;
            obj.LMxx = Ml.*omegaTh;
            obj.TMxy = Mtl.*omegaTh;
            obj.LMxy = - Mlt.*omegaL;
            obj.TNxx =- Nt.*thetadot;
            obj.LNxx = - Nl.*vlong;
            obj.TNxy =- Ntl.*vlong;
            obj.LNxy =  - Nlt.*thetadot;
            obj.qt = 1/2*real(fullqtheta);
            obj.ql = 1/2*real(fullqlong);
            obj.qtimag = 1/2*imag(fullqtheta);
            obj.qlimag = 1/2*imag(fullqlong);
        end
        function plot_cyl_power_flow_part(obj, strng)
            if isempty(obj.qt)
                obj.calculate_cyl_power_flow;
            end
            [Thetas, rs, Zs] = obj.get_all_node_coordinates_in_cylindrical_coors();
            rad = mean(rs);
            switch strng
                case 'Q'
                    v1 = obj.TQ;
                    v2 = obj.LQ;
                case 'Mxx'                    
                    v1 = obj.TMxx;
                    v2 = obj.LMxx;
                case 'Mxy'                    
                    v1 = obj.TMxy;
                    v2 = obj.LMxy;
                case 'Nxx'                    
                    v1 = obj.TNxx;
                    v2 = obj.LNxx;
                case 'Nxy'
                    v1 = obj.TNxy;
                    v2 = obj.LNxy;
            end
            quiver(Thetas'*rad, Zs', v1, v2)
        end
        function a = plot_cyl_power_flow(obj)
            a = figure;
            if isempty(obj.qt)
                obj.calculate_cyl_power_flow;
            end
            [Thetas, Zs, rad] = obj.theta_longitudinal_coordinates_with_fixed_radius();
            quiver(gca, Thetas'*rad, Zs', obj.qt, obj.ql)
%             disps = obj.get_all_node_value_by_type('rad');
            hold on;
            obj.plot_flat_nodal_result('rad')
        end
        function plot_3D_cyl_power_flow(obj)
            if isempty(obj.qt)
                obj.calculate_cyl_power_flow();
            end
            [Thetas, rs, ~] = obj.get_all_node_coordinates_in_cylindrical_coors();
            rad = mean(rs);
            [X,Y,Z] = obj.get_all_node_coordinates();
            center = [mean(X), mean(Y)];
            height = max(Z)-min(Z);
            qx = obj.qt'.*sin(Thetas);
            qy = obj.qt'.*-cos(Thetas);
            qz = obj.ql';
            quiver3(X,Y,Z, qx, qy, qz)
            hold on;
            [X, Y, Z] = cylinder(rad*.99, 30);
            X = X-center(1);
            Y = Y-center(2);
            Z = height*Z;
            surf(X,Y,Z);
            alpha 0.2
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
            if ~obj.extrapolatedToNodes
                for r = 1:length(obj.elementResultList)
                    obj.extrapolate_single_resultant_to_nodes(obj.elementResultList{r});
                end
                obj.extrapolatedToNodes = true;
            end
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
            Thetas = wrapTo2Pi(Thetas);
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
        function [ths, longs, radius] = theta_longitudinal_coordinates_with_fixed_radius(obj)
            [ths, rads, longs] = obj.get_all_node_coordinates_in_cylindrical_coors();
            radius = mean(mean(rads));
        end
        function calculate_angular_rotation_velocities(obj)
            [ths, longs, rad] = obj.theta_longitudinal_coordinates_with_fixed_radius;
            %radial velocity
            rdot = obj.get_all_node_value_by_type('rdot');
%             OLdot = drad^2/(dzdt) = drdotdz
%             Othdot = dradial^2/dsdt = dradial^2/(rad&*dthdt) =
%             drdot/dth *1/a;
            % real and imaginary for Olong then real and imag for Otheta
            data = [real(rdot)', imag(rdot)', real(rdot/rad)', imag(rdot/rad)'];
            obj.Fitter = quinticBSplineSurfaceFitter([ths, longs], data, {"closed", "open"}, [15,15]);
            obj.Fitter.fit_spline_surfaces();
            obj.DotSpline = obj.Fitter.output_solved_spline_evaluator();
            %take derivatives wrt th and long vars
            %for each node evaluate derivatives and save value
            for n = 1:length(obj.nodes)
                
                                                                                   %param, param, fitting data, derive order of 1st and 2nd param
                BetaDotThetaReal= obj.DotSpline.evaluate_spline_number_at_parameters(ths(n), longs(n), 1, [0,1]); %e.g. dw/dl real
                BetaDotThetaImag= obj.DotSpline.evaluate_spline_number_at_parameters(ths(n), longs(n), 2, [0,1]); % dw/dl imag
                betaDTh =  complex(BetaDotThetaReal, BetaDotThetaImag);
                                                                                     %param, param, fitting data, derive order of 1st and 2nd param
                BetaDotLongReal= obj.DotSpline.evaluate_spline_number_at_parameters(ths(n), longs(n), 3, [1,0]); % d(w/a)/dth real
                BetaDotLongImag= obj.DotSpline.evaluate_spline_number_at_parameters(ths(n), longs(n), 4, [1,0]); % d(w/a)/dth real
                betaDL = complex(BetaDotLongReal, BetaDotLongImag);
                % save values to node
                obj.nodes(n).add_disp_or_vel('betaDotTh', betaDTh);
                obj.nodes(n).add_disp_or_vel('betaDotLong',betaDL);
            end
            
        end
    end
end