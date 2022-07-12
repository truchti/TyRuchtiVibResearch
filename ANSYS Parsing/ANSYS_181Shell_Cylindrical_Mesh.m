classdef ANSYS_181Shell_Cylindrical_Mesh < handle
    properties
        elements
        nodes
        splineSegmentNumbers = [15,15]
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
        Splines
        qt; ql; qtimag; qlimag;
        LMxx; LMxy; LNxx; LNxy; LQ;
        TMxx; TMxy; TNxx; TNxy; TQ;
        nodeIndxMap
        thetaShift = 0
        triData
        valid_types = {'w', 'wdot', 'M11', 'M22', 'M12', 'N11', 'N22', 'N12', 'Q1', 'Q2', 'betaDotX', 'betaDotY'}
    end
    methods
        function obj = ANSYS_181Shell_Cylindrical_Mesh(nodes, elements, nodeResultTypes, IndxMap)
            obj.nodes = nodes;
            obj.elements = elements;
            if nargin > 2
                obj.nodeResultList = nodeResultTypes;
                obj.nodeIndxMap = IndxMap;
            end
            obj.extrapolate_all_resultants_to_nodes();
        end
        %% plotting
        function fast_plot(obj, compnt, realPart)
            if isempty(obj.triData)
                obj.fast_plot_data_setup();
            end
            if nargin < 3
                realPart = true;
            end
            ntype = obj.convert_type_to_node_type(compnt);
            X = obj.triData.coors(:,1);
            Y = obj.triData.coors(:,2);
            V = obj.triData.(ntype);
            if realPart
                V = real(V);
            else
                V = imag(V);
            end
            trisurf(obj.triData.connects, X,Y,V);
            shading interp
            colormap jet; colorbar;
            axis equal; view(0,90);
        end
        function plot_component(obj, compnt, realPart, forcePatches)
            % plot the patches that represents each element
            if nargin < 3 
                realPart = true; 
            end
            if nargin < 4
                forcePatches = false;
            end
            if forcePatches
                obj.plot_flat(compnt, realPart);
            else
                obj.fast_plot(compnt, realPart);
            end
        end
        function vMax = plot_RI_component(obj, compnt)
            subplot(2,1,1);
            obj.plot_component(compnt, true)
%             title(['Real ', compnt]);
            
            h = subplot(2,1,2);
            obj.plot_component(compnt, false)
%             title(['Imag ', compnt]);
            vMax = h.CLim;
        end
        function plot_cyl_component(obj, compnt)
            obj.plot_nodal_result(compnt)
        end
        function plot_element_result(obj, type)
            for ele = 1:length(obj.elements)
                resultantValue = obj.elements(ele).get_resultant(type);
                [Xs, Ys, Zs] = obj.get_nodal_coordinates_from_element(ele);
                patch(Xs, Ys, Zs, real(resultantValue))
            end
            colormap(gca, jet); colorbar
            view(0,35); axis equal
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
        function plot_flat(obj, type, isReal)
            for e = 1:length(obj.elements)
                [Xs, Ys, Zs] = obj.get_nodal_coordinates_from_element(e);
                Cir = obj.convert_xy_to_circumference(Xs,Ys);
                Hei = Zs;
                nodeValues = obj.get_nodal_values_for_an_element(type, e);
                if isReal
                    Values = real(nodeValues);
                else
                    Values = imag(nodeValues);
                end
                patch(Cir, Hei, 0*ones(size(Cir)), Values, 'EdgeColor', 'none');
            end
            xlabel('Tangential'); ylabel('Height')
            if isReal
                title(['Real ', type]);
            else
                title(['Imag ', type]);
            end
            colormap(gca, jet); axis equal; colorbar
        end
        function power_flow_component_magnitudes(obj)
            clc
            fprintf('Magnitude of Power Contributions by component in Longitudinal Direction:\n')
            ShearMag = max(max(abs(obj.LQ))); NxxMag = max(max(abs(obj.LNxx)));
            NxyMag = max(max(abs(obj.LNxy))); MxxMag = max(max(abs(obj.LMxx)));
            MxyMag = max(max(abs(obj.LMxy)));
            fprintf('Shear:\t\t\t %d\n', ShearMag); fprintf('Normal:\t\t\t %d\n', NxxMag)
            fprintf('CrossNormal:\t %d\n', NxyMag); fprintf('Moment:\t\t\t %d\n', MxxMag)
            fprintf('CrossMoment:\t %d\n\n', MxyMag)
            ShearMag = max(max(abs(obj.TQ))); NxxMag = max(max(abs(obj.TNxx)));
            NxyMag = max(max(abs(obj.TNxy))); MxxMag = max(max(abs(obj.TMxx)));
            MxyMag = max(max(abs(obj.TMxy)));
            fprintf('Magnitude of Power Contributions by component in Tangential Direction:\n')
            fprintf('Shear:\t\t\t %d\n', ShearMag); fprintf('Normal:\t\t\t %d\n', NxxMag); 
            fprintf('CrossNormal:\t %d\n', NxyMag); fprintf('Moment:\t\t\t %d\n', MxxMag);
            fprintf('CrossMoment:\t %d\n\n', MxyMag)
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
            colormap(gca, jet); colorbar
            view(0,90); axis equal
        end
        %% Power flow
        function calculate_cyl_power_flow(obj)
            %resultants
            [Qt, Ql, Mt, Ml, Mtl, Nt, Nl, Ntl] = obj.get_all_resultants();
            [~, ~, rad] = obj.theta_longitudinal_coordinates_with_fixed_radius();
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
            fullqlong = -Ql.*wdot - Nlt.*vdot - Nl.*vlong + Ml.*omegaTh - Mlt.*omegaL;
            obj.TQ =-Qt.*wdot;
            obj.LQ =  -Ql.*wdot;
            obj.TMxx = - Mt.*omegaL;
            obj.LMxx = Ml.*omegaTh;
            obj.TMxy = Mtl.*omegaTh;
            obj.LMxy = - Mlt.*omegaL;
            obj.TNxx =- Nt.*thetadot*rad;
            obj.LNxx = - Nl.*vlong;
            obj.TNxy =- Ntl.*vlong;
            obj.LNxy =  - Nlt.*thetadot*rad;
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
                    ttl = 'Shear Contribution';
                case 'Mxx'                    
                    v1 = obj.TMxx;
                    v2 = obj.LMxx;
                    ttl = 'Bending Moment Contribution';
                case 'Mxy'                    
                    v1 = obj.TMxy;
                    v2 = obj.LMxy;
                    ttl = 'CrossBending Moment Contribution';
                case 'Nxx'                    
                    v1 = obj.TNxx;
                    v2 = obj.LNxx;
                    ttl = 'Normal Contribution';
                case 'Nxy'
                    v1 = obj.TNxy;
                    v2 = obj.LNxy;
                    ttl = 'Cross Normal Contribution';
            end
            quiver(Thetas'*rad, Zs', real(v1), real(v2))
            title(ttl);
        end
        function a = plot_cyl_power_flow(obj, backPlot)
            if nargin< 2  
                backPlot = false;
            end
            a = figure;
            if isempty(obj.qt)
                obj.calculate_cyl_power_flow;
            end
            [Thetas, Zs, rad] = obj.theta_longitudinal_coordinates_with_fixed_radius();
            quiver(gca, Thetas'*rad, Zs', obj.qt, obj.ql)
            if backPlot
                hold on;
                obj.plot_flat_nodal_result('rad')
                hold off;
            end
            hold on
%             scatter([0, pi*.02], [.135 .025],100, 'xr')
            hold off;
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
        %% accessing node/element values
        function [Xs, Ys, Zs] = get_nodal_coordinates_from_element(obj, elementNumber)
            ele = obj.elements(elementNumber);
            allCoors = [];
            for j = 1:length(ele.nodeNums)
                indx = obj.nodeIndxMap(ele.nodeNums(j));
                allCoors = [allCoors; obj.nodes(indx).coors];
            end
            if nargout == 1
                Xs = allCoors;
            else
                Xs = allCoors(:,1);
                Ys = allCoors(:,2);
                Zs = allCoors(:,3);
            end
        end
        function values = get_nodal_values_for_an_element(obj, type, eleNum)
            ele = obj.elements(eleNum);
            nodeType = obj.convert_type_to_node_type(type);
            values = [];
            for i = 1:length(ele.nodeNums)
                indx = obj.nodeIndxMap(ele.nodeNums(i));
                values = [values; obj.nodes(indx).(nodeType)]; %#ok<*AGROW>
            end
        end
        function extrapolate_all_resultants_to_nodes(obj)
            if ~obj.extrapolatedToNodes
                fprintf('Extrapolating Resultants: ');
                prev = '';
                for r = 1:length(obj.elementResultList)
                    s = obj.elementResultList{r};
                    fprintf(repmat('\b', 1, length(prev))); fprintf(s);
                    obj.extrapolate_single_resultant_to_nodes(s);
                    prev = s;
                end
                fprintf('\n')
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
            Thetas = wrapTo2Pi(Thetas+obj.thetaShift);
        end
        function values = get_all_node_value_by_type(obj, type)
            values = zeros(size(obj.nodes));
            for n = 1:length(obj.nodes)
                values(n) = obj.nodes(n).get_value_by_type(type);
            end
        end
        %% setup plotter
        function fast_plot_data_setup(obj)
            wb = waitbar(0, 'Configuring Data');
            [T, R, Z] = obj.get_all_node_coordinates_in_cylindrical_coors();
            obj.triData.coors = [T*mean(R), Z];
            dt = delaunayTriangulation(T,Z);
            obj.triData.connects = dt.ConnectivityList;
            waitbar(1/11, wb);
            for i = 1:length(obj.valid_types)
                vtype = obj.valid_types{i};
                typ = obj.convert_type_to_node_type(vtype);
                obj.triData.(typ) = obj.get_all_node_value_by_type(typ);
                waitbar((i+1)/11, wb);
            end
            close(wb)
        end 
        %% get resultants
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
        %% angular velocities
        function [ths, longs, radius] = theta_longitudinal_coordinates_with_fixed_radius(obj)
            [ths, rads, longs] = obj.get_all_node_coordinates_in_cylindrical_coors();
            radius = mean(mean(rads));
        end
        function calculate_angular_rotation_velocities(obj)
            [ths, longs, rad] = obj.theta_longitudinal_coordinates_with_fixed_radius;
            rdot = obj.get_all_node_value_by_type('rdot');
            % real and imaginary for Olong then real and imag for Otheta
            data = [real(rdot)', imag(rdot)', real(rdot/rad)', imag(rdot/rad)'];
            obj.Fitter = quinticBSplineSurfaceFitter([ths, longs], data, {"closed", "open"}, obj.splineSegmentNumbers);
            obj.Fitter.fit_spline_surfaces();
            obj.DotSpline = obj.Fitter.output_solved_spline_evaluator();
            %take derivatives wrt th and long vars
            %for each node evaluate derivatives and save value
            for n = 1:length(obj.nodes)         %param, param, fitting data, derive order of 1st and 2nd param
                BetaDotThetaReal= obj.DotSpline.evaluate_spline_number_at_parameters(ths(n), longs(n), 1, [0,1]); %e.g. dw/dl real
                BetaDotThetaImag= obj.DotSpline.evaluate_spline_number_at_parameters(ths(n), longs(n), 2, [0,1]); % dw/dl imag
                betaDTh =  complex(BetaDotThetaReal, BetaDotThetaImag);
                                                                                     %param, param, fitting data, derive order of 1st and 2nd param
                BetaDotLongReal= obj.DotSpline.evaluate_spline_number_at_parameters(ths(n), longs(n), 1, [1,0]); % d(w/a)/dth real
                BetaDotLongImag= obj.DotSpline.evaluate_spline_number_at_parameters(ths(n), longs(n), 2, [1,0]); % d(w/a)/dth real
                betaDL = complex(BetaDotLongReal, BetaDotLongImag);
                % save values to node
                obj.nodes(n).add_disp_or_vel('betaDotTh', betaDTh);
                obj.nodes(n).add_disp_or_vel('betaDotLong',betaDL);
            end            
        end
        function calc_derivatives_with_splines(obj)
            [ths, longs, radius] = obj.theta_longitudinal_coordinates_with_fixed_radius;
            longs = longs - min(longs); % ensure knot vector starts at zero
            %radial velocity
            r = obj.get_all_node_value_by_type('u'); % the radial component is stored in u in each node
            rdot = obj.get_all_node_value_by_type('udot'); % the radial component is stored in u in each node
            % real and imaginary for Olong then real and imag for Otheta
            data = [real(r)', imag(r)', real(rdot)', imag(rdot)']; 
            obj.Fitter = quinticBSplineSurfaceFitter([ths, longs], data, {"closed", "open"}, obj.splineSegmentNumbers);
            obj.Fitter.fit_spline_surfaces();
            obj.Splines = obj.Fitter.output_solved_spline_evaluator();
        end
    end
    methods (Static = true)
        function [Cir] = convert_xy_to_circumference(Xs, Ys)
            thetas = atan2(Ys,Xs);
            Thetas = wrapTo2Pi(thetas);
            if max(Thetas) - min(Thetas) > 3
                Thetas(Thetas<3) = Thetas(Thetas<3)+ 2*pi;
            end
            rad = mean((Xs.^2 + Ys.^2).^(1/2));
            Cir = Thetas* rad;
        end
            
        function nodeType = convert_type_to_node_type(resType)
            switch resType
                case {'N1', 'N11', 'Nx', 'Nl', 'NL'}
                    nodeType = 'N11';
                case {'N2', 'N22', 'Ntheta', 'N0', 'NO', 'Nth'}
                    nodeType = 'N22';
                case {'N12', 'N21','N0L', 'Nthx', 'NL0', 'NLO', 'NOL' 'Nxth'}
                    nodeType = 'N12';
                case {'M1', 'M11', 'Mx', 'Ml', 'ML'}
                    nodeType = 'M11';
                case {'M2', 'M22', 'Mtheta', 'M0', 'MO', 'Mth'}
                    nodeType = 'M22';
                case {'M12', 'M21','M0L', 'Mthx', 'ML0', 'MLO', 'MOL' 'Mxth'}
                    nodeType = 'M12';
                case {'Q1', 'Ql', 'Qx', 'QL'}
                    nodeType = 'Q1';
                case {'Q2', 'QO', 'Q0', 'Qth', 'Qtheta'}
                    nodeType = 'Q2';
                case {'0dotL', '0dotx', 'OdotL', 'Odotx', 'ThetaDotL', 'thetaDotL', 'BL',}
                    nodeType = 'betaDotX';
                case {'0dot0', '0doty', 'OdotO', 'Odoty', 'ThetaDotO', 'thetaDotO', 'ThetaDot0', 'thetaDot0', 'B0'}
                    nodeType = 'betaDotY';
                case {'r', 'rad', 'w'}
                    nodeType = 'u';
                case {'rdot', 'raddot', 'wdot'}
                    nodeType = 'udot';
                case {'th', 'theta', '0', 'O'}
                    nodeType = 'v';
                case {'l', 'long', 'x'}
                    nodeType ='w';
                otherwise
                    nodeType = resType;                    
            end
        end
    end
end