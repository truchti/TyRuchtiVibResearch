classdef Cylinder_Power_Flow <handle
    % This class can take in a SLDV_parser, a SLDV_data or no input in the
    % no input case it creates a default SLDV_FastScan parser. If a parser
    % is input then the data is extracted using the parser. If the raw data
    % is processed separately the processed data can be passed instead of a
    % parser
    
    % in all cases 'w' stands for radial deflection 'v' for theta
    % deflection and 'u' for longitudinal deflection
    properties 
        laserData
        cylindricalData
        rectToCylConverter
        cylinderProperties
        numberOfEvaluationPts = [50 55];
        orders = [9,9]
        qTh
        qL
    end
    properties (Access = private)
        thPtsCurrent
        lPtsCurrent
        surfaces
        dataList = {'w', 'u', 'v', 'wdot', 'vdot', 'udot'};
        
    end
    properties (Dependent)
        D
        K
    end
    methods
        function obj = Cylinder_Power_Flow(laserData, cylinderProps)
            obj.laserData = laserData;
            obj.cylinderProperties = cylinderProps;
        end
        function [qt, ql] = calculate_power_flow(obj)
            obj.convert_to_cylindrical_coordinates();
            obj.fit_surfaces();
            [qt, ql] = obj.compute_cylindrical_power_flow;
        end
    end
    methods
        function plot(obj)
            [TH, L] = ndgrid(obj.thPtsCurrent, obj.thPtsCurrent);
            quiver(TH, L,obj.qTh, obj.qL)
        end
        function  convert_to_cylindrical_coordinates(obj)
                obj.rectToCylConverter = CylCoorParameterizer(obj.laserData);
                obj.rectToCylConverter.calculate_cylindrical_data();
                obj.cylindricalData = obj.rectToCylConverter.export_cylindrical_data;
        end
        function val =  plot_component(obj, component, isReal, thPts, lPts, trim)
            if nargin < 6 
                trim = false;
            end
            if nargin < 4
                [thPts, lPts] = obj.default_plot_points();
            end
            if nargin < 3
                isReal = true;
            end
            switch component
                case {'w','wdot', 'u', 'udot', 'v', 'vdot'}
                    val = obj.plot_disp_or_vel(component, isReal, thPts, lPts, trim);
                case {'Qt', 'Ql', 'Mt', 'Ml', 'Mlt', 'Mtl', 'Nt', 'Nl', 'Ntl', 'Nlt'}
                    val = obj.plot_resultant(component, isReal, thPts, lPts, trim);
                otherwise
                    val = obj.plot_derivative(component, isReal, thPts, lPts, trim);
            end
        end
        function [qt, ql] = compute_cylindrical_power_flow(obj, thPts, lPts)
            [qts, qls] = obj.calculate_shear_power_flow(thPts, lPts);
            [qln, qtn] = obj.calculate_normal_power_flow(thPts, lPts);
            [qtb, qlb] = obj.calculate_bend_power_flow(thPts, lPts);
            [qtt, qlt] = obj.calculate_twist_power_flow(thPts, lPts);
            qt = qts + qtn + qtb + qtt;
            ql = qls + qln + qlb + qlt;
            obj.qTh = qt;
            obj.qL = ql;
        end
        %% set parameters
        function set_order_of_fit(obj, thOrder, lOrder)
            if nargin < 3
                lOrder = thOrder;
            end
            if length(thOrder) >1
                lOrder = [];
            end
            obj.orders = [thOrder, lOrder];
            obj.fit_surfaces();
        end
        function set_number_evaluation_points(obj, nTh, nL)
            obj.numberOfEvaluationPnts = [nTh, nL];
        end
        %% resultants
        function [Nt, Nl, Ntl, Nlt] = calculate_normal_resultants(obj, thPts, lPts)
            nu = obj.cylinderProperties.poissons;
            r = obj.cylinderProperties.radius;
            DD = obj.D;   KK = obj.K;
            [Vt, Vl] = obj.calculate_1st_order_th_disp_derivs(thPts, lPts);
            [Ut, Ul] = obj.calculate_1st_order_lng_disp_derivs(thPts, lPts);
            [Wtt, Wtl, Wll] = obj.calculate_2nd_order_rad_disp_derivs(thPts, lPts);
            W = obj.surfaces.evaluate;
            Nl = DD * (Ul + nu/r*(Vt + W)) - KK/r*Wll;
            Nt = DD * (1/r*Vt + nu*Ul + W/r) + KK/r^3 *(W + Wtt);
            Ntl = DD*(1-nu)/2*(Vl+1/r*Ut) + KK*(1-nu)/2 *(Wtl/r^2 + Ut/r^3);
            Nlt = DD*(1-nu)/2*(Vl+1/r*Ut) + KK*(1-nu)/2*(Vl - Wtl)/r^2;
        end
        function [Mt, Ml, Mtl, Mlt] = calculate_moment_resultants(obj, thPts, lPts)
            r = obj.cylinderProperties.radius;
            nu = obj.cylinderProperties.poissons;
            KK = obj.K;
            W = obj.surfaces.evaluate_complex_derivative_at_points(thPts, lPts, 'w');
            [Vt, Vl] = obj.calculate_1st_order_th_disp_derivs(thPts, lPts);
            [Ut, Ul] = obj.calculate_1st_order_lng_disp_derivs(thPts, lPts);
            [Wtt, Wtl, Wll] = obj.calculate_2nd_order_rad_disp_derivs(thPts, lPts);
            Mt= KK *(W/r^2 +Wtt/r^2 + nu*Wll);
            Ml = KK*(nu/r^2*Wtt + Wll -1/r * Ul - nu/r^2*Vt);
            Mlt = KK*(1-nu)/r*(Wtl-Vl);
            Mtl = KK*(1-nu)/2*(2*Wtl/r+ Vl/r + Ut/r^2);
        end
        function [Qt, Ql] = calculate_shear_resultants(obj, thPts, lPts)
            r = obj.cylinderProperties.radius;
            nu = obj.cylinderProperties.poissons;
            [Wttt, Wttl, Wtll, Wlll] = obj.calculate_3rd_order_rad_disp_derivs(thPts, lPts);
            dMt_dth_a = obj.K *(obj.drad_dth/r^3 +Wttt/r^3 + nu/r*Wtll);
            dMtl_dth_a =  obj.K*(1-nu)/2*(2*Wttl/r^2+ obj.d2tang_dthdz/r^2 + obj.d2long_dth2/r^3);
            dMl_dz =  obj.K*(nu/r^2*Wttl + Wlll -1/r * obj.d2long_dz2 - nu/r^2*obj.d2tang_dthdz);
            dMlt_dz = obj.K*(1-nu)/r*(Wtll-obj.d2tang_dz2);
           
            Ql = dMl_dz + dMtl_dth_a;
            Qt = dMt_dth_a + dMlt_dz;
        end
        %% resultant deriviatives
        function [Vd, Vp] = calculate_1st_order_th_disp_derivs(obj, thPts, lPts)
            Vd = obj.surfaces.evaluate_complex_displacement_derivative_at_points(thPts, lPts, 'v', [1,0]);
            Vp = obj.surfaces.evaluate_complex_displacement_derivative_at_points(thPts, lPts, 'v', [0,1]);
        end
        function [Ud, Up] = calculate_1st_order_lng_disp_derivs(obj, thPts, lPts)
            Ud = obj.surfaces.evaluate_complex_derivative_at_points(thPts, lPts, 'u', [1,0]);
            Up = obj.surfaces.evaluate_complex_derivative_at_points(thPts, lPts, 'u', [0,1]);
        end
        function [Wdd, Wdp, Wpp] = calculate_2nd_order_rad_disp_derivs(obj, thPts, lPts)
            Wdd = obj.surfaces.evaluate_complex_derivative_at_points(thPts, lPts, 'w', [2,0]);
            Wdp = obj.surfaces.evaluate_complex_derivative_at_points(thPts, lPts, 'w', [1,1]);
            Wpp = obj.surfaces.evaluate_complex_derivative_at_points(thPts, lPts, 'w', [0,2]);
        end
        function [Wddd, Wddp, Wdpp, Wppp] = calculate_3rd_order_rad_disp_derivs(obj, thPts, lPts)
            Wddd = obj.surfaces.evaluate_complex_derivative_at_points(thPts, lPts, 'w', [3,0]);
            Wddp = obj.surfaces.evaluate_complex_derivative_at_points(thPts, lPts, 'w', [2,1]);
            Wdpp = obj.surfaces.evaluate_complex_derivative_at_points(thPts, lPts, 'w', [1,2]);
            Wppp = obj.surfaces.evaluate_complex_derivative_at_points(thPts, lPts, 'w', [0,3]);
        end
        %% velocities
        function [vTh, vR, vL] = calculate_velocities(obj, thPts, lPts)
            % leaving derivative input blank just returns complex value of
            % surface at pts
            vTh = obj.surfaces.evaluate_complex_derivative_at_points(thPts, lPts, 'vdot');
            vR = obj.surfaces.evaluate_complex_derivative_at_points(thPts, lPts, 'wdot');
            vL = obj.surfaces.evaluate_complex_derivative_at_points(thPts, lPts, 'udot');
        end
        function [vThS, vRS, vLS] = calculate_conjugate_velocities(obj, thPts, lPts)
            [vTh, vR, vL] = obj.calculate_velocities(thPts, lPts);
            vThS = conj(vTh);
            vRS = conj(vR);
            vLS = conj(vL);            
        end
        function [OthDot, OlDot] = calculate_ang_velocities(obj, thPts, lPts)
            OthDot = obj.surfaces.evaluate_complex_derivative_at_points(thPts, lPts, 'wdot', [0,1]);
            OlDot = obj.surfaces.evaluate_complex_derivative_at_points(thPts, lPts, 'wdot', [1,0]);
        end
        function [OthDotS, OlDotS] = calculate_conjugate_ang_velocities(obj, thPts, lPts)
            [OthDot, OlDot] = obj.calculate_ang_velocities(thPts, lPts);
            OthDotS = conj(OthDot);
            OlDotS = conj(OlDot);
        end
        %% power flow parts
        function [qts, qls] = calculate_shear_power_flow(obj, thPts, lPts)
            [Qt, Ql] = obj.calculate_shear_resultants(thPts, lPts);
            [~, vRStar, ~] = obj.calculate_conjugate_velocities(thPts, lPts);
            qls = Ql.*vRStar;
            qts = Qt.*vRStar;
        end
        function [qtn, qln] = calculate_normal_power_flow(obj, thPts, lPts)
            [Nt, Nl, Ntl, Nlt] = obj.calculate_normal_resultants(thPts, lPts);
            [vThStar, ~, vLStar] = obj.calculate_conjugate_velocities(thPts, lPts);
            qln = Nlt.*vThStar + Nl.*vLStar;
            qtn = Nt.*vThStar + Ntl.*vLStar;
        end
        function [qtb, qlb] = calculate_bend_power_flow(obj, thPts, lPts)
            [Mt, Ml, ~, ~] = obj.calculate_moment_resultants(thPts, lPts);
            [omegaThS, omegaLS] = obj.calculate_conjugate_ang_velocities(thPts, lPts);
            qlb = Ml.*omegaThS;
            qtb = Mt.*omegaLS;
        end
        function [qtt, qlt] = calculate_twist_power_flow(obj, thPts, lPts)
            [~, ~, Mtl, Mlt] = obj.calculate_moment_resultants(thPts, lPts);
            
            [omegaThS, omegaLS] = obj.calculate_conjugate_ang_velocities(thPts, lPts);
            qlt = Mlt.*omegaLS;
            qtt =  Mtl.*omegaThS;
        end
        %% fitting surfaces
        function fit_surfaces(obj)
            [data_th, data_lg] = obj.get_x_and_y_parameters_for_fitting();
            fitter = ChebyshevSurfaceFitter(obj.orders, data_th, data_lg, []);
            obj.surfaces = CylinderPowerFlowSurfStruct();
            for i = 1:length(obj.dataList) % for all displacements and velocities
                rp = [obj.dataList{i} 'R']; ip = [obj.dataList{i} 'I'];
                vals = obj.get_measured_values_to_fit(rp); % get data for specific part
                surf = fitter.fit_surface_to_new_values(vals); % fit surface
                obj.surfaces.set_surface(rp, surf);% store surface
                vals = obj.get_measured_values_to_fit(ip);
                surf = fitter.fit_surface_to_new_values(vals);
                obj.surfaces.set_surface(ip, surf)
            end
        end
        function [xs, ys] = get_x_and_y_parameters_for_fitting(obj)
            xs = obj.laserData.coordinates(:,1);
            ys = obj.laserData.coordinates(:,2);
        end
        function vals = get_measured_values_to_fit(obj, type)
            if contains(type, 'dot') 
                if contains(type, 'R')
                    var = 'realVel';
                else
                    var = 'imagVel';
                end
            else
                if contains(type, 'R')
                    var = 'realDisp';
                else
                    var = 'imagDisp';
                end
            end
            if contains(type, 'w')
                indx = 2;
            elseif contains(type, 'v')
                indx = 1;
            else
                indx = 3;
            end            
            vals = obj.laserData.(var)(:,indx);
        end
        %% getters
        function val = get.D(obj)
            val = obj.cylinderProperties.E*obj.cylinderProperties.thickness/(1-obj.cylinderProperties.poissons^2);
        end
        function val = get.K(obj)
            val = obj.cylinderProperties.E*obj.cylinderProperties.thickness^3/(12*(1-obj.cylinderProperties.poissons^2));
        end
    end
end