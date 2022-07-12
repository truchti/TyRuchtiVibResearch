classdef Distributed_Force_Plate_Simulator < handle
    properties
        material
        forces
        geometry
        mesh
        saveFolder
    end
    properties(Dependent)
        fileName
        naturalFrequencies
    end
    properties (Hidden = true, Dependent = true)
        heights
        widths
        D
        lambda
        time
    end
    properties (Hidden = true)
        timesteps = 20
        nOmega= [];
        complexDampingPhase = [];
        WtotalD, WtotalF
        data
        dWdx, dWdy, dWdt
        dWdx2, dWdxy, dWdy2
        dWdx3, dWdx2y, dWdxy2, dWdy3
        dWdyt, dWdxt
    end
    properties (Access = private)
        timeV, Phi
        FW, FWt
        PNwidth, PNheight
        x1t, x2t, y1t, y2t
        C, C_m, C_mn
        cosXdiff, cosYdiff
        fx, dfx, dfx2, dfx3
        fy, dfy, dfy2, dfy3
        fxy, dfxydx, dfxydy
        dfxydx2, dfxydxy, dfxydy2, 
        dfxydx3, dfxydx2y, dfxydxy2, dfxydy3
        ft, dft
    end
    methods
        %% constructor and switch functions
        function obj = Distributed_Force_Plate_Simulator(geometry, forces, material, mesh)
            if nargin < 4
                warning('Not enough info');
            else
                obj.geometry = geometry;
                obj.forces = forces;
                obj.material = material;
                obj.mesh = mesh;
            end
        end
        function disable_force(obj, force_number)
            switch force_number
                case 1
                    obj.force1.enabled = false;
                case 2
                    obj.force2.enabled = false;
            end
        end
        function enable_force(obj, force_number)
            switch force_number
                case 1
                    obj.force1.enabled = true;
                case 2
                    obj.force2.enabled = true;
            end
        end
        %% data generation functions
        function simulate_and_save(obj)
            obj.simulate_data();
            obj.save_data();
        end
        function simulate_data(obj)
            obj.calculate_mode_frequencies_and_phases();
            obj.calculate_displacement_over_time();
            obj.calculate_spectral_terms();
            obj.collectData();
        end
        function save_data(obj)
            if ~isempty(obj.data)
                obj.write_velocity_file();
                obj.write_displacement_file();
            else
                warning("data has not been simulated")
            end
        end
        function props = export_plate_power_flow_properties(obj)
            props = Plate_Properties(obj.geometry.width, obj.geometry.height, obj.geometry.thickness, obj.material.E, obj.material.poisson);
        end
        %% plotting functions
        function plot_RI_component(obj)
        end
%         function [qx, qy] = plot_power_flow(obj)
%             [Wid, Hei] = ndgrid(obj.widths, obj.heights);
%             
%             quiver(Wid, Hei, qx, qy);
%             axis equal
%         end
        
        function [qxp, qyp] = plot_power_flow(obj, type)
            if nargin < 2
                type = 'Power Flow';
            end
            switch type
                case {'bending', 'b', 'bend'}
                    [qxp, qyp] = obj.calculate_bending_power_flow();
                case {'twisting', 't', 'twist'}
                    [qxp, qyp] = obj.calculate_twisting_power_flow();
                case {'shear', 's'}
                    [qxp, qyp] = obj.calculate_shear_power_flow();
                otherwise
                    [qxp, qyp] = obj.calculate_power_flow();
            end
            [Wid, Hei] = ndgrid(obj.widths, obj.heights);
            quiver(Wid, Hei, qxp, qyp, 'color', [1 0 .2]);
            title(type)
        end
        function [X,Y,Z] = show_displacement(obj, isReal)
            if nargin< 2
                isReal = true;
            end
            value = obj.WtotalF;
            if isReal
                value = real(value);
                st = 'real ';
            else
                value = imag(value);
                st = 'imag ';
            end
            [H, W] = ndgrid(obj.heights, obj.widths);
            surf(H, W, value', value')
%             view(0,90)
            X = H; Y = W; Z = value;
            colorbar
            colormap jet
%             axis equal
            xlabel('Y')
            ylabel('X')
            title(strcat(st, 'displacement' ))
        end
        function [X,Y,Z] = show_velocity(obj, type, isReal)
            switch type
                case 'wdot'
                    value = obj.dWdt;
                case {'Oxdot', 'dWdyt'}
                    value = obj.dWdyt;
                case {'Oydot', 'dWdxt'}
                    value = obj.dWdxt;
            end
            if isReal
                value = real(value);
                st = 'real ';
            else
                value = imag(value);
                st = 'imag ';
            end
            [H, W] = ndgrid(obj.heights, obj.widths);
            surf(H, W, value', value')
%             view(0,90)
            X = H; Y = W; Z = value;
            colorbar
            colormap jet
%             axis equal
            xlabel('Y')
            ylabel('X')
            title(strcat(st, type))
        end
        function [X,Y,Z] = show_resultant(obj, type, isReal)
            if nargin < 3
                isReal = true;
            end
            switch type
                case 'Mx'
                    [value, ~, ~] = obj.calculate_moment();
                case 'My'
                    [~, value, ~] = obj.calculate_moment();
                case 'Mxy'
                    [~, ~, value] = obj.calculate_moment();
                case 'Qx'
                    [value,~] = obj.calculate_shear();
                case 'Qy'
                    [~, value] = obj.calculate_shear();
            end
            if isReal
                value = real(value);
                st = "real ";
            else
                value = imag(value);
                st = "imag ";
            end
            [H, W] = ndgrid(obj.heights, obj.widths);
            surf(H, W, value', value')
            X = H; Y = W; Z = value;
            xlabel('Y')
            ylabel('X')
%             view(0,90)
            colorbar
            colormap jet
            title(strcat(st, type))
        end
        function [X,Y,Z] = show_derivative(obj, type, isReal)
            if nargin < 3
                isReal = true;
            end
            switch type
                case 'dwdx'
                    value = obj.dWdx;
                case 'dwdy'
                    value = obj.dWdy;
                case 'dwdxy'
                    value = obj.dWdxy;
                case 'dwdx2'
                    value = obj.dWdx2;
                case 'dwdy2'
                    value = obj.dWdy2;
                case 'dwdxy2'
                    value = obj.dWdxy2;
                case 'dwdx2y'
                    value = obj.dWdx2y;
                case 'dwdx3'
                    value = obj.dWdx3;
                case 'dwdy3'
                    value = obj.dWdy3;
            end
            if isReal
                value = real(value);
                st = "real ";
            else
                value = imag(value);
                st = "imag ";
            end
            [H, W] = ndgrid(obj.heights, obj.widths);
            surf(H, W, value', value')
            X = H; Y = W; Z = value;
            xlabel('Y')
            ylabel('X')
%             view(0,90)
            colorbar
            colormap jet
            title(strcat(st, type))
        end
        function [X,Y,Z] = show_power_component(obj,type)
            [qx, qy] = obj.calculate_shear_power_flow();
            switch type
                case {'qx'}
                    value = qx;
                case {'qy'}
                    value = qy;
            end
            st = '';
            [H, W] = ndgrid(obj.heights, obj.widths);
            surf(H, W, value', value')
            X = H; Y = W; Z = value;
            xlabel('Y')
            ylabel('X')
%             view(0,90)
            colorbar
            colormap jet
            title(strcat(st, type))            
        end
        function [X,Y,Z] = show_3D(obj, type, isReal)
            if nargin < 3
                isReal = true;
            end
            switch type
                case {'qx', 'qy'}
                    [X,Y,Z] = obj.show_power_component(type);
                case 'w'
                    [X,Y,Z] = obj.show_displacement(isReal);
                case {'Mx', 'My', 'Mxy', 'Qx', 'Qy'}
                    [X,Y,Z] = obj.show_resultant(type, isReal);
                case{'wdot', 'Oxdot', 'dWdyt', 'Oydot', 'dWdxt'}
                    [X,Y,Z] = obj.show_velocity(type, isReal);
                otherwise
                    [X,Y,Z] = obj.show_derivative(type, isReal);
            end
        end
        function animate_displacement_in_time(obj, imagData, loops)
            if nargin < 3
                loops = 3;
            end
            if nargin < 2
                imagData = false;
            end
            [Hei, Wid] = meshgrid(obj.heights, obj.widths);
            if imagData
                values = imag(obj.WtotalD);
            else
                values = real(obj.WtotalD);
            end
            figure(1)
            s = surf(Wid,Hei,values(:,:,1),values(:,:,1));
            s.EdgeAlpha = 0.0;
            s.FaceColor = 'interp';
            colorbar
            clim = max(max(max(values)));
            zlim([-clim, clim]);
            caxis([-clim, clim])
            sz = size(values,3);
            for i = 1:loops*sz
                s.ZData = values(:,:,mod(i,sz)+1);
                s.CData = values(:,:,mod(i,sz)+1);
                pause(0.1);
            end
        end
    end
    methods (Access = private) %helper functions
        %% calculation steps
        function calculate_mode_frequencies_and_phases(obj)
            scale = pi^2*sqrt(obj.D/(obj.material.density*obj.geometry.thickness));
            obj.nOmega = zeros(obj.mesh.width_modes,obj.mesh.height_modes);
            obj.complexDampingPhase = zeros(size(obj.nOmega));
            for m = 1:obj.mesh.width_modes
                for n = 1:obj.mesh.height_modes
                    obj.nOmega(m,n) = ((m/obj.geometry.width)^2+(n/obj.geometry.height)^2)*scale;
                    obj.complexDampingPhase(m,n) = atan(2*obj.lambda*(obj.forces(1).omega/obj.nOmega(m,n))/(1-(obj.forces(1).omega/obj.nOmega(m,n))^2));
                end
            end
        end
        function calculate_displacement_over_time(obj)
            %precalculate values that will be used a lot
            obj.preCalculate_oft_used_terms();
            obj.preallocate_zero_arrays(); % make sure displacment arrays start at 0
            for i = 1:length(obj.forces)
                obj.preCalculate_force_terms(i);
                for M = 1:obj.mesh.width_modes 
                    obj.calculate_force_effect_m_mode(M);
                    obj.calculate_horizontal_disp_and_deriv_vectors(M);
                    for N = 1:obj.mesh.height_modes
                        obj.calculate_force_effect_n_mode(N);
                        obj.calculate_vertical_disp_and_deriv_vectors(N);
                        obj.calc_mode_disp_as_func_of_x_y_t_and_derivs(M,N);
                    end
                end
            end
        end
        function calculate_spectral_terms(obj)
            timeDim = 3;% time is the third dimension of the data
            props = {'dWdt', 'dWdx', 'dWdy', 'dWdx2', 'dWdxy', 'dWdy2', 'dWdx3', 'dWdx2y', 'dWdxy2', 'dWdy3', 'dWdyt', 'dWdxt'};
            fullFFT = 1/(size(obj.WtotalD,3))*fft(obj.WtotalD, [], timeDim); % take the FFT over time dimension
            [~, indx] = max(squeeze(fullFFT(2,2,:)));
            obj.WtotalF = fullFFT(:,:,indx);
            for i = 1:length(props)
                FFT = 1/(size(obj.(props{i}),3))*fft(obj.(props{i}), [], timeDim);
                obj.(props{i}) = FFT(:,:,indx);
            end
        end
        function [qx, qy] = calculate_power_flow(obj)
            [sx, sy] = obj.calculate_shear_power_flow();
            [bx, by] = obj.calculate_bending_power_flow();
            [tx, ty] = obj.calculate_twisting_power_flow();
            qx = sx + bx + tx;
            qy = sy + by + ty;
        end
        %% power flow helpers
        function [qxs, qys] = calculate_shear_power_flow(obj)
            [Qx, Qy] = obj.calculate_shear();
            qxs = 0.5 * real(Qx .* conj(obj.dWdt));
            qys = 0.5 * real(Qy .* conj(obj.dWdt));
        end
        function [qxb, qyb] = calculate_bending_power_flow(obj)
            [Mx, My, ~] = obj.calculate_moment();
            qxb = -0.5 * real(Mx .* conj(obj.dWdxt));
            qyb = -0.5 * real(My .* conj(obj.dWdyt));
        end
        function [qxt, qyt] = calculate_twisting_power_flow(obj)
            [~, ~, Mxy] = obj.calculate_moment();
            qxt = -0.5 * real(Mxy .* conj(obj.dWdyt));
            qyt = -0.5 * real(Mxy .* conj(obj.dWdxt));
        end        
        function [Mx, My, Mxy] = calculate_moment(obj)
            nu = obj.material.poisson;
            Mx = obj.D * (obj.dWdx2 + nu * obj.dWdy2);
            My = obj.D * (obj.dWdy2 + nu * obj.dWdx2);
            Mxy = obj.D * (1-nu) * obj.dWdxy;
        end
        function [Qx, Qy] = calculate_shear(obj)
            Qx = obj.D * (obj.dWdx3 + obj.dWdxy2);
            Qy = obj.D * (obj.dWdx2y + obj.dWdy3);
        end
        %% calculate f(x,y,t) helpers
        function calc_mode_disp_as_func_of_x_y_t_and_derivs(obj, M,N) 
            TS = obj.calculate_total_scalar(M,N);
            obj.calculate_spatial_function_and_derivatives(TS);
            obj.calcualte_time_function_and_derivative(M,N);
            obj.add_modal_contribution_of_temporal_spatial_func_deriv();
        end 
        function calculate_spatial_function_and_derivatives(obj, TS)
            % because fx terms are columns and fy terms are rows vector multiplying 
            %gives a matrix that is length(fx) by length(fy)
            obj.fxy = TS * obj.fx *obj.fy;
            obj.dfxydx = TS * obj.dfx * obj.fy;
            obj.dfxydy = TS * obj.fx * obj.dfy;
            obj.dfxydx2 = TS * obj.dfx2 * obj.fy;
            obj.dfxydxy = TS * obj.dfx * obj.dfy;
            obj.dfxydy2 = TS * obj.fx * obj.dfy2;
            obj.dfxydx3 = TS * obj.dfx3 * obj.fy;
            obj.dfxydx2y = TS * obj.dfx2 * obj.dfy;
            obj.dfxydxy2 = TS * obj.dfx * obj.dfy2;
            obj.dfxydy3 =  TS * obj.fx * obj.dfy3;                         
        end
        function calcualte_time_function_and_derivative(obj, M, N)
            obj.ft(1,1,:) = cos(obj.FWt-obj.Phi(M,N))+1i*sin(obj.FWt-obj.Phi(M,N));% create time vector with phase shift as a 1 by 1 by t vector
            omg = obj.forces(1).omega;
            obj.dft(1,1,:) = 1i*omg*obj.ft(1,1,:);
        end
        function add_modal_contribution_of_temporal_spatial_func_deriv(obj)
            %because the first two dims of times are 1 the .* causes the 
            %entire matrix amp to be multiplied and by each value of times 
            %and gives an M by N by t array
            P1 = {'WtotalD', 'dWdx', 'dWdy', 'dWdx2', 'dWdxy', 'dWdy2', 'dWdx3', 'dWdx2y', 'dWdxy2', 'dWdy3'};
            P2 = {'fxy', 'dfxydx', 'dfxydy', 'dfxydx2', 'dfxydxy', 'dfxydy2', 'dfxydx3', 'dfxydx2y', 'dfxydxy2', 'dfxydy3'};
            for i = 1:length(P1)
                obj.(P1{i}) = obj.(P1{i})+obj.(P2{i}) .* obj.ft;
            end
            obj.dWdt = obj.dWdt  + obj.fxy .* obj.dft;
            obj.dWdyt = obj.dWdyt + obj.dfxydy .* obj.dft;
            obj.dWdxt = obj.dWdxt + obj.dfxydx .* obj.dft;
        end
        function calculate_force_effect_m_mode(obj, num)
            obj.cosXdiff = cos(num*obj.x1t)-cos(num*obj.x2t);
            obj.C_m = obj.C/num;
        end
        function calculate_force_effect_n_mode(obj, num)
            obj.cosYdiff = cos(num*obj.y1t)-cos(num*obj.y2t);
            obj.C_mn = obj.C_m/num;
        end
        function calculate_horizontal_disp_and_deriv_vectors(obj, num)
            din = num*pi/obj.geometry.width;
            obj.fx = sin(num*obj.PNwidth)';
            obj.dfx = din*(cos(num*obj.PNwidth))';
            obj.dfx2 = din^2 * -(sin(num*obj.PNwidth))';
            obj.dfx3 = din^3 * -(cos(num*obj.PNwidth))';
        end
        function calculate_vertical_disp_and_deriv_vectors(obj, num)
            din = num*pi/obj.geometry.height;
            obj.fy = sin(num*obj.PNheight);
            obj.dfy = din*(cos(num*obj.PNheight));
            obj.dfy2 = din^2 * -(sin(num*obj.PNheight));
            obj.dfy3 = din^3 * -(cos(num*obj.PNheight));
        end
        function preallocate_zero_arrays(obj)
            obj.WtotalD = zeros(length(obj.widths), length(obj.heights), length(obj.time));
            obj.dWdt = zeros(length(obj.widths), length(obj.heights), length(obj.time));
            obj.dWdx = zeros(length(obj.widths), length(obj.heights), length(obj.time));
            obj.dWdy = zeros(length(obj.widths), length(obj.heights), length(obj.time));
            obj.dWdx2 = zeros(length(obj.widths), length(obj.heights), length(obj.time));
            obj.dWdxy = zeros(length(obj.widths), length(obj.heights), length(obj.time));
            obj.dWdy2 = zeros(length(obj.widths), length(obj.heights), length(obj.time));
            obj.dWdx3 = zeros(length(obj.widths), length(obj.heights), length(obj.time));
            obj.dWdx2y = zeros(length(obj.widths), length(obj.heights), length(obj.time));
            obj.dWdxy2 = zeros(length(obj.widths), length(obj.heights), length(obj.time));
            obj.dWdy3 = zeros(length(obj.widths), length(obj.heights), length(obj.time));
            obj.dWdyt = zeros(length(obj.widths), length(obj.heights), length(obj.time));
            obj.dWdxt = zeros(length(obj.widths), length(obj.heights), length(obj.time));
        end
        function preCalculate_oft_used_terms(obj) %% check for height and width mesh
            obj.timeV = obj.time;
            obj.PNwidth = pi*obj.widths/obj.geometry.width;
            obj.PNheight = pi*obj.heights/obj.geometry.height;
        end
        function preCalculate_force_terms(obj, fNum)
            obj.FW = obj.forces(fNum).omega;
            obj.FWt = obj.timeV*obj.forces(fNum).omega;
            obj.x1t = obj.forces(fNum).xbounds(1)*pi/obj.geometry.width;
            obj.x2t = obj.forces(fNum).xbounds(2)*pi/obj.geometry.width;
            obj.y1t = obj.forces(fNum).ybounds(1)*pi/obj.geometry.height;
            obj.y2t = obj.forces(fNum).ybounds(2)*pi/obj.geometry.height;
            obj.C = 4*obj.forces(fNum).magnitude/(pi^2*obj.material.density*obj.geometry.thickness);
            obj.Phi =  obj.complexDampingPhase + deg2rad(obj.forces(fNum).phase);
        end
        function denom = calculate_denominator(obj,M,N)
            denom =(obj.nOmega(M,N)^2-obj.forces(1).omega^2)^2+obj.material.eta^2*obj.nOmega(M,N)^4;
        end
        function imagCoef = calculate_complex_damping_coefficient(obj, M,N)
            imagCoef = (1-1i*obj.material.eta)*obj.nOmega(M,N)^2-obj.forces(1).omega^2;
        end
        function totalScalar = calculate_total_scalar(obj, M, N)
            denom = obj.calculate_denominator(M,N);
            cmplxDamp= obj.calculate_complex_damping_coefficient(M,N);
            totalScalar = obj.C_mn* obj.cosXdiff * obj.cosYdiff * cmplxDamp / denom;
        end
        %% writing files helpers
        function collectData(obj)
            % put into column vectors
            [Wid, Hei] = ndgrid(obj.widths, obj.heights); % height and width mesh
            tempData.xCoor =reshape(Wid, numel(Wid),1);
            tempData.yCoor =reshape(Hei, numel(Hei), 1);
            tempData.zCoor = ones(size(tempData.yCoor));
            % no displacement in X or Y
            tempData.rX = zeros(size(tempData.xCoor));
            tempData.iX = zeros(size(tempData.xCoor));
            tempData.rY = zeros(size(tempData.yCoor));
            tempData.iY = zeros(size(tempData.yCoor));
            tempData.rXv = zeros(size(tempData.xCoor));
            tempData.iXv = zeros(size(tempData.xCoor));
            tempData.rYv = zeros(size(tempData.yCoor));
            tempData.iYv = zeros(size(tempData.yCoor));
            %
            fullDisp = reshape(obj.WtotalF, numel(obj.WtotalF), 1);
            fullVel = 1i*obj.forces(1).omega*fullDisp;
            tempData.rZ = real(fullDisp);
            tempData.iZ = imag(fullDisp);
            tempData.rZv = real(fullVel);
            tempData.iZv = imag(fullVel);
            obj.data = tempData;
        end
        function write_velocity_file(obj)
            obj.write_file('vel')
        end
        function write_displacement_file(obj)
            obj.write_file('disp')
        end
        function write_file(obj, type)
            switch type
                case 'vel'
                    tstr = 'Velocity.txt';
                case 'disp'
                    tstr = 'Displacement.txt';
            end
            old = cd(obj.saveFolder);
            DataFileName = strcat(obj.fileName, "_Simulated_", tstr);
            fid = fopen(DataFileName,'wt');
            obj.write_header(fid, type, obj.forces(1).frequency)
            obj.write_data(fid, obj.data, type);
            fclose(fid);
            cd(old);
        end
        
    end
    methods % getters
        function value = get.fileName(obj)
            date = datestr(floor(now));
            date = strrep(date,'-','_');
            value = strcat('plate_simulation_2_forces_', date, '_', num2str(obj.forces(1).frequency), '_Hz');
        end
        function value = get.D(obj)
            value = obj.material.E*obj.geometry.thickness^3/(12*(1-obj.material.poisson^2));
        end
        function value = get.heights(obj)
            value = linspace(0,obj.geometry.height,obj.mesh.height_divisions+1);
        end
        function value = get.widths(obj)
            value = linspace(0, obj.geometry.width, obj.mesh.width_divisions+1);
        end
        function value = get.time(obj)
            t =  linspace(0, 1/obj.forces(1).frequency, obj.timesteps+1);
            value = t(1:end-1);
        end
        function value = get.lambda(obj)
            value = obj.material.density*obj.geometry.thickness*obj.material.eta/obj.forces(1).omega;
        end
        function value = get.naturalFrequencies(obj)
            value = obj.nOmega/(2*pi);
        end
    end
    methods (Static) % for writing laser files
        function write_header(fid, type, frequency)
            switch type
                case 'vel'
                    Type = 'Velocity';
                    unit = 'm/s';
                case 'disp'
                    Type = 'Displacement';
                    unit = 'm';
            end
            fprintf(fid, strcat("Source File Name:\t", "Simulated Data\n"));
            fprintf(fid, strcat("Signal: FFT - Vib 3D ", Type, " - Real & Imag.\n\n"));
            fprintf(fid, strcat("Band No.:\t 1\nFrequency:\t", num2str(frequency), " Hz\n\n\n"));
            fprintf(fid, "Interpolated:\tYes\nFiltered:\tYes\n\n");
            fprintf(fid, strcat("Index\tX\tY\tZ\tReal X [", unit,"]\tReal Y [", unit,"]\tReal Z [", unit,"]\tImaginary X [", unit,"]\t Imaginary Y [", unit,"]\tImaginary Z [", unit,"]\n"));
        end
        function write_data(fid, data, type)
            dataFormat = strcat('%d',repmat('\t%g', 1, 9), '\n');
            X = data.xCoor;  Y = data.yCoor;   Z = data.zCoor;
            switch type
                case 'vel'
                    rX = data.rXv;  rY = data.rYv;   rZ = data.rZv;
                    iX = data.iXv;  iY = data.iYv;   iZ = data.iZv;
                case 'disp'
                    rX = data.rX;  rY = data.rY;   rZ = data.rZ;
                    iX = data.iX;  iY = data.iY;   iZ = data.iZ;
            end
            for point = 1:length(X)
                fprintf(fid, dataFormat, point, X(point), Y(point), Z(point), rX(point), rY(point), rZ(point), iX(point), iY(point), iZ(point));
            end
        end
    end
end
