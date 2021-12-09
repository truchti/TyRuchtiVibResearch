classdef SLDV_data 
    properties
        frequency
        coordinates
        realDisp
        imagDisp
        realVel
        imagVel
    end
    properties(Hidden = true)
        validPoints
    end
    methods
        function obj = SLDV_data(coors, reDisp, imDisp, reVel, imVel, freq)
            if nargin == 1
                if isa(coors, 'SLDV3D_FastScan_parser')
                    obj = obj.process_fastscan_parser(coors);
                elseif isa(coors, 'SLDV3D_FFT_parser')
%                     obj.process_full_fft_parser(coors);
                end
            else
                obj.frequency = freq;
                obj.coordinates = coors;
                obj.realDisp = reDisp;
                obj.imagDisp = imDisp;
                obj.realVel = reVel;
                obj.imagVel = imVel;
            end
            obj.validPoints = true(length(coors),1);
        end
        function plot_coordinates(obj)
            [X, Y, Z] = obj.separate_coordinates;
            scatter3(X,Y,Z);
        end
        function plot_displacement(obj, imaginary)
            if nargin < 2 
                imaginary = false;
            end
            [X,Y,Z]= obj.separate_coordinates;
            [U,V,W] = obj.separate_displacements(imaginary);
            quiver3(X,Y,Z,U,V,W)
            xlabel('X')
            ylabel('Y')
            zlabel('Z')
        end
        function plot_velocity(obj, imaginary)
            if nargin < 2 
                imaginary = false;
            end
            [X,Y,Z]= obj.separate_coordinates;
            [U,V,W] = obj.separate_velocities(imaginary);
            quiver3(X,Y,Z,U,V,W)
        end
        function identify_outliers_visually(obj)
            obj.plot_coordinates;
            view(0,0)
            title('Select All Points that are invalid')
            brsh = brush;
            set(brsh,'Enable','on','ActionPostCallback',@brushedDataCallback);
        end
        function obj = remove_points_by_index(obj, indexes)
            if isempty(obj.validPoints)
                invalid = zeros(length(obj.coordinates));
            else
                invalid = ~obj.validPoints;
            end
            for indx = 1:length(indexes)
                invalid(indexes(indx)) = 1;
            end
            obj.validPoints = ~logical(invalid);
        end
        function obj = disable_points_based_on_error_function(obj, errorFunc, tol)
            coors = obj.coordinates;
            for i = 1:length(coors)
                valid = errorFunc(coors(i,:)) < tol;
                obj.validPoints(i) = valid;
            end
            obj.validPoints = logical(obj.validPoints);
        end
        function plot_valid_points(obj)
            [X,Y,Z] = obj.separate_coordinates;
            scatter3(X(obj.validPoints), Y(obj.validPoints), Z(obj.validPoints));
        end
        function validObj = create_data_structure_with_only_valid_points(obj)
            coors = obj.coordinates(obj.validPoints,:);
            reDisp = obj.realDisp(obj.validPoints,:);
            imDisp = obj.imagDisp(obj.validPoints,:);
            reVel = obj.realVel(obj.validPoints,:);
            imVel = obj.imagVel(obj.validPoints,:);
            freq = obj.frequency;
            validObj = SLDV_data(coors, reDisp, imDisp, reVel, imVel, freq);
        end
        function groupData = get_array_of_complex_displacement_and_velocity_values(obj)
            groupData = [obj.realDisp, obj.imagDisp, obj.realVel, obj.imagVel];
        end
    end
    methods (Hidden = true)
        function obj = process_fastscan_parser(obj, parser)
            obj.frequency = parser.frequency;
            obj.coordinates = parser.coordinates;
            obj.realDisp = parser.realDisp;
            obj.imagDisp = parser.imagDisp;
            obj.realVel = parser.realVel;
            obj.imagVel = parser.imagVel;
        end
        function [rotObj] = translate_and_rotate_data(obj, center, axis)%, forceLocations)
            points = obj.coordinates-center;
            M = obj.get_rotation_matrix(axis);
            newCoordinates = points*M;
            newRealDisp = obj.realDisp*M;
            newImagDisp = obj.imagDisp*M;
            newRealVel = obj.realVel*M;
            newImagVel = obj.imagVel*M;
            rotObj = SLDV_data(newCoordinates, newRealDisp, newImagDisp, newRealVel, newImagVel, obj.frequency);
        end
        function [X, Y, Z] = separate_coordinates(obj)
            X = obj.coordinates(:,1);
            Y = obj.coordinates(:,2);
            Z = obj.coordinates(:,3);
        end
        function [U, V, W] = separate_displacements(obj, imaginary)
            if nargin<2 || ~imaginary
                U = obj.realDisp(:,1);
                V = obj.realDisp(:,2);
                W = obj.realDisp(:,3);
            else
                U = obj.imagDisp(:,1);
                V = obj.imagDisp(:,2);
                W = obj.imagDisp(:,3);
            end                
        end
        function [U, V, W] = separate_velocities(obj, imaginary)
            if nargin<2 || ~imaginary
                U = obj.realVel(:,1);
                V = obj.realVel(:,2);
                W = obj.realVel(:,3);
            else
                U = obj.imagVel(:,1);
                V = obj.imagVel(:,2);
                W = obj.imagVel(:,3);
            end       
        end
    end
    methods (Hidden = true, Static = true)
        function M = get_rotation_matrix( axis)
            nx = axis(1); ny = axis(2); nz = axis(3);
            d = 1+ny;
            M = [1-nx^2/d, nx,  -(nx*nz)/d;
                -nx, ny, -nz;
                -(nx*nz)/d, nz, 1-nz^2/d];
        end
    end
end