classdef SLDV_FFT_Data < handle
    properties
        coordinates
        displacements
        velocities
        frequencies
    end
    properties (Hidden)
        points
    end
    methods
        function obj = SLDV_FFT_Data(frequencies, coords, displacements, velocities)
            if isa(coords,'GeoPoint')
                obj.points = coords;
                obj.coordinates = obj.convert_geo_points_to_coordinates(coords);
            else
                obj.coordinates = coords;
            end
            obj.frequencies = frequencies;
            obj.displacements = displacements;
            obj.velocities = velocities;
        end
        function value = get_type_rOrI_at_point_dim_freq(obj, velOrDisp, realOrImag, point, dim, freq)
            vdri = obj.process_type_and_reality_bools(velOrDisp, realOrImag);
            switch vdri
                case 'VR'
                    allvalues = real(obj.velocities);
                case 'DR'
                    allvalues = real(obj.displacements);
                case 'VI'
                    allvalues = imag(obj.velocities);
                case 'DI'
                    allvalues = imag(obj.displacements);
                otherwise
                    allvalues = [];
            end
            value = obj.get_value_at_point_dim_freq(allvalues,point,dim, freq);
        end
        function value = get_value_at_point_dim_freq(obj, allVals, point,dim,freq)
            freq_index = obj.frequencies == freq;
            value = (allVals(point,dim,freq_index));
        end
        function value = x_d_r(obj, point, frequency)
            value = obj.get_type_rOrI_at_point_dim_freq('D', 'R', point, 1, frequency);
        end
        function value = y_d_r(obj, point, frequency)
            value = obj.get_type_rOrI_at_point_dim_freq('D', 'R', point, 2, frequency);
        end
        function value = z_d_r(obj, point, frequency)
            value = obj.get_type_rOrI_at_point_dim_freq('D', 'R', point, 3, frequency);
        end
        function value = x_d_i(obj, point, frequency)
            value = obj.get_type_rOrI_at_point_dim_freq('D', 'I', point, 1, frequency);
        end
        function value = y_d_i(obj, point, frequency)
            value = obj.get_type_rOrI_at_point_dim_freq('D', 'I', point, 2, frequency);
        end
        function value = z_d_i(obj, point, frequency)
            value = obj.get_type_rOrI_at_point_dim_freq('D', 'I', point, 3, frequency);
        end
        function value = x_v_r(obj, point, frequency)
            value = obj.get_type_rOrI_at_point_dim_freq('V', 'R', point, 1, frequency);
        end
        function value = y_v_r(obj, point, frequency)
            value = obj.get_type_rOrI_at_point_dim_freq('V', 'R', point, 2, frequency);
        end
        function value = z_v_r(obj, point, frequency)
            value = obj.get_type_rOrI_at_point_dim_freq('V', 'R', point, 3, frequency);
        end
        function value = x_v_i(obj, point, frequency)
            value = obj.get_type_rOrI_at_point_dim_freq('V', 'I', point, 1, frequency);
        end
        function value = y_v_i(obj, point, frequency)
            value = obj.get_type_rOrI_at_point_dim_freq('V', 'I', point, 2, frequency);
        end
        function value = z_v_i(obj, point, frequency)
            value = obj.get_type_rOrI_at_point_dim_freq('V', 'I', point, 3, frequency);
        end
        function plot_point_cloud(obj)
            scatter3(obj.coordinates(:,1), obj.coordinates(:,2), obj.coordinates(:,3))
        end
    end
    methods (Hidden = true, Static = true)
        function chrs = process_type_and_reality_bools(vOrD, rOrI)
            if vOrD
                vd = 'V';
            else
                vd = 'D';
            end
            if rOrI
                ri = 'R';
            else
                ri = 'I';
            end
            chrs = strcat(vd,ri);
        end
        function coors = convert_geo_points_to_coordinates(geoPoints)
            coors = zeros(length(geoPoints), 3);
            for pt = 1:length(geoPoints)
                coors(pt,:) = geoPoints(pt).get_coordinate;
            end                
        end 
    end
end