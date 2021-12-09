classdef GeoPoint
    %UNTITLED21 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        number
        xCoordinate
        yCoordinate
        zCoordinate
    end
    methods
        function obj = GeoPoint(pointDataString)
            if nargin ==1
                pointDataString(strfind(pointDataString, '=')) = [];
                obj.number = sscanf(pointDataString, '%f');
                obj.xCoordinate = parse_for_x(pointDataString);
                obj.yCoordinate = parse_for_y(pointDataString);
                obj.zCoordinate = parse_for_z(pointDataString);
            end
        end
        function value = get_coordinate(obj)
            value = [obj.xCoordinate, obj.yCoordinate, obj.zCoordinate];
        end
    end
end
function number = parse_for_number_after_key(pointDataString, Key)
    Index = strfind(pointDataString, Key);
    number = sscanf(pointDataString(Index(1) + length(Key):end), '%g', 1);
end
function value = parse_for_x(pointDataString)
    value = parse_for_number_after_key(pointDataString, 'x');
end
function value = parse_for_y(pointDataString)
    value = parse_for_number_after_key(pointDataString, 'y');
end
function value = parse_for_z(pointDataString)
    value = parse_for_number_after_key(pointDataString, 'z');
end
