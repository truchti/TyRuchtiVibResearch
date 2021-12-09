classdef FFT_file_data < handle
    
    properties
        fileAddress
        nameOfDataSet
        point
        scanType
        outputType
        componentType
        measurementType
        frequencies
        realData
        imaginaryData
        dataHeaders
        units
        rawData
        
    end
    properties (Hidden)
        fileID
        
    end
    methods
        function obj = FFT_file_data(fileToBeRead)
            if nargin ==1
                obj.fileAddress = fileToBeRead;
            end            
        end 
        function parse(obj)
            obj.fileID = fopen(obj.fileAddress);
            obj.parse_header();
            obj.parse_data();
            fclose(obj.fileID);
        end
        function parse_header(obj)
            for count = 1:7
                obj.parse_header_line(fgetl(obj.fileID));
            end                
        end
        function parse_data(obj)
            obj.rawData = textscan(obj.fileID, '%f%f%f');
            obj.parse_raw_data;
        end  
        function parse_header_line(obj, fileLine)
            splitStrings = strsplit(fileLine, ':');
            if(length(splitStrings) == 2)
                obj.parse_file_info(splitStrings)
            else
                splitStrings = strsplit(fileLine, '\t');
                if ~isempty(splitStrings{1})
                    obj.parse_data_units(splitStrings)
                end
            end
        end
        function parse_file_info(obj, splitStrings)
            switch splitStrings{1}
                case 'Source File Name'
                    obj.parse_file_name(splitStrings{2});
                case 'Point Index'
                    obj.parse_point_data(splitStrings{2});
                case 'Component'
                    obj.parse_component(splitStrings{2});
                case 'Signal'
                    obj.parse_signal(splitStrings{2});
                otherwise
                    warning('not expected data format')
            end
        end
        function parse_data_units(obj, splitStrings)
            if contains(splitStrings{1}, '[')
                obj.units = splitStrings;
            else
                obj.dataHeaders = splitStrings;
                obj.dataHeaders{2} = strcat(splitStrings{2},' -', obj.measurementType(4:end));
                obj.dataHeaders{3} = strcat(splitStrings{3},' -', obj.measurementType(4:end));
            end            
        end        
        function parse_file_name(obj, string)
            obj.nameOfDataSet = strtrim(string);
        end
        function parse_point_data(obj, string)
            obj.point = GeoPoint(strtrim(string));
        end
        function parse_component(obj, string)
            obj.componentType = strtrim(string);
        end
        function parse_signal(obj, string)
            signal = strsplit(string, ' - ');
            obj.scanType = strtrim(signal{1});
            obj.measurementType = strtrim(signal{2});
            obj.outputType = strtrim(signal{3});
        end
        function parse_raw_data(obj)
            obj.parse_frequencies();
            obj.parse_real_data();
            obj.parse_imaginary_data(); 
        end
        function parse_frequencies(obj)
            col = find(strcmp(obj.dataHeaders(:), 'Frequency'));
            obj.frequencies = obj.rawData{col};
        end
        function parse_real_data(obj)
            col = find(contains(obj.dataHeaders(:), 'Real'));
            obj.realData = obj.rawData{col};
        end
        function parse_imaginary_data(obj)
            col = find(contains(obj.dataHeaders(:), 'Imaginary'));
            obj.imaginaryData = obj.rawData{col};
        end
        function choose_file_to_parse(obj)
            obj.fileAddress = uigetfile(pwd, 'Select a file to parse');
        end
    end
end

