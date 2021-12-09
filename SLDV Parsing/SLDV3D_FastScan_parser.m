classdef SLDV3D_FastScan_parser < handle
    properties
        velFile
        dispFile
        nameOfData
        frequency
        pointNumbers
        coordinates
        realVel
        imagVel
        realDisp
        imagDisp
        scanType
        outputType
        measurementType = {}
    end
    properties (Hidden = true)
        dispFileID
        velFileID
        rawData
        isParsed = false;
        band
        filtered 
        interpolated
        dataHeaders = {}
        SLDVdata
    end
    methods 
        function obj = SLDV3D_FastScan_parser(dispFile, velFile)
            if nargin ==2
                obj.dispFile = dispFile;
                obj.velFile = velFile;
            end
        end
        function parse(obj)
            % this function takes the ASCII displacement and velocity files created by the SLDV when 
            % running a fast scan and compiles a data structure that
            % condenses all the info
            
            % if the files aren't given open a selection menu
            if isempty(obj.dispFile)
                obj.choose_displacement_file();
            end
            if isempty(obj.velFile)
                obj.choose_velocity_file();
            end
            % open both the displacement and velocity files
            obj.dispFileID = fopen(obj.dispFile);
            obj.velFileID = fopen(obj.velFile);
            % parse the headers of both files
            obj.parse_header();
            % parse the data of both files
            obj.parse_data();
            % close files
            fclose(obj.dispFileID);
            fclose(obj.velFileID);
            % mark that the data has been parsed
            obj.isParsed = true;            
        end
        function SLDVdata = export_data(obj)
            % export the condensed matlab data structure of the data
            if isempty(obj.SLDVdata)
                obj.create_SLDV_data_structure
            end
            SLDVdata = obj.SLDVdata;
        end
    end
    methods (Access = private)
        % set up functions
        function create_SLDV_data_structure(obj)
            coors = obj.coordinates;
            reDisp = obj.realDisp;
            imDisp = obj.imagDisp;
            reVel = obj.realVel;
            imVel = obj.imagVel;
            obj.SLDVdata = SLDV_data(coors, reDisp, imDisp, reVel, imVel, obj.frequency);
        end
        function choose_displacement_file(obj)
            [fileName, pathName] = uigetfile('*.txt', 'Select the displacement file to parse');
            obj.dispFile = strcat(pathName, fileName);
        end
        function choose_velocity_file(obj)
            [fileName, pathName] = uigetfile('*.txt', 'Select the velocity file to parse');
            obj.velFile = strcat(pathName, fileName);
        end
        % Header Parsing function and helpers
        function parse_header(obj)
            % This function reads the header lines of the file
            for count = 1:11
                obj.parse_header_line(fgetl(obj.dispFileID));
                obj.parse_header_line(fgetl(obj.velFileID));
            end                
        end
        function parse_header_line(obj, fileLine)
            splitStrings = strsplit(fileLine, ':');
            if(length(splitStrings) == 2)
                obj.parse_file_info(splitStrings)
            else
                splitStrings = strsplit(fileLine, '\t');
                if ~isempty(splitStrings{1})
                    obj.dataHeaders = [obj.dataHeaders; splitStrings];
                end
            end
        end
        function parse_file_info(obj, splitStrings)
            switch splitStrings{1}
                case 'Source File Name'
                    obj.parse_file_name(splitStrings{2});
                case 'Signal'
                    obj.parse_signal(splitStrings{2});
                case 'Band No.'
                    obj.parse_band(splitStrings{2});
                case 'Frequency'
                    obj.parse_frequency(splitStrings{2});
                case 'Interpolated'
                    obj.parse_interpolated(splitStrings{2});
                case 'Filtered'
                    obj.parse_filtered(splitStrings{2});
                otherwise
                    warning('not expected data format')
            end
        end
        function parse_file_name(obj, string)
            if ~isempty(obj.nameOfData)
                checkstring = strtrim(string);
                if ~strcmp(checkstring, obj.nameOfData)
                    warning("Velocity and Displacement file names don't match")
                end
            else
                obj.nameOfData = strtrim(string);
            end
        end
        function parse_signal(obj, string)
            signal = strsplit(string, ' - ');
            obj.scanType = strtrim(signal{1});
            obj.measurementType = [obj.measurementType,strtrim(signal{2})];
            obj.outputType = strtrim(signal{3});
        end
        function parse_band(obj, string)
            obj.band = str2double(string);
        end
        function parse_frequency(obj, string)
            splitString = strsplit(string, ' ');
            if strcmp(splitString{2}, 'Hz')
                obj.frequency = str2double(splitString{1});
            end
        end
        function parse_interpolated(obj, string)
            if strcmp(string,'Yes')
                obj.interpolated = true;
            else
                obj.interpolated = false;
            end
        end
        function parse_filtered(obj, string)
            if strcmp(string,'Yes')
                obj.filtered = true;
            else
                obj.filtered = false;
            end
        end
        % Data Parsing function and helpers
        function parse_data(obj)
            % this function takes all of the raw data from a the
            % displacement and velocity files and joins them into one
            % data structure
            obj.read_in_raw_data(obj.dispFileID)
            obj.parse_point_numbers();
            obj.parse_point_coordinates();
            obj.parse_displacement
            obj.read_in_raw_data(obj.velFileID)
            obj.parse_velocity();
        end
        function read_in_raw_data(obj, fileID)
            % This function reads in all of the data from a single file starting
            % after the header lines
            scanString = repmat('%f', 1, length(obj.dataHeaders));
            obj.rawData = textscan(fileID, scanString);
        end
        function parse_point_numbers(obj)
            indx = obj.findCellIndex(obj.dataHeaders, 'Index');
            obj.pointNumbers = obj.rawData{indx};
        end
        function parse_point_coordinates(obj)
            obj.coordinates = [obj.rawData{2}, obj.rawData{3}, obj.rawData{4}];
        end
        function parse_displacement(obj)
            obj.realDisp = [obj.rawData{5}, obj.rawData{6}, obj.rawData{7}];
            obj.imagDisp = [obj.rawData{8}, obj.rawData{9}, obj.rawData{10}];   
        end
        function parse_velocity(obj)
            obj.realVel = [obj.rawData{5}, obj.rawData{6}, obj.rawData{7}];
            obj.imagVel = [obj.rawData{8}, obj.rawData{9}, obj.rawData{10}];   
        end
    end
    methods (Static = true, Hidden = true)
        function indx = findCellIndex(cellArray, str2find)
            searchresults = strfind(cellArray, str2find);
            indx = find(~(cellfun('isempty',searchresults))); %#ok<STRCL1>
        end
    end
end
