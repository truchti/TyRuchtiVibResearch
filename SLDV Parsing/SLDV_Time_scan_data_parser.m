classdef SLDV_Time_scan_data_parser   
    properties
        mainDirectory = '';
        fileName = '';
        scanType  = '';
        outputType  = '';
        componentType = '';
        measurmentTypes  = '';
        points =[];
        times = [];
        displacement = [];
        velocity = [];
        acceleration = [];
        headers = []
        units = [];
    end
    
    methods
        function obj = SLDV_Time_scan_data_parser(directory)
            if nargin == 1
                obj.mainDirectory = directory;
            end
        end
        function obj = choose_main_directory(obj)
            obj.mainDirectory = uigetdir(pwd, 'Select a folder');
        end
        function obj = parse(obj)
        end
        function obj = read_first_file(obj, file_name)
            %% read in header information
            fileID = fopen(file_name);
            line = fgetl(fileID);
            C = strsplit(line, ':');
            while(~isempty(C))
                switch C{1}
                    case 'Source File Name'
                       obj.fileName = strtrim(C{2});
                    case 'Point Index'
                        obj.points = GeoPoint(strtrim(C{2}));
                    case 'Component'
                        obj.componentType = strtrim(C{2});
                    case 'Signal'
                        signal = strsplit(C{2}, ' - ');
                        obj.scanType = strtrim(signal{1});
                         measureType = strtrim(signal{2});
                         obj.measurmentTypes = measureType;
                        obj.outputType = strtrim(signal{3});
                    otherwise
                end
                line = fgetl(fileID);
                if(isempty(line))
                    C = {};
                else
                    C = strsplit(line, ':');
                end
            end
            %% read in data header and units
            obj.headers = strsplit(fgetl(fileID), '\t');
            obj.units = strsplit(fgetl(fileID), '\t');
            data = textscan(fileID, '%f%f');
            obj.times = data{1};
            switch measureType
                case 'Vib Displacement'
                    obj.displacement = data{2};
                case 'Vib Velocity'
                    obj.velocity = data{2};
                case 'Vib Acceleration'
                    obj.acceleration = data{2};
                otherwise
                    warning("unknown data type")
            end
            fclose(fileID);

        end
        function obj = read_file(obj, file_name)
           fileID = fopen(file_name);
           line = fgetl(fileID);
           C = strsplit(line, ':');
            while(~isempty(C))
                switch C{1}
                    case 'Source File Name'
                       nameCheck = strtrim(C{2});
                       check_for_match(obj, 'FileName' , nameCheck);
                    case 'Point Index'
                        pointCheck = GeoPoint(strtrim(C{2}));
                        check_for_match(obj, 'Point' , pointCheck)
                    case 'Component'
                        obj.componentType = strtrim(C{2});
                        check_for_match(obj, 'Component', nameCheck);
                    case 'Signal'
                        signal = strsplit(C{2}, ' - ');
                        obj.scanType = strtrim(signal{1});
                        check_for_match(obj, property, nameCheck);
                        measureType = strtrim(signal{2});
                        check_for_match(obj, property, nameCheck);
                        obj.outputType = strtrim(signal{3});
                        check_for_match(obj, property, nameCheck);
                    otherwise
                end
                line = fgetl(fileID);
                if(isempty(line))
                    C = {};
                else
                    C = strsplit(line, ':');
                end
            end
            %% read in data header and units
            obj.headers = strsplit(fgetl(fileID), '\t');
            obj.units = strsplit(fgetl(fileID), '\t');
            data = textscan(fileID, '%f%f');
            obj.times = data{1};
            switch measureType
                case 'Vib Displacement'
                    obj.displacement = data{2};
                case 'Vib Velocity'
                    obj.velocity = data{2};
                case 'Vib Acceleration'
                    obj.acceleration = data{2};
                otherwise
                    warning("unknown data type")
            end
            fclose(fileID);
        end
        function obj = read_sub_directory(obj)
            obj.mainDirectory = directory;
        end
        function obj = read_all_subs(obj)

        end
        %% setters
        function obj = set.mainDirectory(obj, directory)
            obj.mainDirectory = directory;
        end
        function obj = set.fileName(obj, fileName)
            obj.fileName = fileName;
        end
        function obj = set.scanType(obj, scanType)
            obj.scanType = scanType;
        end
        function obj = set.outputType(obj, outputType)
            obj.outputType = outputType;
        end
        function obj = set.componentType(obj, componentType)
            obj.componentType = componentType;
        end
        function obj = set.measurmentTypes(obj, measurementTypes)
            obj.measurmentTypes = measurementTypes;
        end
        function obj = set.times(obj, timeVector)
            obj.times = timeVector;
        end
        %% getters
        function directory = get.mainDirectory(obj)
            directory = obj.mainDirectory;
        end
        function fileName = get.fileName(obj)
            fileName = obj.fileName;
        end
        function scanType = get.scanType(obj)
            scanType = obj.scanType;
        end
        function outputType = get.outputType(obj)
            outputType = obj.outputType;
        end
        function componentType = get.componentType(obj)
            componentType = obj.componentType;
        end
        function measureTypes = get.measurmentTypes(obj)
            measureTypes = obj.measurmentTypes;
        end
        function times = get.times(obj)
            times = obj.times;
        end
%         function obj = set_position(obj)
%         end
%         function obj = set_velocity(obj)
%         end
%         function obj = set_acceleration(obj)
%         end
        %%

        function outputArg = method1(obj,inputArg)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
        end
        
        
    end
end





