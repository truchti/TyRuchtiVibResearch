classdef SLDV_FFT_parser < handle

    properties
        directory = '';
        nameOfDataSet
        scanType
        outputType
        componentType
        measurementTypes = {};
    end
    properties (Hidden)
        tempFileData
        extension = '.txt';
        mapPointNumberToIndex = containers.Map('KeyType','int32','ValueType','int32');
        points = []
        dataHeaders
        units
        frequencies
        realVelocity = [];
        realDisplacement = [];
        imaginaryVelocity = [];
        imaginaryDisplacement = [];
        data = powerflowDataStructure();
    end
    methods
        function obj = SLDV_FFT_parser(folder)
            if nargin == 1
                obj.directory = folder;
            end
        end
        function choose_directory(obj)
            obj.directory = uigetdir(pwd, 'Select a folder');
        end
        function parse(obj)
            obj.clean_parser;
            if ~exist(obj.directory, 'dir')
                obj.choose_directory;
            end
            oldFolder = cd(obj.directory);
            if(obj.check_for_sub_directories)
                subDirectories = obj.get_sub_directories;
                for i = 1:length(subDirectories)
                    fprintf("Now parsing files in %s subfolder\n", subDirectories{i})
                    obj.parse_directory(subDirectories{i})
                end
            else
                obj.parse_directory
            cd(oldFolder);
            end
        end
%         function export_data_as_powerflow_structure(obj)
%             obj.data.directoryToSaveIn = obj.directory;
%             obj.data.nameOfDataSet = obj.nameOfDataSet;
%             obj.data.frequencies = obj.frequencies;
%             obj.data.points = obj.points;
%             obj.data.displacements = obj.realDisplacement + sqrt(-1)*obj.imaginaryDisplacement;
%             obj.data.velocities = obj.realVelocity + sqrt(-1)*obj.imaginaryVelocity;
%             obj.data.save_data();
%         end
%         function parse_and_export_powerflow_structure(obj)
%             obj.parse;
%             obj.export_data_as_powerflow_structure;
%         end
%         function compute_and_save_partials_in_frequency_range(obj, frequencyRange)
%             obj.data.calculate_and_save_partials(frequencyRange)
%         end
    end
    methods (Hidden)
        function trueOrFalse = check_for_sub_directories(obj)
            trueOrFalse = false;
            if ~isempty(obj.get_sub_directories)
                trueOrFalse = true;
            end
        end
        function subDirectories = get_sub_directories(obj)
            allfiles = dir(obj.directory);
            dirFlags = [allfiles.isdir];
            subs = allfiles(dirFlags);
            subs = subs(3:end);
            subDirectories = {};
            for i = 1:length(subs)
                subDirectories{i} = subs(i).name; %#ok<AGROW>
            end
        end
        function parse_directory(obj, directory)
            if(nargin == 1)
                directory = obj.directory;
            end
            old = cd(directory);
            filesInDirectory = dir(strcat('*', obj.extension));
            for file = 1:length(filesInDirectory)
                obj.parse_file(filesInDirectory(file).name);
                if mod(file, 25) == 0
                    fprintf("Just finished parsing file %d of %d in this directory\n", file, length(filesInDirectory))
                end
            end
            cd(old);
        end
        function parse_file(obj, file)
            obj.tempFileData = FFT_file_data(file);
            obj.tempFileData.parse;
            obj.parse_check_and_merge_temp_data();
        end
        function parse_check_and_merge_temp_data(obj)
            obj.tempFileData.parse;
            if(obj.check_compatability_of_new_data)
                merge_new_data_into_parser(obj)
            else
                warning("Data from the following file was not compatable so it was skipped. \n File:\n")
                obj.tempFileData.fileAddress
            end
        end
        function trueOrFalse = check_compatability_of_new_data(obj)
            trueOrFalse = false;
            if  isempty(strcat(obj.nameOfDataSet, obj.scanType, obj.componentType, obj.outputType))
                obj.fill_file_information();
                trueOrFalse = true;
            elseif (strcmp(obj.nameOfDataSet, obj.tempFileData.nameOfDataSet) && ...
                    strcmp(obj.scanType, obj.tempFileData.scanType) && ...
                    strcmp(obj.componentType, obj.tempFileData.componentType) && ...
                    strcmp(obj.outputType, obj.tempFileData.outputType))
                trueOrFalse = true;
            end
        end
        function fill_file_information(obj)
           obj.nameOfDataSet =  obj.tempFileData.nameOfDataSet;
           obj.scanType =  obj.tempFileData.scanType;
           obj.outputType =  obj.tempFileData.outputType;
           obj.componentType =  obj.tempFileData.componentType;
        end
        function merge_new_data_into_parser(obj)
            pointIndex = obj.find_or_add_point_index();
            obj.add_point_using_index(pointIndex);
            obj.check_for_and_add_new_measurement_types()
            obj.add_data_based_on_types(pointIndex);
        end
        function check_for_and_add_new_measurement_types(obj)
            if isempty(find(strcmp(obj.measurementTypes(:), obj.tempFileData.measurementType), 1))
                obj.measurementTypes{end+1} = obj.tempFileData.measurementType;
            end
            for header = 1:length(obj.tempFileData.dataHeaders)
                if isempty(find(strcmp(obj.dataHeaders(:), obj.tempFileData.dataHeaders{header}), 1))
                    obj.dataHeaders{end+1} = obj.tempFileData.dataHeaders{header};
                    obj.units{end+1} = obj.tempFileData.units{header};
                end
            end
        end
        function index = find_or_add_point_index(obj)
            pointNumber = obj.tempFileData.point.number;
            if(isKey(obj.mapPointNumberToIndex, pointNumber))
                index = obj.mapPointNumberToIndex(pointNumber);
            else
                index = obj.create_new_point_index_in_map();
            end
        end
        function index = create_new_point_index_in_map(obj)
%             obj.points = [obj.points, obj.tempFileData.point];
            index = length(obj.points) + 1;
            obj.mapPointNumberToIndex(obj.tempFileData.point.number) = index; 
        end
        function index = add_point_using_index(obj, index)
            obj.points(index) = obj.tempFileData.point;
        end
        function add_data_based_on_types(obj, pointIndex)
            
            if isempty(obj.frequencies)
                obj.frequencies = obj.tempFileData.frequencies;
            else
                obj.check_new_frequencies();
            end
            if strcmp(obj.tempFileData.measurementType, "Vib Velocity")
                obj.add_velocity(pointIndex)
            elseif strcmp(obj.tempFileData.measurementType, "Vib Displacement")
                obj.add_displacement(pointIndex)
            end
        end
        function add_velocity(obj, pointIndex)
            obj.realVelocity(:, pointIndex) = obj.tempFileData.realData;
            obj.imaginaryVelocity(:, pointIndex) = obj.tempFileData.imaginaryData;
        end
        function add_displacement(obj, pointIndex)
            obj.realDisplacement(:, pointIndex) = obj.tempFileData.realData;
            obj.imaginaryDisplacement(:, pointIndex) = obj.tempFileData.imaginaryData;
        end
        function check_new_frequencies(obj)
            if( ~isequal(obj.frequencies, obj.tempFileData.frequencies))
                warning('The frequencies to not line up for the following point\n')
                disp(obj.tempFileData.point.number)
            end            
        end
        function clean_parser(obj)
            obj.points = GeoPoint.empty;
            obj.mapPointNumberToIndex = containers.Map('KeyType','int32','ValueType','int32');
            obj.dataHeaders = {};
            obj.units = {};
            obj.frequencies
            obj.realVelocity = [];
            obj.realDisplacement = [];
            obj.imaginaryVelocity = [];
            obj.imaginaryDisplacement = [];
            obj.nameOfDataSet = '';
            obj.scanType = '';
            obj.outputType = '';
            obj.componentType = '';
            obj.measurementTypes = {};
        end
    end
end

