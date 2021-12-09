classdef SLDV3D_FFT_parser < handle
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
        isParsed = false;
        saveDataStructure;
    end
    methods 
        function obj = SLDV3D_FFT_parser(folder)
            if nargin == 1 && exist(folder, 'dir')
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
            end
            cd(oldFolder);
            obj.isParsed = true;
            if isempty(obj.realVelocity) && ~isempty(obj.realDisplacement)
                obj.calculate_velocites_from_displacements_and_frequency();
            end
        end
        function save_as_data_structure(obj, frequency)
            if nargin == 1
                frequency = obj.frequencies;
            end
            obj.create_save_data_structure(frequency);
            defaultName = strcat(obj.nameOfDataSet);
            filter = '*.mat';
            windowTitle = 'Save data file';
            [file, path] = uiputfile(filter, windowTitle, defaultName);
            
            if ~ischar(file) || ~ischar(path)
                warning("File Not Saved")
            else
                old = cd(path);
                fileName = fullfile(path,file);
                dataStructure = obj.saveDataStructure;
                save(fileName, 'dataStructure')
                cd(old)
            end
        end
    end
    methods (Hidden)
        function calculate_velocites_from_displacements_and_frequency(obj)
            obj.realVelocity = zeros(size(obj.realDisplacement));
            obj.imaginaryVelocity = zeros(size(obj.imaginaryDisplacement));
            omegas = 2*pi*obj.frequencies;
            for frq_num = 1:length(obj.frequencies)
                obj.realVelocity(frq_num, :, :) = -omegas(frq_num)*obj.imaginaryDisplacement(frq_num,:,:);
                obj.imaginaryVelocity(frq_num, :,:) = -omegas(frq_num)*obj.realDisplacement(frq_num,:,:);
            end
        end
        function create_save_data_structure(obj, frequency)
            indexes = find(obj.frequencies == frequency);
            obj.saveDataStructure = SLDV_FFT_Data(frequency, obj.points, permute(complex(obj.realDisplacement(indexes,:,:), obj.imaginaryDisplacement(indexes,:,:)), [3 2 1]),  ...
                                             permute(complex(obj.realVelocity(indexes,:,:), obj.imaginaryVelocity(indexes,:,:)), [3 2 1]));
        end
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
            wb = waitbar(0, 'Parsing');
            for file = 1:length(filesInDirectory)
                obj.parse_file(filesInDirectory(file).name);
                waitbar(file/length(filesInDirectory), wb);
            end
            close(wb)
            cd(old);
        end
        function parse_file(obj, file)
            obj.tempFileData = FFT_3D_file_data(file);
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
            elseif (strcmp(obj.scanType, obj.tempFileData.scanType) && ...
                    strcmp(obj.componentType, obj.tempFileData.componentType) && ...
                    strcmp(obj.outputType, obj.tempFileData.outputType))
%                 strcmp(obj.nameOfDataSet, obj.tempFileData.nameOfDataSet) && ...
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
            if contains(obj.tempFileData.measurementType, 'Velocity')
                obj.add_velocity(pointIndex)
            elseif contains(obj.tempFileData.measurementType, 'Displacement')
                obj.add_displacement(pointIndex)
            else 
                warning("There is no match for the measurementType %s\n", obj.tempFileData.measurementType)
            end
        end
        function add_velocity(obj, pointIndex)
            obj.realVelocity(:, :, pointIndex) = obj.tempFileData.realData;
            obj.imaginaryVelocity(:, :, pointIndex) = obj.tempFileData.imaginaryData;
        end
        function add_displacement(obj, pointIndex)
            obj.realDisplacement(:, :, pointIndex) = obj.tempFileData.realData;
            obj.imaginaryDisplacement(:, :, pointIndex) = obj.tempFileData.imaginaryData;
        end
        function check_new_frequencies(obj)
            if( ~isequal(obj.frequencies, obj.tempFileData.frequencies))
                warning('The frequencies to not line up for point %d\n', obj.tempFileData.point.number)
            end            
        end
        function clean_parser(obj)
            obj.points = GeoPoint.empty;
            obj.mapPointNumberToIndex = containers.Map('KeyType','int32','ValueType','int32');
            obj.dataHeaders = {};
            obj.units = {};
            obj.frequencies = [];
            obj.realVelocity = [];
            obj.realDisplacement = [];
            obj.imaginaryVelocity = [];
            obj.imaginaryDisplacement = [];
            obj.nameOfDataSet = '';
            obj.scanType = '';
            obj.outputType = '';
            obj.componentType = '';
            obj.measurementTypes = {};
            obj.isParsed = false;
        end
    end
end