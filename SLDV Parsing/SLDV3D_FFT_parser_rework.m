classdef SLDV3D_FFT_parser_rework < handle
    properties
        directory = '';
        nameOfDataSet
        pointNumbers
        coordinates
        frequencies
        scanType
        outputType
        componentType
        measurementTypes = {};
    end
    properties (Hidden)
        tempFileData
        extension = '.txt';
        dataHeaders
        units
        realVelocity = [];
        realDisplacement = [];
        imaginaryVelocity = [];
        imaginaryDisplacement = [];
        isParsed = false;
    end
    methods 
        function obj = SLDV3D_FFT_parser_rework(folder)
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
                warning("Calculating velocity from displacement")
                obj.calculate_velocities_from_displacements_and_frequency();
            end
            if isempty(obj.realDisplacement) && ~isempty(obj.realVelocity)
                warning("Calculating displacement from velocity")
                obj.calculate_displacements_from_velocities_and_frequency()
            end
        end
        function data = export_SLDV_data_at_frequency(obj, frequency)
            if nargin == 1
                frequency = obj.frequencies(1);
            end
            data = obj.create_SLDV_data_structure(frequency);          
        end
    end
    methods (Hidden)
        function calculate_velocities_from_displacements_and_frequency(obj)
            obj.realVelocity = zeros(size(obj.realDisplacement));
            obj.imaginaryVelocity = zeros(size(obj.imaginaryDisplacement));
            omegas = 2*pi*obj.frequencies;
            for frq_num = 1:length(obj.frequencies)
                obj.realVelocity(frq_num, :, :) = -omegas(frq_num)*obj.imaginaryDisplacement(frq_num,:,:);
                obj.imaginaryVelocity(frq_num, :,:) = -omegas(frq_num)*obj.realDisplacement(frq_num,:,:);
            end
        end
        function calculate_displacements_from_velocities_and_frequency(obj)
            obj.realDisplacement = zeros(size(obj.realVelocity));
            obj.imaginaryDisplacement = zeros(size(obj.imaginaryVelocity));
            omegas = 2*pi*obj.frequencies;
            for frq_num = 1:length(obj.frequencies)
                obj.realDisplacement(frq_num, :, :) = obj.imaginaryVelocity(frq_num,:,:)/-omegas(frq_num);
                obj.imaginaryDisplacement(frq_num, :,:) = obj.realVelocity(frq_num,:,:)/-omegas(frq_num);
            end
        end
        function data = create_SLDV_data_structure(obj, frequency)
            coors = obj.coordinates;
            reDisp = obj.get_displacements_at_frequency(frequency);
            imDisp = obj.get_displacements_at_frequency(frequency, true);
            reVel  = obj.get_velocities_at_frequency(frequency);
            imVel  = obj.get_velocities_at_frequency(frequency,true);
            data = SLDV_data(coors, reDisp, imDisp, reVel, imVel, frequency);
        end
        function displacements = get_displacements_at_frequency(obj, frequency, imaginary)
            if nargin < 3
                imaginary = false;
            end
            index = obj.get_index_of_frequency(frequency);
            if imaginary
                displacements = permute(obj.imaginaryDisplacement(index,:,:), [3 2 1]);
            else
                displacements = permute(obj.realDisplacement(index,:,:), [3 2 1]);
            end
        end
        function velocities = get_velocities_at_frequency(obj, frequency, imaginary)
            if nargin < 3
                imaginary = false;
            end
            index = obj.get_index_of_frequency(frequency);
            if imaginary
                velocities = permute(obj.imaginaryVelocity(index,:,:), [3 2 1]);
            else
                velocities = permute(obj.realVelocity(index,:,:), [3 2 1]);
            end
        end
        function index = get_index_of_frequency(obj, frequency)
            index = find(obj.frequencies == frequency);
            if isempty(index)
                warning("frequency requested not in data set. Using nearest frequency")
                t = find(obj.frequencies > frequency, 1, 'first');
                if abs(obj.frequencies(t)-frequency) < abs(freqency-obj.frequencies(t-1))
                    index = t;
                else
                    index = t-1;
                end
            end
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
%             pointIndex = obj.find_or_add_point_index();
            obj.add_data_point();
            obj.check_for_and_add_new_measurement_types()
            obj.add_data_based_on_types();
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
        
        function add_data_point(obj)
            if isempty(find( obj.pointNumbers == obj.tempFileData.pointNumber, 1))
                obj.pointNumbers = [obj.pointNumbers; obj.tempFileData.pointNumber];
                obj.coordinates = [obj.coordinates; obj.tempFileData.coordinate];
            end
        end
        function add_data_based_on_types(obj)
            if isempty(obj.frequencies)
                obj.frequencies = obj.tempFileData.frequencies;
            else
                obj.check_new_frequencies();
            end
            if contains(obj.tempFileData.measurementType, 'Velocity')
                obj.add_velocity()
            elseif contains(obj.tempFileData.measurementType, 'Displacement')
                obj.add_displacement()
            else 
                warning("There is no match for the measurementType %s\n", obj.tempFileData.measurementType)
            end
        end
        function add_velocity(obj)
            if isempty(obj.realVelocity)
                i = 0;
            else
                i = size(obj.realVelocity,3);
            end
            obj.realVelocity(:, :, i+1) = obj.tempFileData.realData;
            obj.imaginaryVelocity(:, :, i+1) = obj.tempFileData.imaginaryData;
        end
        function add_displacement(obj)
            if isempty(obj.realDisplacement)
                i = 0;
            else
                i = size(obj.realDisplacement,3);
            end
            obj.realDisplacement(:, :, i+1) = obj.tempFileData.realData;
            obj.imaginaryDisplacement(:, :, i+1) = obj.tempFileData.imaginaryData;
        end
        function check_new_frequencies(obj)
            if( ~isequal(obj.frequencies, obj.tempFileData.frequencies))
                warning('The frequencies to not line up for point %d\n', obj.tempFileData.point.number)
            end            
        end
        function clean_parser(obj)
%             obj.points = GeoPoint.empty;
%             obj.mapPointNumberToIndex = containers.Map('KeyType','int32','ValueType','int32');
            obj.pointNumbers = [];
            obj.coordinates = [];
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