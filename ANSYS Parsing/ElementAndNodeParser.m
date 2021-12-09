classdef ElementAndNodeParser < handle
    % this class takes in the ds.dat file from ansys work bench results and
    % pulls out the node locations element numbers and the nodes associated
    % with each element
    % it can also take in exported ansys nodal or elemental results and
    % store them in the cooresponding data structure,
    % it then links the nodes and elements and exports a mesh object
    properties
        nodes = ANSYSnode.empty
        elements = ANSYSshell181.empty
        datFile
        dataDirectory
    end
    properties(Hidden = true)
        foundNodes = false
        foundElements = false
        fid
        nodeData
        eleData
        realElementData
        imagElementData
        realNodeData
        imagNodeData
        projectName = 'Model'
        listOfElemResults = {'N11', 'N22', 'N12', 'M11', 'M22', 'M12', 'Q1', 'Q2'};
        listOfNodalResults = {'u', 'v', 'w', 'udot', 'vdot', 'wdot'};
    end
    methods
        function obj = ElementAndNodeParser(file, dataDirectory)
            if nargin < 1
                obj.choose_ds_dat_file();
            else
                obj.datFile = file;
            end
            if nargin < 2
                obj.choose_data_directory();
            else
            	obj.dataDirectory = dataDirectory;
            end
        end
        function parse(obj)
            % open file
            obj.fid = fopen(obj.datFile);
            %read in line until I find either nblock or eblock
            while(~(obj.foundNodes && obj.foundElements))
                string = fgetl(obj.fid);
                if ischar(string) && contains(string, 'eblock')
                    obj.read_in_elements();
                end
                if ischar(string) && contains(string, 'nblock')
                    obj.read_in_nodes();
                end
            end
            fclose(obj.fid);
            obj.create_nodes_from_data();
            obj.create_elements_from_data();
            obj.read_in_all_element_results();
            obj.read_in_all_nodal_results();
            obj.add_element_relation_to_nodes();
        end        
        function mesh = export_mesh_object(obj)
            mesh = ANSYS_181Shell_Mesh(obj.nodes, obj.elements);
            mesh.extrapolate_all_resultants_to_nodes();
            mesh.calculate_theta_dots();
        end
        end
    methods (Hidden = true)
        function set_as_cylinder(obj)
            obj.listOfNodalResults = {'t','r', 'z', 'tdot', 'rdot', 'zdot'};
            obj.isCylinder = true;
        end
        function choose_ds_dat_file(obj)
            [fileName, pathName] = uigetfile('*.dat', 'Select the ds.dat file for the ANSYS Results');
            obj.datFile = strcat(pathName, fileName);
        end
        function choose_data_directory(obj)
            obj.dataDirectory = uigetdir(pwd, 'Select a folder containing resultant files');
        end
        function create_nodes_from_data(obj)
            numbers = obj.nodeData{1};
            xs = obj.nodeData{2};
            ys = obj.nodeData{3};
            zs = obj.nodeData{4};
            for node = length(numbers)-1:-1:1
                obj.nodes(numbers(node)) = ANSYSnode(numbers(node), xs(node), ys(node), zs(node));
            end
        end
        function create_elements_from_data(obj)
            numbers = obj.eleData{1};
            node1Nums = obj.eleData{2};
            node2Nums = obj.eleData{3};
            node3Nums = obj.eleData{4};
            node4Nums = obj.eleData{5};
            for ele = length(numbers)-1:-1:1
                obj.elements(ele) = ANSYSshell181(numbers(ele), [node1Nums(ele), node2Nums(ele), node3Nums(ele), node4Nums(ele)]);
            end
        end
        function add_element_relation_to_nodes(obj)
            for i = 1:length(obj.elements)
                el = obj.elements(i);
                for j = 1:length(el.nodeNums)
                    obj.nodes(el.nodeNums(j)).link_to_an_element(el.number);
                end
            end
        end
        function read_in_all_element_results(obj)
            if isempty(obj.dataDirectory)
                obj.choose_data_directory();
            end
            for i = 1:length(obj.listOfElemResults)
                realFileName = [obj.dataDirectory, '\', obj.projectName, '_', obj.listOfElemResults{i}, ' Real.csv'];
                imagFileName  = [obj.dataDirectory,  '\', obj.projectName, '_', obj.listOfElemResults{i}, ' Imag.csv'];
                obj.parse_complex_element_resultants(realFileName, imagFileName)
                obj.add_resultant_data_to_respective_element_object(obj.listOfElemResults{i});
            end
        end
        function read_in_all_nodal_results(obj)
            if isempty(obj.dataDirectory)
                obj.choose_data_directory();
            end
            for i = 1:length(obj.listOfNodalResults)
                realFileName = [obj.dataDirectory, '\', obj.projectName, '_', obj.listOfNodalResults{i}, ' Real.csv'];
                imagFileName = [obj.dataDirectory, '\', obj.projectName, '_', obj.listOfNodalResults{i}, ' Imag.csv'];
                obj.parse_complex_node_results(realFileName, imagFileName)
                obj.add_data_to_respective_node_object(obj.listOfNodalResults{i});
            end
        end
        function read_in_nodes(obj)
            if ~obj.foundNodes
                % read past the () line
                fgetl(obj.fid);
                % auto stops when it reached end because -1 line doesn't match
                % format
                obj.nodeData = textscan(obj.fid, '%f%f%f%f');
                obj.foundNodes = true;
            else
                warning("There were remote Point Nodes and they are not included")
            end
        end
        function read_in_elements(obj)
            % read past the () line
            fgetl(obj.fid);
            data = textscan(obj.fid, repmat('%d', [1 15]));
            obj.eleData = data(11:15);
            obj.foundElements = true;
        end
        function add_resultant_data_to_respective_element_object(obj, type)
            for i = 1:length(obj.elements)
                 if ~isempty(obj.realElementData)
                    realPart = obj.realElementData(i);
                else
                    realPart = 0;
                end
                if ~isempty(obj.imagElementData)
                    imagPart = obj.imagElementData(i);
                else
                    imagPart = 0;
                end
                obj.elements(i).add_resultant(type, realPart+1i*imagPart)
            end
        end
        function add_data_to_respective_node_object(obj, type)
            for i = 1:length(obj.nodes)
                if ~isempty(obj.realNodeData)
                    realPart = obj.realNodeData(i);
                else
                    realPart = 0;
                end
                if ~isempty(obj.imagNodeData)
                    imagPart = obj.imagNodeData(i);
                else
                    imagPart = 0;
                end
                obj.nodes(i).add_disp_or_vel(type, realPart + 1i*imagPart)
            end
            obj.realNodeData = [];
            obj.imagNodeData = [];
        end
        function parse_complex_element_resultants(obj, realFile, imagFile)
            if nargin < 2
                [fileName, pathName] = uigetfile('.txt', 'Select the txt file with the real element results');
                realFile = strcat(pathName, fileName);
                [fileName, pathName] = uigetfile('.txt', 'Select the txt file with the imaginary element results');
                imagFile = strcat(pathName, fileName);
            end
            realID = fopen(realFile);
            % first line skip
            fgetl(realID);
            rawData = textscan(realID, '%d%f');
            fclose(realID);
            if (~isempty(rawData{1}))
                obj.realElementData = containers.Map(rawData{1}, rawData{2});
            end
            imagID = fopen(imagFile);
            % first line skip
            fgetl(imagID);
            rawData = textscan(imagID, '%d%f');
            fclose(imagID);
            if(~isempty(rawData{1}))
                obj.imagElementData = containers.Map(rawData{1}, rawData{2});
            end
        end
        function parse_complex_node_results(obj, realFile, imagFile)
            if nargin < 2
                [fileName, pathName] = uigetfile('.txt', 'Select the txt file with the real nodal displacement');
                realFile = strcat(pathName, fileName);
                [fileName, pathName] = uigetfile('.txt', 'Select the txt file with the imaginary nodal displacement');
                imagFile = strcat(pathName, fileName);
            end
            realID = fopen(realFile);
            if realID > 1
            % first line skip
                fgetl(realID);
                rawData = textscan(realID, '%d%f');
                fclose(realID);
                obj.realNodeData = containers.Map(rawData{1}, rawData{2});
            else
                wrn = ['Bad File: ', realFile];
                warning(wrn)
            end           
            imagID = fopen(imagFile);
            % first line skip
            if imagID >1
                fgetl(imagID);
                rawData = textscan(imagID, '%d%f');
                fclose(imagID);
                obj.imagNodeData = containers.Map(rawData{1}, rawData{2});
            else
                wrn = ['Bad File: ', imagFile];
                warning(wrn)
            end 
        end
    end
end