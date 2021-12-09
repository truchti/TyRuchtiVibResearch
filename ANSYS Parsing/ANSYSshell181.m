classdef ANSYSshell181 < handle
    properties
        number
        nodeNums
    end
    properties (Hidden = true)
        N11
        N22
        N12
        M11
        M22
        M12
        Q1
        Q2
    end
    methods
        function obj =ANSYSshell181(number, nodeNumbers)
            if nargin == 2
                obj.number = number ;
                obj.nodeNums = nodeNumbers;
            else
                if nargin ~= 0
                    warning("Incorrect inputs")
                end
            end            
        end
    end
    methods (Hidden = true)
        function add_resultant(obj, type, value)
            switch type
                case 'N11'
                    obj.N11 = value;
                case 'N22'
                    obj.N22 = value;
                case 'N12' 
                    obj.N12 = value;
                case 'M11'
                    obj.M11 = value;
                case 'M22'
                    obj.M22 = value;
                case 'M12'
                    obj.M12 = value;
                case 'Q1'
                    obj.Q1 = value;
                case 'Q2'
                    obj.Q2 = value;
            end
        end
        function value = get_resultant(obj, type)
            switch type
                case 'N11'
                    value = obj.N11;
                case 'N22'
                    value = obj.N22;
                case 'N12' 
                    value = obj.N12;
                case 'M11'
                    value = obj.M11;
                case 'M22'
                    value = obj.M22;
                case 'M12'
                    value = obj.M12;
                case 'Q1'
                    value = obj.Q1;
                case 'Q2'
                    value = obj.Q2;
            end
        end
    end
end