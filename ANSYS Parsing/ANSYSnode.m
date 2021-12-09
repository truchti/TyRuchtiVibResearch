classdef ANSYSnode < handle
    properties
        number
        x
        y
        z
    end
    properties(Dependent)
        coors
    end
    properties (Hidden = true)
        u
        v
        w
        udot
        vdot
        wdot
        thetaDotX
        thetaDotY
        associatedElements
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
        function obj = ANSYSnode(number, x,y,z)
            if nargin < 4
            else
                obj.number = number;
                obj.x = x;
                obj.y = y;
                obj.z = z;
            end
        end
        function value = get.coors(obj)
            value = [obj.x, obj.y, obj.z];
        end
    end
    methods(Hidden = true)
        function value = get_disp_or_vel(obj, type)
            switch type
                case 'u'
                    value = obj.get_displacement('u');
                case 'v'
                    value = obj.get_displacement('v');
                case 'w'
                    value = obj.get_displacement('w');
                case {'udot', 'du'}
                    value = obj.get_velocity('u');
                case {'vdot', 'dv'}
                    value = obj.get_velocity('v');
                case {'wdot', 'dw'}
                    value = obj.get_velocity('w');
                case 'thetaDotX'
                    value = obj.get_velocity('thetaX');
                case 'thetaDotY'
                    value = obj.get_velocity('thetaY');
            end
        end
        function add_disp_or_vel(obj, type, value)
            switch type
                case 'u'
                    obj.u = value;
                case 'v'
                    obj.v = value;
                case 'w'
                    obj.w = value;
                case {'udot', 'du'}
                    obj.udot = value;
                case {'vdot', 'dv'}
                    obj.vdot = value;
                case {'wdot', 'dw'}
                    obj.wdot = value;
                case 'thetaDotX'
                    obj.thetaDotX = value;
                case 'thetaDotY'
                    obj.thetaDotY = value;
                otherwise
                    warning("Tried to set a displacement or velocity of a non-supported type")
            end
        end
        function value =get_value_by_type(obj, type)
            switch type
                case {'N11', 'N22', 'N12', 'M11', 'M22', 'M12', 'Q1', 'Q2'}
                    value = obj.get_resultant(type);
                case {'u', 'v', 'w', 'du', 'dv', 'dw', 'udot', 'vdot', 'wdot', 'thetaDotX', 'thetaDotY'}
                    value = obj.get_disp_or_vel(type);
            end
        end
        function value = get_displacement(obj, type)
            switch type
                case 'u'
                    value = obj.u;
                case 'v' 
                    value = obj.v;
                case 'w'
                    value = obj.w;
            end
        end
        function value = get_velocity(obj, type)
            switch type
                case 'u'
                    value = obj.udot;
                case 'v' 
                    value = obj.vdot;
                case 'w'
                    value = obj.wdot;
                case 'thetaX'
                    if isempty(obj.thetaDotX)
                        value = NaN;
                    else
                        value = obj.thetaDotX;
                    end
                case 'thetaY'
                    if isempty(obj.thetaDotY)
                        value =NaN;
                    else
                        value = obj.thetaDotY;
                    end
            end
        end
        function set_extrapolated_resultant(obj, type, value)
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
            if isempty(value)
                warning(['resultant ', type, ' has not been calculated'])
            end
        end
        function link_to_an_element(obj, elementNumber)
            obj.associatedElements = unique([obj.associatedElements, elementNumber]);            
        end
    end
end