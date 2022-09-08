classdef PowerFlowResultantsCylindrical < handle
    properties
        radDisp
        thetaDisp
        zDisp
        radVel
        thetaVel
        zVel
        cylProps PFCylProperties
    end
    properties (Hidden = true)
        D;  K;
        Nt; Nl; Ntl; Nlt;
        Mt; Ml; Mtl; Mlt;
        Qt; Ql;
    end
    methods 
        function obj = PowerFlowResultantsCylindrical(Surfaces, cylinderProperties)
            %Surfaces is a 12 by 1 array of chebPoly Surfaces 
            %cylinderProperties describe cylinder geometery and material
            obj.radDisp = Surfaces([1,4]);
            obj.thetaDisp = Surfaces([2,5]);
            obj.zDisp = Surfaces([3,6]);
            obj.radVel = Surface([7,10]);
            obj.thetaVel = Surfaces([8, 11]);
            obj.zVel = Surfaces([9, 12]);
            obj.cylProps = cylinderProperties;
            obj.calculate_coefficients();
        end
    end
    methods (Hidden = true)
        function calculate_coefficients(obj)
            mu = obj.cylProps.poissons;
            % D = E*h/(1-v^2)
            obj.D = obj.cylProps.E*obj.cylProps.thickness/(1-mu^2);
            % K = Eh^3/(12*(1-v^2))
            obj.K = obj.cylProps.E*obj.cylProps.thickness^3/(12*(1-mu^2));
        end
    end
end