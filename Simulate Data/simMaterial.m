classdef simMaterial
    properties
        E
        poisson
        eta
        density
    end
      methods
        function obj = simMaterial(E, poi, eta, den)
            if nargin < 1
                obj = obj.create_default;
            else
                obj.E = E;
                obj.poisson = poi;
                obj.eta = eta;
                obj.density = den;
            end
        end
        function obj = create_default(obj)
            obj.E = 2.1e11;
            obj.poisson = .28;
            obj.eta = .001;
            obj.density = 7800;
        end
        
    end
end

