classdef PlatePowerFlowDataPlotter
    properties
        xi
        eta
        qx
        qy
    end
    methods
        function obj = PlatePowerFlowDataPlotter(xi,eta,qx,qy)
            obj.xi = xi;
            obj.eta = eta;
            obj.qx = qx;
            obj.qy = qy;
        end
        function plot_power_flow(obj, radius, forceLocs)
            [Z, Theta] = meshgrid(obj.eta, obj.xi);
            if radius ~= 0
                Theta = radius*Theta;
            end
            quiver(Theta, Z, obj.qx, obj.ql,'color', [1 0 .2]);
            if nargin>2 && ~isempty(forceLocs)
                hold on;
                forceLocs(:,1) = wrapTo2Pi(forceLocs(:,1));
                
                scatter(forceLocs(:,1), forceLocs(:,2),'r');
                hold off;
            end
        end
    end
end