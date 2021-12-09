classdef PowerFlowDataPlotter
    properties
        xi
        eta
        qt
        ql
    end
    methods
        function obj = PowerFlowDataPlotter(xi,eta,qt,ql)
            obj.xi = xi;
            obj.eta = eta;
            obj.qt = qt;
            obj.ql = ql;
        end
        function plot_flat(obj, radius, forceLocs)
            [Z, Theta] = meshgrid(obj.eta, obj.xi);
            if radius ~= 0
                    Theta = radius*Theta;
            end
            quiver(Theta, Z, obj.qt, obj.ql,'color', [1 0 .2]);
            if nargin>2 && ~isempty(forceLocs)
                hold on;
                forceLocs(:,1) = wrapTo2Pi(forceLocs(:,1));
                
                scatter(forceLocs(:,1), forceLocs(:,2),'r');
                hold off;
            end
        end
        function plot_as_cylinder(obj, radius, height)
            [Z, Theta] = meshgrid(obj.eta, obj.xi);
            X = radius*cos(Theta);
            Y = radius*sin(Theta);
            qz = obj.ql;
            qx = -radius*obj.qt.*sin(Theta);
            qy = radius*obj.qt.*cos(Theta);
            [cx, cy,cz] = cylinder(radius, 40);
            cz = cz*height;
            s = surf(cx,cy,cz);
            s.EdgeColor = 'none';
            s.FaceAlpha = 0;
            hold on;
            quiver3(X,Y,Z,qx,qy,qz,'color', [1,0,0]);
            axis equal
            hold off;
        end
    end
end