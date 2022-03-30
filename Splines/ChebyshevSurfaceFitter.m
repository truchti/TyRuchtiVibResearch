classdef ChebyshevSurfaceFitter < handle
    properties
        fitSurface
    end
    properties (Access = private)
        data_x
        data_y
        values
        order
    end
    properties (Dependent = true) 
        ChebOrderX
        ChebOrderY
    end
    properties (Hidden = true, Dependent = true)
        rngX
        rngY
    end
    properties (Hidden = true, Access = private)
        coeffNet
        dummySurf
    end
    methods
        function obj = ChebyshevSurfaceFitter(order, data_x, data_y, values)
            if length(order) == 2
                obj.order = order;
            end
            if length(order) ==1
                obj.order = [order, order];
            end
            if nargin> 2
                obj.data_x = data_x;
                obj.data_y = data_y;
                obj.values = values;
            end
            obj.dummySurf = ChebyshevSurface(zeros(obj.order(1)+1, obj.order(2)+1),obj.rngX, obj.rngY);
            obj.fit_data();
            obj.create_fit_surface_object();
        end
        function fit_data(obj)
            [~, Txim, Tetan] = obj.dummySurf.eval(obj.data_x,obj.data_y);
            obj.coeffNet = numel(obj.values)*(Txim\diag(obj.values)/(Tetan'));            
        end
        function srf = create_fit_surface_object(obj)
            obj.fitSurface = ChebyshevSurface(obj.coeffNet, obj.rngX, obj.rngY);
            if nargout >0
                srf = obj.fitSurface();
            end
        end
        function plot_data_and_fit(obj)
            scatter3(obj.data_x, obj.data_y, obj.values);
            hold on;
            obj.plot_fit();
        end
        function plot_fit(obj)
            x =linspace(obj.rngX(1), obj.rngX(2), 37);
            y = linspace(obj.rngY(1), obj.rngY(2), 39);
%             [X,Y] = ndgrid(x,y);
            obj.fitSurface.plot(x,y);            
        end
        function value = get.rngX(obj)
            value = [min(obj.data_x), max(obj.data_x)];
        end
        function value = get.rngY(obj)
            value = [min(obj.data_y), max(obj.data_y)];
        end
        function value = get.ChebOrderX(obj)
            value = obj.order(1);
        end
        function value = get.ChebOrderY(obj)
            value = obj.order(2);
        end
    end
end