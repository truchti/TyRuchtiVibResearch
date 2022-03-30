classdef ChebyshevSurface
    properties
        coeffNet
        Xrange
        Yrange
    end
    properties (Dependent = true)
        orderX
        orderY
    end
    properties (Hidden = true)
        XbasisEvaluator
        YbasisEvaluator
        mapXInput
        mapYInput
    end
    methods
        function obj = ChebyshevSurface(coeffsNet, rangeX, rangeY)
            obj.coeffNet = coeffsNet;
            obj.Xrange = rangeX;
            obj.Yrange = rangeY;
            obj = obj.create_mapping_functions();
            obj.XbasisEvaluator = ChebyshevBasisEvaluator(obj.orderX);
            obj.YbasisEvaluator = ChebyshevBasisEvaluator(obj.orderY);
        end
        function [fx, Txim, Tetan] = eval(obj, xpts, ypts)
            xis = obj.mapXInput(xpts);
            Txim = obj.XbasisEvaluator.eval(xis);
            etas = obj.mapYInput(ypts);
            Tetan = obj.YbasisEvaluator.eval(etas);
             fx = Txim*obj.coeffNet*Tetan';
             
        end
        function dfdx = deriv_eval(obj, xpts, ypts, deriv)
            xis = obj.mapXInput(xpts);
            dTxim = obj.XbasisEvaluator.deriv_eval(xis, deriv(1));
            dxidx = mean(diff(xis)./diff(xpts)); % assumes linear mapping
            etas = obj.mapYInput(ypts);
            dTetan = obj.YbasisEvaluator.deriv_eval(etas, deriv(2));
            detady = mean(diff(etas)./diff(ypts)); % assumes linear mapping
            dfdxieta = dTxim *obj.coeffNet * dTetan';
            dfdx = dfdxieta * (dxidx^(deriv(1)) * detady^(deriv(2)));
        end
        function plot(obj, xpts, ypts)
            fxy = obj.eval(xpts, ypts);
            [X,Y] = ndgrid(xpts, ypts);
            surf(X, Y, fxy)
        end
        function plot_deriv(obj, xpts, ypts, deriv)
            dfxy = obj.deriv_eval(xpts, ypts, deriv);
            [X,Y] = ndgrid(xpts, ypts);
            surf(X, Y, dfxy)
        end
        function value = get.orderX(obj)
            value = size(obj.coeffNet,1)-1;
        end
        function value = get.orderY(obj)
            value = size(obj.coeffNet,2)-1;
        end
        function obj = create_mapping_functions(obj)
            ax = min(obj.Xrange);
            bx = max(obj.Xrange);
            ay = min(obj.Yrange);
            by = max(obj.Yrange);
            obj.mapXInput = @(x) ((x-ax)/(bx-ax)-.5)*2;
            obj.mapYInput = @(y) ((y-ay)/(by-ay)-.5)*2;
        end
    end
end
