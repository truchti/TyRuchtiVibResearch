classdef ChebyshevPolyFit < handle
    properties
        fitPoly
    end
    properties (Hidden = true)
        order
        data_x
        values
        fitCoeffs
    end
    properties (Access = private)
        dummyPoly
    end
    methods
        function obj = ChebyshevPolyFit(order, data_x, values)
            if nargin > 1
                obj.order = order;
                obj.data_x = data_x;
                obj.values = values;
            end
            rng = [min(data_x), max(data_x)];
            obj.dummyPoly = ChebyshevPolynomial(zeros(order+1,1), rng);
            obj.find_coefficients();
        end
        function chebPoly = output_cheb_poly_obj(obj)
            chebPoly = ChebyshevPolynomial(obj.fitCoeffs,  [min(obj.data_x), max(obj.data_x)]);
        end
    end
    methods %(Access = private)
        function find_coefficients(obj)
            [~, T] = obj.dummyPoly.eval(obj.data_x);
            obj.fitCoeffs = T\obj.values;
            obj.fitPoly = ChebyshevPolynomial(obj.fitCoeffs, obj.dummyPoly.inputRange);
        end
        function plot_data(obj)
            scatter(obj.data_x, obj.values);
        end
        function plot_fit(obj,points)
            if nargin < 2
                points = linspace(min(obj.data_x), max(obj.data_x), 1000);
            end
            fx = obj.fitPoly.eval(points);
            plot(points,fx)
        end
        function plot_data_and_fit(obj)
            obj.plot_data()
            hold on;
            obj.plot_fit()
        end
    end
end

            