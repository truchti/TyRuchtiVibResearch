classdef SurfacePlotter < handle
    properties
        Surface
        radius
        realPlot = true;
    end
    properties (Hidden = true)
        Ths
        Zs
        Values

    end
    methods
        function obj = SurfacePlotter(Surf, Ths, Zs, rad)
            if nargin < 4
                rad = Inf;
            end
            obj.Surface = Surf;
            obj.radius = rad;
            if nargin < 2
                obj.Ths = linspace(0, 2*pi, 87);
                obj.Zs = linspace(Surf.heightRange(1), Surf.heightRange(2), 93);
            else
                obj.Ths = Ths;
                obj.Zs = Zs;
            end
                obj.get_plot_values();

        end
        function plot(obj, Ths, Zs)
            if nargin > 2 
                obj.get_plot_values(Ths, Zs);
            end
            [Thetas, Heights] = ndgrid(obj.Ths, obj.Zs);
            if nargin < 4
                useReal = obj.realPlot;
            end
            if useReal 
                values = real(obj.Values);
            else
                values = imag(obj.Values);
            end
            s = surf(Thetas, Heights, values);
            s.EdgeAlpha = .05;
        end
        function plot_as_cylinder(obj, Ths, Zs)
            if isinf(obj.radius)
                obj.plot(Ths,Zs);% if no radius then don't plot as a cylinder;
                return
            end
            [Thetas, Heights] = ndgrid(obj.Ths, obj.Zs);
            Xs = obj.radius*cos(Thetas);
            Ys = obj.radius*sin(Thetas);
            s =surf(Xs, Ys, Heights, obj.Values); 
            s.EdgeAlpha = 0;
        end
        function values = get_plot_values(obj, Ths, Zs)
            if nargin <2
                Ths = obj.Ths;
                Zs = obj.Zs;
            end
            obj.Ths = Ths;            obj.Zs = Zs;
            obj.Values = obj.Surface.evaluate(Ths,Zs);
        end
    end
%     methods (Static)
%         function [Ths, Zs] = grid_and_vect_parameters(Ths, Zs)
%             [T, Z] = ndgrid(Ths, Zs);
%             Ths = reshape(T, numel(T),1);
%             Zs  = reshape(Z, numel(T),1);
%         end
%     end
end