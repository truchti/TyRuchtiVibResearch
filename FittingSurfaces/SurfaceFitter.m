classdef SurfaceFitter < handle
    % this class find the surface with the least squares error using the
    % given spline evaluators. By default this is a Quintic B-spline tensor
    % mesh. xi and eta parameters should be N x 1 vectors where xiParams(i) and
    % etaParam(i) represent the coordinates of the ith point
    % value should be an N by 1 vector with the surface value at the ith
    % coordinate pair. values can also be an N by M matrix in this case a
    % separate surface will be created for each column of values.
    properties
        xiParams
        etaParams
        values
    end
    methods 
        function obj = SurfaceFitter(xiParam, etaParam, values)
        end
        function fit_surface(obj)
        end
    end
end
