classdef ClosedQuinticSplineEvaluator < SplineEvaluator
    properties
         order = 6;
    end
    properties (Access = private)
        N = [ 1/120, -1/24,  1/12, -1/12,  1/24, -1/120;
                0,  1/24,  -1/6,   1/4,  -1/6,   1/24;
                0,  1/12,  -1/6,     0,   1/6,  -1/12;
                0,  1/12,   1/6,  -1/2,   1/6,   1/12;
                0,  1/24,  5/12,     0, -5/12,  -1/24;
                0, 1/120, 13/60, 11/20, 13/60,  1/120];
        D = [0     0     0     0     0     0;
             5     0     0     0     0     0;
             0     4     0     0     0     0;
             0     0     3     0     0     0;
             0     0     0     2     0     0;
             0     0     0     0     1     0];
       
    end
    methods
    end
end