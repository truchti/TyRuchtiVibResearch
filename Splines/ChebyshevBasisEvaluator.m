classdef ChebyshevBasisEvaluator < handle
    %% This class evaluates the quintic Spline Basis functions for a given knot Vector and parameter
    properties
        order
    end
    properties (Access = private)
        currentXi
        Txi
    end
    methods
        function obj = ChebyshevBasisEvaluator(order)
            obj.order = order;
        end
        function fx = eval(obj,xi)
            fx = obj.evaluate_T_basis_at_xi(xi);
        end
        function dfx = deriv_eval(obj,xi, deriv)
            obj.set_current_xi(xi);
            if deriv == 0
                dfx = obj.Txi;
            else
                dfx = obj.derivative_of_all_basis_functions_at_current_xi(deriv);
            end
        end
        function plot_cheb_polys(obj,points)
            if nargin < 2
                points = linspace(-1, 1, 1000);
            end
            fx = obj.evaluate_T_basis_at_xi(points);
            plot(points,fx)
        end
        function plot_cheb_n_and_derivatives(obj,n, pts)
            if nargin < 3
                pts = linspace(-1, 1, 1000);
            end
            f = obj.eval(pts);
            df = obj.deriv_eval(pts,1);
            d2f = obj.deriv_eval(pts,2);
            d3f = obj.deriv_eval(pts,3);
            
            plot(pts, f(:,n+1), 'k', pts, df(:,n+1), 'r', pts, d2f(:,n+1), 'b', pts, d3f(:,n+1), 'g')
            legend('Tn', 'dTn', 'd^2Tn', 'd^3Tn')
        end
    end
    methods (Hidden = true)
        function T = evaluate_T_basis_at_xi(obj, xi)
            xi = reshape(xi, numel(xi), 1);
            T = nan(length(xi), obj.order+1);
            T(:,1) = ones(size(xi));
            T(:,2) = xi;
            for i = 3:obj.order+1
                T(:,i) = 2*xi.*T(:,i-1)-T(:,i-2);
            end
        end
        function set_current_xi(obj, xi)
            assert(min(xi) >= -1); assert(max(xi) <= 1);
            obj.Txi = obj.evaluate_T_basis_at_xi(xi);
            obj.currentXi = xi;
        end
        function dTn = derivative_of_all_basis_functions_at_current_xi(obj, deriv)
            p = deriv;
            dTn = zeros(size(obj.Txi));
            for n= 0:obj.order
                i = n+1; % Tn is zero based matlab is 1 based so index equals n +1
                dT = zeros(size(obj.Txi(:,i)));
                for k = 0:(n-p)
                    if mod(n-p,2) == mod(k,2)
                        if k == 0
                            nCk = nchoosek((n+p-k)/2-1, (n-p-k)/2);
                            fq = a_fact_o_b_fact((n+p+k)/2-1, (n-p+k)/2);
                            dT = 1/2*(nCk*fq* obj.Txi(:,k+1));
                        else
                            nCk = nchoosek((n+p-k)/2-1, (n-p-k)/2);
                            fq = a_fact_o_b_fact((n+p+k)/2-1, (n-p+k)/2);
                            dT = dT+(nCk*fq* obj.Txi(:,k+1));
                        end
                    end
                end
                dTn(:,i) = 2^p*n*dT;
            end
        end
    end
end
function val = a_fact_o_b_fact(a,b)    
        val = factorial(a)/ factorial(b);
end