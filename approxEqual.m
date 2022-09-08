function tf = approxEqual(a, b, tol)
    if nargin < 2
        warning("requires 2 inputs");
        tf = nan;
        return
    end
    if nargin < 3
        tol = 1e-6;
    end
    [ra, ca] = size(a);
    [rb, cb] = size(b);
    if rb == 1 && cb == 1
        b = b*ones(ra, ca);
        [rb, cb] = size(b);
    end
    if (ra ~= rb || ca ~= cb)
        warning("Sizes don't match");
        tf = nan;
        return
    end
    tol = ones(rb,cb)*tol;
    tf = (a<= b+tol)&(a>= b-tol);
end 