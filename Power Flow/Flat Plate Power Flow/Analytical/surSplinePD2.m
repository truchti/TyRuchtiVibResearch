%This function takes in a structure of splines and interpolation points and
%returns the derivatives of the splines as well as the
function [dVd, der] = surSplinePD2(splines, points,dimen)
%dSur/dy
if length(splines(:,1)) == 1
    for k = 1:length(splines(1,:))
        der(k) = fnder(splines(1,k));
        dVd(:,k) = fnval(points,der(k));
    end
else
    for i = 1:length(splines(:,1))
        for k = 1:length(splines(1,:))
            der(i,k) = fnder(splines(i,k));
            dVd(:,i,k)= fnval(points,der(i,k));
        end
    end
end
if ndims(dVd == 3)
    if (dimen == 1)
        dVd = permute(dVd,[2,1,3]);
    end
else
    dVd = dVd';
end




