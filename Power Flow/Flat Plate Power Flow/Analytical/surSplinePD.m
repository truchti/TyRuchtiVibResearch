%This function takes in a structure of splines and interpolation points and
%returns the derivatives of the splines as well as the 
function [dVd, der] = surSplinePD(splines, points,dimen)
%dSur/dy

for i = 1:length(splines(:,1))
    for k = 1:length(splines(1,:))
        der(i,k) = fnder(splines(i,k));
        dVd(:,i,k)= fnval(points,der(i,k));
    end
end
if (dimen ==1)
    if ndims(dVd == 3)
        dVd = permute(dVd,[2,1,3]);
    end
else
    dVd = dVd'
end

    

