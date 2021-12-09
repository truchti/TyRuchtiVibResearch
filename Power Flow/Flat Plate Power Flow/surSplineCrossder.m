function ddxy = surSplineCrossder(splines, Xq, Yq)

xf = Xq(1,:); % x evaluation points
yf = Yq(:,1)'; % y evaluation points
der = splines; %splines that hold y constant
if length(splines(:,1)) == 1
    dVd = zeros(size(Xq));
    ddxy = zeros(size(Xq));
    %for each spline in the 
    for k = 1:length(splines)
        %take the derivative of the spline wrt x
        der(1,k) = fnder(splines(1,k));
        % evaluate derivative at the x points
        dVd(k,:) = fnval(xf,der(k));
    end
    % create a set of splines where x is constant
    for i = length(xf):-1:1
        altsplines(i) = spapi(5,Yq(:,i),dVd(:,i));
        altder(i) = fnder(altsplines(1,i));
        ddxy(:,i) = fnval(yf,altder(i));
    end
    
else
    tmax = length(splines(1,:));
    dVd = zeros([size(Xq) tmax]);
    ddxy = zeros([size(Xq) tmax]);
    for i = 1:length(splines(:,1))
        for k = 1:length(splines(1,:))
            der(i,k) = fnder(splines(i,k));
            dVd(i,:,k)= fnval(xf,der(i,k));
        end
    end
    for k = 1:tmax
        for i = 1:length(xf)
            altsplines(i,k) = spapi(5,Yq(:,i),dVd(:,i,k));
            ddxy(:,i,k) = fnval(yf,altsplines(i));
        end
    end
end
