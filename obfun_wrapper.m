function [y,dydx]=obfun_wrapper(x,triAreas,func)
    if nargin < 3
        func = @symmetricDirichlet;
    end

    mid = (size(x,1)/2);
    nT = (size(x,1)/4);
    Js = zeros(2,2,nT);
    Js(:,1,:) = reshape(x(1:mid),2,1,[]);
    Js(:,2,:) = reshape(x((mid+1):end),2,1,[]);
    dlJs = dlarray(Js);
    [energy, energyGrad] = dlfeval(func,dlJs,triAreas);
    dydx = extractdata(reshape(permute(energyGrad,[1 3 2]),[],1));
    y = extractdata(energy);
end