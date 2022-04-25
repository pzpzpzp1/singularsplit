% [c,ceq,gradc,gradceq] = ellipseparabola(x)
function [y,ceq,dydx,gradceq]=jacdets(x)

    func=@jacdet;
    Js = permute(reshape(x,2,[],2),[1 3 2]);
    dlJs = dlarray(Js);
    [energy, energyGrad] = dlfeval(func,dlJs);
    
    xr = reshape(x,2,[],2); nT = size(xr,2);
    ii = repelem([1:nT],1,4);
    jj=reshape(permute(reshape(1:numel(x),2,nT,2),[1 3 2]),[],1);
    kk = reshape(extractdata(energyGrad),[],1);
    dydx = sparse(ii,jj,kk)';

    y = extractdata(energy);
    gradceq=[];
    ceq=[];
end

% NEGATIVE jacdet
function [y,dydx]=jacdet(J)
    a = J(1,1,:);
    b = J(1,2,:);
    c = J(2,1,:);
    d = J(2,2,:);
    epsilon = .001;
    y = epsilon - reshape(a.*d-b.*c,[],1); % bound jacdet above epsilon
    yout = sum(y);
    dydx = dlgradient(yout,J);
end

function verify
N=randi(100)*4;
x = randn(N,1);
dx = randn(N,1);
eps = 1e-3;
[y,ceq,dydx,gradceq]=jacdets(x);
yp=jacdets(x+dx*eps);
ym=jacdets(x-dx*eps);
fdiff = (yp-ym)/(2*eps);
adiff = dydx'*dx;
norm(fdiff-adiff)

end