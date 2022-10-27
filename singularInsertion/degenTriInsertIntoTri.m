clear all; close all; 
sig = .1;
area = @(v1,v2,v3) abs(cross([v1-v2, 0],[v3-v2, 0])*[0;0;1])/2;
x0 = [1,1] + (rand(1,2)-.5)*sig;
x1 = [-.5,0] + (rand(1,2)-.5)*sig;
x2 = [1,-1] + (rand(1,2)-.5)*sig;
x3 = [.25,0] + (rand(1,2)-.5)*sig;
e = [.2,0] + (rand(1,2)-.5)*sig; e=e/norm(e); % normalized edge extend direction
ts = randn(2,5)*.4; % target frames
T = [0,2,4; 3,2,1; 3,1,0; 3,4,0; 3,4,2;]+1; 

%% alpha 0 
X0=[x0;x1;x2;x3]; T0 = [0,3,2; 3,1,2;  0,3,1] + 1;
e03 = x3-x0; 
e23 = x3-x2; 
e31 = x1-x3;
for i=1:3; area0s(i)=area(X0(T0(i,1),:),X0(T0(i,2),:),X0(T0(i,3),:)); end;
bc0s = (X0(T0(:,1),:)+X0(T0(:,2),:)+X0(T0(:,3),:))/3;
cvx_begin
%     cvx_solver mosek
    cvx_precision best
    variable fs0(2,3); f1=fs0(:,1);f2=fs0(:,2);f3=fs0(:,3);
    dual variables lam1 lam3 lam6
    objval0 = dot(sum((fs0-ts(:,[1 2 3])).^2,1), area0s(:)');
    minimize objval0
    subject to
        lam1 : (f2-f1)'*e03'==0
        lam3 : (f3-f1)'*e23'==0
        lam6 : (f2-f3)'*e31'==0
cvx_end
lam0 = [lam1 lam1 lam3 lam3 0 lam6]; % augmented dual vars

%% solve for missing frames f4 f5
z=[0,0]; 
e04=e03; e24=e23; e02 = x2-x0;
C0 = [z,e03,z,-e03,z;...
    -e04,z,z,e04,z;...
    z,z,e23,z,-e23;...
    -e24,z,z,z,e24;...
    z,z,z,e,-e;...
    z,e31,-e31,z,z;];
D = [z,z,z,z,z;...
    -e,z,z,e,z;...
    z,z,z,z,z;...
    -e,z,z,z,e;...
    z,z,z,z,z;...
    z,z,z,z,z;];
dAda = [cross([e02,0],[e,0]); 0 0 0; 0 0 0; cross([e03,0],[e,0]); cross([e23,0],[e,0]);]/2; dAda = abs(dAda(:,3)).*[-1,0,0,1,1]';

cvx_begin
    cvx_precision high
    variable f45(4,1)
    fs = [fs0(:);f45];    
    variable ldots(1,3)
    BB1 = [2*dAda(4)*(f45(1:2)-ts(:,4)); 2*dAda(5)*(f45(3:4)-ts(:,5))]; %fterm
    BB2 = -(lam0*D(:,7:10))'; % lterm
    BB = -(BB1+BB2);
    minimize norm(f45-randn(4,1)) % arbitrary convex obj. system is fully determined anyway.
    subject to
        C0*fs == 0
        ldots*C0([1 3 5],7:10) == BB'
cvx_end
fs0Aug = [fs0, reshape(f45,2,2)];

dualpart = lam0*D*fs(:);
l2345 = sum(reshape(fs-ts(:),2,[]).^2,1);
primalpart = l2345*dAda;
adiff = primalpart-dualpart; % dual variables are on LHS of : so equality constraints are subtracted from obj in lagrangian

%% get fdiff
alphas = linspace(0,.0025,10); alphas = alphas(2:end);
for i=1:numel(alphas)
    alpha = alphas(i);
    
    x4 = x3 + e*alpha; X=[x0;x1;x2;x3;x4];
    Xs{i}=X;
    e03 = x3-x0; e04 = x4-x0;
    e23 = x3-x2; e24 = x4-x2;
    e34 = x4-x3; e31 = x1-x3;
    e01 = x1-x0; e21 = x1-x2;
    for j=1:5; areas(j,i)=area(X(T(j,1),:),X(T(j,2),:),X(T(j,3),:)); end;
    bcs(:,:,i) = (X(T(:,1),:)+X(T(:,2),:)+X(T(:,3),:))/3;

    cvx_begin
        cvx_precision best
        variable fs(2,5); f1=fs(:,1);f2=fs(:,2);f3=fs(:,3);f4=fs(:,4);f5=fs(:,5);
        dual variables lam1 lam2 lam3 lam4 lam5 lam6
        objval = dot(sum((fs-ts).^2,1), areas(:,i)');
        minimize objval
        subject to
            lam1 : (f2-f4)'*e03'==0
            lam2 : (f4-f1)'*e04'==0
            lam3 : (f3-f5)'*e23'==0
            lam4 : (f5-f1)'*e24'==0
            lam5 : (f4-f5)'*e'==0
            lam6 : (f2-f3)'*e31'==0
    cvx_end
    fs_alpha(:,:,i)=fs;
    objval_alpha(:,i)=objval;
    lams_alpha(:,i)=[lam1';lam2';lam3';lam4';lam5';lam6';];
end


%% viz
figure; title('primal obj'); plot(alphas,objval_alpha,'.-'); xline(0); yline(objval0,'green');
ab = polyfit(alphas,objval_alpha,1); fdiff = ab(1);

[adiff fdiff primalpart dualpart]
[fs0Aug(:)'; reshape(fs_alpha(:,:,1),[],1)']
