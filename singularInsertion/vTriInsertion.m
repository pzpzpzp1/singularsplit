clear all; close all; 
sig = .00;
area = @(v1,v2,v3) abs(cross([v1-v2, 0],[v3-v2, 0])*[0;0;1])/2;
x0 = [-.05,.65] + randn(1,2)*sig;
x1 = [.95,.07] + randn(1,2)*sig;
x2 = [.05,-.45] + randn(1,2)*sig;
x3 = [.3,0] + randn(1,2)*sig;
e = [.8,-.06] + randn(1,2)*sig; e=e/norm(e); % normalized edge extend direction
% ts = repmat([1,0]',1,5) + randn(2,5)*.2; % target frames
ts = randn(2,5)*.4; % target frames
% ts(:,[4 5]) = abs(randn(2)*2); % target frames
T = [0,3,2;0,4,1;4,1,2;0,3,4;2,3,4]+1; 
load ts.mat; % save('ts.mat','ts')


alphas = linspace(0,.0025,10); alphas = alphas(2:end);
for i=1:numel(alphas)
    alpha = alphas(i);
    
    x4 = x3 + e*alpha; X=[x0;x1;x2;x3;x4];
    Xs{i}=X;
    e03 = x3-x0; e04 = x4-x0;
    e23 = x3-x2; e24 = x4-x2;
    e34 = x4-x3; e41 = x1-x4;
    e01 = x1-x0; e21 = x1-x2;
    areas(1,i)=area(x0,x3,x2);
    areas(2,i)=area(x0,x4,x1);
    areas(3,i)=area(x4,x1,x2);
    areas(4,i)=area(x0,x3,x4);
    areas(5,i)=area(x2,x3,x4);
    bcs(1,:,i)=mean([x0;x3;x2]);
    bcs(2,:,i)=mean([x0;x4;x1]);
    bcs(3,:,i)=mean([x4;x1;x2]);
    bcs(4,:,i)=mean([x0;x3;x4]);
    bcs(5,:,i)=mean([x2;x3;x4]);
    
    cvx_begin
        cvx_precision best
        variable fs(2,5); f1=fs(:,1);f2=fs(:,2);f3=fs(:,3);f4=fs(:,4);f5=fs(:,5);
        dual variables lam1 lam2 lam3 lam4 lam5 lam6
        objval = dot(sum((fs-ts).^2,1), areas(:,i)');
        minimize objval
        subject to
            lam1 : (f1-f4)'*e03'==0
            lam2 : (f4-f2)'*e04'==0
            lam3 : (f1-f5)'*e23'==0
            lam4 : (f5-f3)'*e24'==0
            lam5 : (f4-f5)'*e'==0
            lam6 : (f2-f3)'*e41'==0
    cvx_end
    fs_alpha(:,:,i)=fs;
    objval_alpha(:,i)=objval;
    lams_alpha(:,i)=[lam1';lam2';lam3';lam4';lam5';lam6';];
end

% get lambda dots by linreg on lam(alpha).
tt=polyfit(alphas,lams_alpha([1],:),1); ldot(1)=tt(1);
tt=polyfit(alphas,lams_alpha([2],:),1); ldot(2)=tt(1);
tt=polyfit(alphas,lams_alpha([3],:),1); ldot(3)=tt(1);
tt=polyfit(alphas,lams_alpha([4],:),1); ldot(4)=tt(1);
tt=polyfit(alphas,lams_alpha([5],:),1); ldot(5)=tt(1);
tt=polyfit(alphas,lams_alpha([6],:),1); ldot(6)=tt(1);

%% alpha 0 
X0=[x0;x1;x2;x3]; T0 = [0,3,2; 0,3,1;3,1,2;] + 1;
e03 = x3-x0; 
e23 = x3-x2; 
e31 = x1-x3;
area0s(1,1,1)=area(x0,x3,x2);
area0s(1,1,2)=area(x0,x3,x1);
area0s(1,1,3)=area(x3,x1,x2);
bc0s(1,:,i)=mean([x0;x3;x2]);
bc0s(2,:,i)=mean([x0;x3;x1]);
bc0s(3,:,i)=mean([x3;x1;x2]);
cvx_begin
%     cvx_solver mosek
    cvx_precision best
    variable fs0(2,3); f1=fs0(:,1);f2=fs0(:,2);f3=fs0(:,3);
    dual variables lam1 lam3 lam6
    objval0 = dot(sum((fs0-ts(:,[1 2 3])).^2,1), area0s(:)');
    minimize objval0
    subject to
        lam1 : (f1-f2)'*e03'==0
        lam3 : (f1-f3)'*e23'==0
        lam6 : (f2-f3)'*e31'==0
cvx_end
lam0 = [lam1 lam1 lam3 lam3 0 lam6]; % augmented dual vars
% 2*area0s(1)*(fs0(:,1)-ts(:,1)) - lam1*e03' - lam3*e23'
% 2*area0s(2)*(fs0(:,2)-ts(:,2)) + lam1*e03' - lam6*e31'
% 2*area0s(3)*(fs0(:,3)-ts(:,3)) + lam3*e23' + lam6*e31'

%% viz
figure; title('primal obj'); plot(alphas,objval_alpha,'.-'); xline(0); yline(objval0,'green');
% figure; plot(alphas,reshape(permute(lams_alpha,[2 1 3]),[],numel(alphas))','.-'); legend; xline(0); yline(0);
% figure; plot((reshape(permute(lams_alpha,[2 1 3]),[],numel(alphas))./alphas)','.-'); legend;

figure; cla; hold all; title('dual abs');
plot(alphas,reshape(permute(lams_alpha([1],:),[2 1 3]),[],numel(alphas))','r.-','color','r');  xline(0); yline(0);
plot(alphas,reshape(permute(lams_alpha([2],:),[2 1 3]),[],numel(alphas))','.-','color',[1,.5,.5]); xline(0); yline(0);
plot(alphas,reshape(permute(lams_alpha([3],:),[2 1 3]),[],numel(alphas))','g.-');  xline(0); yline(0);
plot(alphas,reshape(permute(lams_alpha([4],:),[2 1 3]),[],numel(alphas))','.-','color',[.5,1,.5]); xline(0); yline(0);
plot(alphas,reshape(permute(lams_alpha([5],:),[2 1 3]),[],numel(alphas))','k.-');  xline(0); yline(0);
plot(alphas,reshape(permute(lams_alpha([6],:),[2 1 3]),[],numel(alphas))','b.-');  xline(0); yline(0);
yline(lam1(1),'r');
yline(lam3(1),'g');
yline(lam6(1),'b');

figure; cla; hold all; title('dual rel');
plot(alphas,reshape(permute(lams_alpha([1],:)-lam1(1),[2 1 3]),[],numel(alphas))','r.-','color','r');  xline(0); yline(0);
plot(alphas,reshape(permute(lams_alpha([2],:)-lam1(1),[2 1 3]),[],numel(alphas))','.-','color',[1,.5,.5]); xline(0); yline(0);
plot(alphas,reshape(permute(lams_alpha([3],:)-lam3(1),[2 1 3]),[],numel(alphas))','g.-');  xline(0); yline(0);
plot(alphas,reshape(permute(lams_alpha([4],:)-lam3(1),[2 1 3]),[],numel(alphas))','.-','color',[.5,1,.5]); xline(0); yline(0);
plot(alphas,reshape(permute(lams_alpha([5],:),[2 1 3]),[],numel(alphas))','k.-');  xline(0); yline(0);
plot(alphas,reshape(permute(lams_alpha([6],:)-lam6(1),[2 1 3]),[],numel(alphas))','b.-');  xline(0); yline(0);
% yline(lam1(1),'r');
% yline(lam3(1),'g');
% yline(lam6(1),'b');


figure; cla; axis equal; hold all;
for i=1:numel(alphas)
    alpha = alphas(i);
    
    cla; patch('faces',T,'vertices',Xs{i},'facecolor','green')
    quiver(bcs(:,1,i),bcs(:,2,i), reshape(fs_alpha(1,:,i),[],1), reshape(fs_alpha(2,:,i),[],1), 'r')
    quiver(bcs(:,1,i),bcs(:,2,i), reshape(ts(1,:),[],1), reshape(ts(2,:),[],1), 'b')
    drawnow; pause(.01)
end

%% env thm
n=3; ab = polyfit(alphas(1:n), objval_alpha(1:n),1); 
fdiff1 = ab(1);
fdiff2 = (objval_alpha(2)-objval_alpha(1))/(alphas(2)-alphas(1));
fdiff3 = (objval_alpha(1)-objval0)/(alphas(1));
dualpart0 = lam1'*(fs_alpha(:,4,1)-fs_alpha(:,2,1))'*e' + ...
            lam3'*(fs_alpha(:,5,1)-fs_alpha(:,3,1))'*e';
dualpart = lams_alpha(1,1)*(fs_alpha(:,4,1)-fs_alpha(:,2,1))'*e' + ...
    lams_alpha(3,1)*(fs_alpha(:,5,1)-fs_alpha(:,3,1))'*e';
l2345 = sum(reshape(fs_alpha(:,[2 3 4 5],1)-ts(:,[2 3 4 5]),2,[]).^2,1);
dAda = [cross([e01,0],[e,0]); cross([e21,0],[e,0]); ...
        cross([e03,0],[e,0]); cross([e23,0],[e,0]); ]/2;
dAda = abs(dAda(:,3)).*[-1,-1,1,1]';
primalpart = l2345*dAda;
adiff = primalpart-dualpart; % dual variables are on LHS of : so equality constraints are subtracted from obj in lagrangian
[fdiff1 fdiff2 fdiff3 adiff primalpart dualpart dualpart0]

%% test adjoint condition on f45
z=[0,0]; 
C = [e03,z,z,-e03,z;...
    z,-e04,z,e04,z;...
    e23,z,z,z,-e23;...
    z,z,-e24,z,e24;...
    z,z,z,e,-e;...
    z,e41,-e41,z,z;];
C0 = [e03,z,z,-e03,z;...
    z,-e03,z,e03,z;...
    e23,z,z,z,-e23;...
    z,z,-e23,z,e23;...
    z,z,z,e,-e;...
    z,e31,-e31,z,z;];
D = [z,z,z,z,z;...
    z,-e,z,e,z;...
    z,z,z,z,z;...
    z,z,-e,z,e;...
    z,z,z,z,z;...
    z,-e,e,z,z;];
f1_on_e03 = dot(fs0(:,1),e03);
f1_on_e23 = dot(fs0(:,1),e23);
x45e = reshape(fs_alpha(:,[4 5],1),[],1);
A=[e03 0 0; 0 0 e23; e -e];
B = [f1_on_e03; f1_on_e23; 0];
A*x45e - B;
cvx_begin
    variable x45p(4,1)
    minimize norm(x45p-x45e)
    subject to
        A*x45p==B
cvx_end
[x45p x45e]
augfs0 = [fs0, reshape(x45p,2,2)];
augarea0 = [area0s(:);0;0];
C0*augfs0(:)
kkt = 2*reshape(augarea0'.*(augfs0-ts),[],1)' - lam0*C0;
augdAda=[0;dAda];
BB1 = [2*augdAda(4)*(augfs0(:,4)-ts(:,4)); 2*augdAda(5)*(augfs0(:,5)-ts(:,5))]; %fterm
BB2 = -[lam0(2)*e'; lam0(4)*e']; % lterm
BB12 = [BB1, BB2]; BB = -(BB1+BB2);
AA=[e03',z',e';z',e23',-e'];
lambda_dot_12_34_5 = AA\BB;
[AA*lambda_dot_12_34_5 BB]
[lambda_dot_12_34_5, [ldot(2)-ldot(1);ldot(4)-ldot(3);ldot(5)]]

%% derive f45 using adjoint condition and integrability

cvx_begin
    variable f45(4,1)
    variable ldots(3,1)
    
    augdAda=[0;dAda];
    BB1 = [2*augdAda(4)*(f45(1:2)-ts(:,4)); 2*augdAda(5)*(f45(3:4)-ts(:,5))]; %fterm
    BB2 = -[lam0(2)*e'; lam0(4)*e']; % lterm
    BB = -(BB1+BB2);
    AA=[e03',z',e';z',e23',-e'];
%     lambda_dot_12_34_5 = AA\BB;
    minimize norm(f45-randn(4,1)) % arbitrary convex obj. system is fully determined anyway.
    subject to
        C0*[fs0(:); f45]==0
        AA*ldots == BB

cvx_end
%% THEY MATCH!!! 
[x45e, f45]






