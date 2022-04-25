
clear all; close all; 

%% load data
meshname = 'data/polar_wedge_no_singularities.obj';

[V,T,UV,TF,N,NF] = readOBJ(meshname);
[u,v] = loadFrameText('data/frames.txt');
x = [reshape(u',[],1); reshape(v',[],1)];
data0 = getMeshData(V,T); data=data0; V0=V; T0=T;
intedgeind = randi(sum(~data.isBoundaryEdge));
% intedgeind = 1079; % wroks ok. depends a lot on primal convergence.
% intedgeind = 6587; % seems to work. for specific ep
% intedgeind = 3044; % works well
% intedgeind = 5360; % diagonal. works well
% intedgeind = 7948; % close to boundary. doesn't work well.
intedgeind = 6304; % circumcentric orientation. doesn't work well/
optol = 1e-7;

v1 = data.edges(intedgeind,1);
v2 = data.edges(intedgeind,2);
pdir = V0(v2,1:2)-V0(v1,1:2);
triinds = data.edges2triangles(intedgeind,:);
t1 = data0.triangles(triinds(1),:);
if numel(strfind([t1 t1(1)],[v1 v2]))~=0
    highertriind = triinds(1); 
    lowertriind = triinds(2); 
else
    highertriind = triinds(2);
    lowertriind = triinds(1); 
end

figure; hold all; axis equal; axis tight manual; title('neighborhood');
TV = data.triangleBarycenters;
patch('faces',T,'Vertices',V,'facecolor','blue','edgecolor','black');
patch('faces',T(highertriind,:),'Vertices',V,'facecolor','green','edgecolor','black','facealpha',.5);
patch('faces',T(lowertriind,:),'Vertices',V,'facecolor','red','edgecolor','black','facealpha',.5);
patch('vertices',V,'faces',data.edges(intedgeind,[1 2 1]),'edgecolor','cyan','linewidth',2);
scatter(TV(:,1),TV(:,2),'k','filled')
scatter(V0(v1,1),V0(v1,2),'g','filled')
scatter(V0(v2,1),V0(v2,2),'r','filled')
sc = 2;
bl = mean(V0(data.edges(intedgeind,:),:)) - [.1 .1 0]*sc;
tr = mean(V0(data.edges(intedgeind,:),:)) + [.1 .1 0]*sc;
xlim([bl(1) tr(1)]); ylim([bl(2) tr(2)])

% perturbfactors = [.003 .001 .0005 .0003 .0001 .00005]; 
perturbfactors = [.003 .001 .0005 .0003 .0001]; 
N = numel(perturbfactors);
for i=1:N
    perturbfactor=perturbfactors(i);
    %% select perturb edge. displace vertex
    V=V0;
    V(v2,1:2) = V(v1,1:2) + pdir * perturbfactor;
    data = getMeshData(V,T);

    %% verify integrability
    intEdgeInds = find(~data.isBoundaryEdge);
    tri2tri = data.edges2triangles(~data.isBoundaryEdge,:);
    intedgevectors = data.vertices(data.edges(~data.isBoundaryEdge,1),1:2)-data.vertices(data.edges(~data.isBoundaryEdge,2),1:2);
    
    nie = size(tri2tri,1);
    ii = repelem([1:nie],2);
    jj1 = ((tri2tri(:,1)-1)*2+[1:2])';
    jj2 = ((tri2tri(:,2)-1)*2+[1:2])';
    kk = intedgevectors';
    plus = sparse(ii(:),jj1(:),kk(:),nie,2*data.numTriangles);
    minus = sparse(ii(:),jj2(:),-kk(:),nie,2*data.numTriangles);
    edgedotMatrixBase = plus + minus;
    edgedotMatrix = blkdiag(edgedotMatrixBase,edgedotMatrixBase);
    
    %% perturb edge in integrability
    nT = data.numTriangles;
    ulower = (lowertriind-1)*2+[1:2];
    uhigher = (highertriind-1)*2+[1:2];
    vlower = (lowertriind-1)*2+[1:2] + 2*nT;
    vhigher = (highertriind-1)*2+[1:2] + 2*nT;
    
    evec = data.vertices(data.edges(intedgeind,2),1:2)-data.vertices(data.edges(intedgeind,1),1:2);
    edgedotMatrix(intedgeind,:)=0;
    edgedotMatrix(intedgeind+nie,:)=0;
    edgedotMatrix(intedgeind,ulower)=evec;
    edgedotMatrix(intedgeind,vhigher)=evec;
    edgedotMatrix(intedgeind+nie,vlower)=evec;
    edgedotMatrix(intedgeind+nie,uhigher)=-evec;
    
    %% verify boundary constraints
    be2btri = data.edges2triangles(data.isBoundaryEdge,1);
    be = data.vertices(data.edges(data.isBoundaryEdge,1),1:2) - data.vertices(data.edges(data.isBoundaryEdge,2),1:2);
    u_norm_ind = abs(dot(u(be2btri,:),be,2)) < 1e-6;
    v_norm_ind = abs(dot(v(be2btri,:),be,2)) < 1e-6;
    assert(all(u_norm_ind|v_norm_ind))
    
    nbe = sum(data.isBoundaryEdge);
    ii = repelem([1:nbe],2);
    jj = ((be2btri(:)-1)*2+[1:2])';
    kk = be';
    boundaryNormalMatrixBase = sparse(ii(:),jj(:),kk(:),nbe,2*data.numTriangles);
    UboundaryNormalMatrix = boundaryNormalMatrixBase(u_norm_ind,:);
    VboundaryNormalMatrix = boundaryNormalMatrixBase(v_norm_ind,:);
    boundaryNormalMatrix = blkdiag(UboundaryNormalMatrix, VboundaryNormalMatrix);
    
    %% reopt
    nareas = data.triangleAreas/sum(data.triangleAreas); assert(~any(nareas < 0))
    fun = @(x) obfun_wrapper(x,nareas);
    constraintGradientMatrix = [edgedotMatrix;boundaryNormalMatrix];
    options = optimoptions('fmincon','Display','Iter','CheckGradients',false,'SpecifyObjectiveGradient',true,'Algorithm','interior-point','HessianApproximation','lbfgs','SpecifyConstraintGradient',true,...
        'OptimalityTolerance',optol);
    [x,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,x,[],[],constraintGradientMatrix,zeros(size(constraintGradientMatrix,1),1),[],[],[],options);
    xs{i} = x;

    J = permute(reshape(x,2,[],2),[1 3 2]);
    u = reshape(J(:,1,:),2,[])';
    v = reshape(J(:,2,:),2,[])';
    us{i}=u;
    vs{i}=v;
    fvals(i)=fval;
    datas{i}=data;
    lambdas{i}=lambda;

    dets = cross([us{i} us{i}(:,1)*0],[vs{i} vs{i}(:,1)*0]); 
    mindets(i) = min(dets(:,3));
    areas(:,i) = data.triangleAreas([highertriind, lowertriind]); 
end
assert(all(mindets>0))

%% get virtual frames
clear uvhigher uvlower;
for i=1:N
    uvlower(:,i) = [us{i}(lowertriind,:)'; vs{i}(lowertriind,:)'];
    uvhigher(:,i) = [us{i}(highertriind,:)'; vs{i}(highertriind,:)'];

end
K = numel(perturbfactors)-1;
figure; hold all; 
plot(perturbfactors(end-K:end), uvlower(:,end-K:end),'r.-')
plot(perturbfactors(end-K:end), uvhigher(:,end-K:end),'r.-')
title('primal var convergence');
for i=1:4
    pfit = polyfit(perturbfactors(end-K:end), uvlower(i,end-K:end), 2);
    uvlowerc(i) = pfit(end);
    pfit = polyfit(perturbfactors(end-K:end), uvhigher(i,end-K:end), 2);
    uvhigherc(i) = pfit(end);
end
uvlowerc = reshape(uvlowerc,2,2);
uvhigherc = reshape(uvhigherc,2,2);
Ehigherc = symmetricDirichlet(uvhigherc,1);
Elowerc = symmetricDirichlet(uvlowerc,1);
Ehlc = [Ehigherc Elowerc];

[Ehlc fval K]
[uvhigherc;uvlowerc]
gt = [uvhigherc(:,1); uvlowerc(:,1); uvhigherc(:,2); uvlowerc(:,2); ];

%% obtain primal virtual frames (this is the key part to test)
ehat = pdir;
v1edges = setdiff(find(any(ismember(data.edges,v1),2)), intedgeind);
v2edges = setdiff(find(any(ismember(data.edges,v2),2)), intedgeind);
e1 = intersect(data.triangles2edges([highertriind],:),v1edges);
e2 = intersect(data.triangles2edges([lowertriind],:),v1edges);
e3 = intersect(data.triangles2edges([highertriind],:),v2edges);
e4 = intersect(data.triangles2edges([lowertriind],:),v2edges);
t1 = setdiff(data.edges2triangles(e1,:),highertriind);
t2 = setdiff(data.edges2triangles(e2,:),lowertriind);
t5 = setdiff(data.edges2triangles(e3,:),highertriind);
t6 = setdiff(data.edges2triangles(e4,:),lowertriind);
edgev = data.vertices(data.edges(:,2),1:2)-data.vertices(data.edges(:,1),1:2);
e1v = edgev(e1,:);
e2v = edgev(e2,:);
if data.edges(e1,2)~=v1; e1v = -e1v; end;
if data.edges(e2,2)~=v1; e2v = -e2v; end;
e3v = edgev(e3,:);
e4v = edgev(e4,:);
u1 = u(t1,:);
vv1 = v(t1,:);
u2 = u(t2,:);
vv2 = v(t2,:);
u5o = u(t5,:); % uv56 not converged
v5o = v(t5,:);
u6o = u(t6,:);
v6o = v(t6,:);

%% get frames 5,6 on convergence
clear uv5 uv6;
for i=1:N
    uv5(:,i) = [us{i}(t5,:)'; vs{i}(t5,:)'];
    uv6(:,i) = [us{i}(t6,:)'; vs{i}(t6,:)'];    
end
K = numel(perturbfactors)-1;
figure; hold all; title('frames 5,6 convergence')
plot(perturbfactors(end-K:end), uv5(:,end-K:end),'r.-')
plot(perturbfactors(end-K:end), uv6(:,end-K:end),'r.-')
for i=1:4
    pfit = polyfit(perturbfactors(end-K:end), uv5(i,end-K:end), 2);
    uv5c(i) = pfit(end);
    pfit = polyfit(perturbfactors(end-K:end), uv6(i,end-K:end), 2);
    uv6c(i) = pfit(end);
end
uv5c = reshape(uv5c,2,2);
uv6c = reshape(uv6c,2,2);
u5 = uv5c(:,1)';
v5 = uv5c(:,2)';
u6 = uv6c(:,1)';
v6 = uv6c(:,2)';

%{
% get duv(t3,t5,t4,t6)
t3 = highertriind; t4 = lowertriind;
for i=1:numel(perturbfactors)
    u3s(:,i) = us{i}(t3,:);
    u5s(:,i) = us{i}(t5,:);
    v3s(:,i) = vs{i}(t3,:);
    v5s(:,i) = vs{i}(t5,:);
    u4s(:,i) = us{i}(t4,:);
    u6s(:,i) = us{i}(t6,:);
    v4s(:,i) = vs{i}(t4,:);
    v6s(:,i) = vs{i}(t6,:);
end
xx = norm(pdir)*perturbfactors;
% figure; plot(xx,[u3s;u5s;v3s;v5s;u4s;u6s;v4s;v6s;],'r.-')
K = 3;
for i=1:numel(perturbfactors)
    for j=1:2
        pfit = polyfit(xx(end-K:end),u3s(j,end-K:end),2);
        du3(j,1) = pfit(2);
        pfit = polyfit(xx(end-K:end),u5s(j,end-K:end),2);
        du5(j,1) = pfit(2);
        pfit = polyfit(xx(end-K:end),v3s(j,end-K:end),2);
        dv3(j,1) = pfit(2);
        pfit = polyfit(xx(end-K:end),v5s(j,end-K:end),2);
        dv5(j,1) = pfit(2);

        pfit = polyfit(xx(end-K:end),u4s(j,end-K:end),2);
        du4(j,1) = pfit(2);
        pfit = polyfit(xx(end-K:end),u6s(j,end-K:end),2);
        du6(j,1) = pfit(2);
        pfit = polyfit(xx(end-K:end),v4s(j,end-K:end),2);
        dv4(j,1) = pfit(2);
        pfit = polyfit(xx(end-K:end),v6s(j,end-K:end),2);
        dv6(j,1) = pfit(2);
    end
end
p1a = [dot(u3(:)-u5(:),ehat); dot(v3(:)-v5(:),ehat)]/norm(ehat);
p1b = [dot(du3-du5,e1v); dot(dv3-dv5,e1v)];
p2a = [dot(u4(:)-u6(:),ehat); dot(v4(:)-v6(:),ehat)]/norm(ehat);
p2b = [dot(du4-du6,e2v); dot(dv4-dv6,e2v)];
[p1a p1b]
[p2a p2b]
p1 = [dot(u3(:)-u1(:),ehat); dot(v3(:)-vv1(:),ehat)];
p2 = [dot(u3(:)-u5(:),ehat); dot(v3(:)-v5(:),ehat)];
p3 = [dot(u4(:)-u6(:),ehat); dot(v4(:)-v6(:),ehat)];
p4 = [dot(u4(:)-u2(:),ehat); dot(v4(:)-vv2(:),ehat)];
[p1-p2 p3-p4]

% check duv with symDir grad
[E5,g5]=obfun_wrapper([u5(:); v5(:)], 1, @symmetricDirichlet);
[E6,g6]=obfun_wrapper([u6(:); v6(:)], 1, @symmetricDirichlet);
g5./[du5; dv5]
g6./[du6; dv6]
%}

z = zeros(1,2);
% nareas = data0.triangleAreas([highertriind, lowertriind]); 
% nareas = data.triangleAreas([highertriind, lowertriind]); 
% nareas = data.edgeLengths([e1,e2]);
nareas = [mean(areas(1,:)./areas(2,:)); 1];
fun = @(x) obfun_wrapper(x,nareas);
% full constr - nonconverged 5,6 tris.
Af = [e1v z z z; z z e1v z; z e2v z z; z z z e2v; z ehat ehat z; -ehat z z ehat;...
    e3v z z z; z z e3v z; z e4v z z; z z z e4v];
Bf = [dot(u1,e1v); dot(vv1,e1v); dot(u2,e2v); dot(vv2,e2v); 0; 0;...
    dot(u5o,e3v); dot(v5o,e3v); dot(u6o,e4v); dot(v6o,e4v); ];
gt - Af\Bf
% redundant edge removed
Am = [e1v z z z; z z e1v z; z e2v z z; z z z e2v; z ehat ehat z; -ehat z z ehat];
Bm = [dot(u1,e1v); dot(vv1,e1v); dot(u2,e2v); dot(vv2,e2v); 0; 0];
% redundant edge approximated
% eps = 1.7e-6;
eps = 1.7e-8;
e3p = e1v + eps * ehat/norm(ehat);
e4p = e2v + eps * ehat/norm(ehat);
% e1p = e1v - eps * ehat/norm(ehat);
% e2p = e2v - eps * ehat/norm(ehat);
e1p = e1v ;
e2p = e2v ;
% ehp = eps*ehat/norm(ehat);
ehp = ehat/norm(ehat);
A = [e1p z z z; z z e1p z; z e2p z z; z z z e2p; z ehp ehp z; -ehp z z ehp;...
    e3p z z z; z z e3p z; z e4p z z; z z z e4p];
B = [dot(u1,e1p); dot(vv1,e1p); dot(u2,e2p); dot(vv2,e2p); 0; 0;...
    dot(u5,e3p); dot(v5,e3p); dot(u6,e4p); dot(v6,e4p); ];
% eps2 = 1e-10; A(7:end,:)=A(7:end,:)*eps2; B(7:end)=B(7:end)*eps2;
ad = A\B; mstd = [mean(ad./gt) std(ad./gt)]
fac = mstd(1); [ad/fac gt]
norm(ad)/norm(gt)

konst = 0.00000083254366589661970791;
fac2 = konst/eps;
fac2 = norm(ad)/norm(gt);
[fac2 norm(gt)]

nsp = null(Am);
x0 = Am\Bm;
sol = [nsp ad]\-x0;
remade = nsp*sol(1:2)+x0;

% options = optimoptions('fmincon','Display','Iter','CheckGradients',false,'SpecifyObjectiveGradient',true,'Algorithm','interior-point','HessianApproximation','lbfgs','SpecifyConstraintGradient',true,...
%     'OptimalityTolerance',optol,'ConstraintTolerance',1e-6, 'EnableFeasibilityMode',true);
% x0 = randn(8,1);
% x0 = [uvhigherc(:,1); uvlowerc(:,1); uvhigherc(:,2); uvlowerc(:,2)];
% x0 = xs{end}([(highertriind-1)*2+[1:2] (lowertriind-1)*2+[1:2] (highertriind-1)*2+[1:2] + data.numTriangles*2 (lowertriind-1)*2+[1:2]+data.numTriangles*2]);
% % [x,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,x0,[],[],A,B,[],[],@jacdets,options);
% [x,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,x0,[],[],A,B,[],[],[],options);
x = ad;
u3 = x(1:2);
u4 = x(3:4);
v3 = x(5:6);
v4 = x(7:8);


uvhigher2 = [u3 v3];
uvlower2 = [u4 v4];
Ehigher2 = symmetricDirichlet(uvhigher2,1);
Elower2 = symmetricDirichlet(uvlower2,1);
Ehlc = [Ehigher2 Elower2];

%% replay alternatively obtained primal virtual frames in small perturbation mesh
% x = xs{end};
% J = permute(reshape(x,2,[],2),[1 3 2]);
% u = reshape(J(:,1,:),2,[])';
% v = reshape(J(:,2,:),2,[])';
% u(highertriind,:) = uvhigher2(1:2);
% v(highertriind,:) = uvhigher2(3:4);
% u(lowertriind,:) = uvlower2(1:2);
% v(lowertriind,:) = uvlower2(3:4);
% x = [reshape(u',[],1);reshape(v',[],1)];
% 
% nareas = data.triangleAreas/sum(data.triangleAreas); assert(~any(nareas < 0))
% fun = @(x) obfun_wrapper(x,nareas);
% constraintGradientMatrix = [edgedotMatrix;boundaryNormalMatrix];
% options = optimoptions('fmincon','Display','Iter','CheckGradients',false,'SpecifyObjectiveGradient',true,'Algorithm','interior-point','HessianApproximation','lbfgs','SpecifyConstraintGradient',true,...
%     'OptimalityTolerance',optol);
% [x,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,x,[],[],constraintGradientMatrix,zeros(size(constraintGradientMatrix,1),1),[],[],[],options);
% J = permute(reshape(x,2,[],2),[1 3 2]);
% u = reshape(J(:,1,:),2,[])';
% v = reshape(J(:,2,:),2,[])';
% [u(highertriind,:) us{end}(highertriind,:);...
% v(highertriind,:) vs{end}(highertriind,:);...
% u(lowertriind,:) us{end}(lowertriind,:);...
% v(lowertriind,:) vs{end}(lowertriind,:)]



%% get neighboring triangles
t1234 = setdiff(find(any(ismember(data.triangles,v2),2)),[highertriind, lowertriind]);
nt = numel(t1234);
for i=1:nt
    ti = t1234(i);
    uvi = [us{end}(ti,:)' vs{end}(ti,:)'];
    Eic(i) = symmetricDirichlet(uvi,1);
end

%% get neighboring edges and dual variables
einds = find(any(data.edges == v2,2));
clear lambdaedges;
for i=1:N
    lambdaedges(:,i) = lambdas{i}.eqlin(einds);
end
figure; largeinds = (vecnorm(lambdaedges,2,2) > .01); 
subplot(2,1,1);plot(perturbfactors,lambdaedges(largeinds,:), 'r.-'); title('nonzero dual vars');
subplot(2,1,2);plot(perturbfactors,lambdaedges(~largeinds,:), 'g.-'); title('zero dual vars');

EA_unnormalized = [Eic Ehlc] .* data.triangleAreas([t1234', lowertriind, highertriind])';

areaFromV123 = @(v123) dot(cross(v123(2,:)-v123(1,:),v123(3,:)-v123(1,:)),[0 0 1])/2;
syms alpha real;
V = datas{end}.vertices;
p2 = V(v1,:) + [pdir/norm(pdir) 0] * alpha;

tnhl = [t1234; highertriind; lowertriind];
for i=1:numel(tnhl)
    ti = tnhl(i);
    vind = find(data.triangles(ti,:) == v2);
    v123 = sym(V(data.triangles(ti,:),:));
    v123(vind,:)=p2;
    Ai(i) = areaFromV123(v123);
end
dAda = [double(diff(Ai,alpha))]/sum(data0.triangleAreas);
adiff = dot(dAda, [Eic Ehlc]);


xx=perturbfactors((end-K):end)*norm(pdir);yy=fvals((end-K):end);
slope = polyfit(xx,yy,2);
figure; hold all; 
plot(perturbfactors*norm(pdir), fvals, '.-')
xlabel('distance');
ylabel('local min sym dir');
title(sprintf('curvature:%f, slope:%f, yinter:%f',slope));
xsp = linspace(0,abs(slope(2)/slope(1)),100);
plot(xsp,polyval(slope, xsp))
fdiff = slope(2);

[adiff fdiff]

%% visualize mesh
%{
figure('units','normalized','outerposition',[0 0 1 1]);
hold all; axis equal; axis tight manual;
xlim([1.0628       1.1844]); ylim([0.39906      0.49494]);
i=1;
    data=datas{i};
    T=data.triangles;V=data.vertices;
    TV = data.triangleBarycenters;
    u=us{i};
    v=vs{i};
    g{1}=patch('faces',T,'Vertices',V,'facecolor','blue','edgecolor','black');
patch('faces',T([t1,t2,t3,t4],:),'Vertices',V,'facecolor','green','edgecolor','black','facealpha',.5);
patch('vertices',V,'faces',data.edges(einds,[1 2 1]),'edgecolor','green','linewidth',2);
    g{2}=patch('vertices',V,'faces',data.edges(intedgeind,[1 2 1]),'edgecolor','cyan','linewidth',2);
    g{3}=quiver([TV(:,1)],[TV(:,2)],[u(:,1)],[u(:,2)],'r','ShowArrowHead','off');
    g{4}=quiver([TV(:,1)],[TV(:,2)],[v(:,1)],[v(:,2)],'g','ShowArrowHead','off');
%}
