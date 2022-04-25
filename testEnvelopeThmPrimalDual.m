
clear all; close all; 

%% load data
meshname = 'data/polar_wedge_no_singularities.obj';

[V,T,UV,TF,N,NF] = readOBJ(meshname);
[u,v] = loadFrameText('data/frames.txt');
x = [reshape(u',[],1); reshape(v',[],1)];
data0 = getMeshData(V,T); data=data0; V0=V; T0=T;
e0i = randi(sum(~data.isBoundaryEdge));
% e0i = 1079; % wroks ok. depends a lot on primal convergence.
% e0i = 6587; % seems to work. for specific ep
e0i = 3044; % works well
% e0i = 5360; % diagonal. works well
% e0i = 7948; % close to boundary. doesn't work well.
% e0i = 6304; % circumcentric orientation. doesn't work well/
optol = 1e-9;

p1 = data.edges(e0i,1);
p2 = data.edges(e0i,2);
e0 = V0(p2,1:2)-V0(p1,1:2);
triinds = data.edges2triangles(e0i,:);
tri1 = data0.triangles(triinds(1),:);
if numel(strfind([tri1 tri1(1)],[p1 p2]))~=0
    th = triinds(1); 
    tl = triinds(2); 
else
    th = triinds(2);
    tl = triinds(1); 
end
ehat = e0/norm(e0);
p1edges = setdiff(find(any(ismember(data.edges,p1),2)), e0i);
p2edges = setdiff(find(any(ismember(data.edges,p2),2)), e0i);
e1i = intersect(data.triangles2edges([th],:),p1edges);
e2i = intersect(data.triangles2edges([tl],:),p1edges);
e3i = intersect(data.triangles2edges([th],:),p2edges);
e4i = intersect(data.triangles2edges([tl],:),p2edges);
t0 = setdiff(data.edges2triangles(e1i,:),th);
t2 = setdiff(data.edges2triangles(e2i,:),tl);
t1 = setdiff(data.edges2triangles(e3i,:),th);
t3 = setdiff(data.edges2triangles(e4i,:),tl);
edgev = data.vertices(data.edges(:,2),1:2)-data.vertices(data.edges(:,1),1:2);
e1 = edgev(e1i,:);
e2 = edgev(e2i,:);
if data.edges(e1i,2)~=p1; e1 = -e1; end;
if data.edges(e2i,2)~=p1; e2 = -e2; end;
e3 = edgev(e3i,:);
e4 = edgev(e4i,:);
p4 = setdiff(data.triangles(th,:),[p1 p2]);
p7 = setdiff(data.triangles(tl,:),[p1 p2]);
p3 = setdiff(data.triangles(t0,:),[p1 p4]);
p5 = setdiff(data.triangles(t1,:),[p2 p4]);
p6 = setdiff(data.triangles(t2,:),[p1 p7]);
p8 = setdiff(data.triangles(t3,:),[p2 p7]);
e5i = find(all(sort(data.edges,2)==sort([p3 p4],2),2));
e6i = find(all(sort(data.edges,2)==sort([p3 p1],2),2));
e7i = find(all(sort(data.edges,2)==sort([p6 p1],2),2));
e8i = find(all(sort(data.edges,2)==sort([p6 p7],2),2));
e9i = find(all(sort(data.edges,2)==sort([p5 p4],2),2));
e10i = find(all(sort(data.edges,2)==sort([p2 p5],2),2));
e11i = find(all(sort(data.edges,2)==sort([p8 p2],2),2));
e12i = find(all(sort(data.edges,2)==sort([p7 p8],2),2));


%% display
figure; clf; hold all; axis tight manual; title('neighborhood');
TV = data.triangleBarycenters;
patch('faces',T,'Vertices',V,'facecolor','blue','edgecolor','black','facealpha','.5');
patch('faces',T(th,:),'Vertices',V,'facecolor','green','edgecolor','black','facealpha',.3);
patch('faces',T(tl,:),'Vertices',V,'facecolor','red','edgecolor','black','facealpha',.3);
scatter(TV(:,1),TV(:,2),'k','filled')
scatter(V0(p1,1),V0(p1,2),'g','filled')
scatter(V0(p2,1),V0(p2,2),'r','filled')
sc = [.4 .9 0];
bl = mean(V0(data.edges(e0i,:),:)) - [.1 .1 0].*sc;
tr = mean(V0(data.edges(e0i,:),:)) + [.1 .1 0].*sc;
xlim([bl(1) tr(1)]); ylim([bl(2) tr(2)])
text(data.vertices(p1,1),data.vertices(p1,2),'<--p1','color','y')
text(data.vertices(p2,1),data.vertices(p2,2),'<--p2','color','y')
text(data.vertices(p3,1),data.vertices(p3,2),'<--p3','color','y')
text(data.vertices(p4,1),data.vertices(p4,2),'<--p4','color','y')
text(data.vertices(p5,1),data.vertices(p5,2),'<--p5','color','y')
text(data.vertices(p6,1),data.vertices(p6,2),'<--p6','color','y')
text(data.vertices(p7,1),data.vertices(p7,2),'<--p7','color','y')
text(data.vertices(p8,1),data.vertices(p8,2),'<--p8','color','y')
text(data.triangleBarycenters(th,1),data.triangleBarycenters(th,2),'<--th','color','c')
text(data.triangleBarycenters(tl,1),data.triangleBarycenters(tl,2),'<--tl','color','c')
text(data.triangleBarycenters(t0,1),data.triangleBarycenters(t0,2),'<--t0','color','c')
text(data.triangleBarycenters(t1,1),data.triangleBarycenters(t1,2),'<--t1','color','c')
text(data.triangleBarycenters(t2,1),data.triangleBarycenters(t2,2),'<--t2','color','c')
text(data.triangleBarycenters(t3,1),data.triangleBarycenters(t3,2),'<--t3','color','c')
edgeBC = (data.vertices(data.edges(:,1),:)+data.vertices(data.edges(:,2),:))/2;
text(edgeBC(e0i,1),edgeBC(e0i,2),'<--e0','color',[.5 .9 .2])
text(edgeBC(e1i,1),edgeBC(e1i,2),'<--e1','color',[.5 .9 .2])
text(edgeBC(e2i,1),edgeBC(e2i,2),'<--e2','color',[.5 .9 .2])
text(edgeBC(e3i,1),edgeBC(e3i,2),'<--e3','color',[.5 .9 .2])
text(edgeBC(e4i,1),edgeBC(e4i,2),'<--e4','color',[.5 .9 .2])
text(edgeBC(e5i,1),edgeBC(e5i,2),'<--e5','color',[.5 .9 .2])
text(edgeBC(e6i,1),edgeBC(e6i,2),'<--e6','color',[.5 .9 .2])
text(edgeBC(e7i,1),edgeBC(e7i,2),'<--e7','color',[.5 .9 .2])
text(edgeBC(e8i,1),edgeBC(e8i,2),'<--e8','color',[.5 .9 .2])
text(edgeBC(e9i,1),edgeBC(e9i,2),'<--e9','color',[.5 .9 .2])
text(edgeBC(e10i,1),edgeBC(e10i,2),'<--e10','color',[.5 .9 .2])
text(edgeBC(e11i,1),edgeBC(e11i,2),'<--e11','color',[.5 .9 .2])
text(edgeBC(e12i,1),edgeBC(e12i,2),'<--e12','color',[.5 .9 .2])

% perturbfactors = [.003 .001 .0005 .0003 .0001 .00005]; 
perturbfactors = [.0005 .0001 .00001]; 
N = numel(perturbfactors);
for i=1:N
    % select perturb edge. displace vertex
    perturbfactor=perturbfactors(i);
    V=V0; V(p2,1:2) = V(p1,1:2) + e0 * perturbfactor;
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
    ulower = (tl-1)*2+[1:2];
    uhigher = (th-1)*2+[1:2];
    vlower = (tl-1)*2+[1:2] + 2*nT;
    vhigher = (th-1)*2+[1:2] + 2*nT;
    
    evec = e0; % data.vertices(data.edges(e0i,2),1:2)-data.vertices(data.edges(e0i,1),1:2);
    edgedotMatrix(e0i,:)=0;
    edgedotMatrix(e0i+nie,:)=0;
    edgedotMatrix(e0i,ulower)=evec;
    edgedotMatrix(e0i,vhigher)=evec;
    edgedotMatrix(e0i+nie,vlower)=evec;
    edgedotMatrix(e0i+nie,uhigher)=-evec;
    
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
    constraintGradientMatrices{i} = constraintGradientMatrix;
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
    areas(:,i) = data.triangleAreas([th, tl]); 
    grads{i}=grad;
    
    opts(i) = norm(grad + (lambdas{i}.eqlin'*constraintGradientMatrices{i})');
end
assert(all(mindets>0))

%% get virtual frames
clear uvhigher uvlower;
for i=1:N
    uvlower(:,i) = [us{i}(tl,:)'; vs{i}(tl,:)'];
    uvhigher(:,i) = [us{i}(th,:)'; vs{i}(th,:)'];
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

%% get neighboring frames 0123 on convergence
clear uv0 uv1 uv2 uv3;
for i=1:N
    uv0(:,i) = [us{i}(t0,:)'; vs{i}(t0,:)'];
    uv1(:,i) = [us{i}(t1,:)'; vs{i}(t1,:)'];    
    uv2(:,i) = [us{i}(t2,:)'; vs{i}(t2,:)'];
    uv3(:,i) = [us{i}(t3,:)'; vs{i}(t3,:)'];    
end
for i=1:4
    pfit = polyfit(perturbfactors(end-K:end), uv0(i,end-K:end), 2);
    uv0c(i) = pfit(end);
    pfit = polyfit(perturbfactors(end-K:end), uv1(i,end-K:end), 2);
    uv1c(i) = pfit(end);
    pfit = polyfit(perturbfactors(end-K:end), uv2(i,end-K:end), 2);
    uv2c(i) = pfit(end);
    pfit = polyfit(perturbfactors(end-K:end), uv3(i,end-K:end), 2);
    uv3c(i) = pfit(end);
end
uv0k=[uv0 uv0c'];
uv1k=[uv1 uv1c'];
uv2k=[uv2 uv2c'];
uv3k=[uv3 uv3c'];
figure; hold all; title('frames 0123 convergence')
plot([perturbfactors((end-K+1):end) 0], uv0k(:,end-K:end),'r.-')
plot([perturbfactors((end-K+1):end) 0], uv1k(:,end-K:end),'g.-')
plot([perturbfactors((end-K+1):end) 0], uv2k(:,end-K:end),'b.-')
plot([perturbfactors((end-K+1):end) 0], uv3k(:,end-K:end),'y.-')

%% get neighboring dual vars on convergence
clear lambs;
einds = [e0i e1i e2i e3i e4i e5i e6i e7i e8i e9i e10i e11i e12i];
for i=1:N
    for j=1:numel(einds)
        lambs(:,i,j) = lambdas{i}.eqlin([einds(j) einds(j)+nie]);
    end
end
for i=1:2
    for j=1:numel(einds)
        pfit = polyfit(perturbfactors(end-K:end), lambs(i,end-K:end,j), 2);
        lambcs(i,j) = pfit(end);
    end
end
figure; hold all; title('lambdas convergence')
for i=1:2
    for j=1:numel(einds)
        xx = [perturbfactors((end-K+1):end) 0];
        yy = [lambs(i,(end-K+1):end,j) lambcs(i,j)];
        plot(xx,yy,'.-')
    end
end

%% get nullspace for virtual frames
u0 = us{end}(t0,:);
v0 = vs{end}(t0,:);
u2 = us{end}(t2,:);
v2 = vs{end}(t2,:);
z = zeros(1,2);
Am = [e1 z z z; z z e1 z; z e2 z z; z z z e2; z ehat ehat z; -ehat z z ehat];
Bm = [dot(u0,e1); dot(v0,e1); dot(u2,e2); dot(v2,e2); 0; 0];


%% verify kkt at each epsilon
for k=1:K
    tti = th;
    areas = datas{k}.triangleAreas;
    uvti = [us{k}(tti,:)' vs{k}(tti,:)'];
    nareas = areas/sum(areas); fun = @(x) obfun_wrapper(x,nareas);
    [eval, egrad] = fun(xs{k});
    constraintGradientMatrix=constraintGradientMatrices{k};
    norm(egrad + (lambdas{k}.eqlin'*constraintGradientMatrix)') - opts(k)
    
    nT = data.numTriangles;
    einds
    xinds = [(tti-1)*2+[1:2] (tti-1)*2+[1:2]+2*nT];
    xs{k}(xinds) - uvti(:);
    fgradi = egrad(xinds);
    lambdaConstrgrad = (lambdas{k}.eqlin'*constraintGradientMatrix)';
    constraintForceI = lambdaConstrgrad(xinds);
    [norm(constraintForceI) norm([constraintForceI + fgradi])]
    [constraintForceI fgradi]

    [Ei, EgradI] = obfun_wrapper(xs{k}(xinds),1);
    lfullinds = data.triangles2edges(tti,:);
    lambdainds = [lfullinds lfullinds+nie];
    llocalinds = find(ismember(einds, lfullinds));
    
    localConstraintGradientMatrix = full(constraintGradientMatrix(lambdainds, xinds));    
    localLambda = lambdas{k}.eqlin(lambdainds);

    localConstraintGradientMatrixConverged = localConstraintGradientMatrix; 
    localConstraintGradientMatrixConverged(3,:)=localConstraintGradientMatrixConverged(2,:);
    localConstraintGradientMatrixConverged(6,:)=localConstraintGradientMatrixConverged(5,:)
    localLambdaConverged = reshape(lambcs(:,llocalinds)',[],1);

    [(localLambdaConverged'*localConstraintGradientMatrixConverged)',...
        (localLambdaConverged'*localConstraintGradientMatrix)',...
        (localLambda'*localConstraintGradientMatrix)',...
        constraintForceI]/nareas(th)
    [localLambdaConverged localConstraintGradientMatrixConverged]
    % found 

end

%{
Results are only that the dual variables do not go to 0. the permuted edge
lambda does go to 0 if that edge's constraint maintains unit length. the
edges adjacent to degenerate triangles have flipped sign dual variables so
that the energy gradient can sum to 0. provides no information about the
degenerate frame.
%}

error('no point continuing');


%% 
z = zeros(1,2);
% full constr - nonconverged 5,6 tris.
Af = [e1 z z z; z z e1 z; z e2 z z; z z z e2; z ehat ehat z; -ehat z z ehat;...
    e3 z z z; z z e3 z; z e4 z z; z z z e4];
Bf = [dot(u1,e1); dot(vv1,e1); dot(u2,e2); dot(vv2,e2); 0; 0;...
    dot(u5o,e3); dot(v5o,e3); dot(u6o,e4); dot(v6o,e4); ];
gt - Af\Bf
% redundant edge removed
Am = [e1 z z z; z z e1 z; z e2 z z; z z z e2; z ehat ehat z; -ehat z z ehat];
Bm = [dot(u1,e1); dot(vv1,e1); dot(u2,e2); dot(vv2,e2); 0; 0];
% redundant edge approximated
% eps = 1.7e-6;
eps = 1.7e-8;
e3p = e1 + eps * ehat/norm(ehat);
e4p = e2 + eps * ehat/norm(ehat);
% e1p = e1v - eps * ehat/norm(ehat);
% e2p = e2v - eps * ehat/norm(ehat);
e1p = e1 ;
e2p = e2 ;
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
t1234 = setdiff(find(any(ismember(data.triangles,p2),2)),[th, tl]);
nt = numel(t1234);
for i=1:nt
    ti = t1234(i);
    uvi = [us{end}(ti,:)' vs{end}(ti,:)'];
    Eic(i) = symmetricDirichlet(uvi,1);
end

%% get neighboring edges and dual variables
einds = find(any(data.edges == p2,2));
clear lambdaedges;
for i=1:N
    lambdaedges(:,i) = lambdas{i}.eqlin(einds);
end
figure; largeinds = (vecnorm(lambdaedges,2,2) > .01); 
subplot(2,1,1);plot(perturbfactors,lambdaedges(largeinds,:), 'r.-'); title('nonzero dual vars');
subplot(2,1,2);plot(perturbfactors,lambdaedges(~largeinds,:), 'g.-'); title('zero dual vars');

EA_unnormalized = [Eic Ehlc] .* data.triangleAreas([t1234', tl, th])';

areaFromV123 = @(v123) dot(cross(v123(2,:)-v123(1,:),v123(3,:)-v123(1,:)),[0 0 1])/2;
syms alpha real;
V = datas{end}.vertices;
p2 = V(p1,:) + [e0/norm(e0) 0] * alpha;

tnhl = [t1234; th; tl];
for i=1:numel(tnhl)
    ti = tnhl(i);
    vind = find(data.triangles(ti,:) == p2);
    v123 = sym(V(data.triangles(ti,:),:));
    v123(vind,:)=p2;
    Ai(i) = areaFromV123(v123);
end
dAda = [double(diff(Ai,alpha))]/sum(data0.triangleAreas);
adiff = dot(dAda, [Eic Ehlc]);


xx=perturbfactors((end-K):end)*norm(e0);yy=fvals((end-K):end);
slope = polyfit(xx,yy,2);
figure; hold all; 
plot(perturbfactors*norm(e0), fvals, '.-')
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
