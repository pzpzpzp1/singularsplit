clear all; close all; 

%% load data
meshname = 'data/polar_wedge_no_singularities.obj';

[V,T,UV,TF,N,NF] = readOBJ(meshname); 
[u,v] = loadFrameText('data/frames.txt');
x = [reshape(u',[],1); reshape(v',[],1)];
data0 = getMeshData(V,T); data=data0; V0=V; T0=T;
intedgeind = randi(sum(~data.isBoundaryEdge));
intedgeind = 6587; % 34.893      -36.281
% intedgeind = 6304; % 18.375      -13.597: [-34.627      -6.1913  35.067       4.9743      -26.593      -17.474      -26.328      -18.432] 
% optol = 1e-9;
optol = 1e-6;

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

perturbfactors = [.0005 .0003 .0001 0]; 
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
    % nareas([highertriind, lowertriind])=0; % TEST IF VIRTUAL FRAME HAS INFLUENCE ON CONVERGENT FRAME VECTOR
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
end
assert(all(mindets>0))

%% get virtual frames
clear uvhigher uvlower;
for i=1:N
    uvlower(:,i) = [us{i}(lowertriind,:)'; vs{i}(lowertriind,:)'];
    uvhigher(:,i) = [us{i}(highertriind,:)'; vs{i}(highertriind,:)'];
end
K = 3;
figure; hold all; 
plot(perturbfactors(end-K:end), uvlower(:,end-K:end),'r.-')
plot(perturbfactors(end-K:end), uvhigher(:,end-K:end),'r.-')
title('primal var convergence');
for i=1:4
    pfit = polyfit(perturbfactors(end-K:end), uvlower(i,end-K:end), 1);
    uvlowerc(i) = pfit(2);
    pfit = polyfit(perturbfactors(end-K:end), uvhigher(i,end-K:end), 1);
    uvhigherc(i) = pfit(2);
end
uvlowerc = reshape(uvlowerc,2,2);
uvhigherc = reshape(uvhigherc,2,2);
Ehigherc = symmetricDirichlet(uvhigherc,1);
Elowerc = symmetricDirichlet(uvlowerc,1);
Ehlc = [Ehigherc Elowerc];

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
