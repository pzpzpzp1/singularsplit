clear all; close all; 

%% load data
meshname = 'data/polar_wedge_no_singularities.obj';

[V,T,UV,TF,N,NF] = readOBJ(meshname);
[u,v] = loadFrameText('data/frames.txt');
x = [reshape(u',[],1); reshape(v',[],1)];
data0 = getMeshData(V,T); data=data0; V0=V;

perturbfactors = [.003 .001 .0005 .0003 .0001]; 
N = numel(perturbfactors);
for i=1:N
    perturbfactor=perturbfactors(i);
    %% select perturb edge. displace vertex
    % [~,v1] = min(vecnorm(data.vertices - [1.09788 .437436 0],2,2))
    % [~,v2] = min(vecnorm(data.vertices - [ 1.12134 .446783 0],2,2))
    % [~,lowertriind] = min(vecnorm(TV - [1.10935 .453734 0],2,2));
    % [~,highertriind] = min(vecnorm(TV - [1.10951 .430588 0],2,2));
    v1 = 1428; v2 = 1429;
    [~, intedgeind] = min(vecnorm(sort(data0.edges(~data0.isBoundaryEdge,:),2) - sort([v1 v2],2),2,2));
    lowertriind = 1595;
    highertriind = 1594;
    pdir = V0(v2,1:2)-V0(v1,1:2);
    V(v2,1:2) = V(v1,1:2) + pdir * perturbfactor;
    data = getMeshData(V,T);

%{
figure; hold all; axis equal; axis tight manual;
TV = data.triangleBarycenters;
g{1}=patch('faces',T,'Vertices',V,'facecolor','blue','edgecolor','black');
g{2}=patch('vertices',V,'faces',data.edges(intedgeind,[1 2 1]),'edgecolor','cyan','linewidth',2);
g{3}=quiver([TV(:,1)],[TV(:,2)],[u(:,1)],[u(:,2)],'r','ShowArrowHead','off');
g{4}=quiver([TV(:,1)],[TV(:,2)],[v(:,1)],[v(:,2)],'g','ShowArrowHead','off');
%}

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
    options = optimoptions('fmincon','Display','Iter','CheckGradients',false,'SpecifyObjectiveGradient',true,'Algorithm','interior-point','HessianApproximation','lbfgs','SpecifyConstraintGradient',true);
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

%% get virtual frames
clear uvhigher uvlower;
for i=1:N
    uvlower(:,i) = [us{i}(lowertriind,:)'; vs{i}(lowertriind,:)'];
    uvhigher(:,i) = [us{i}(highertriind,:)'; vs{i}(highertriind,:)'];
end
K = 4;
figure; hold all; 
plot(perturbfactors(end-K:end), uvlower(:,end-K:end),'r.-')
plot(perturbfactors(end-K:end), uvhigher(:,end-K:end),'r.-')
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
TV = data.triangleBarycenters;
[d1,t1] = min(vecnorm(TV-[1.11294 .467137 0],2,2));
[d2,t2] = min(vecnorm(TV-[1.12493 .460186 0],2,2));
[d3,t3] = min(vecnorm(TV-[1.12526 .436616 0],2,2));
[d4,t4] = min(vecnorm(TV-[1.11343 .420421 0],2,2));
t1234 = [t1 t2 t3 t4];
uv1 = [us{end}(t1,:)' vs{end}(t1,:)'];
uv2 = [us{end}(t2,:)' vs{end}(t2,:)'];
uv3 = [us{end}(t3,:)' vs{end}(t3,:)'];
uv4 = [us{end}(t4,:)' vs{end}(t4,:)'];
E1c = symmetricDirichlet(uv1,1);
E2c = symmetricDirichlet(uv2,1);
E3c = symmetricDirichlet(uv3,1);
E4c = symmetricDirichlet(uv4,1);
E1234c = [E1c E2c E3c E4c];

%% get neighboring edges and dual variables
einds = find(any(data.edges == v2,2));
clear lambdaedges;
for i=1:N
    lambdaedges(:,i) = lambdas{i}.eqlin(einds);
end
figure; plot(perturbfactors,lambdaedges, 'r.-')

[E1234c Ehlc] .* data.triangleAreas([t1234, lowertriind, highertriind])'

areaFromV123 = @(v1,v2,v3) dot(cross(v2-v1,v3-v1),[0 0 1])/2;

syms alpha real;
V = datas{end}.vertices;
p2 = V(v1,:) + [pdir/norm(pdir) 0] * alpha;
A1 = areaFromV123(p2, data.vertices(data.triangles(t1,2),:), data.vertices(data.triangles(t1,3),:));
A2 = areaFromV123(p2, data.vertices(data.triangles(t2,2),:), data.vertices(data.triangles(t2,3),:));
A3 = areaFromV123(data.vertices(data.triangles(t3,1),:), data.vertices(data.triangles(t3,2),:), p2);
A4 = areaFromV123(data.vertices(data.triangles(t4,1),:), data.vertices(data.triangles(t4,2),:), p2);
Al = areaFromV123(data.vertices(data.triangles(lowertriind,1),:), p2, data.vertices(data.triangles(lowertriind,3),:));
Ah = areaFromV123(data.vertices(data.triangles(highertriind,1),:), p2, data.vertices(data.triangles(highertriind,3),:));
dAda = [double(diff(A1,alpha)),...
double(diff(A2,alpha)),...
double(diff(A3,alpha)),...
double(diff(A4,alpha)),...
double(diff(Ah,alpha)),...
double(diff(Al,alpha))]/sum(data0.triangleAreas);
adiff = dot(dAda, [E1234c Ehlc]);


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

