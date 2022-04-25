close all;

%% load data
meshname = 'data/polar_wedge_no_singularities.obj';
% removedEdges = [4509 4455 4452 4449 4446 4443 4557 4410 4539];
removedEdges = [];

[V,T,UV,TF,N,NF] = readOBJ(meshname);
[u,v] = loadFrameText('frames.txt'); invT = false;

if invT
    J = reshape([u v]',2,2,[]);
    a = J(1,1,:);
    b = J(1,2,:);
    c = J(2,1,:);
    d = J(2,2,:);
    detJ = a.*d-b.*c;
    invJ = [d -b; -c a]./detJ;
    J = permute(invJ,[2 1 3]);
    u = reshape(J(:,1,:),2,[])';
    v = reshape(J(:,2,:),2,[])';
end
u0=u;v0=v;
x0 = [reshape(u',[],1); reshape(v',[],1)];
J0 = permute(reshape(x0,2,[],2),[1 3 2]);
data = getMeshData(V,T);

%% verify integrability
intEdgeInds = find(~data.isBoundaryEdge);
tri2tri = data.edges2triangles(~data.isBoundaryEdge,:);
intedgevectors = data.vertices(data.edges(~data.isBoundaryEdge,1),1:2)-data.vertices(data.edges(~data.isBoundaryEdge,2),1:2);
intedgevectors = intedgevectors./vecnorm(intedgevectors,2,2);
intedgevectors = intedgevectors.*(abs(randn(size(intedgevectors,1),1))*4+1);

nie = size(tri2tri,1);
ii = repelem([1:nie],2);
jj1 = ((tri2tri(:,1)-1)*2+[1:2])';
jj2 = ((tri2tri(:,2)-1)*2+[1:2])';
kk = intedgevectors';
plus = sparse(ii(:),jj1(:),kk(:),nie,2*data.numTriangles);
minus = sparse(ii(:),jj2(:),-kk(:),nie,2*data.numTriangles);
edgedotMatrixBase = plus + minus;
edgedotMatrix = blkdiag(edgedotMatrixBase,edgedotMatrixBase);

uv_int = edgedotMatrix*[reshape(u',[],1);reshape(v',[],1)];
fprintf('max integrability violation: %f\n',max(abs([uv_int])))

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
boundaryViolation = boundaryNormalMatrix*[reshape(u',[],1); reshape(v',[],1)];
fprintf('max boundary tangency violation: %f\n',max(abs(boundaryViolation)))

%% get objective energy and gradient
Js = zeros(2,2,data.numTriangles);
Js(:,1,:) = u';
Js(:,2,:) = v';
dlJs = dlarray(Js);
nareas = data.triangleAreas/sum(data.triangleAreas);
[energy, energyGrad] = dlfeval(@symmetricDirichlet,dlJs,nareas);

energyGradReshaped = extractdata(reshape(permute(energyGrad,[1 3 2]),[],1));
edgedotMatrix([removedEdges(:); removedEdges(:)+nie],:)=0;

%% get dual variables
constraintGradientMatrix = [edgedotMatrix;boundaryNormalMatrix];
dualvars = constraintGradientMatrix'\-energyGradReshaped;
kktopt = norm(dualvars'*constraintGradientMatrix + energyGradReshaped');
fprintf('kkt condition optmality at %f\n',kktopt)
% dualvar_nullspace = eigs(constraintGradientMatrix*constraintGradientMatrix',5,'sm')';

%% re-solve for more optimal point
options = optimoptions('fmincon','Display','Iter','CheckGradients',false,'SpecifyObjectiveGradient',true,'Algorithm','interior-point','HessianApproximation','lbfgs',...
    'ConstraintTolerance',1e-6,'OptimalityTolerance',1e-8);
%     func = @symmetrizedSymmetricDirichlet;
%     func = @symmetrizedDirichlet;
func = @symmetricDirichlet;
fun = @(x) obfun_wrapper(x,nareas,func);
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,x0,[],[],constraintGradientMatrix,zeros(size(constraintGradientMatrix,1),1),[],[],[],options);
% xp = x0 + randn(size(x0))*.01;
% [x,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,xp,[],[],constraintGradientMatrix,zeros(size(constraintGradientMatrix,1),1),[],[],[],options);

%% verify kkt
[~, egrad] = fun(x);
kktopt = norm(lambda.eqlin'*constraintGradientMatrix + egrad');
integLambdasU = lambda.eqlin(1:nie);
integLambdasV = lambda.eqlin((nie+1):2*nie);
feasibility = abs(constraintGradientMatrix*x);

figure; hold all; rotate3d on; axis equal; title('constraint satisfaction')
patch('faces',T,'Vertices',V,'facecolor','blue','edgecolor','none')
edgescore = (max(reshape(feasibility(1:2*nie),[],2)')');
xx = reshape(data.vertices(data.edges(~data.isBoundaryEdge,:)',1),2,[]);
yy = reshape(data.vertices(data.edges(~data.isBoundaryEdge,:)',2),2,[]);
xx(3,:)=nan;yy(3,:)=nan;xx = xx(:);yy = yy(:);
patch(xx,yy,yy*0,repelem(edgescore,3,1),'edgecolor','interp','linewidth',2)
colorbar;colormap(inferno);

figure; hold all; rotate3d on; axis equal; title('integrability: dual var')
patch('faces',T,'Vertices',V,'facecolor','blue','edgecolor','none')
edgescore = vecnorm([integLambdasU integLambdasV],2,2);
xx = reshape(data.vertices(data.edges(~data.isBoundaryEdge,:)',1),2,[]);
yy = reshape(data.vertices(data.edges(~data.isBoundaryEdge,:)',2),2,[]);
xx(3,:)=nan;yy(3,:)=nan;xx = xx(:);yy = yy(:);
patch(xx,yy,yy*0,repelem(edgescore,3,1),'edgecolor','interp','linewidth',2)
colorbar;colormap(inferno);

figure; 
subplot(1,2,1); hold all; rotate3d on; axis equal; title('integrability: split lambdas')
patch('faces',T,'Vertices',V,'facecolor','blue','edgecolor','none')
edgescore = vecnorm([integLambdasU],2,2);
xx = reshape(data.vertices(data.edges(~data.isBoundaryEdge,:)',1),2,[]);
yy = reshape(data.vertices(data.edges(~data.isBoundaryEdge,:)',2),2,[]);
xx(3,:)=nan;yy(3,:)=nan;xx = xx(:);yy = yy(:);
patch(xx,yy,yy*0,repelem(edgescore,3,1),'edgecolor','interp','linewidth',2)
colorbar;colormap(inferno);
subplot(1,2,2); hold all; rotate3d on; axis equal; title('integrability: split lambdas')
patch('faces',T,'Vertices',V,'facecolor','blue','edgecolor','none')
edgescore = vecnorm([integLambdasV],2,2);
xx = reshape(data.vertices(data.edges(~data.isBoundaryEdge,:)',1),2,[]);
yy = reshape(data.vertices(data.edges(~data.isBoundaryEdge,:)',2),2,[]);
xx(3,:)=nan;yy(3,:)=nan;xx = xx(:);yy = yy(:);
patch(xx,yy,yy*0,repelem(edgescore,3,1),'edgecolor','interp','linewidth',2)
colorbar;colormap(inferno);

lambdadgdx = full(abs(lambda.eqlin.*vecnorm(constraintGradientMatrix,2,2)));
figure; hold all; rotate3d on; axis equal; title('dual var * |dgdx|')
patch('faces',T,'Vertices',V,'facecolor','blue','edgecolor','none')
edgescore = vecnorm(reshape(lambdadgdx(1:2*nie),[],2),2,2);
xx = reshape(data.vertices(data.edges(~data.isBoundaryEdge,:)',1),2,[]);
yy = reshape(data.vertices(data.edges(~data.isBoundaryEdge,:)',2),2,[]);
xx(3,:)=nan;yy(3,:)=nan;xx = xx(:);yy = yy(:);
patch(xx,yy,yy*0,repelem(edgescore,3,1),'edgecolor','interp','linewidth',2)
colorbar;colormap(inferno);

lambdadgdx = full(abs(lambda.eqlin.*vecnorm(constraintGradientMatrix,2,2)));
figure; hold all; rotate3d on; title('dual var * |dgdx| / edgelength'); axis equal;
patch('faces',T,'Vertices',V,'facecolor','blue','edgecolor','none')
edgescore = vecnorm(reshape(lambdadgdx(1:2*nie),[],2),2,2);
edgescore = edgescore./data.edgeLengths(1:nie);
xx = reshape(data.vertices(data.edges(~data.isBoundaryEdge,:)',1),2,[]);
yy = reshape(data.vertices(data.edges(~data.isBoundaryEdge,:)',2),2,[]);
xx(3,:)=nan;yy(3,:)=nan;xx = xx(:);yy = yy(:);
patch(xx,yy,yy*0,repelem(edgescore,3,1),'edgecolor','interp','linewidth',2)
colorbar;colormap(inferno);
bestscore = edgescore;
[~,maxind]=max(bestscore);
patch('vertices',data.vertices,'faces',data.edges(maxind,[1 2 1]),'edgecolor','green','linewidth',3)

faceScore = vecnorm(reshape(permute(reshape(lambda.eqlin'*constraintGradientMatrix,2,[],2),[2 1 3]),[],4),2,2)./data.triangleAreas;
figure; hold all; rotate3d on; axis equal; title('energygradient per face')
patch('faces',T,'Vertices',V,'facecolor','flat','edgecolor','none','FaceVertexCdata',faceScore)
colorbar;colormap(inferno)

%% update u,v
J = permute(reshape(x,2,[],2),[1 3 2]);
u = reshape(J(:,1,:),2,[])';
v = reshape(J(:,2,:),2,[])';

%% vis
% close all;
TV = data.triangleBarycenters;
figure; hold all; axis equal; rotate3d on;
patch('faces',T,'Vertices',V,'facecolor','blue','edgecolor','black')
quiver([TV(:,1);TV(:,1)],[TV(:,2);TV(:,2)],[u0(:,1);v0(:,1)],[u0(:,2);v0(:,2)],'r','ShowArrowHead','off')
quiver([TV(:,1);TV(:,1)],[TV(:,2);TV(:,2)],[u(:,1);v(:,1)],[u(:,2);v(:,2)],'g','ShowArrowHead','off')
patch('vertices',data.vertices,'faces',data.edges(removedEdges,[1 2 1]),'edgecolor','cyan','linewidth',3)
% patch('vertices',data.vertices,'faces',data.edges(4455,[1 2 1]),'edgecolor','cyan','linewidth',3)

iu = integrateVecfield(data,u);
iv = integrateVecfield(data,v);
newE = [iu iv];
newE = newE-mean(newE);

figure; hold all; axis equal; rotate3d on;
patch('faces',T,'Vertices',newE,'facecolor','blue','edgecolor','black')
scatter(newE(1,1),newE(1,2),'r','filled')
scatter(newE(end,1),newE(end,2),'g','filled')


