clear all; close all; 

%% toggles
nonInversionConstraint = true; % use jacdet > 0 or not.
flipSingularPermutation = true; % true means increasing energy. false means decrease

%% load data
meshname = 'data/polar_wedge_no_singularities.obj';

[V,T,UV,TF,N,NF] = readOBJ(meshname);
[u,v] = loadFrameText('data/frames.txt');
x = [reshape(u',[],1); reshape(v',[],1)];

%% select perturb edge. displace vertex
% [~,intedgeind]=min(vecnorm([data.vertices(data.edges(~data.isBoundaryEdge,1),:) data.vertices(data.edges(~data.isBoundaryEdge,2),:)] - [1.63744 .652415 0 1.6609 .661762 0],2,2));
% [~,lowertriind] = min(vecnorm(TV - [1.65095 .640671 0],2,2));
% [~,highertriind] = min(vecnorm(TV - [1.6469 .673558 0],2,2));
intedgeind = 4198;
lowertriind = 2928;
highertriind = 2929;
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

%% get old
oldConstraintGradientMatrix = [edgedotMatrix];
oldViolation = max(abs(oldConstraintGradientMatrix*x));

%% perturb edge in integrability
nT = data.numTriangles;
ulower = (lowertriind-1)*2+[1:2];
uhigher = (highertriind-1)*2+[1:2];
vlower = (lowertriind-1)*2+[1:2] + 2*nT;
vhigher = (highertriind-1)*2+[1:2] + 2*nT;

evec = data.vertices(data.edges(intedgeind,2),1:2)-data.vertices(data.edges(intedgeind,1),1:2);
edgedotMatrix(intedgeind,:)=0;
edgedotMatrix(intedgeind+nie,:)=0;
if ~flipSingularPermutation
    % correct permutation. decrease energy
    edgedotMatrix(intedgeind,ulower)=evec;
    edgedotMatrix(intedgeind,vhigher)=evec;
    edgedotMatrix(intedgeind+nie,vlower)=evec;
    edgedotMatrix(intedgeind+nie,uhigher)=-evec;
else
    % wrong permute. increase energy
    edgedotMatrix(intedgeind,ulower)=evec;
    edgedotMatrix(intedgeind,vhigher)=-evec;
    edgedotMatrix(intedgeind+nie,vlower)=evec;
    edgedotMatrix(intedgeind+nie,uhigher)=evec;
end


%% reopt
freeTriangles = [lowertriind highertriind]';
for i=1:1
    vinds = unique(data.triangles(freeTriangles,:));
    freeTriangles = find(any(ismember(data.triangles,vinds),2));
end

constraintGradientMatrix = [edgedotMatrix];
freeTrianglesXindsU = (freeTriangles-1)*2+[1:2];
freeTrianglesXindsV = (freeTriangles-1)*2+[1:2] + data.numTriangles*2;
freeXinds = [freeTrianglesXindsU(:); freeTrianglesXindsV(:)];
fixedXinds = setdiff(1:numel(x),freeXinds);
fixedX = x; fixedX(freeXinds)=0;
RHS = -constraintGradientMatrix*fixedX;
A = constraintGradientMatrix(:,freeXinds);
freeX = A\RHS;
newViolation = max(abs(A*freeX-RHS));

nareas = data.triangleAreas/sum(data.triangleAreas); nareas = nareas(freeTriangles);
fun = @(x) obfun_wrapper(x,nareas);
options = optimoptions('fmincon','Display','Iter','CheckGradients',false,'SpecifyObjectiveGradient',true,'Algorithm','interior-point','HessianApproximation','lbfgs','SpecifyConstraintGradient',true);
[freeX,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,x(freeXinds),[],[],A,zeros(size(RHS,1),1),[],[],@jacdets,options);
xs{i} = x;


[oldViolation, newViolation]

x(freeXinds) = freeX;
J = permute(reshape(x,2,[],2),[1 3 2]);
u = reshape(J(:,1,:),2,[])';
v = reshape(J(:,2,:),2,[])';
dets = cross([u u(:,1)*0],[v v(:,1)*0]); 
mindets = min(dets(:,3));


%% visss
close all;
figure('units','normalized','outerposition',[0 0 1 1]);
hold all; axis equal; axis tight manual;
xlim([1.5196       1.8296]); ylim([0.53584      0.78032]);
T=data.triangles;V=data.vertices;
TV = data.triangleBarycenters;
g{1}=patch('faces',T,'Vertices',V,'facecolor','blue','edgecolor','black');
g{2}=patch('vertices',V,'faces',data.edges(intedgeind,[1 2 1]),'edgecolor','cyan','linewidth',2);
g{3}=quiver([TV(:,1)],[TV(:,2)],[u(:,1)],[u(:,2)],'r','ShowArrowHead','off');
g{4}=quiver([TV(:,1)],[TV(:,2)],[v(:,1)],[v(:,2)],'g','ShowArrowHead','off');

