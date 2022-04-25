clear all; close all; 

%% toggles
nonInversionConstraint = true; % use jacdet > 0 or not.
flipSingularPermutation = true; % true means increasing energy. false means decrease

%% load data
meshname = 'data/polar_wedge_no_singularities.obj';

[V,T,UV,TF,N,NF] = readOBJ(meshname);
[u,v] = loadFrameText('data/frames.txt');
x = [reshape(u',[],1); reshape(v',[],1)];

% run analyze polar to get a good initialization
perturbfactors = fliplr([0.05 0.03 0.01]); % high res sample
% perturbfactors = fliplr([.001 .0025 .005 .01 .05 .1 .2 .3 .4 .44 .43 .2]); % back and forth 
% perturbfactors = fliplr([.01 .03 .05 .1 .2 .3 .4 .43 .44 .46 .45 .44 .43 .4 .3 .2]); % overshoot discontinuity

N = numel(perturbfactors);
for i=1:N
    perturbfactor=perturbfactors(i);
    %% select perturb edge. displace vertex
    % [~,intedgeind]=min(vecnorm([data.vertices(data.edges(~data.isBoundaryEdge,1),:) data.vertices(data.edges(~data.isBoundaryEdge,2),:)] - [1.63744 .652415 0 1.6609 .661762 0],2,2));
    % [~,lowertriind] = min(vecnorm(TV - [1.65095 .640671 0],2,2));
    % [~,highertriind] = min(vecnorm(TV - [1.6469 .673558 0],2,2));
    intedgeind = 4198;
    lowertriind = 2928;
    highertriind = 2929;
    pdir = [0.0234600000000000	0.00934699999999999];
    V(1452,1:2) = V(1451,1:2) + pdir * perturbfactor;
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
    nareas = data.triangleAreas/sum(data.triangleAreas);
    fun = @(x) obfun_wrapper(x,nareas);
    constraintGradientMatrix = [edgedotMatrix;boundaryNormalMatrix];
    options = optimoptions('fmincon','Display','Iter','CheckGradients',false,'SpecifyObjectiveGradient',true,'Algorithm','interior-point','HessianApproximation','lbfgs','SpecifyConstraintGradient',true);
    if nonInversionConstraint
        [x,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,x,[],[],constraintGradientMatrix,zeros(size(constraintGradientMatrix,1),1),[],[],@jacdets,options);
    else
        [x,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,x,[],[],constraintGradientMatrix,zeros(size(constraintGradientMatrix,1),1),[],[],[],options);
    end
    xs{i} = x;

    J = permute(reshape(x,2,[],2),[1 3 2]);
    u = reshape(J(:,1,:),2,[])';
    v = reshape(J(:,2,:),2,[])';
    us{i}=u;
    vs{i}=v;
    fvals(i)=fval;
    datas{i}=data;

    dets = cross([us{i} us{i}(:,1)*0],[vs{i} vs{i}(:,1)*0]); 
    mindets(i) = min(dets(:,3));
end

close all;
figure('units','normalized','outerposition',[0 0 1 1]);
hold all; axis equal; axis tight manual;
xlim([1.5196       1.8296]); ylim([0.53584      0.78032]);
vname='singularpair';
vid1 = VideoWriter([vname '_uncomp.avi'],'Uncompressed AVI'); 
vid2 = VideoWriter([vname '.avi'],'Motion JPEG AVI'); vid2.Quality=100;
open(vid1); open(vid2);
for i=1:N
    try; for j=1:numel(g);delete(g{j}); end; catch; end;
    data=datas{i};
    T=data.triangles;V=data.vertices;
    TV = data.triangleBarycenters;
    u=us{i};
    v=vs{i};
    title(sprintf('symdir: %f',fvals(i)));
    g{1}=patch('faces',T,'Vertices',V,'facecolor','blue','edgecolor','black');
    g{2}=patch('vertices',V,'faces',data.edges(intedgeind,[1 2 1]),'edgecolor','cyan','linewidth',2);
    g{3}=quiver([TV(:,1)],[TV(:,2)],[u(:,1)],[u(:,2)],'r','ShowArrowHead','off');
    g{4}=quiver([TV(:,1)],[TV(:,2)],[v(:,1)],[v(:,2)],'g','ShowArrowHead','off');
    pause(.2)
    frame = getframe(gcf);
    writeVideo(vid1,frame);
    writeVideo(vid2,frame);
end
close(vid1)
close(vid2)


xx=perturbfactors(1:5);yy=fvals(1:5);
slope = polyfit(xx*norm(pdir),yy,1);
figure; hold all; 
plot(perturbfactors*norm(pdir), fvals, '.-')
xlabel('distance');
ylabel('local min sym dir');
title(sprintf('slope:%f, yinter:%f',slope));



% xx=perturbfactors(1:5);yy=fvals(1:5);
% slope = polyfit(xx*norm(pdir),yy,1);
% figure; hold all; 
% plot(perturbfactors(1:3)*norm(pdir), fvals(1:3), '.-')
% xlabel('distance');
% ylabel('local min sym dir');
% title(sprintf('slope:%f, yinter:%f',slope));


