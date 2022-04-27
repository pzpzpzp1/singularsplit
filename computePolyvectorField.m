
function [x,u,v,V,T,data,lambda,fval,C,exitflag] = computePolyvectorField(V,T,optol,x0,permuteEdgeVinds)
    if nargin==0
        %% load data
        meshname = 'data/polar_wedge_no_singularities.obj';
        
        [V,T,UV,TF,N,NF] = readOBJ(meshname);
        
        [u,v] = loadFrameText('frames.txt'); 
        x0 = [reshape(u',[],1); reshape(v',[],1)];
        
        optol = 1e-7;
        permuteEdgeVinds = [1061 1062];
    end
    data = getMeshData(V,T);
    constol = 1e-14;
    xr = reshape(x0,2,[],2);
    u = xr(:,:,1)';
    v = xr(:,:,2)';

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
    if numel(permuteEdgeVinds) ~= 0
        eind = find(all(sort(data.edges,2) == sort(permuteEdgeVinds,2),2));
        triinds = data.edges2triangles(eind,:);
        t1 = data.triangles(triinds(1),:);
        t2 = data.triangles(triinds(2),:);
    
        intedgeind = eind;
        if numel(strfind([t1 t1(1)],permuteEdgeVinds))~=0
            highertriind = triinds(1); 
            lowertriind = triinds(2); 
        else
            highertriind = triinds(2);
            lowertriind = triinds(1); 
        end
    
        nT = data.numTriangles;
        ulower = (lowertriind-1)*2+[1:2];
        uhigher = (highertriind-1)*2+[1:2];
        vlower = (lowertriind-1)*2+[1:2] + 2*nT;
        vhigher = (highertriind-1)*2+[1:2] + 2*nT;
        
        evec = data.vertices(permuteEdgeVinds(2),1:2)-data.vertices(permuteEdgeVinds(1),1:2);
        evec=evec/norm(evec);
        edgedotMatrix(intedgeind,:)=0;
        edgedotMatrix(intedgeind+nie,:)=0;
        % correct permutation. decrease energy
        edgedotMatrix(intedgeind,ulower)=evec;
        edgedotMatrix(intedgeind,vhigher)=evec;
        edgedotMatrix(intedgeind+nie,vlower)=evec;
        edgedotMatrix(intedgeind+nie,uhigher)=-evec;
    end

    uv_int = edgedotMatrix*x0;
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
    
    constraintGradientMatrix = [edgedotMatrix;boundaryNormalMatrix];
    C = constraintGradientMatrix;
    options = optimoptions('fmincon','Display','Iter','CheckGradients',false,'SpecifyObjectiveGradient',true,'Algorithm','interior-point','HessianApproximation','lbfgs',...
        'OptimalityTolerance',optol,'ConstraintTolerance',constol);
    func = @symmetricDirichlet;
    nareas = data.triangleAreas/sum(data.triangleAreas);
    fun = @(x) obfun_wrapper(x,nareas,func);
    [x,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,x0,[],[],constraintGradientMatrix,zeros(size(constraintGradientMatrix,1),1),[],[],[],options);
    
    %% vis
    %{

    % update u,v
    J = permute(reshape(x,2,[],2),[1 3 2]);
    u = reshape(J(:,1,:),2,[])';
    v = reshape(J(:,2,:),2,[])';
    
    close all;
    TV = data.triangleBarycenters;
    figure; hold all; axis equal; rotate3d on;
    patch('faces',T,'Vertices',V,'facecolor','blue','edgecolor','black')
    patch('faces',data.edges(intedgeind,[1 2 1]),'Vertices',V,'facecolor','blue','edgecolor','cyan','linewidth',3)
    quiver([TV(:,1);TV(:,1)],[TV(:,2);TV(:,2)],[u(:,1);v(:,1)],[u(:,2);v(:,2)],'g','ShowArrowHead','off')

    iu = integrateVecfield(data,u);
    iv = integrateVecfield(data,v);
    newE = [iu iv];
    newE = newE-mean(newE);

    figure; hold all; axis equal; rotate3d on;
    patch('faces',T,'Vertices',newE,'facecolor','blue','edgecolor','black')
    scatter(newE(1,1),newE(1,2),'r','filled')
    scatter(newE(end,1),newE(end,2),'g','filled')
    %}
end

