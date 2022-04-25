clear all; close all; 

%% load data
meshname = 'data/polar_wedge_no_singularities.obj';
[x0,u0,v0,V0,T0,data0] = analyzePolar();
close all;
u=u0;v=v0;x=x0;
data=data0;V=V0;T=T0;

boundaryVertices = unique(data.edges(data.isBoundaryEdge,:));
fullyIntEdges = all(~ismember(data.edges, boundaryVertices),2);
% patch('faces',T,'Vertices',V,'facecolor','blue','edgecolor','black');
% patch('vertices',V,'faces',data.edges(fullyIntEdges,[1 2 1]),'edgecolor','cyan','linewidth',2);
halfedges = [data0.edges(fullyIntEdges,:);fliplr(data0.edges(fullyIntEdges,:))];
    
for heind=1:size(halfedges,1)
    he = halfedges(heind,:);
    eind = find(all(data0.edges == he,2) | all(data0.edges == [he(2) he(1)],2));
    triinds = data.edges2triangles(eind,:);
    t1 = data0.triangles(triinds(1),:);
    t2 = data0.triangles(triinds(2),:);

%     perturbfactors = fliplr([0.05 0.04 0.03 0.02 0.01]);
    perturbfactors = fliplr([.03 0.02 0.01]);
    N = numel(perturbfactors);
    x = x0;
    for i=1:N
        perturbfactor=perturbfactors(i);
        
        %% select perturb edge. displace vertex
        intedgeind = eind;
        if numel(strfind([t1 t1(1)],he))~=0
            highertriind = triinds(1); 
            lowertriind = triinds(2); 
        else
            highertriind = triinds(2);
            lowertriind = triinds(1); 
        end
        
        
        v1 = he(1); v2 = he(2);
        V=V0;
        pdir = V(v2,1:2) - V(v1,1:2);
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
        u_norm_ind = abs(dot(u0(be2btri,:),be,2)) < 1e-6;
        v_norm_ind = abs(dot(v0(be2btri,:),be,2)) < 1e-6;
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
        options = optimoptions('fmincon','Display','Iter','CheckGradients',false,'SpecifyObjectiveGradient',true,'Algorithm','interior-point','HessianApproximation','lbfgs');
        [x,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,x,[],[],constraintGradientMatrix,zeros(size(constraintGradientMatrix,1),1),[],[],[],options);
        xs{heind,i} = x;
    
        J = permute(reshape(x,2,[],2),[1 3 2]);
        u = reshape(J(:,1,:),2,[])';
        v = reshape(J(:,2,:),2,[])';
        us{heind,i}=u;
        vs{heind,i}=v;
        fvals(heind,i)=fval;
        dets = cross([u u(:,1)*0],[v v(:,1)*0]); 
        mindets(heind,i) = min(dets(:,3));
    end

    xx=perturbfactors*norm(pdir);yy=fvals(heind,:);
    pfit{heind} = polyfit(xx,yy,2);
    xxs(heind,:)=xx;
    yys(heind,:)=yy;
    %{
    figure; plot(xx,yy,'.-');    title(num2str(mindets(heind,:)))
    %}
end

%% visualize. 
figure; hold all;
plot(xxs',yys','r.-')
xline(0,'k')
t = linspace(0,max(perturbfactors)*max(data.edgeLengths),100);
for i=1:size(xxs,1)
    yy2 = polyval(pfit{i},t);
    plot(t,yy2,'g-')

    k = 2;    xx = xxs(i,1:k);    yy = yys(i,1:k);
    pfitk = polyfit(xx,yy,2);
    yyk = polyval(pfitk,t);
    plot(t,yyk,'c-')

    slpk(i) = pfitk(2);
end

colors = slpk; colors = colors - min(colors) + .01; colors = colors/max(colors);
for ii=1:numel(pfit); slp(ii)=pfit{ii}(2); end
figure; hold all; axis equal;
patch('faces',T,'vertices',V,'facecolor','green','edgecolor','none')
fullyIntEdges2 = data0.edges(fullyIntEdges,:);
for ii=1:numel(slp)
    ves = data0.vertices(fullyIntEdges2(ii,:),:);
    plot(ves(:,1),ves(:,2),'color',[colors(ii) 0 0],'linewidth',2)
%     patch('faces',data0.edges(),'vertices',V,'facecolor','green')
end


