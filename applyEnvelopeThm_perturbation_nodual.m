
clear all; close all; 

%% load data
[x,u,v,V,T,data] = analyzePolar(1e-7,1e-7); close all;
data0 = getMeshData(V,T); data=data0; V0=V; T0=T;
areaFromV123 = @(v123) dot(cross(v123(2,:)-v123(1,:),v123(3,:)-v123(1,:)),[0 0 1])/2;
syms alpha real;

boundaryVerts = unique(data.edges(data.isBoundaryEdge,:));
innerEdges = find(all(~ismember(data.edges, boundaryVerts),2))'; nie = numel(innerEdges);
adiffs = nan(nie,1);
for j = 1:nie
    intedgeind  = innerEdges(j);
    % intedgeind = randi(sum(~data.isBoundaryEdge));
    % intedgeind = 1079; % wroks ok. depends a lot on primal convergence.
    % intedgeind = 6587; % seems to work. for specific ep
    % intedgeind = 3044; % works well
    % intedgeind = 5360; % diagonal. works well
    % intedgeind = 7948; % close to boundary. doesn't work well.
    % intedgeind = 6304; % circumcentric orientation. doesn't work well/
    
    % flip edges
%     v1 = data.edges(intedgeind,2);
%     v2 = data.edges(intedgeind,1);
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
    
%     figure; hold all; axis equal; axis tight manual; title('neighborhood');
%     TV = data.triangleBarycenters;
%     patch('faces',T,'Vertices',V,'facecolor','blue','edgecolor','black');
%     patch('faces',T(highertriind,:),'Vertices',V,'facecolor','green','edgecolor','black','facealpha',.5);
%     patch('faces',T(lowertriind,:),'Vertices',V,'facecolor','red','edgecolor','black','facealpha',.5);
%     patch('vertices',V,'faces',data.edges(intedgeind,[1 2 1]),'edgecolor','cyan','linewidth',2);
%     scatter(TV(:,1),TV(:,2),'k','filled')
%     scatter(V0(v1,1),V0(v1,2),'g','filled')
%     scatter(V0(v2,1),V0(v2,2),'r','filled')
%     sc = 2;
%     bl = mean(V0(data.edges(intedgeind,:),:)) - [.1 .1 0]*sc;
%     tr = mean(V0(data.edges(intedgeind,:),:)) + [.1 .1 0]*sc;
%     xlim([bl(1) tr(1)]); ylim([bl(2) tr(2)])
    
    % J = permute(reshape(x,2,[],2),[1 3 2]);
    % uvlowerc = J(:,:,lowertriind);
    % uvhigherc = J(:,:,highertriind);
    % Ehigherc = symmetricDirichlet(uvhigherc,1);
    % Elowerc = symmetricDirichlet(uvlowerc,1);
    % Ehlc = [Ehigherc Elowerc];
    % gt = [uvhigherc(:,1); uvlowerc(:,1); uvhigherc(:,2); uvlowerc(:,2); ];
    
    %% obtain primal virtual frames (this is the key part to test)
    ehat = pdir/norm(pdir);
    v1edges = find(any(ismember(data.edges,v1),2));
    t5 = highertriind; t6 = lowertriind;
    t5edges = find(ismember(sort(data.edges,2), sort(reshape(data.triangles(t5,[1 2 2 3 3 1]),2,3)',2), 'rows'));
    t6edges = find(ismember(sort(data.edges,2), sort(reshape(data.triangles(t6,[1 2 2 3 3 1]),2,3)',2), 'rows'));
    e1 = setdiff(intersect(t5edges, v1edges), intedgeind);
    e2 = setdiff(intersect(t6edges, v1edges), intedgeind);
    t1 = setdiff(data.edges2triangles(e1,:),t5);
    t2 = setdiff(data.edges2triangles(e2,:),t6);
    edgev = data.vertices(data.edges(:,2),1:2)-data.vertices(data.edges(:,1),1:2);
    e1v = edgev(e1,:);
    e2v = edgev(e2,:);
    if data.edges(e1,2)~=v1; e1v = -e1v; end;
    if data.edges(e2,2)~=v1; e2v = -e2v; end;
    u1 = u(t1,:);
    vv1 = v(t1,:);
    u2 = u(t2,:);
    vv2 = v(t2,:);
    u5 = u(t5,:);
    v5 = v(t5,:);
    u6 = u(t6,:);
    v6 = v(t6,:);
    p3 = setdiff(data.edges(e1,:),v1);
    p4 = setdiff(data.edges(e2,:),v1);
    
    z = zeros(1,2);
    % redundant edge removed
    Am = [e1v z z z; z z e1v z; z e2v z z; z z z e2v; z ehat ehat z; -ehat z z ehat];
    Bm = [dot(u1,e1v); dot(vv1,e1v); dot(u2,e2v); dot(vv2,e2v); 0; 0];
    % redundant edge approximated
    clear ad;
    epss = 1e-8;
%     epss = 1e-7;
    for i=1:numel(epss)
        eps = epss(i);
%         e3p = e1v;
%         e4p = e2v;
%         e1p = e1v - eps * ehat;
%         e2p = e2v - eps * ehat;
        e3p = e1v + eps * ehat;
        e4p = e2v + eps * ehat;
        e1p = e1v;
        e2p = e2v;
        ehp = ehat;
        A = [e1p z z z; z z e1p z; z e2p z z; z z z e2p; z ehp ehp z; -ehp z z ehp;...
            e3p z z z; z z e3p z; z e4p z z; z z z e4p];
        B = [dot(u1,e1p); dot(vv1,e1p); dot(u2,e2p); dot(vv2,e2p); 0; 0;...
            dot(u5,e3p); dot(v5,e3p); dot(u6,e4p); dot(v6,e4p); ];
        eps2 = 1e-3; A(7:end,:)=A(7:end,:)*eps2; B(7:end)=B(7:end)*eps2;
        ad(:,i) = A\B; 
    end
%     [ttt ad]
    ads{j} = ad;
    %{
    figure; 
    for i=1:8
        subplot(8,1,i);
        loglog(epss,abs(ad(i,:)-ad(i,end)),'r.-')
    end
    figure; plot(epss,ad','r.-')
    %}

    x = ad(:,end);
    u3 = x(1:2);
    u4 = x(3:4);
    v3 = x(5:6);
    v4 = x(7:8);
    
    uvhigher2 = [u3 v3];
    uvlower2 = [u4 v4];
    Ehigher2 = symmetricDirichlet(uvhigher2,1);
    Elower2 = symmetricDirichlet(uvlower2,1);
    Ehlc = [Ehigher2 Elower2];
    
    %% get neighboring triangles
    t12inds = [highertriind, lowertriind];
    nt = numel(t12inds);
    for i=1:nt
        ti = t12inds(i);
        uvi = [u(ti,:)' v(ti,:)'];
        Eic(i) = symmetricDirichlet(uvi,1);
    end
    
    % EA_unnormalized = dot(Eic, data.triangleAreas(t12inds));
    % [Eic Ehlc] .* data.triangleAreas([t1234', lowertriind, highertriind])';
    V = data.vertices;
    newp = data.numVertices+1;
    V(newp,:)=V(v1,:);
    vx = V(v1,:) + [pdir/norm(pdir) 0] * alpha;
    t12 = data.triangles(t12inds,:);
    t12(t12==v1)=newp;
    tvis = [t12 ; p3 v1 newp; v1 p4 newp];
    for i=1:size(tvis,1)
        t123 = tvis(i,:);
        vind = find(t123 == newp);
        v123 = sym(V(t123,:));
        v123(vind,:)=vx;
        Ai(i) = areaFromV123(v123);
    end
    dAda = [double(diff(Ai,alpha))]/sum(data0.triangleAreas);
    adiff = dot(dAda, [Eic Ehlc]);
    adiffs(j) = adiff;
    j/nie
end


% lambdadgdx = full(abs(lambda.eqlin.*vecnorm(constraintGradientMatrix,2,2)));
figure; hold all; rotate3d on; title('score'); axis equal;
patch('faces',T,'Vertices',V,'facecolor','blue','edgecolor','none')
edgescore = adiffs;
% edgescore(edgescore<0)=0;
xx = reshape(data.vertices(data.edges(innerEdges,:)',1),2,[]);
yy = reshape(data.vertices(data.edges(innerEdges,:)',2),2,[]);
xx(3,:)=nan;yy(3,:)=nan;xx = xx(:);yy = yy(:);
patch(xx,yy,yy*0,repelem(edgescore,3,1),'edgecolor','interp','linewidth',2)
colorbar;colormap(inferno);
[~,maxind]=max(edgescore);
patch('vertices',data.vertices,'faces',data.edges(innerEdges(maxind),[1 2 1]),'edgecolor','green','linewidth',3)


