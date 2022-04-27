
clear all; close all; 

%% load data
optol = 1e-6;
[x,u,v,V,T,data,lambda,fval0] = analyzePolar(1e-9,1e-9); close all;
data0 = getMeshData(V,T); data=data0; V0=V; T0=T;
areaFromV123 = @(v123) dot(cross(v123(2,:)-v123(1,:),v123(3,:)-v123(1,:)),[0 0 1])/2;
syms alpha real;

boundaryVerts = unique(data.edges(data.isBoundaryEdge,:));
innerEdges = find(all(~ismember(data.edges, boundaryVerts),2))'; niie = numel(innerEdges);
nie=sum(~data.isBoundaryEdge);
adiff_part1 = nan(niie,1);
adiff_part2 = nan(niie,1);
fdiffs = nan(niie,1);
order = randperm(niie); order = 1:niie;
order(order==2807)=[]; order = [2807 order];
for j = 1:niie
    
    intedgeind  = innerEdges(order(j));
    
    % flip edges
%     p1 = data.edges(intedgeind,2);
%     p2 = data.edges(intedgeind,1);
    p1 = data.edges(intedgeind,1);
    p2 = data.edges(intedgeind,2);
    pdir = (V0(p2,1:2)-V0(p1,1:2))';
    triinds = data.edges2triangles(intedgeind,:);
    tri1 = data0.triangles(triinds(1),:);
    if numel(strfind([tri1 tri1(1)],[p1 p2]))~=0
        t1 = triinds(1); 
        t3 = triinds(2); 
    else
        t1 = triinds(2);
        t3 = triinds(1); 
    end

    ehat = pdir/norm(pdir);
    p1edges = find(any(ismember(data.edges,p1),2));
    t1edges = data.triangles2edges(t1,:);
    t3edges = data.triangles2edges(t3,:);
    e1 = setdiff(intersect(t1edges, p1edges), intedgeind);
    e2 = setdiff(intersect(t3edges, p1edges), intedgeind);
    t0 = setdiff(data.edges2triangles(e1,:),t1);
    t2 = setdiff(data.edges2triangles(e2,:),t3);
    edgev = data.vertices(data.edges(:,2),1:2)-data.vertices(data.edges(:,1),1:2);
    e1v = edgev(e1,:)';
    e2v = edgev(e2,:)';
    l1flip=false;l2flip=false;
    if data.edges(e1,2)~=p1; e1v = -e1v; l1flip=~l1flip; end;
    if data.edges(e2,2)~=p1; e2v = -e2v; l2flip=~l1flip; end;
    p4 = setdiff(data.edges(e1,:),p1);
    p6 = setdiff(data.edges(e2,:),p1);
    u0 = u(t0,:)';
    v0 = v(t0,:)';
    u1 = u(t1,:)';
    v1 = v(t1,:)';
    u2 = u(t2,:)';
    v2 = v(t2,:)';
    u3 = u(t3,:)';
    v3 = v(t3,:)';

    %% build perturbation frames
    eps = 1e-4; % 1e-5, 1e-4 seem ok. too small and hit optimality tolerance too quick obtaining inaccurate answer. too large and it never finishes solving.
    z = zeros(1,2);
    Am = [e1v' z z z; z z e1v' z; z e2v' z z; z z z e2v'; z ehat' ehat' z; -ehat' z z ehat'];
    Bm = [dot(u0,e1v); dot(v0,e1v); dot(u2,e2v); dot(v2,e2v); 0; 0];
    e3p = e1v + eps * ehat;
    e4p = e2v + eps * ehat;
    e1p = e1v;    e2p = e2v;
    ehp = ehat;
    A = [e1p' z z z; z z e1p' z; z e2p' z z; z z z e2p'; z ehp' ehp' z; -ehp' z z ehp';...
        e3p' z z z; z z e3p' z; z e4p' z z; z z z e4p'];
    B = [dot(u0,e1p); dot(v0,e1p); dot(u2,e2p); dot(v2,e2p); 0; 0;...
        dot(u1,e3p); dot(v1,e3p); dot(u3,e4p); dot(v3,e4p); ];
    xperturb = A\B;

    %% build new mesh topology
    newv = data.vertices(p1,:) + eps*[ehat' 0];
    newvind = data.numVertices+1;
    Vnew = [V0; newv]; 
    Tnew = T;     Tnew(t1,:)=[p2 p4  newvind];    Tnew(t3,:)=[p6 p2 newvind];    Tnew = [Tnew; p4 p1 newvind ; p1 p6 newvind ];
    permuteEdgeVinds = [p1 newvind];
    x0 = reshape([reshape(x,[],2); reshape(xperturb,[],2)],[],1);
    [newX,newU,newV,~,~,newData,newLambda,fval1,C,exitflag(j)] = computePolyvectorField(Vnew,Tnew,optol,x0,permuteEdgeVinds);
    fdiffs(j) = (fval1-fval0)/eps; 
    % [intedgeind=201: -.17 | gt =  -0.47667     -0.86683     -0.50003   -0.0069196      0.61278      -4.1552     -0.40546      -3.4958] 
    % [intedgeind=3044: -34.112 | gt = -31.896       4.2846      -31.552      -8.3002       35.544      -6.0746      -28.372      -8.4041] 

    xr = reshape(newX,[],2);
    gt = reshape(xr(end-3:end,:),[],1);
    
    %TODO: CHECK LAMBDAS. TRY OPTIMIZING dL/dalpha w.r.t. virtual frames

    %% set primal virtual frames (this is the key part to test)
    % lopital's rule area ratio
    %{
    A4_div_A5 = abs(norm(cross([e1v' 0],[ehat' 0]))/norm(cross([e2v' 0],[ehat' 0])));
    fun = @(x)obfun_wrapper(x,[A4_div_A5 1]);
    options = optimoptions('fmincon','Display','Iter','CheckGradients',false,'SpecifyObjectiveGradient',true,'Algorithm','interior-point','HessianApproximation','lbfgs','SpecifyConstraintGradient',true,...
        'OptimalityTolerance',optol);
    xl = fmincon(fun,x0,[],[],Am,Bm,[],[],[],options); % use min energy frames
    %}
    xl = gt; % use gt frames
%     xl = x0; % use perturbation frames
    u4 = xl(1:2);
    u5 = xl(3:4);
    v4 = xl(5:6);
    v5 = xl(7:8);
    
    uv4 = [u4 v4];
    uv5 = [u5 v5];
    E4 = symmetricDirichlet(uv4,1);
    E5 = symmetricDirichlet(uv5,1);
    E45c = [E4 E5];
    
    %% Energy envelope gradient
    t12inds = [t1, t3];
    nt = numel(t12inds);
    for i=1:nt
        ti = t12inds(i);
        uvi = [u(ti,:)' v(ti,:)'];
        Eic(i) = symmetricDirichlet(uvi,1);
    end
    V = data.vertices;
    newp = data.numVertices+1;
    V(newp,:)=V(p1,:);
    vx = V(p1,:) + [ehat' 0] * alpha;
    t12 = data.triangles(t12inds,:);
    t12(t12==p1)=newp;
    tvis = [t12 ; p4 p1 newp; p1 p6 newp];
    for i=1:size(tvis,1)
        t123 = tvis(i,:);
        vind = find(t123 == newp);
        v123 = sym(V(t123,:));
        v123(vind,:)=vx;
        Ai(i) = areaFromV123(v123);
    end
    dAda = [double(diff(Ai,alpha))]/sum(data0.triangleAreas);
    adiff_part1(j) = dot(dAda, [Eic E45c]);

    %% get constraint envelope gradient
    if data.edges2triangles(e1,1)~=t1
        l1flip=~l1flip;
    end
    if data.edges2triangles(e2,1)~=t3
        l2flip=~l2flip;
    end
    l1u = lambda.eqlin(e1);
    l1v = lambda.eqlin(e1+nie);
    l2u = lambda.eqlin(e2);
    l2v = lambda.eqlin(e2+nie);
    if l1flip;
        l1u = -l1u;
        l1v = -l1v;
    end
    if l2flip;
        l2u = -l2u;
        l2v = -l2v;
    end
    adiff_part2(j) = l1u*dot(u1-u4,ehat) + l1v*dot(v1-v4,ehat) + l2u*dot(u3-u5,ehat) + l2v*dot(v3-v5,ehat);

    %% get xperturb that minimizes env grad: dLda
    %{
    fun = @(x)obfun_wrapper(x, 1);
    xo=Am\Bm;N=null(Am); dt=1e-1;
    uv45 = gt;
    for t=1:500
        % project onto nullspace
        n=N\(uv45-xo);
        uv45 = xo+N*n;
        u4 = uv45(1:2);
        u5 = uv45(3:4);
        v4 = uv45(5:6);
        v5 = uv45(7:8);

        % get d_dLdalpha_d_uv45
        [~,dEduv4] = fun([u4;v4]);
        [~,dEduv5] = fun([u5;v5]);
        duv4 = dEduv4*dAda(3) - (l1u*[ehat; z'] + l1v*[z'; ehat]);
        duv5 = dEduv5*dAda(4) - (l2u*[ehat; z'] + l2v*[z'; ehat]);
        duv45 = reshape([reshape(duv4,2,2); reshape(duv5,2,2); ],[],1);
        
        % project onto tangent space
        duv45 = N*(N\duv45);

        % step
        uv45 = uv45 - dt * duv45;
        norm(uv45-gt)
    end
    uv45_min_dLda = uv45;
    %}

    %% break at set percentage
%     if j/niie > .03758
%         break;
%     end
end
adiffs = adiff_part1 + adiff_part2;
adiffs2 = adiff_part1 - adiff_part2;

edgescore = -fdiffs;
%{
edgescore = -fdiffs;
edgescore = -adiffs;
edgescore = adiff_part1;
edgescore = -adiff_part2;
edgescore = adiffs2;
%}
edgescore(abs(edgescore)>prctile(abs(edgescore),99.99))=nan;
edgescore(edgescore<prctile(edgescore,.01))=nan;
figure; hold all; rotate3d on; title('score'); axis equal;
patch('faces',T,'Vertices',V,'facecolor','blue','edgecolor','none')
xx = reshape(data.vertices(data.edges(innerEdges(order),:)',1),2,[]);
yy = reshape(data.vertices(data.edges(innerEdges(order),:)',2),2,[]);
xx(3,:)=nan;yy(3,:)=nan;xx = xx(:);yy = yy(:);
patch(xx,yy,yy*0,repelem(edgescore,3,1),'edgecolor','interp','linewidth',2)
colorbar;colormap(inferno);
[~,maxind]=max(edgescore);
patch('vertices',data.vertices,'faces',data.edges(innerEdges(order(maxind)),[1 2 1]),'edgecolor','green','linewidth',3)

[~,perm] = sort(fdiffs);
figure; 
clf; 
subplot(2,1,1); hold all;
plot(adiffs(perm)*130,'r')
plot(-adiffs2(perm)*100,'b')
plot(fdiffs(perm),'g','linewidth',2)
subplot(2,1,2); hold all;
plot(adiff_part1(perm)*500,'b')
plot(adiff_part2(perm)*-80,'g')
