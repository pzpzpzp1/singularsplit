
clear all; 

% variables
syms v0lx v0ly v0rx v0ry v1x v1y v2x v2y v3x v3y v4x v4y v5x v5y v6x v6y real;
v0l = [v0lx v0ly]';v0r = [v0rx v0ry]';v1 = [v1x v1y]';v2 = [v2x v2y]';v3 = [v3x v3y]';v4 = [v4x v4y]';v5 = [v5x v5y]';v6 = [v6x v6y]';
e1l = v1-v0l;e1r = v1-v0r;e4l = v4-v0l;e4r = v4-v0r;e2 = v2-v0l;e3 = v3-v0l;e5 = v5-v0r;e6 = v6-v0r;elr = v0l-v0r;
syms u1x u1y w1x w1y u4x u4y w4x w4y real;
u1 = [u1x u1y]';w1 = [w1x w1y]';u4 = [u4x u4y]';w4 = [w4x w4y]';
syms u12x u12y u23x u23y u34x u34y real;syms u45x u45y u56x u56y u61x u61y real;syms w12x w12y w23x w23y w34x w34y real;syms w45x w45y w56x w56y w61x w61y real;
u12 = [u12x u12y]';w12 = [w12x w12y]';u23 = [u23x u23y]';w23 = [w23x w23y]';u34 = [u34x u34y]';w34 = [w34x w34y]';u45 = [u45x u45y]';w45 = [w45x w45y]';u56 = [u56x u56y]';w56 = [w56x w56y]';u61 = [u61x u61y]';w61 = [w61x w61y]';
syms v0rx2 v0lx2 v0ry2 v0ly2 real; % augment elr variables

vertvars = [v0l; v0r];
framevars = [u1;w1;u4;w4;u12;w12;u23;w23;u34;w34;u45;w45;u56;w56;u61;w61];
realframevars = [u12;w12;u23;w23;u34;w34;u45;w45;u56;w56;u61;w61];

% contraints
g2u = dot(u12-u23,e2); g2w = dot(w12-w23,e2);
g3u = dot(u34-u23,e3); g3w = dot(w34-w23,e3);
g5u = dot(u45-u56,e5); g5w = dot(w45-w56,e5);
g6u = dot(u61-u56,e6); g6w = dot(w61-w56,e6);

g1u = dot(u61-u12,e1r); g1w = dot(w61-w12,e1r);
g4u = dot(u34-u45,e4r); g4w = dot(w34-w45,e4r);

g1lu = dot(e1l,u1-u12);g1lw = dot(e1l,w1-w12);
g1ru = dot(e1r,u1-u61);g1rw = dot(e1r,w1-w61);
g4lu = dot(e4l,u4-u34);g4lw = dot(e4l,w4-w34);
g4ru = dot(e4r,u4-u45);g4rw = dot(e4r,w4-w45);
g41uw = dot(elr,w4+u1); g41wu = dot(elr,w1-u4); % permutation on 0 len edge.
g = [g2u g2w g3u g3w g5u g5w g6u g6w g1lu g1lw g1ru g1rw g4lu g4lw g4ru g4rw g41uw g41wu]';
g_collapsed = [g2u g2w g3u g3w g5u g5w g6u g6w g1u g1w g4u g4w];

% constraint gradients
for i=1:numel(vertvars)
    dgdv(:,i)=diff(g,vertvars(i));
end
for i=1:numel(framevars)
    dgdf(:,i)=diff(g,framevars(i));
end
for i=1:numel(realframevars)
    dgcdf(:,i)=diff(g_collapsed,realframevars(i));
end

% constraints involving virtual frames gradients. 
vframes = [u1x u1y w1x w1y u4x u4y w4x w4y];
vframeCons = g(has(g,vframes));
% constraint augment: decide vertex split orientation
vframeCons(end+1,:) = [(v0lx2 - v0rx2)*(u1x + w4x) + (v0ly2 - v0ry2)*(u1y + w4y)];
vframeCons(end+1,:) = [- (v0lx2 - v0rx2)*(u4x - w1x) - (v0ly2 - v0ry2)*(u4y - w1y)];
for i=1:numel(vframes)
    dgdv_vframes(:,i)=diff(vframeCons,vframes(i));
end


%% pick feasible point
% non generic initial triangulation. inner 2 verts are coincident at 0.
% frame field is constant. virtual edge is unspecified. 
v0lx=0;v0ly=0;v0rx=0;v0ry=0;
v1x=cos(0*2*pi/6+pi/2)+randn*.08; v1y=sin(0*2*pi/6+pi/2)+randn*.08;
v2x=cos(1*2*pi/6+pi/2)+randn*.08; v2y=sin(1*2*pi/6+pi/2)+randn*.08;
v3x=cos(2*2*pi/6+pi/2)+randn*.08; v3y=sin(2*2*pi/6+pi/2)+randn*.08;
v4x=cos(3*2*pi/6+pi/2)+randn*.08; v4y=sin(3*2*pi/6+pi/2)+randn*.08;
v5x=cos(4*2*pi/6+pi/2)+randn*.08; v5y=sin(4*2*pi/6+pi/2)+randn*.08;
v6x=cos(5*2*pi/6+pi/2)+randn*.08; v6y=sin(5*2*pi/6+pi/2)+randn*.08;

% choose constant generic field, then make it a non-constant generic, integrable field.
a=1+randn*.08;b=randn*.08;c=randn*.08;d=1+randn*.08;
u12x=a; u12y=b; u23x=a; u23y=b; u34x=a; u34y=b;
u45x=a; u45y=b; u56x=a; u56y=b; u61x=a; u61y=b;
w12x=c; w12y=d; w23x=c; w23y=d; w34x=c; w34y=d;
w45x=c; w45y=d; w56x=c; w56y=d; w61x=c; w61y=d;

nsp=null(double(subs(dgcdf)));
realframevars_c = double(subs(realframevars)) + nsp*randn(size(nsp,2),1)/20;
u12x=realframevars_c(1);u12y=realframevars_c(2);w12x=realframevars_c(3);w12y=realframevars_c(4);
u23x=realframevars_c(5);u23y=realframevars_c(6);w23x=realframevars_c(7);w23y=realframevars_c(8);
u34x=realframevars_c(9);u34y=realframevars_c(10);w34x=realframevars_c(11);w34y=realframevars_c(12);
u45x=realframevars_c(13);u45y=realframevars_c(14);w45x=realframevars_c(15);w45y=realframevars_c(16);
u56x=realframevars_c(17);u56y=realframevars_c(18);w56x=realframevars_c(19);w56y=realframevars_c(20);
u61x=realframevars_c(21);u61y=realframevars_c(22);w61x=realframevars_c(23);w61y=realframevars_c(24);
double(subs(g_collapsed))

% augment elr to pick a free direction
randedge = abs(randn(2,1))*.1;
v0rx2 = randedge(1);
v0lx2 = -randedge(1);
v0ry2 = randedge(2);
v0ly2 = -randedge(2);

dgdv_vframes_c = double(subs(dgdv_vframes));
rhs = double(subs(simplify(vframeCons - dgdv_vframes*vframes')));
fullnsp = null(dgdv_vframes_c(1:end-2,:))
nsp = null(dgdv_vframes_c);
fprintf('Addition of virtual edge constraint decreased virtual frame nullspace by %d dimension.',rank(fullnsp)-rank(nsp));

%% choose virtual frame
% vframe with augment edge with nsp
vframesc = (dgdv_vframes_c\-rhs) + nsp*randn(size(nsp,2),1); 
u1x = vframesc(1);u1y = vframesc(2);w1x = vframesc(3);w1y = vframesc(4);
u4x = vframesc(5);u4y = vframesc(6);w4x = vframesc(7);w4y = vframesc(8);
[double(subs(g)); double(subs(vframeCons((end-1):end)))]'

% check amips. non 0 dets means amips energy wont blow up. 
J1 = (double(subs([u1 w1])));
J4 = (double(subs([u4 w4])));
amips = @(J) trace(J'*J)/det(J) + det(J) + 1/det(J);
amipsE = [amips(J1) amips(J4)]

%% check existence of tangent space
dgdvc = double(subs(dgdv));
dgdfc = double(subs(dgdf));
x = dgdfc\dgdvc;

r_resid = rank(dgdfc*x-dgdvc);
r_dgdv = rank(dgdvc);
r_dgdf = rank(dgdfc);
nc = size(dgdvc,1);
nf = size(dgdfc,2);
nv = size(dgdvc,2);

fprintf('There is a %d dimensional space of vertex motions that have no associated integrable frame field perturbation.\n',r_resid);
fprintf('The full space of vertex movements is %d dimensional.\n',nv);
fprintf('There is a %d dim space of valid vertex perturbations.\n',nv-r_resid);

% manual vertex perturbations
% first two columns are if verts move together. obviously that one is ok.
% latter two are if verts move apart in x or y axis. these seem ok depending on if the
% virtual frames are parallel aligned to x or y. 
vpert = [1 0 1 0; 0 1 0 1; 1 0 -1 0; 0 1 0 -1; randedge' -randedge']'; 
fpert = dgdfc\(dgdvc*vpert);
vecnorm(dgdfc*fpert - (dgdvc*vpert),2,1)

%{
Takeaway!!
The virtual frames of augmented edge produce valid virtual vertex pull-apart movmenets in that augmented edge direction! 
augmented edge constraint removes two dims from virtual frames: generically 4 -> 2.
%}

