close all;
rad1 = .05;
rad2 = 1;
t2 = linspace(0,2*pi,4); t2=t2(1:3);
t1 = t2 + t2(2)/2;
V1 = rad1*[cos(t1') sin(t1') t1'*0];
V2 = rad2*[cos(t2') sin(t2') t2'*0];
V = [V1; V2];
T = [1 2 5; 3 1 4; 2 3 6; 2 6 5; 1 5 4; 3 4 6; 1 3 2];

V0=V;T0=T;
writeOBJ('triSeam0.obj',V0,T0)

for i=1:3
    E = edges(T);
    newV = (V(E(:,1),:)+V(E(:,2),:))/2;
    
    [Vt,T,SS,J] = loop(V,T,1);
    V=[V;newV];
        
    figure; axis equal; patch('vertices',V,'faces',T,'facecolor','g');
    writeOBJ(sprintf('triSeam%d.obj',i),V,T);
end


