function [u,v] = loadFrameText(filename)
    if nargin==0
        filename='frames.txt';
    end

    fid = fopen(filename,'r');
    A = fscanf(fid,'%f %f %f\n');
    A = reshape(A,3,2,[]);
    fclose(fid);

    f = A(1:2,:,:);
    u=permute(f(:,1,:),[1 3 2])';
    v=permute(f(:,2,:),[1 3 2])';
end