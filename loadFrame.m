function f = loadFrame(filename)

    fid = fopen(filename,'r');
    l1 = fgetl ( fid );
    A = fscanf(fid,'∥(%f, %f, %f)∥ = %f\n');
    A = reshape(A,4,[])';
    % vecnorm(A(:,1:3),2,2) - A(:,4)
    f = A(:,1:3) .* A(:,4);
    fclose(fid);

end