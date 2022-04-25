function [U,d,V,S] = compute_signed_SVD_batch(F)
%This function computes the signed SVD according to the paper
%Invertible Finite Elements for Robust Simulation of Large Deformation, by
%Irving et al.
%Inputs:    
%       F: nxnxN Jacobian matrix. n=3 for tetrahedra, n=2 for triangle.
%Outputs:
%       U: nxnxN rotation matrix (detU = 1)
%       d: nx1xN singular value matrix
%       V: nxnxN rotation amtrix (detU = 1)
%       S: nxnxN diagonal matrix of signed singular values. only returned
%       if nargout = 4.
n = size(F,3);
dim = size(F,1);
%[U,S,V] = batchSVD3x3Eigen(F);
% not sure why, but have to do this.
% when calling batchop svd from a CPU array, 
% the output V is actually V'. So, to get the true V,
% need to run this below.
% when calling from a GPU, then we get true V.
[U,S,V] = batchop('svd',F,3);
if(isa(V,'double'))
     V = permute(V,[2 1 3]);
end
diags = reshape(bsxfun(@plus,(0:n-1)*9,(1:4:9).'),3,1,[]);
d = squeeze(S(diags));
% assume positive eigenvalues...
%S = abs(S);
% sort 
% [d,ind] = sort(S,'descend');
% 
% % how to reorder the columns of each V, per page. Would also need to do
% this for U?
% n = size(V,3);
% p = [ones(n,1) (1:n)';2*ones(n,1) (1:n)';3*ones(n,1) (1:n)']; 
% prhs = [reshape(squeeze(ind)',[],1), p(:,2)];
% k = [ones(3*n,1);2*ones(3*n,1);3*ones(3*n,1)];
% indlhs = sub2ind(size(V),k,repmat(p(:,1),3,1),repmat(p(:,2),3,1));
% indrhs = sub2ind(size(V),k,repmat(prhs(:,1),3,1),repmat(prhs(:,2),3,1));
% V(indlhs) = V(indrhs);

%reorder the singular values
% n = size(S,2);
% p = [(1:n)';(1:n)';(1:n)']; 
% prhs = reshape(squeeze(ind)',[],1);
% k = [ones(n,1);2*ones(n,1);3*ones(n,1)];
% indlhs = sub2ind(size(S),k,p);
% indrhs = sub2ind(size(S),prhs,p);
% S(indlhs) = S(indrhs);
% turn S into a 3x3xN tensor



% do proper negations here. Ensure detU=detV=1
detV = det3x3(V);
ind = detV < 0;
iF = find(ind);
V(:,3,iF) = -1 * V(:,3,iF);
d(3,iF) = -1*d(3,iF);
% negate the singular value if detU < 0.
detU = det3x3(U);
ind = detU < 0;
iF = find(ind);
U(:,3,iF) = -1 * U(:,3,iF);
d(3,iF) = -1*d(3,iF);
if(nargout == 4)
    S(3,3,iF) = d(3,iF);
end
end

