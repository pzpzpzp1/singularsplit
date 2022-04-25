function [E, grad] = compute_f(J,f,varargin)
%Function computes...
%Inputs:
%       J: Jacobian marix, 3x3xN
%       f: pointer to function that operates on a matrix. Input to f is a
%       set of singular values, 3xN. f outputs a vector 1xN, and the
%       gradient with respect to the singular values, 3xN.
%       varargin:
%           signed_svd: flag if wanting to use the signed SVD. Boolean,
%           true or false.
%           alpha: scalar, what to weight the energy and gradients by. Default
%               is 1
%           weights: weight for the contribution of each tetrahedron. By
%           default all 1. size is 1xN
%           use_GPU: whether to use the GPU or not
%Outputs:
        % E: scalar, energy
        % grad: 3x3xN: gradient with respect to the Jacobian matrix.
N = size(J,3);
dim = size(J,1);
use_GPU = false;
if(nargin == 2)
    alpha = 1;
    weights = ones(1,N);
    signed_SVD = true;
    use_GPU = false;
elseif( nargin == 3)
    signed_SVD = varargin{1};
    alpha = 1;
    weights = ones(1,N);
elseif(nargin == 4)
    signed_SVD = varargin{1};
    alpha = varargin{2};
    weights = ones(1,N);
elseif(nargin == 5)
    signed_SVD = varargin{1};
    alpha = varargin{2};
    weights = varargin{3};
elseif(nargin == 6)
    signed_SVD = varargin{1};
    alpha = varargin{2};
    weights = varargin{3};
    use_GPU = varargin{4};

end


% S = zeros(dim,N);
% U = zeros(dim,dim,N);
% V = zeros(dim,dim,N);
% if(isa(J,'gpuArray'))
%    S = gpuArray(single(S));
%    U = gpuArray(single(U));
%    V = gpuArray(single(V));
% end
%J = double(J);
[U,S,V] = compute_signed_SVD_batch(J);
% for i = 1 : size(J,3)
%     if(signed_SVD)
%         [u,s,v] = compute_signed_SVD(J(:,:,i));
%     else
%         [u,s,v] = svd(J(:,:,i));
%     end
%     s = diag(s);
%     S(:,i) = s;
%     U(:,:,i) = u;
%     V(:,:,i) = v;
%  end
%J = single(J);
% call to f
[E,grad_f] = f(S);
E = alpha*E.*weights;grad_f = alpha*weights.*grad_f;

if nargout < 2
    return;
end
% compute the gradient. use CPU or GPU
if(~use_GPU)
%     grad = zeros(dim, dim, N);
%     for i = 1: dim
%         for j = 1: dim
%             for k = 1 : dim
%                 grad(i,j,:) = squeeze(grad(i,j,:)) + squeeze(U(i,k,:)).*squeeze(grad_f(k,:))'.*squeeze(V(j,k,:));
%             end
%         end
%     end
    grad = zeros(dim, dim, N);
    diags = reshape(bsxfun(@plus,(0:N-1)*9,(1:4:9).'),3,1,[]);
    grad(diags) = grad_f;
    grad = reshape(grad,size(J,1),size(J,2),[]);
    grad = pagemtimes(U,grad);
    grad = pagemtimes(grad,'none',V,'transpose');
else
    grad = gpuArray(zeros(dim, dim, N));
    diags = reshape(bsxfun(@plus,(0:N-1)*9,(1:4:9).'),3,1,[]);
    grad(diags) = grad_f;
    grad = reshape(grad,size(J,1),size(J,2),[]);
    Vt = pagefun(@transpose, V);
    grad = pagefun(@mtimes, U, grad);
    grad = pagefun(@mtimes, grad, Vt);
end
end


