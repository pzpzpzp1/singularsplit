function [E, grad_f] = AMIPS_energy(S, varargin)
% Function computes the ARAP energy and gradient with respect to the
% singular values.
%Input:
%       S: DxN matrix of singular values. D is 2 for surface and 3 for
%       volume. N is the number of tetrahedra or triangles.
%       varargin: weights. 1xN scalar to weight each SVD contribution. If 
%                    none passed, then assume equal weight.
%Output:
%       E: 1xN ARAP energy per tetrahedron
%       grad_f: DxN: gradient w.r.t. singular values
% AMIPS: 

E = sum(log(S).^2,1);
grad_f = 2*log(S)./S;

end