function [Q] = GetQ(M)
% Return the Q factor of the polar decomposition of M
[U,S,V ] = svd(M);
Q = U*V';
return;
