function [Q, A] = GetPolar(M)
% Polar decomposition of M
[U,S,V] = svd(M,0);
Q = U*V';
A = V*S*V';
return
