function sp = EstimSprd(G, GradGX, Xs)
% function EstimSprd estimates the spectral radius of the linear operator 
% corresponding the aligned NEPv at Xs.
% INPUT:
%   G(X)        - coefficient matrix of the aligned NEPv;
%   GradGX(E)   - directional derivative of G at Xs, i.e., GradGX(E) = gradG(Xs, E);
%   Xs          - target solution of the NEPv.

n = size(Xs,1);
k = size(Xs,2);
% eigenvalue decomposition
Gx = G(Xs);  % smallest to largest eigenvalues
[Vx, Dx] = eig(Gx); 
[~, idx] = sort(real(diag(Dx)), 'descend');
Vx = Vx(:,idx); Dx = Dx(idx,idx);

% construct the linear operator 
lam = real(diag(Dx)); % may remove real
lam1 = lam(1:k); 
lam2 = lam(k+1:end);
Ss = -1./(lam2 - lam1');
Vs = Vx(:,1:k); Vst = Vx(:,k+1:end);
Op = @(X) Ss.*( Vst'*GradGX(Vs, Vst*X)*Vs );
    
% compute spectral radius of Op by eigs
Afun = @(x) reshape(Op(reshape(x, n-k, k)),[],1); 
lammax = eigs(Afun, (n-k)*k, 1, 'largestabs');
sp = abs(lammax);

return 
% END OF EstimSprd


