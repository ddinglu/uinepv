function sp = EstimSprds(G, GradGXs, Xs)
% function EstimSprds is a more efficient implementation of EstimSprd
% which estimates the spectral radius of the linear operator 
% corresponding the aligned NEPv at Xs.
% INPUT:
%   G(X)        - coefficient matrix of the aligned NEPv;
%   GradGX(E)   - directional derivative of G at Xs, i.e., GradGX(E) = gradG(Xs, E);
%   Xs          - target solution of the NEPv.

n = size(Xs,1);
k = size(Xs,2);
% eigenvalue decomposition
Gx = G(Xs);  % smallest to largest eigenvalues
Vs = Xs;
lam1 = real(diag(Vs'*Gx*Vs));

Vst0 = null(Vs');
Gst = Vst0'*Gx*Vst0;
[Vx, Dx] = eig(Gst); 
[~, idx] = sort(real(diag(Dx)), 'descend');
Vx = Vx(:,idx); Dx = Dx(idx,idx);
Vst = Vst0 * Vx;
lam2 = diag(Dx);

% construct the linear operator 
Ss = -1./(lam2 - lam1');
Op = @(X) Ss.*( Vst'*GradGXs(Vst*X)*Vs );
    
% compute spectral radius of Op by eigs
Afun = @(x) reshape(Op(reshape(x, n-k, k)),[],1); 
lammax = eigs(Afun, (n-k)*k, 1);
sp = abs(lammax);

return 
% END OF EstimSprds


